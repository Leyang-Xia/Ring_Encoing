#include "xd/search.hpp"
#include <chrono>
#include <climits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace xd {
namespace {

void check_search_deadline(const BoundedScheduleOptions &options) {
    if (options.deadline != std::chrono::steady_clock::time_point::max() &&
        std::chrono::steady_clock::now() >= options.deadline)
        throw SearchDeadlineExceeded();
}

std::uint64_t residual_pair_key(int left, int right) {
    return (std::uint64_t(std::uint32_t(left)) << 32) | std::uint32_t(right);
}
}

int combine_depth(const std::vector<int> &depths) {
    if (depths.empty()) return 0;
    std::vector<int> work = depths;
    for (int d : work) require(d >= 0, "negative arrival depth");
    int level = 0;
    while (work.size() > 1) {
        std::sort(work.begin(), work.end());
        int ready = 1 + std::max(work[0], work[1]);
        work.erase(work.begin(), work.begin() + 2);
        work.push_back(ready);
        level = std::max(level, ready);
    }
    return level;
}

std::vector<int> normalize_terms(const std::vector<int> &terms) {
    std::vector<int> out = terms;
    std::sort(out.begin(), out.end());
    std::vector<int> result;
    for (int term : out) {
        require(term >= 0, "negative term");
        if (!result.empty() && result.back() == term) result.pop_back();
        else result.push_back(term);
    }
    return result;
}

bool dominates(const std::vector<int> &a, const std::vector<int> &b) {
    if (a.size() != b.size()) return false;
    bool strict = false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] > b[i]) return false;
        strict = strict || a[i] < b[i];
    }
    return strict;
}

bool dominates(const Representation &a, const Representation &b) {
    if (a.future_support > b.future_support) return false;
    return a.terms.size() <= b.terms.size() && a.arrival_depth <= b.arrival_depth &&
           (a.terms.size() < b.terms.size() || a.arrival_depth < b.arrival_depth);
}

void pareto_insert(ParetoFrontier &frontier, Representation candidate, std::size_t max_representations) {
    candidate.terms = normalize_terms(candidate.terms);
    for (const auto &entry : frontier.entries)
        if (dominates(entry, candidate) ||
            (entry.arrival_depth == candidate.arrival_depth && entry.terms == candidate.terms)) return;
    frontier.entries.erase(std::remove_if(frontier.entries.begin(), frontier.entries.end(),
                                          [&](const Representation &entry) { return dominates(candidate, entry); }),
                           frontier.entries.end());
    frontier.entries.push_back(std::move(candidate));
    std::sort(frontier.entries.begin(), frontier.entries.end(), [](const Representation &x, const Representation &y) {
        if (x.future_support != y.future_support)
            return x.future_support > y.future_support;
        if (x.terms.size() != y.terms.size()) return x.terms.size() < y.terms.size();
        if (x.arrival_depth != y.arrival_depth) return x.arrival_depth < y.arrival_depth;
        if (x.mask != y.mask) return x.mask < y.mask;
        return x.terms < y.terms;
    });
    if (max_representations && frontier.entries.size() > max_representations)
        frontier.entries.resize(max_representations);
}

int candidate_depth(int left_depth, int right_depth) {
    require(left_depth >= 0 && right_depth >= 0, "negative candidate depth");
    return 1 + std::max(left_depth, right_depth);
}

bool candidate_allowed(int left_depth, int right_depth, int max_depth) {
    require(max_depth >= 0, "negative maximum depth");
    return candidate_depth(left_depth, right_depth) <= max_depth;
}

int lookup_node(const std::unordered_map<BitMask, int, BitMaskHash> &nodes,
                const BitMask &mask) {
    auto it = nodes.find(mask);
    return it == nodes.end() ? -1 : it->second;
}

void validate_target_width(const BitRows &targets, int inputs) {
    require(inputs > 0, "invalid input count");
    for (const auto &target : targets)
        require(target.bits() <= static_cast<std::size_t>(inputs), "target mask exceeds input width");
}

std::string canonical_representation(const Representation &representation) {
    auto terms = normalize_terms(representation.terms);
    std::ostringstream out;
    out << representation.arrival_depth << ':';
    for (int term : terms) out << term << ',';
    out << ':' << representation.mask.decimal();
    return out.str();
}

SearchState SearchState::initial(const BitRows &targets, int inputs,
                                 const BoundedScheduleOptions &options) {
    check_search_deadline(options);
    require(inputs > 0, "invalid input count");
    require(options.max_depth >= 0 && options.beam_width > 0 && options.branch_width > 0 &&
                options.max_representations > 0 && options.restarts > 0 &&
                options.timeout_seconds > 0 && options.lookahead_width >= 0 &&
                options.lookahead_depth >= 0 && options.residual_pairs >= 0 &&
                options.residual_targets >= 0 &&
                (options.residual_pairs == 0 || options.residual_targets > 0),
            "invalid bounded search options");
    validate_target_width(targets, inputs);
    SearchState state;
    state.inputs_ = inputs;
    state.options_ = options;
    state.nodes_.reserve(inputs);
    for (int i = 0; i < inputs; ++i) {
        check_search_deadline(options);
        state.nodes_.push_back(NodeRecord{BitMask::unit(i), 0, -1, -1});
        state.node_by_mask_.emplace(state.nodes_.back().mask, i);
    }
    state.frontiers_.resize(targets.size());
    for (std::size_t target = 0; target < targets.size(); ++target) {
        check_search_deadline(options);
        std::vector<int> terms;
        for (int bit = 0; bit < inputs; ++bit) {
            check_search_deadline(options);
            if (targets[target].bit(bit)) terms.push_back(bit);
        }
        auto representation = state.make_representation(terms);
        if (representation.arrival_depth > options.max_depth) {
            state.valid_ = false;
            return state;
        }
        pareto_insert(state.frontiers_[target], std::move(representation),
                      static_cast<std::size_t>(options.max_representations));
    }
    state.rebuild_pair_index();
    state.rebuild_residual_pair_index();
    return state;
}

int SearchState::node_depth(int node) const {
    require(node >= 0 && static_cast<std::size_t>(node) < nodes_.size(), "invalid search node");
    return nodes_[node].depth;
}

Representation SearchState::make_representation(const std::vector<int> &terms) const {
    check_search_deadline(options_);
    Representation result;
    result.terms = normalize_terms(terms);
    for (int node : result.terms) {
        check_search_deadline(options_);
        require(node >= 0 && static_cast<std::size_t>(node) < nodes_.size(),
                "invalid representation node");
        result.mask ^= nodes_[node].mask;
    }
    std::vector<int> depths;
    depths.reserve(result.terms.size());
    for (int node : result.terms) {
        check_search_deadline(options_);
        depths.push_back(node_depth(node));
    }
    result.arrival_depth = combine_depth(depths);
    return result;
}

void SearchState::rebuild_pair_index() {
    check_search_deadline(options_);
    pair_index_.clear();
    for (std::size_t target = 0; target < frontiers_.size(); ++target) {
        check_search_deadline(options_);
        for (std::size_t representation = 0; representation < frontiers_[target].entries.size();
             ++representation) {
            check_search_deadline(options_);
            const auto &entry = frontiers_[target].entries[representation];
            for (std::size_t i = 0; i < entry.terms.size(); ++i) {
                for (std::size_t j = i + 1; j < entry.terms.size(); ++j) {
                    check_search_deadline(options_);
                    auto edge = nodes_[entry.terms[i]].mask ^ nodes_[entry.terms[j]].mask;
                    pair_index_[edge].push_back(PairRef{static_cast<int>(target),
                                                        static_cast<int>(representation),
                                                        entry.terms[i], entry.terms[j]});
                }
            }
        }
    }
}

void SearchState::adjust_residual_pair_index(const ParetoFrontier &frontier,
                                              int node_count, int delta) {
    if (options_.residual_pairs == 0 || options_.residual_targets == 0) return;
    for (const auto &entry : frontier.entries) {
        std::set<std::uint64_t> seen;
        if (options_.allow_empty_toggle) {
            for (int left = 0; left < node_count; ++left)
                for (int right = left + 1; right < node_count; ++right)
                    seen.insert(residual_pair_key(left, right));
        } else {
            for (int present : entry.terms) {
                for (int other = 0; other < node_count; ++other) {
                    if (present == other) continue;
                    seen.insert(residual_pair_key(std::min(present, other),
                                                  std::max(present, other)));
                }
            }
        }
        for (auto key : seen) {
            auto &value = residual_pair_support_[key];
            value += delta;
            if (value < 0)
                throw std::runtime_error("invalid residual pair support key=" +
                                         std::to_string(key) + " delta=" +
                                         std::to_string(delta));
            if (!value) residual_pair_support_.erase(key);
        }
    }
}

void SearchState::rebuild_residual_pair_index() {
    residual_pair_support_.clear();
    for (const auto &frontier : frontiers_)
        adjust_residual_pair_index(frontier, static_cast<int>(nodes_.size()), 1);
}

std::vector<Candidate> SearchState::candidates() const {
    check_search_deadline(options_);
    std::map<std::pair<int, int>, Candidate> by_pair;
    for (const auto &[edge, refs] : pair_index_) {
        check_search_deadline(options_);
        for (const auto &ref : refs) {
            check_search_deadline(options_);
            if (!candidate_allowed(node_depth(ref.left), node_depth(ref.right),
                                   options_.max_depth))
                continue;
            auto key = std::make_pair(ref.left, ref.right);
            auto it = by_pair.find(key);
            if (it == by_pair.end()) {
                by_pair.emplace(key, Candidate{ref.left, ref.right,
                                               candidate_depth(node_depth(ref.left),
                                                              node_depth(ref.right)),
                                               static_cast<int>(refs.size())});
            }
        }
    }
    std::vector<Candidate> result;
    for (const auto &[key, candidate] : by_pair) result.push_back(candidate);
    check_search_deadline(options_);
    std::sort(result.begin(), result.end(), [](const Candidate &a, const Candidate &b) {
        if (a.support != b.support) return a.support > b.support;
        if (a.arrival_depth != b.arrival_depth) return a.arrival_depth < b.arrival_depth;
        if (a.left != b.left) return a.left < b.left;
        return a.right < b.right;
    });
    return result;
}

std::vector<Candidate> SearchState::toggle_candidates() const {
    check_search_deadline(options_);
    if (options_.residual_pairs == 0 || options_.residual_targets == 0) return {};
    struct PairSummary {
        int left = -1;
        int right = -1;
        int support = 0;
        int depth = 0;
    };
    std::vector<PairSummary> pairs;
    pairs.reserve(residual_pair_support_.size());
    for (const auto &[key, support] : residual_pair_support_) {
        int left = static_cast<int>(key >> 32);
        int right = static_cast<int>(std::uint32_t(key));
        if (support > 0 && candidate_allowed(node_depth(left), node_depth(right),
                                             options_.max_depth))
            pairs.push_back(PairSummary{left, right, support,
                                         candidate_depth(node_depth(left), node_depth(right))});
    }
    std::sort(pairs.begin(), pairs.end(), [](const PairSummary &a, const PairSummary &b) {
        if (a.support != b.support) return a.support > b.support;
        if (a.depth != b.depth) return a.depth < b.depth;
        if (a.left != b.left) return a.left < b.left;
        return a.right < b.right;
    });
    if (pairs.size() > static_cast<std::size_t>(options_.residual_pairs))
        pairs.resize(static_cast<std::size_t>(options_.residual_pairs));
    std::vector<Candidate> result;
    for (const auto &pair : pairs) {
        std::vector<Candidate> pair_candidates;
        for (std::size_t target = 0; target < frontiers_.size(); ++target) {
            for (std::size_t representation = 0;
                 representation < frontiers_[target].entries.size(); ++representation) {
                check_search_deadline(options_);
                const auto &entry = frontiers_[target].entries[representation];
                const bool has_left = std::find(entry.terms.begin(), entry.terms.end(), pair.left) != entry.terms.end();
                const bool has_right = std::find(entry.terms.begin(), entry.terms.end(), pair.right) != entry.terms.end();
                if (!has_left && !has_right && !options_.allow_empty_toggle) continue;
                int node = lookup_node(node_by_mask_, nodes_[pair.left].mask ^ nodes_[pair.right].mask);
                if (node < 0) node = static_cast<int>(nodes_.size());
                std::vector<int> terms = entry.terms;
                auto toggle_term = [&](int term) {
                    auto it = std::find(terms.begin(), terms.end(), term);
                    if (it == terms.end()) terms.push_back(term);
                    else terms.erase(it);
                };
                toggle_term(pair.left);
                toggle_term(pair.right);
                toggle_term(node);
                const int new_node_depth = node < static_cast<int>(nodes_.size())
                                               ? node_depth(node)
                                               : pair.depth;
                std::vector<int> arrival_depths;
                arrival_depths.reserve(terms.size());
                for (int term : terms)
                    arrival_depths.push_back(term == node ? new_node_depth : node_depth(term));
                const int new_arrival = combine_depth(arrival_depths);
                if (new_arrival > options_.max_depth) continue;
                const int term_reduction = static_cast<int>(entry.terms.size()) -
                                           static_cast<int>(terms.size());
                int future_support = 0;
                const auto new_mask = nodes_[pair.left].mask ^ nodes_[pair.right].mask;
                for (int term : terms) {
                    if (term == node && node >= static_cast<int>(nodes_.size())) continue;
                    const auto edge = new_mask ^ nodes_[term].mask;
                    auto edge_it = pair_index_.find(edge);
                    if (edge_it != pair_index_.end())
                        future_support = std::max(future_support,
                                                  static_cast<int>(edge_it->second.size()));
                }
                if (term_reduction <= 0 && new_arrival >= entry.arrival_depth &&
                    pair.support <= 1 && future_support == 0)
                    continue;
                Candidate candidate{pair.left, pair.right, pair.depth, pair.support,
                                    static_cast<int>(target), static_cast<int>(representation), true};
                candidate.affected_targets = pair.support;
                candidate.term_reduction = term_reduction;
                candidate.depth_slack = options_.max_depth - new_arrival;
                candidate.future_support = future_support;
                candidate.toggle_penalty = std::max(0, -term_reduction);
                pair_candidates.push_back(candidate);
            }
        }
        std::sort(pair_candidates.begin(), pair_candidates.end(), [](const Candidate &a, const Candidate &b) {
            if (a.term_reduction != b.term_reduction) return a.term_reduction > b.term_reduction;
            if (a.future_support != b.future_support) return a.future_support > b.future_support;
            if (a.affected_targets != b.affected_targets)
                return a.affected_targets > b.affected_targets;
            if (a.depth_slack != b.depth_slack) return a.depth_slack > b.depth_slack;
            if (a.toggle_penalty != b.toggle_penalty)
                return a.toggle_penalty < b.toggle_penalty;
            if (a.target != b.target) return a.target < b.target;
            return a.representation < b.representation;
        });
        if (pair_candidates.size() > static_cast<std::size_t>(options_.residual_targets))
            pair_candidates.resize(static_cast<std::size_t>(options_.residual_targets));
        result.insert(result.end(), pair_candidates.begin(), pair_candidates.end());
    }
    std::sort(result.begin(), result.end(), [](const Candidate &a, const Candidate &b) {
        if (a.term_reduction != b.term_reduction) return a.term_reduction > b.term_reduction;
        if (a.future_support != b.future_support) return a.future_support > b.future_support;
        if (a.affected_targets != b.affected_targets)
            return a.affected_targets > b.affected_targets;
        if (a.depth_slack != b.depth_slack) return a.depth_slack > b.depth_slack;
        if (a.toggle_penalty != b.toggle_penalty)
            return a.toggle_penalty < b.toggle_penalty;
        if (a.support != b.support) return a.support > b.support;
        if (a.target != b.target) return a.target < b.target;
        if (a.representation != b.representation) return a.representation < b.representation;
        if (a.left != b.left) return a.left < b.left;
        return a.right < b.right;
    });
    return result;
}

std::vector<Candidate> SearchState::bridge_candidates() const {
    check_search_deadline(options_);
    if (options_.residual_pairs == 0) return {};
    std::vector<Candidate> result;
    for (std::size_t ta = 0; ta < frontiers_.size(); ++ta) {
        for (std::size_t ra = 0; ra < frontiers_[ta].entries.size(); ++ra) {
            const auto &a = frontiers_[ta].entries[ra];
            for (std::size_t tb = ta + 1; tb < frontiers_.size(); ++tb) {
                std::size_t rb_start = 0;
                for (std::size_t rb = rb_start; rb < frontiers_[tb].entries.size(); ++rb) {
                    const auto &b = frontiers_[tb].entries[rb];
                    for (int left : a.terms) {
                        for (int right : b.terms) {
                            if (left == right) continue;
                            int lo = std::min(left, right), hi = std::max(left, right);
                            if (!candidate_allowed(node_depth(lo), node_depth(hi), options_.max_depth))
                                continue;
                            Candidate candidate{lo, hi,
                                                candidate_depth(node_depth(lo), node_depth(hi)), 1,
                                                static_cast<int>(ta), static_cast<int>(ra), false,
                                                true, static_cast<int>(tb), static_cast<int>(rb)};
                            candidate.affected_targets = ta == tb ? 1 : 2;
                            candidate.term_reduction = 2;
                            candidate.depth_slack = options_.max_depth - candidate.arrival_depth;
                            result.push_back(candidate);
                        }
                    }
                }
            }
        }
    }
    std::sort(result.begin(), result.end(), [](const Candidate &a, const Candidate &b) {
        if (a.affected_targets != b.affected_targets)
            return a.affected_targets > b.affected_targets;
        if (a.depth_slack != b.depth_slack) return a.depth_slack > b.depth_slack;
        if (a.left != b.left) return a.left < b.left;
        if (a.right != b.right) return a.right < b.right;
        if (a.target != b.target) return a.target < b.target;
        return a.target2 < b.target2;
    });
    if (result.size() > static_cast<std::size_t>(options_.residual_pairs))
        result.resize(static_cast<std::size_t>(options_.residual_pairs));
    return result;
}

SearchState SearchState::apply(const Candidate &candidate) const {
    check_search_deadline(options_);
    SearchState child = *this;
    const int parent_node_count = static_cast<int>(nodes_.size());
    child.rebuild_residual_pair_index();
    if (!valid_ || candidate.left < 0 || candidate.right < 0 ||
        candidate.left >= candidate.right ||
        static_cast<std::size_t>(candidate.left) >= nodes_.size() ||
        static_cast<std::size_t>(candidate.right) >= nodes_.size() ||
        !candidate_allowed(node_depth(candidate.left), node_depth(candidate.right),
                           options_.max_depth) || pair_support(candidate) == 0) {
        child.valid_ = false;
        return child;
    }
    auto result_mask = nodes_[candidate.left].mask ^ nodes_[candidate.right].mask;
    int new_depth = candidate_depth(node_depth(candidate.left), node_depth(candidate.right));
    int node = lookup_node(child.node_by_mask_, result_mask);
    if (node < 0) {
        node = static_cast<int>(child.nodes_.size());
        child.nodes_.push_back(NodeRecord{result_mask, new_depth,
                                          candidate.left, candidate.right});
        child.node_by_mask_.emplace(result_mask, node);
    }

    for (std::size_t target = 0; target < child.frontiers_.size(); ++target) {
        check_search_deadline(options_);
        child.adjust_residual_pair_index(frontiers_[target], parent_node_count, -1);
        ParetoFrontier next;
        for (const auto &entry : frontiers_[target].entries) {
            check_search_deadline(options_);
            pareto_insert(next, entry, static_cast<std::size_t>(options_.max_representations));
            auto left = std::find(entry.terms.begin(), entry.terms.end(), candidate.left);
            auto right = std::find(entry.terms.begin(), entry.terms.end(), candidate.right);
            if (left == entry.terms.end() || right == entry.terms.end()) continue;
            std::vector<int> terms = entry.terms;
            terms.erase(std::remove(terms.begin(), terms.end(), candidate.left), terms.end());
            terms.erase(std::remove(terms.begin(), terms.end(), candidate.right), terms.end());
            terms.push_back(node);
            auto replacement = child.make_representation(terms);
            if (replacement.arrival_depth <= options_.max_depth)
                pareto_insert(next, std::move(replacement),
                              static_cast<std::size_t>(options_.max_representations));
        }
        if (next.entries.empty()) {
            child.valid_ = false;
            return child;
        }
        child.frontiers_[target] = std::move(next);
        child.adjust_residual_pair_index(child.frontiers_[target],
                                         static_cast<int>(child.nodes_.size()), 1);
    }
    child.rebuild_pair_index();
    return child;
}

SearchState SearchState::apply_toggle(int left, int right, std::size_t target,
                                      std::size_t representation,
                                      int future_support) const {
    check_search_deadline(options_);
    SearchState child = *this;
    const int parent_node_count = static_cast<int>(nodes_.size());
    child.rebuild_residual_pair_index();
    if (!valid_ || left < 0 || right < 0 || left >= right ||
        static_cast<std::size_t>(left) >= nodes_.size() ||
        static_cast<std::size_t>(right) >= nodes_.size() || target >= frontiers_.size() ||
        representation >= frontiers_[target].entries.size()) {
        child.valid_ = false;
        return child;
    }
    const auto &entry = frontiers_[target].entries[representation];
    const bool has_left = std::find(entry.terms.begin(), entry.terms.end(), left) != entry.terms.end();
    const bool has_right = std::find(entry.terms.begin(), entry.terms.end(), right) != entry.terms.end();
    if (!has_left && !has_right && !options_.allow_empty_toggle) {
        child.valid_ = false;
        return child;
    }
    const auto result_mask = nodes_[left].mask ^ nodes_[right].mask;
    const int proposed_depth = candidate_depth(node_depth(left), node_depth(right));
    if (!candidate_allowed(node_depth(left), node_depth(right), options_.max_depth)) {
        child.valid_ = false;
        return child;
    }
    int node = lookup_node(child.node_by_mask_, result_mask);
    if (node < 0) {
        node = static_cast<int>(child.nodes_.size());
        child.nodes_.push_back(NodeRecord{result_mask, proposed_depth, left, right});
        child.node_by_mask_.emplace(result_mask, node);
    }
    std::vector<int> terms = entry.terms;
    auto toggle_term = [&](int term) {
        auto it = std::find(terms.begin(), terms.end(), term);
        if (it == terms.end()) terms.push_back(term);
        else terms.erase(it);
    };
    toggle_term(left);
    toggle_term(right);
    toggle_term(node);
    auto replacement = child.make_representation(terms);
    replacement.future_support = future_support;
    if (replacement.arrival_depth > options_.max_depth) {
        child.valid_ = false;
        return child;
    }
    ParetoFrontier next = frontiers_[target];
    child.adjust_residual_pair_index(frontiers_[target], parent_node_count, -1);
    pareto_insert(next, std::move(replacement),
                  static_cast<std::size_t>(options_.max_representations));
    if (next.entries.empty()) {
        child.valid_ = false;
        return child;
    }
    child.frontiers_[target] = std::move(next);
    child.adjust_residual_pair_index(child.frontiers_[target],
                                     static_cast<int>(child.nodes_.size()), 1);
    child.rebuild_pair_index();
    return child;
}

SearchState SearchState::apply_bridge(int left, int right, std::size_t left_target,
                                      std::size_t left_representation,
                                      std::size_t right_target,
                                      std::size_t right_representation) const {
    check_search_deadline(options_);
    SearchState child = *this;
    if (!valid_ || left_target >= frontiers_.size() || right_target >= frontiers_.size() ||
        left_representation >= frontiers_[left_target].entries.size() ||
        right_representation >= frontiers_[right_target].entries.size() || left == right)
        return child.valid_ = false, child;
    if (!candidate_allowed(node_depth(std::min(left, right)), node_depth(std::max(left, right)),
                           options_.max_depth))
        return child.valid_ = false, child;
    const auto mask = nodes_[left].mask ^ nodes_[right].mask;
    int node = lookup_node(child.node_by_mask_, mask);
    if (node < 0) {
        node = static_cast<int>(child.nodes_.size());
        child.nodes_.push_back(NodeRecord{mask, candidate_depth(node_depth(left), node_depth(right)), left, right});
        child.node_by_mask_.emplace(mask, node);
    }
    auto update = [&](std::size_t target, std::size_t representation) {
        const auto &entry = frontiers_[target].entries[representation];
        std::vector<int> terms = entry.terms;
        auto toggle = [&](int term) {
            auto it = std::find(terms.begin(), terms.end(), term);
            if (it == terms.end()) terms.push_back(term); else terms.erase(it);
        };
        toggle(left); toggle(right); toggle(node);
        auto replacement = child.make_representation(terms);
        if (replacement.arrival_depth > options_.max_depth) return false;
        ParetoFrontier next = frontiers_[target];
        pareto_insert(next, std::move(replacement), static_cast<std::size_t>(options_.max_representations));
        child.frontiers_[target] = std::move(next);
        return true;
    };
    if (!update(left_target, left_representation)) { child.valid_ = false; return child; }
    if (!update(right_target, right_representation)) { child.valid_ = false; return child; }
    child.rebuild_pair_index();
    child.rebuild_residual_pair_index();
    return child;
}

bool SearchState::complete() const {
    check_search_deadline(options_);
    if (!valid_) return false;
    for (const auto &frontier : frontiers_) {
        check_search_deadline(options_);
        bool done = false;
        for (const auto &entry : frontier.entries) {
            check_search_deadline(options_);
            done = done || entry.terms.size() <= 1;
        }
        if (!done) return false;
    }
    return true;
}

int SearchState::best_cost_lower_bound() const {
    check_search_deadline(options_);
    int bound = 0;
    for (const auto &frontier : frontiers_) {
        int target_bound = INT_MAX;
        for (const auto &entry : frontier.entries) {
            std::set<int> live;
            std::vector<int> pending = entry.terms;
            while (!pending.empty()) {
                check_search_deadline(options_);
                int node = pending.back();
                pending.pop_back();
                if (node < inputs_ || !live.insert(node).second) continue;
                require(static_cast<std::size_t>(node) < nodes_.size(), "invalid closure node");
                require(nodes_[node].left >= 0 && nodes_[node].right >= 0,
                        "closure gate is missing a producer");
                pending.push_back(nodes_[node].left);
                pending.push_back(nodes_[node].right);
            }
            target_bound = std::min(target_bound, static_cast<int>(live.size()));
        }
        if (target_bound != INT_MAX) bound = std::max(bound, target_bound);
    }
    // Pair replacement preserves disjoint primary-input supports within each
    // representation, starting from unit inputs. Producers are cancellation
    // free, so each closure above uses at most weight(target)-1 gates. Any XOR
    // DAG computing that target needs at least weight(target)-1 gates, even
    // with cancellation. This remains safe if later reuse drops old producers.
    return bound;
}

int SearchState::cost_priority() const {
    check_search_deadline(options_);
    // This union is useful for ordering states because it estimates the
    // amount of materialized structure. It includes mutually exclusive
    // Pareto alternatives, so it is deliberately not a lower bound.
    // Nodes can be created by a transition whose replacement is later removed
    // by Pareto pruning.  Such nodes are not part of the exported circuit and
    // must not inflate the lower bound used for branch-and-bound pruning.
    std::set<int> live;
    std::vector<int> pending;
    for (const auto &frontier : frontiers_) {
        check_search_deadline(options_);
        for (const auto &entry : frontier.entries) {
            check_search_deadline(options_);
            pending.insert(pending.end(), entry.terms.begin(), entry.terms.end());
        }
    }
    while (!pending.empty()) {
        check_search_deadline(options_);
        const int node = pending.back();
        pending.pop_back();
        if (node < inputs_ || !live.insert(node).second) continue;
        require(static_cast<std::size_t>(node) < nodes_.size(), "invalid live search node");
        require(nodes_[node].left >= 0 && nodes_[node].right >= 0,
                "live gate is missing a producer");
        pending.push_back(nodes_[node].left);
        pending.push_back(nodes_[node].right);
    }
    int bound = static_cast<int>(live.size());
    for (const auto &frontier : frontiers_) {
        check_search_deadline(options_);
        int remaining = INT_MAX;
        for (const auto &entry : frontier.entries) {
            check_search_deadline(options_);
            remaining = std::min(remaining, static_cast<int>(entry.terms.size()));
        }
        if (remaining != INT_MAX)
            bound = std::max(bound, static_cast<int>(live.size()) +
                                      std::max(0, remaining - 1));
    }
    return bound;
}

std::string SearchState::fingerprint() const {
    check_search_deadline(options_);
    std::ostringstream out;
    out << "targets=" << frontiers_.size() << ';';
    for (const auto &frontier : frontiers_) {
        int min_terms = INT_MAX;
        int min_depth = INT_MAX;
        std::vector<int> depths;
        for (const auto &entry : frontier.entries) {
            min_terms = std::min(min_terms, static_cast<int>(entry.terms.size()));
            min_depth = std::min(min_depth, entry.arrival_depth);
            depths.push_back(entry.arrival_depth);
        }
        std::sort(depths.begin(), depths.end());
        out << min_terms << ':' << min_depth << ':' << frontier.entries.size() << ':';
        for (int depth : depths) out << depth << ',';
        out << ';';
    }
    std::vector<int> supports;
    supports.reserve(pair_index_.size());
    for (const auto &[edge, refs] : pair_index_) supports.push_back(static_cast<int>(refs.size()));
    std::sort(supports.begin(), supports.end(), std::greater<int>());
    if (supports.size() > 8) supports.resize(8);
    out << "support=";
    for (int support : supports) out << support << ',';
    int toggle_representations = 0;
    for (const auto &frontier : frontiers_)
        if (frontier.entries.size() > 1) ++toggle_representations;
    out << "toggle=" << toggle_representations
        << ";closure=" << cost_priority()
        << ";nodes=" << nodes_.size() << ";pairs=" << pair_index_.size();
    return out.str();
}

const ParetoFrontier &SearchState::frontier(std::size_t index) const {
    require(index < frontiers_.size(), "frontier index out of range");
    return frontiers_[index];
}

int SearchState::pair_support(const Candidate &candidate) const {
    check_search_deadline(options_);
    if (candidate.left < 0 || candidate.right < 0 ||
        static_cast<std::size_t>(candidate.left) >= nodes_.size() ||
        static_cast<std::size_t>(candidate.right) >= nodes_.size())
        return 0;
    auto edge = nodes_[candidate.left].mask ^ nodes_[candidate.right].mask;
    auto it = pair_index_.find(edge);
    if (it == pair_index_.end()) return 0;
    int support = 0;
    for (const auto &ref : it->second)
        if (ref.left == candidate.left && ref.right == candidate.right) ++support;
    return support;
}

std::string SearchState::signature() const {
    check_search_deadline(options_);
    std::vector<int> order(nodes_.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (nodes_[a].mask != nodes_[b].mask) return nodes_[a].mask < nodes_[b].mask;
        if (nodes_[a].depth != nodes_[b].depth) return nodes_[a].depth < nodes_[b].depth;
        return a < b;
    });
    std::vector<int> canonical(nodes_.size());
    std::ostringstream out;
    for (std::size_t i = 0; i < order.size(); ++i) {
        check_search_deadline(options_);
        canonical[order[i]] = static_cast<int>(i);
        out << nodes_[order[i]].mask.decimal() << '@' << nodes_[order[i]].depth << ';';
    }
    out << '|';
    for (const auto &frontier : frontiers_) {
        check_search_deadline(options_);
        out << '[';
        for (const auto &entry : frontier.entries) {
            check_search_deadline(options_);
            out << entry.arrival_depth << ':';
            for (int node : entry.terms) out << canonical[node] << ',';
            out << entry.mask.decimal() << ';';
        }
        out << ']';
    }
    return out.str();
}

Circuit SearchState::export_circuit() const {
    check_search_deadline(options_);
    require(valid_ && complete(), "cannot export incomplete search state");
    GreedyNetwork circuit(inputs_);
    std::vector<int> mapped(nodes_.size(), -1);
    for (int i = 0; i < inputs_; ++i) mapped[i] = i;
    for (std::size_t i = static_cast<std::size_t>(inputs_); i < nodes_.size(); ++i) {
        check_search_deadline(options_);
        require(nodes_[i].left >= 0 && nodes_[i].right >= 0, "missing search producer");
        mapped[i] = circuit.xor_nodes(mapped[nodes_[i].left], mapped[nodes_[i].right]);
    }
    std::vector<int> outputs;
    outputs.reserve(frontiers_.size());
    for (const auto &frontier : frontiers_) {
        check_search_deadline(options_);
        const Representation *done = nullptr;
        for (const auto &entry : frontier.entries) {
            check_search_deadline(options_);
            if (entry.terms.size() <= 1) {
                done = &entry;
                break;
            }
        }
        require(done != nullptr, "incomplete search frontier");
        outputs.push_back(done->terms.empty() ? inputs_ : mapped[done->terms.front()]);
    }
    auto result = circuit.export_circuit(outputs);
    return result;
}

BoundedScheduleResult synthesize_bounded(const BitRows &targets, int inputs,
                                         const BoundedScheduleOptions &options) {
    validate_target_width(targets, inputs);
    require(options.max_depth >= 0 && options.beam_width > 0 && options.branch_width > 0 &&
                options.max_representations > 0 && options.restarts > 0 &&
                options.timeout_seconds > 0 && options.lookahead_width >= 0 &&
                options.lookahead_depth >= 0 && options.residual_pairs >= 0 &&
                options.residual_targets >= 0 &&
                (options.residual_pairs == 0 || options.residual_targets > 0),
            "invalid bounded search options");
    BoundedScheduleResult result;
    result.requested_depth = options.max_depth;
    result.depth_lower_bound = depth_lower_bound(targets);
    if (options.max_depth < result.depth_lower_bound) {
        result.status = "infeasible_by_lower_bound";
        return result;
    }

    const auto started = std::chrono::steady_clock::now();
    const auto deadline = started + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                                      std::chrono::duration<double>(options.timeout_seconds));
    auto search_options = options;
    search_options.deadline = deadline;
    SearchState initial;
    try {
        initial = SearchState::initial(targets, inputs, search_options);
    } catch (const SearchDeadlineExceeded &) {
        if (result.status != "verified") result.status = "timeout";
        return result;
    }
    if (!initial.valid()) {
        if (result.status != "verified") result.status = "search_failed";
        return result;
    }
    if (std::chrono::steady_clock::now() >= deadline) {
        if (result.status != "verified") result.status = "timeout";
        return result;
    }
    std::vector<SearchState> beam{initial};
    std::unordered_map<std::string, int> transposition;
    try {
        transposition.emplace(initial.signature(), initial.best_cost_lower_bound());
    } catch (const SearchDeadlineExceeded &) {
        result.status = "timeout";
        return result;
    }
    auto consider_candidate = [&](Circuit circuit) {
        if (std::chrono::steady_clock::now() >= deadline) throw SearchDeadlineExceeded();
        circuit.verify(targets);
        auto profile = depth_profile(circuit);
        if (profile.depth > options.max_depth) return;
        if (result.status != "verified" || circuit.cost() < result.circuit.cost()) {
            result.circuit = std::move(circuit);
            result.status = "verified";
        }
    };
    const auto deadline_reached = [&]() { return std::chrono::steady_clock::now() >= deadline; };
    ScheduleOptions baseline_options;
    baseline_options.share_threshold = 2;
    baseline_options.depth_aware = true;
    try {
        if (!deadline_reached())
            consider_candidate(synthesize_schedule(targets, inputs, baseline_options, 0, deadline));
        if (!deadline_reached()) {
            auto baseline = synthesize_schedule(targets, inputs, baseline_options, 0, deadline);
            consider_candidate(rebalanced(baseline, deadline));
        }
    } catch (const SearchDeadlineExceeded &) {
        if (result.status != "verified") result.status = "timeout";
        return result;
    }
    if (deadline_reached()) {
        if (result.status != "verified") result.status = "timeout";
        return result;
    }
    bool timed_out = false;
    try {
    for (int restart = 0; restart < options.restarts && !timed_out; ++restart) {
        if (restart > 0) beam = {initial};
        while (!beam.empty()) {
            if (deadline_reached()) {
                timed_out = true;
                break;
            }
            std::vector<SearchState> next;
            for (const auto &state : beam) {
                ++result.explored_states;
                if (state.complete()) {
                    if (deadline_reached()) throw SearchDeadlineExceeded();
                    auto circuit = state.export_circuit();
                    auto profile = depth_profile(circuit);
                    circuit.verify(targets);
                    if (deadline_reached()) throw SearchDeadlineExceeded();
                    require(profile.depth <= options.max_depth, "bounded result exceeds depth limit");
                    if (result.status != "verified" || circuit.cost() < result.circuit.cost()) {
                        result.circuit = std::move(circuit);
                        result.status = "verified";
                    }
                    continue;
                }
                auto candidates = state.candidates();
                if (options.residual_pairs > 0) {
                    auto toggles = state.toggle_candidates();
                    result.toggle_candidates += toggles.size();
                    candidates.insert(candidates.end(), toggles.begin(), toggles.end());
                    auto bridges = state.bridge_candidates();
                    result.bridge_candidates += bridges.size();
                    candidates.insert(candidates.end(), bridges.begin(), bridges.end());
                }
                SeedRandom random(options.seed + static_cast<std::uint64_t>(restart));
                if (restart > 0 || options.seed != 0)
                    for (std::size_t i = candidates.size(); i > 1; --i)
                        std::swap(candidates[i - 1], candidates[random.below(i)]);
                if (options.lookahead_width > 0 && options.lookahead_depth > 0 &&
                    candidates.size() > 1) {
                    struct Ranked {
                        Candidate candidate;
                        int bound;
                        std::size_t original;
                    };
                    const std::size_t width = std::min<std::size_t>(
                        candidates.size(), static_cast<std::size_t>(options.lookahead_width));
                    std::vector<Ranked> ranked;
                    ranked.reserve(width);
                    std::vector<bool> scored(candidates.size(), false);
                    for (std::size_t i = 0; i < width; ++i) {
                        auto probe = candidates[i].toggle
                                         ? state.apply_toggle(
                                               candidates[i].left, candidates[i].right,
                                               static_cast<std::size_t>(candidates[i].target),
                                               static_cast<std::size_t>(candidates[i].representation),
                                               candidates[i].future_support)
                                         : state.apply(candidates[i]);
                        if (!probe.valid()) continue;
                        int bound = probe.cost_priority();
                        for (int step = 1; step < options.lookahead_depth && !probe.complete();
                             ++step) {
                            auto next = probe.candidates();
                            if (options.residual_pairs > 0) {
                                auto toggles = probe.toggle_candidates();
                                next.insert(next.end(), toggles.begin(), toggles.end());
                                auto bridges = probe.bridge_candidates();
                                result.bridge_candidates += bridges.size();
                                next.insert(next.end(), bridges.begin(), bridges.end());
                            }
                            if (next.empty()) break;
                            const auto &future = next.front();
                            probe = future.toggle
                                        ? probe.apply_toggle(
                                              future.left, future.right,
                                              static_cast<std::size_t>(future.target),
                                              static_cast<std::size_t>(future.representation),
                                              future.future_support)
                                        : future.bridge
                                            ? probe.apply_bridge(
                                                  future.left, future.right,
                                                  static_cast<std::size_t>(future.target),
                                                  static_cast<std::size_t>(future.representation),
                                                  static_cast<std::size_t>(future.target2),
                                                  static_cast<std::size_t>(future.representation2))
                                            : probe.apply(future);
                            if (!probe.valid()) break;
                            bound = probe.cost_priority();
                        }
                        ranked.push_back(Ranked{candidates[i], bound, i});
                        scored[i] = true;
                    }
                    std::sort(ranked.begin(), ranked.end(), [](const Ranked &a, const Ranked &b) {
                        if (a.bound != b.bound) return a.bound < b.bound;
                        if (a.candidate.support != b.candidate.support)
                            return a.candidate.support > b.candidate.support;
                        if (a.candidate.arrival_depth != b.candidate.arrival_depth)
                            return a.candidate.arrival_depth < b.candidate.arrival_depth;
                        return a.original < b.original;
                    });
                    std::vector<Candidate> ordered;
                    ordered.reserve(candidates.size());
                    for (const auto &entry : ranked) ordered.push_back(entry.candidate);
                    for (std::size_t i = 0; i < candidates.size(); ++i)
                        if (!scored[i]) ordered.push_back(candidates[i]);
                    candidates = std::move(ordered);
                }
                if (candidates.size() > static_cast<std::size_t>(options.branch_width))
                    candidates.resize(static_cast<std::size_t>(options.branch_width));
                for (const auto &candidate : candidates) {
                    auto child = candidate.toggle
                                     ? state.apply_toggle(candidate.left, candidate.right,
                                                          static_cast<std::size_t>(candidate.target),
                                                          static_cast<std::size_t>(candidate.representation),
                                                          candidate.future_support)
                                     : candidate.bridge
                                         ? state.apply_bridge(candidate.left, candidate.right,
                                                              static_cast<std::size_t>(candidate.target),
                                                              static_cast<std::size_t>(candidate.representation),
                                                              static_cast<std::size_t>(candidate.target2),
                                                              static_cast<std::size_t>(candidate.representation2))
                                         : state.apply(candidate);
                    if (!child.valid()) {
                        ++result.pruned_depth;
                        continue;
                    }
                    if (result.status == "verified" &&
                        child.best_cost_lower_bound() >= result.circuit.cost()) {
                        ++result.pruned_bound;
                        continue;
                    }
                    auto key = child.signature();
                    auto bound = child.best_cost_lower_bound();
                    auto seen = transposition.find(key);
                    if (seen != transposition.end() && seen->second <= bound) {
                        ++result.pruned_bound;
                        continue;
                    }
                    transposition[key] = bound;
                    next.push_back(std::move(child));
                }
            }
            std::map<std::string, std::vector<SearchState>> buckets;
            for (auto &candidate : next)
                buckets[candidate.fingerprint()].push_back(std::move(candidate));
            result.diversity_buckets = std::max(result.diversity_buckets, buckets.size());
            for (auto &bucket : buckets) {
                std::sort(bucket.second.begin(), bucket.second.end(),
                          [](const SearchState &a, const SearchState &b) {
                              if (a.cost_priority() != b.cost_priority())
                                  return a.cost_priority() < b.cost_priority();
                              return a.signature() < b.signature();
                          });
                if (bucket.second.size() > 1)
                    bucket.second.resize(1);
            }
            std::vector<std::size_t> cursors(buckets.size(), 0);
            std::vector<SearchState> diverse;
            diverse.reserve(std::min<std::size_t>(buckets.size(), options.beam_width));
            while (diverse.size() < static_cast<std::size_t>(options.beam_width)) {
                bool added = false;
                std::size_t bucket_index = 0;
                for (auto &bucket : buckets) {
                    if (cursors[bucket_index] < bucket.second.size() &&
                        diverse.size() < static_cast<std::size_t>(options.beam_width)) {
                        diverse.push_back(std::move(bucket.second[cursors[bucket_index]++]));
                        added = true;
                    }
                    ++bucket_index;
                }
                if (!added) break;
            }
            beam = std::move(diverse);
        }
    }
    } catch (const SearchDeadlineExceeded &) {
        if (result.status != "verified") result.status = "timeout";
        return result;
    }
    if (result.status == "verified") return result;
    result.status = timed_out ? "timeout" : "search_failed";
    result.circuit = Circuit();
    return result;
}

} // namespace xd
