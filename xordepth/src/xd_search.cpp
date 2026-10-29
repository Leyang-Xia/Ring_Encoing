#include "xd/search.hpp"
#include <chrono>
#include <climits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace xd {
namespace {
using ringlab::require;
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

int lookup_node(const std::unordered_map<ringlab::Mask, int, ringlab::MaskHash> &nodes,
                const ringlab::Mask &mask) {
    auto it = nodes.find(mask);
    return it == nodes.end() ? -1 : it->second;
}

void validate_target_width(const ringlab::Rows &targets, int inputs) {
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

SearchState SearchState::initial(const ringlab::Rows &targets, int inputs,
                                 const BoundedScheduleOptions &options) {
    require(inputs > 0, "invalid input count");
    require(options.max_depth >= 0 && options.beam_width > 0 && options.branch_width > 0 &&
                options.max_representations > 0 && options.restarts > 0 &&
                options.timeout_seconds > 0,
            "invalid bounded search options");
    validate_target_width(targets, inputs);
    SearchState state;
    state.inputs_ = inputs;
    state.options_ = options;
    state.nodes_.reserve(inputs);
    for (int i = 0; i < inputs; ++i) {
        state.nodes_.push_back(NodeRecord{ringlab::Mask::unit(i), 0, -1, -1});
        state.node_by_mask_.emplace(state.nodes_.back().mask, i);
    }
    state.frontiers_.resize(targets.size());
    for (std::size_t target = 0; target < targets.size(); ++target) {
        std::vector<int> terms;
        for (int bit = 0; bit < inputs; ++bit)
            if (targets[target].bit(bit)) terms.push_back(bit);
        auto representation = state.make_representation(terms);
        if (representation.arrival_depth > options.max_depth) {
            state.valid_ = false;
            return state;
        }
        pareto_insert(state.frontiers_[target], std::move(representation),
                      static_cast<std::size_t>(options.max_representations));
    }
    state.rebuild_pair_index();
    return state;
}

int SearchState::node_depth(int node) const {
    require(node >= 0 && static_cast<std::size_t>(node) < nodes_.size(), "invalid search node");
    return nodes_[node].depth;
}

Representation SearchState::make_representation(const std::vector<int> &terms) const {
    Representation result;
    result.terms = normalize_terms(terms);
    for (int node : result.terms) {
        require(node >= 0 && static_cast<std::size_t>(node) < nodes_.size(),
                "invalid representation node");
        result.mask ^= nodes_[node].mask;
    }
    std::vector<int> depths;
    depths.reserve(result.terms.size());
    for (int node : result.terms) depths.push_back(node_depth(node));
    result.arrival_depth = combine_depth(depths);
    return result;
}

void SearchState::rebuild_pair_index() {
    pair_index_.clear();
    for (std::size_t target = 0; target < frontiers_.size(); ++target) {
        for (std::size_t representation = 0; representation < frontiers_[target].entries.size();
             ++representation) {
            const auto &entry = frontiers_[target].entries[representation];
            for (std::size_t i = 0; i < entry.terms.size(); ++i) {
                for (std::size_t j = i + 1; j < entry.terms.size(); ++j) {
                    auto edge = nodes_[entry.terms[i]].mask ^ nodes_[entry.terms[j]].mask;
                    pair_index_[edge].push_back(PairRef{static_cast<int>(target),
                                                        static_cast<int>(representation),
                                                        entry.terms[i], entry.terms[j]});
                }
            }
        }
    }
}

std::vector<Candidate> SearchState::candidates() const {
    std::map<std::pair<int, int>, Candidate> by_pair;
    for (const auto &[edge, refs] : pair_index_) {
        for (const auto &ref : refs) {
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
    std::sort(result.begin(), result.end(), [](const Candidate &a, const Candidate &b) {
        if (a.support != b.support) return a.support > b.support;
        if (a.arrival_depth != b.arrival_depth) return a.arrival_depth < b.arrival_depth;
        if (a.left != b.left) return a.left < b.left;
        return a.right < b.right;
    });
    return result;
}

SearchState SearchState::apply(const Candidate &candidate) const {
    SearchState child = *this;
    if (!valid_ || candidate.left < 0 || candidate.right < 0 ||
        candidate.left == candidate.right ||
        static_cast<std::size_t>(candidate.left) >= nodes_.size() ||
        static_cast<std::size_t>(candidate.right) >= nodes_.size() ||
        !candidate_allowed(node_depth(candidate.left), node_depth(candidate.right),
                           options_.max_depth)) {
        child.valid_ = false;
        return child;
    }
    auto result_mask = nodes_[candidate.left].mask ^ nodes_[candidate.right].mask;
    int node = lookup_node(child.node_by_mask_, result_mask);
    if (node < 0) {
        node = static_cast<int>(child.nodes_.size());
        child.nodes_.push_back(NodeRecord{result_mask, candidate.arrival_depth,
                                          candidate.left, candidate.right});
        child.node_by_mask_.emplace(result_mask, node);
    }

    for (std::size_t target = 0; target < child.frontiers_.size(); ++target) {
        ParetoFrontier next;
        for (const auto &entry : frontiers_[target].entries) {
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
    }
    child.rebuild_pair_index();
    return child;
}

bool SearchState::complete() const {
    if (!valid_) return false;
    for (const auto &frontier : frontiers_) {
        bool done = false;
        for (const auto &entry : frontier.entries)
            done = done || entry.terms.size() <= 1;
        if (!done) return false;
    }
    return true;
}

int SearchState::best_cost_lower_bound() const {
    int bound = static_cast<int>(nodes_.size()) - inputs_;
    for (const auto &frontier : frontiers_) {
        int remaining = INT_MAX;
        for (const auto &entry : frontier.entries)
            remaining = std::min(remaining, static_cast<int>(entry.terms.size()));
        if (remaining != INT_MAX) bound = std::max(bound, static_cast<int>(nodes_.size()) - inputs_ +
                                                    std::max(0, remaining - 1));
    }
    return bound;
}

const ParetoFrontier &SearchState::frontier(std::size_t index) const {
    require(index < frontiers_.size(), "frontier index out of range");
    return frontiers_[index];
}

int SearchState::pair_support(const Candidate &candidate) const {
    if (candidate.left < 0 || candidate.right < 0 ||
        static_cast<std::size_t>(candidate.left) >= nodes_.size() ||
        static_cast<std::size_t>(candidate.right) >= nodes_.size())
        return 0;
    auto edge = nodes_[candidate.left].mask ^ nodes_[candidate.right].mask;
    auto it = pair_index_.find(edge);
    return it == pair_index_.end() ? 0 : static_cast<int>(it->second.size());
}

std::string SearchState::signature() const {
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
        canonical[order[i]] = static_cast<int>(i);
        out << nodes_[order[i]].mask.decimal() << '@' << nodes_[order[i]].depth << ';';
    }
    out << '|';
    for (const auto &frontier : frontiers_) {
        out << '[';
        for (const auto &entry : frontier.entries) {
            out << entry.arrival_depth << ':';
            for (int node : entry.terms) out << canonical[node] << ',';
            out << entry.mask.decimal() << ';';
        }
        out << ']';
    }
    return out.str();
}

ringlab::Circuit SearchState::export_circuit() const {
    require(valid_ && complete(), "cannot export incomplete search state");
    ringlab::GreedyCircuit circuit(inputs_);
    std::vector<int> mapped(nodes_.size(), -1);
    for (int i = 0; i < inputs_; ++i) mapped[i] = i;
    for (std::size_t i = static_cast<std::size_t>(inputs_); i < nodes_.size(); ++i) {
        require(nodes_[i].left >= 0 && nodes_[i].right >= 0, "missing search producer");
        mapped[i] = circuit.xor_nodes(mapped[nodes_[i].left], mapped[nodes_[i].right]);
    }
    std::vector<int> outputs;
    outputs.reserve(frontiers_.size());
    for (const auto &frontier : frontiers_) {
        const Representation *done = nullptr;
        for (const auto &entry : frontier.entries)
            if (entry.terms.size() <= 1) {
                done = &entry;
                break;
            }
        require(done != nullptr, "incomplete search frontier");
        outputs.push_back(done->terms.empty() ? inputs_ : mapped[done->terms.front()]);
    }
    auto result = circuit.export_circuit(outputs);
    return result;
}

BoundedScheduleResult synthesize_bounded(const ringlab::Rows &targets, int inputs,
                                         const BoundedScheduleOptions &options) {
    validate_target_width(targets, inputs);
    require(options.max_depth >= 0 && options.beam_width > 0 && options.branch_width > 0 &&
                options.max_representations > 0 && options.restarts > 0 &&
                options.timeout_seconds > 0,
            "invalid bounded search options");
    BoundedScheduleResult result;
    result.requested_depth = options.max_depth;
    result.depth_lower_bound = depth_lower_bound(targets);
    if (options.max_depth < result.depth_lower_bound) {
        result.status = "infeasible_by_lower_bound";
        return result;
    }

    const auto started = std::chrono::steady_clock::now();
    SearchState initial = SearchState::initial(targets, inputs, options);
    if (!initial.valid()) {
        result.status = "search_failed";
        return result;
    }
    std::vector<SearchState> beam{initial};
    std::unordered_map<std::string, int> transposition;
    transposition.emplace(initial.signature(), initial.best_cost_lower_bound());
    auto consider_incumbent = [&](ringlab::Circuit circuit) {
        circuit.verify(targets);
        auto profile = depth_profile(circuit);
        if (profile.depth > options.max_depth) return;
        if (result.status != "verified" || circuit.cost() < result.circuit.cost()) {
            result.circuit = std::move(circuit);
            result.status = "verified";
        }
    };
    ScheduleOptions baseline_options;
    baseline_options.share_threshold = 2;
    baseline_options.depth_aware = true;
    consider_incumbent(synthesize_schedule(targets, inputs, baseline_options));
    consider_incumbent(rebalanced(synthesize_schedule(targets, inputs, baseline_options)));
    bool timed_out = false;
    for (int restart = 0; restart < options.restarts && !timed_out; ++restart) {
        if (restart > 0) beam = {initial};
        while (!beam.empty()) {
            if (std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count() >=
                options.timeout_seconds) {
                timed_out = true;
                break;
            }
            std::vector<SearchState> next;
            for (const auto &state : beam) {
                ++result.explored_states;
                if (state.complete()) {
                    auto circuit = state.export_circuit();
                    auto profile = depth_profile(circuit);
                    circuit.verify(targets);
                    require(profile.depth <= options.max_depth, "bounded result exceeds depth limit");
                    if (result.status != "verified" || circuit.cost() < result.circuit.cost()) {
                        result.circuit = std::move(circuit);
                        result.status = "verified";
                    }
                    continue;
                }
                auto candidates = state.candidates();
                ringlab::SeedRandom random(options.seed + static_cast<std::uint64_t>(restart));
                if (restart > 0 || options.seed != 0)
                    for (std::size_t i = candidates.size(); i > 1; --i)
                        std::swap(candidates[i - 1], candidates[random.below(i)]);
                if (candidates.size() > static_cast<std::size_t>(options.branch_width))
                    candidates.resize(static_cast<std::size_t>(options.branch_width));
                for (const auto &candidate : candidates) {
                    auto child = state.apply(candidate);
                    if (!child.valid()) {
                        ++result.pruned_depth;
                        continue;
                    }
                    if (result.status == "verified" &&
                        child.best_cost_lower_bound() > result.circuit.cost()) {
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
            std::sort(next.begin(), next.end(), [](const SearchState &a, const SearchState &b) {
                if (a.best_cost_lower_bound() != b.best_cost_lower_bound())
                    return a.best_cost_lower_bound() < b.best_cost_lower_bound();
                return a.signature() < b.signature();
            });
            if (next.size() > static_cast<std::size_t>(options.beam_width))
                next.resize(static_cast<std::size_t>(options.beam_width));
            beam = std::move(next);
        }
    }
    if (result.status == "verified") return result;
    result.status = timed_out ? "timeout" : "search_failed";
    result.circuit = ringlab::Circuit();
    return result;
}

} // namespace xd
