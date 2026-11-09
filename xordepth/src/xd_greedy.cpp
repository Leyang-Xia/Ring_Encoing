#include "xd/search.hpp"

#include <algorithm>
#include <chrono>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <tuple>
#include <utility>

namespace xd {

namespace {

// Encode an unordered node-id pair as one sortable 64-bit key. Node ids are
// non-negative 32-bit values; the smaller id goes to the high half so key
// equality is exactly pair equality.
std::uint64_t mwsq_pair_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (std::uint64_t(static_cast<std::uint32_t>(a)) << 32) |
           std::uint64_t(static_cast<std::uint32_t>(b));
}

// (max support, ties at max) over a flat multiset of pair keys. Equivalent to
// counting a map from pair to multiplicity: a run of equal keys in the sorted
// vector is exactly the map entry, so max/ties are identical by construction.
std::pair<std::size_t, std::size_t> mwsq_summary_from_keys(
    std::vector<std::uint64_t> &keys) {
    if (keys.empty()) return {0, 0};
    std::sort(keys.begin(), keys.end());
    std::size_t max_support = 0;
    std::size_t ties_at_max = 0;
    std::size_t i = 0;
    while (i < keys.size()) {
        std::size_t run = 1;
        while (i + run < keys.size() && keys[i + run] == keys[i]) ++run;
        if (run > max_support) {
            max_support = run;
            ties_at_max = 1;
        } else if (run == max_support) {
            ++ties_at_max;
        }
        i += run;
    }
    return {max_support, ties_at_max};
}

} // namespace

std::pair<std::size_t, std::size_t> mwsq_support_summary(
    const std::vector<std::vector<int>> &rows) {
    // Flat-key counting: the std::map node per distinct pair dominated the
    // greedy runtime; the sorted-run scan computes the identical summary.
    std::vector<std::uint64_t> keys;
    std::size_t total = 0;
    for (const auto &row : rows) total += row.size() * (row.size() - 1) / 2;
    keys.reserve(total);
    for (const auto &row : rows) {
        for (std::size_t i = 0; i + 1 < row.size(); ++i) {
            for (std::size_t j = i + 1; j < row.size(); ++j) {
                keys.push_back(mwsq_pair_key(row[i], row[j]));
            }
        }
    }
    return mwsq_summary_from_keys(keys);
}

} // namespace xd

namespace xd {
namespace {

using Deadline = std::chrono::steady_clock::time_point;
using Expression = std::set<int>;
using Pair = std::pair<int, int>;

void check_deadline(Deadline deadline) {
    if (deadline != Deadline::max() && std::chrono::steady_clock::now() >= deadline)
        throw SearchDeadlineExceeded();
}

int next_depth(int left, int right) {
    const int parent = std::max(left, right);
    require(parent < std::numeric_limits<int>::max(), "greedy depth overflow");
    return parent + 1;
}

struct PairInfo {
    int node = -1;
    int depth = 0;
    BitMask mask;
};

struct NodeStore {
    explicit NodeStore(int input_count) : inputs(input_count) {
        circuit.inputs = input_count;
        circuit.zero_node = input_count;
        masks.reserve(static_cast<std::size_t>(input_count) + 1);
        depths.reserve(static_cast<std::size_t>(input_count) + 1);
        for (int i = 0; i < input_count; ++i) {
            masks.push_back(BitMask::unit(i));
            depths.push_back(0);
            by_mask.emplace(masks.back(), i);
        }
        masks.emplace_back();
        depths.push_back(0);
        by_mask.emplace(BitMask(), input_count);
    }

    int zero() const { return inputs; }

    PairInfo pair_info(int left, int right) const {
        PairInfo result;
        result.mask = masks.at(left) ^ masks.at(right);
        result.depth = next_depth(depths.at(left), depths.at(right));
        auto found = by_mask.find(result.mask);
        if (found != by_mask.end() && depths.at(found->second) <= result.depth) {
            result.node = found->second;
            result.depth = depths.at(found->second);
        }
        return result;
    }

    int get_node(int left, int right, Deadline deadline) {
        check_deadline(deadline);
        auto info = pair_info(left, right);
        if (info.node >= 0) return info.node;
        require(masks.size() < static_cast<std::size_t>(std::numeric_limits<int>::max()),
                "too many greedy nodes");
        const int node = static_cast<int>(masks.size());
        masks.push_back(std::move(info.mask));
        depths.push_back(info.depth);
        circuit.ops.push_back({node, left, right});
        auto found = by_mask.find(masks.back());
        if (found == by_mask.end())
            by_mask.emplace(masks.back(), node);
        else if (info.depth < depths.at(found->second))
            found->second = node;
        check_deadline(deadline);
        return node;
    }

    int inputs;
    Circuit circuit;
    std::vector<BitMask> masks;
    std::vector<int> depths;
    std::unordered_map<BitMask, int, BitMaskHash> by_mask;
};

struct Candidate {
    int left = -1;
    int right = -1;
    int node = -1;
    int depth = 0;
    int gain = 0;
    int score = 0;
    long long delta = 0;
    BitMask mask;
    std::vector<int> uses;
};

Expression replaced(const Expression &expression, int left, int right, int node, int zero) {
    Expression result = expression;
    result.erase(left);
    result.erase(right);
    if (node != zero) {
        auto existing = result.find(node);
        if (existing == result.end())
            result.insert(node);
        else
            result.erase(existing);
    }
    result.erase(zero);
    return result;
}

std::size_t expression_size(const std::vector<Expression> &expressions) {
    std::size_t total = 0;
    for (const auto &expression : expressions) {
        if (total > std::numeric_limits<std::size_t>::max() - expression.size())
            return std::numeric_limits<std::size_t>::max();
        total += expression.size();
    }
    return total;
}

// For D<63, W and 2^D fit in uint64_t with headroom for saturated failure.
// Larger depths use the exactly equivalent earliest-ready/Kraft check instead.
std::uint64_t row_weight(const Expression &expression, const NodeStore &nodes, int max_depth) {
    require(max_depth >= 0 && max_depth < 63, "invalid small completion budget");
    const std::uint64_t budget = std::uint64_t{1} << max_depth;
    std::uint64_t weight = 0;
    for (int term : expression) {
        const int depth = nodes.depths.at(term);
        if (depth > max_depth) return budget + 1;
        const std::uint64_t add = std::uint64_t{1} << depth;
        if (weight > budget - add) return budget + 1;
        weight += add;
    }
    return weight;
}

bool feasible_update(const Expression &expression, std::uint64_t current_weight,
                     int left, int right, int node, int node_depth,
                     const NodeStore &nodes, int max_depth, std::uint64_t &next_weight,
                     int override_node = -1, int override_depth = -1) {
    if (node_depth > max_depth) return false;
    auto arrival = [&](int term) {
        return term == override_node ? override_depth : nodes.depths.at(term);
    };
    const bool node_already_present =
        node >= 0 && node != left && node != right && expression.count(node) != 0;

    if (max_depth >= 63) {
        std::vector<int> arrivals;
        arrivals.reserve(expression.size() + 1);
        for (int term : expression) {
            if (term == left || term == right || (node_already_present && term == node))
                continue;
            arrivals.push_back(nodes.depths.at(term));
        }
        if (!node_already_present && node != nodes.zero()) arrivals.push_back(node_depth);
        next_weight = 0;
        return combine_depth(arrivals) <= max_depth;
    }

    const std::uint64_t budget = std::uint64_t{1} << max_depth;
    if (current_weight > budget) return false;
    const std::uint64_t left_weight = std::uint64_t{1} << arrival(left);
    const std::uint64_t right_weight = std::uint64_t{1} << arrival(right);
    if (current_weight < left_weight) return false;
    std::uint64_t weight = current_weight - left_weight;
    if (weight < right_weight) return false;
    weight -= right_weight;
    const std::uint64_t replacement_weight = std::uint64_t{1} << node_depth;
    if (node_already_present) {
        if (weight < replacement_weight) return false;
        weight -= replacement_weight;
    } else if (node != nodes.zero()) {
        if (weight > budget - replacement_weight) return false;
        weight += replacement_weight;
    }
    next_weight = weight;
    return true;
}

long long candidate_delta(int left, int right, const Candidate &candidate,
                          const NodeStore &nodes, int max_depth) {
    if (max_depth >= 63)
        return candidate.depth - std::max(nodes.depths.at(left), nodes.depths.at(right));
    const auto left_weight = std::uint64_t{1} << nodes.depths.at(left);
    const auto right_weight = std::uint64_t{1} << nodes.depths.at(right);
    const auto replacement_weight = std::uint64_t{1} << candidate.depth;
    return static_cast<long long>(replacement_weight) - static_cast<long long>(left_weight) -
           static_cast<long long>(right_weight);
}

std::vector<int> earliest_ready(const std::vector<Expression> &expressions, NodeStore &nodes,
                                Deadline deadline) {
    std::vector<int> outputs;
    outputs.reserve(expressions.size());
    for (const auto &expression : expressions) {
        check_deadline(deadline);
        std::vector<int> pool(expression.begin(), expression.end());
        while (pool.size() > 1) {
            check_deadline(deadline);
            std::sort(pool.begin(), pool.end(), [&](int a, int b) {
                return nodes.depths.at(a) != nodes.depths.at(b)
                           ? nodes.depths.at(a) < nodes.depths.at(b)
                           : a < b;
            });
            const int left = pool[0], right = pool[1];
            pool.erase(pool.begin(), pool.begin() + 2);
            const int node = nodes.get_node(left, right, deadline);
            if (node == nodes.zero()) continue;
            auto existing = std::find(pool.begin(), pool.end(), node);
            if (existing == pool.end())
                pool.push_back(node);
            else
                pool.erase(existing);
        }
        outputs.push_back(pool.empty() ? nodes.zero() : pool.front());
    }
    return outputs;
}

Circuit greedy_attempt(const BitRows &targets, int inputs, int max_depth,
                       bool lookahead, unsigned seed, Deadline deadline,
                       bool mwsq_tiebreak) {
    NodeStore nodes(inputs);
    std::mt19937 random(seed);
    std::vector<Expression> expressions;
    expressions.reserve(targets.size());
    for (const auto &target : targets) {
        check_deadline(deadline);
        Expression expression;
        for (int bit = 0; bit < inputs; ++bit) {
            check_deadline(deadline);
            if (target.bit(bit)) expression.insert(bit);
        }
        std::vector<int> arrivals(expression.size(), 0);
        require(combine_depth(arrivals) <= max_depth,
                "initial expression exceeds greedy depth limit");
        expressions.push_back(std::move(expression));
    }

    const std::size_t iteration_limit = expression_size(expressions);
    for (std::size_t iteration = 0; iteration < iteration_limit; ++iteration) {
        check_deadline(deadline);
        std::vector<std::uint64_t> weights(expressions.size(), 0);
        if (max_depth < 63) {
            const std::uint64_t budget = std::uint64_t{1} << max_depth;
            for (std::size_t i = 0; i < expressions.size(); ++i) {
                check_deadline(deadline);
                weights[i] = row_weight(expressions[i], nodes, max_depth);
                require(weights[i] <= budget, "greedy lost completion feasibility");
            }
        }

        // Feasible (pair, row) occurrences are collected flat and sorted once
        // per iteration; the previous std::map did an O(log n) lookup plus a
        // node allocation per distinct pair. Groups of equal keys reproduce
        // the map's sorted (left, right) iteration order exactly, and rows
        // inside a group stay in ascending order like the map's uses list.
        struct Occurrence {
            std::uint64_t key;
            int row;
        };
        std::vector<Occurrence> occurrences;
        // NodeStore is immutable until the selected candidate materializes its
        // node after this enumeration, so pair_info results are cached per
        // iteration: a pair occurring in many rows used to recompute the same
        // BitMask xor + hash lookup for every occurrence.
        std::unordered_map<std::uint64_t, PairInfo> pair_cache;
        for (std::size_t i = 0; i < expressions.size(); ++i) {
            const auto &expression = expressions[i];
            for (auto left = expression.begin(); left != expression.end(); ++left) {
                for (auto right = std::next(left); right != expression.end(); ++right) {
                    check_deadline(deadline);
                    const std::uint64_t cache_key = mwsq_pair_key(*left, *right);
                    auto cached = pair_cache.find(cache_key);
                    if (cached == pair_cache.end())
                        cached = pair_cache.emplace(cache_key, nodes.pair_info(*left, *right))
                                     .first;
                    const PairInfo &info = cached->second;
                    std::uint64_t next_weight = 0;
                    if (!feasible_update(expression, weights[i], *left, *right, info.node,
                                         info.depth, nodes, max_depth, next_weight))
                        continue;
                    occurrences.push_back(Occurrence{cache_key, static_cast<int>(i)});
                }
            }
        }
        std::sort(occurrences.begin(), occurrences.end(),
                  [](const Occurrence &a, const Occurrence &b) {
                      if (a.key != b.key) return a.key < b.key;
                      return a.row < b.row;
                  });

        std::vector<Candidate> candidates;
        candidates.reserve(occurrences.size());
        int best_gain = 0;
        for (std::size_t g = 0; g < occurrences.size();) {
            check_deadline(deadline);
            std::size_t group_end = g;
            while (group_end < occurrences.size() &&
                   occurrences[group_end].key == occurrences[g].key)
                ++group_end;
            const PairInfo &info = pair_cache.at(occurrences[g].key);
            Candidate candidate;
            // The set iterates terms in ascending order, so left < right and
            // the packed key decodes back to (left, right).
            candidate.left = static_cast<int>(occurrences[g].key >> 32);
            candidate.right = static_cast<int>(occurrences[g].key & 0xffffffffu);
            candidate.node = info.node;
            candidate.depth = info.depth;
            candidate.mask = info.mask;
            for (std::size_t k = g; k < group_end; ++k)
                candidate.uses.push_back(occurrences[k].row);
            g = group_end;
            candidate.gain = static_cast<int>(candidate.uses.size()) -
                             (candidate.node < 0 ? 1 : 0);
            if (candidate.gain <= 0) continue;
            candidate.delta = candidate_delta(candidate.left, candidate.right, candidate,
                                              nodes, max_depth);
            candidate.score = candidate.gain;
            best_gain = std::max(best_gain, candidate.gain);
            candidates.push_back(std::move(candidate));
        }
        if (candidates.empty()) break;
        check_deadline(deadline);
        std::sort(candidates.begin(), candidates.end(), [](const Candidate &a, const Candidate &b) {
            return std::make_tuple(-a.gain, a.delta, a.depth, a.left, a.right) <
                   std::make_tuple(-b.gain, b.delta, b.depth, b.left, b.right);
        });
        check_deadline(deadline);

        if (lookahead) {
            std::shuffle(candidates.begin(), candidates.end(), random);
            std::stable_sort(candidates.begin(), candidates.end(),
                             [](const Candidate &a, const Candidate &b) {
                                 return std::make_tuple(-a.gain, a.delta, a.depth) <
                                        std::make_tuple(-b.gain, b.delta, b.depth);
                             });
            candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                                            [&](const Candidate &candidate) {
                                                return candidate.gain < best_gain - 1;
                                            }),
                             candidates.end());
            if (candidates.size() > 8) candidates.resize(8);
            for (auto &candidate : candidates) {
                check_deadline(deadline);
                std::map<int, int> future;
                const int candidate_node = candidate.node >= 0
                                               ? candidate.node
                                               : static_cast<int>(nodes.masks.size());
                for (int row_index : candidate.uses) {
                    check_deadline(deadline);
                    const auto next = replaced(expressions[row_index], candidate.left,
                                               candidate.right, candidate_node, nodes.zero());
                    std::uint64_t next_weight = 0;
                    require(feasible_update(expressions[row_index], weights[row_index],
                                            candidate.left, candidate.right, candidate.node,
                                            candidate.depth, nodes, max_depth, next_weight),
                            "selected greedy candidate lost feasibility");
                    for (int other : next) {
                        if (other == candidate_node) continue;
                        check_deadline(deadline);
                        const auto future_mask = candidate.mask ^ nodes.masks.at(other);
                        int future_node = -1;
                        int future_depth = next_depth(candidate.depth, nodes.depths.at(other));
                        auto cached = nodes.by_mask.find(future_mask);
                        if (cached != nodes.by_mask.end() &&
                            nodes.depths.at(cached->second) <= future_depth) {
                            future_node = cached->second;
                            future_depth = nodes.depths.at(future_node);
                        }
                        std::uint64_t ignored = 0;
                        if (feasible_update(next, next_weight, candidate_node, other,
                                            future_node, future_depth, nodes, max_depth, ignored,
                                            candidate_node, candidate.depth))
                            ++future[other];
                    }
                }
                int height = 0;
                for (const auto &entry : future) height = std::max(height, entry.second);
                candidate.score = candidate.gain + std::max(0, height - 1);
            }
        }

        auto rank = [](const Candidate &candidate) {
            return std::make_tuple(candidate.score, candidate.gain, -candidate.delta,
                                   -candidate.depth);
        };
        auto best = rank(candidates.front());
        for (const auto &candidate : candidates) best = std::max(best, rank(candidate));
        std::vector<std::size_t> tied;
        for (std::size_t i = 0; i < candidates.size(); ++i)
            if (rank(candidates[i]) == best) tied.push_back(i);
        check_deadline(deadline);
        std::size_t choice = tied[random() % tied.size()];
        if (mwsq_tiebreak && tied.size() > 1) {
            // MW_SQ lookahead: for every tied candidate, apply the replacement
            // to all rows (an absent node gets the virtual id
            // nodes.masks.size(), exactly like the lookahead pass, so nothing
            // is materialized) and measure the remaining shared structure.
            // The lexicographically largest (max_support, ties_at_max) wins;
            // a surviving tie keeps the seeded random choice.
            const int virtual_node = static_cast<int>(nodes.masks.size());
            std::vector<char> in_uses(expressions.size(), 0);
            std::pair<std::size_t, std::size_t> best_score{0, 0};
            std::vector<std::size_t> winners;
            std::vector<std::uint64_t> keys;
            std::vector<int> replaced_terms;
            for (std::size_t index : tied) {
                check_deadline(deadline);
                const auto &candidate = candidates[index];
                const int candidate_node =
                    candidate.node >= 0 ? candidate.node : virtual_node;
                keys.clear();
                for (int use : candidate.uses) in_uses[use] = 1;
                for (std::size_t r = 0; r < expressions.size(); ++r) {
                    check_deadline(deadline);
                    if (in_uses[r]) {
                        // Terms of expressions[r] after replacing
                        // (left, right) with candidate_node: identical to
                        // replaced() — drop left/right and the zero node,
                        // then toggle candidate_node in or out.
                        replaced_terms.clear();
                        for (int term : expressions[r]) {
                            if (term == candidate.left || term == candidate.right ||
                                term == nodes.zero())
                                continue;
                            replaced_terms.push_back(term);
                        }
                        if (candidate_node != nodes.zero()) {
                            auto it = std::find(replaced_terms.begin(),
                                                replaced_terms.end(), candidate_node);
                            if (it == replaced_terms.end())
                                replaced_terms.push_back(candidate_node);
                            else
                                replaced_terms.erase(it);
                        }
                        for (std::size_t i = 0; i + 1 < replaced_terms.size(); ++i)
                            for (std::size_t j = i + 1; j < replaced_terms.size(); ++j)
                                keys.push_back(mwsq_pair_key(replaced_terms[i],
                                                             replaced_terms[j]));
                    } else {
                        const auto &expr = expressions[r];
                        for (auto a = expr.begin(); a != expr.end(); ++a)
                            for (auto b = std::next(a); b != expr.end(); ++b)
                                keys.push_back(mwsq_pair_key(*a, *b));
                    }
                }
                for (int use : candidate.uses) in_uses[use] = 0;
                check_deadline(deadline);
                const auto score = mwsq_summary_from_keys(keys);
                if (winners.empty() || score > best_score) {
                    best_score = score;
                    winners.clear();
                    winners.push_back(index);
                } else if (score == best_score) {
                    winners.push_back(index);
                }
            }
            choice = winners[random() % winners.size()];
        }
        const auto &selected = candidates[choice];
        const int node = nodes.get_node(selected.left, selected.right, deadline);
        const std::size_t before = expression_size(expressions);
        for (int row_index : selected.uses) {
            check_deadline(deadline);
            expressions[row_index] = replaced(expressions[row_index], selected.left,
                                              selected.right, node, nodes.zero());
        }
        const std::size_t after = expression_size(expressions);
        require(after < before, "greedy replacement did not reduce expression size");
    }

    nodes.circuit.outputs = earliest_ready(expressions, nodes, deadline);
    auto result = nodes.circuit.pruned();
    result.verify(targets);
    require(depth_profile(result).depth <= max_depth, "greedy completion exceeds depth limit");
    return result;
}

} // namespace

Circuit synthesize_completable(const BitRows &targets, int inputs, int max_depth,
                               int restarts, std::uint64_t seed, Deadline deadline,
                               CircuitObserver observer, bool mwsq_tiebreak) {
    check_deadline(deadline);
    require(inputs > 0 && inputs < 1000000, "invalid greedy input count");
    require(max_depth >= 0 && restarts > 0, "invalid greedy options");
    validate_target_width(targets, inputs);
    require(max_depth >= depth_lower_bound(targets), "greedy depth below row-weight bound");

    Circuit best;
    bool found = false;
    std::pair<int, int> best_metric{0, 0};
    for (int restart = 0; restart < restarts; ++restart) {
        try {
            check_deadline(deadline);
            auto candidate = greedy_attempt(targets, inputs, max_depth, restart % 2 == 1,
                                            static_cast<unsigned>(seed + static_cast<std::uint64_t>(restart / 2)),
                                            deadline, mwsq_tiebreak);
            const std::pair<int, int> metric{candidate.cost(), depth_profile(candidate).depth};
            if (!found || metric < best_metric) {
                best = candidate;
                best_metric = metric;
                found = true;
            }
            // Publish each successful attempt immediately; a downstream rewrite
            // may raise a deadline, which the enclosing try handles by keeping
            // the already-found best.
            if (observer) observer(candidate);
        } catch (const SearchDeadlineExceeded &) {
            if (!found) throw;
            break;
        }
    }
    best.verify(targets);
    require(depth_profile(best).depth <= max_depth, "greedy result exceeds depth limit");
    return best;
}

} // namespace xd
