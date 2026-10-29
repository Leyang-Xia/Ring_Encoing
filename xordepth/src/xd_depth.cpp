#include "xd/depth.hpp"
#include <algorithm>
#include <array>
#include <climits>
#include <set>

namespace xd {
namespace {
using ringlab::GreedyCircuit;
using ringlab::MaskHash;
using ringlab::require;

// GreedyCircuit deduplicates by mask but exposes only its mask vector, so keep a
// local node index in step with it instead of reaching into the class.
void index_new(const GreedyCircuit &c, std::unordered_map<Mask, int, MaskHash> &by_mask,
               std::size_t &seen) {
    for (; seen < c.masks.size(); ++seen)
        by_mask.emplace(c.masks[seen], int(seen));
}
int xor_node(GreedyCircuit &c, std::unordered_map<Mask, int, MaskHash> &by_mask,
             std::size_t &seen, int a, int b) {
    int node = c.xor_nodes(a, b);
    index_new(c, by_mask, seen);
    return node;
}
void update_pairs(std::unordered_map<std::uint64_t, int> &freq, const std::vector<int> &row,
                  int delta) {
    for (std::size_t i = 0; i < row.size(); ++i)
        for (std::size_t j = i + 1; j < row.size(); ++j) {
            auto key = (std::uint64_t(row[i]) << 32) | unsigned(row[j]);
            if (delta > 0) {
                ++freq[key];
            } else {
                auto it = freq.find(key);
                require(it != freq.end() && it->second > 0, "invalid pair frequency");
                if (--it->second == 0)
                    freq.erase(it);
            }
        }
}
void record_depth(std::unordered_map<int, int> &depth, int node, int value) {
    auto it = depth.find(node);
    if (it == depth.end())
        depth.emplace(node, value);
    else
        // A reused node keeps the depth it was created with; never report better.
        it->second = std::max(it->second, value);
}
// Huffman accumulation on arrival depth: repeatedly combine the two earliest
// ready signals. This minimizes the depth of one sum for a fixed term set and
// uses exactly (#terms - 1) gates before sharing is accounted for.
int accumulate(GreedyCircuit &c, std::unordered_map<Mask, int, MaskHash> &by_mask,
               std::size_t &seen, std::unordered_map<int, int> &depth, std::vector<int> pool) {
    while (pool.size() > 1) {
        std::sort(pool.begin(), pool.end(), [&](int a, int b) {
            return depth.at(a) != depth.at(b) ? depth.at(a) < depth.at(b) : a < b;
        });
        int a = pool[0], b = pool[1];
        pool.erase(pool.begin(), pool.begin() + 2);
        int node = xor_node(c, by_mask, seen, a, b);
        record_depth(depth, node, 1 + std::max(depth.at(a), depth.at(b)));
        // Equal masks cancel: drop the copy already in the row instead of XORing
        // a node with itself.
        auto duplicate = std::find(pool.begin(), pool.end(), node);
        if (duplicate == pool.end())
            pool.push_back(node);
        else
            pool.erase(duplicate);
    }
    return pool.empty() ? c.inputs : pool.front();
}
} // namespace

DepthProfile depth_profile(const Circuit &c) {
    require(c.inputs > 0 && c.inputs < 1000000, "invalid circuit input count");
    DepthProfile p;
    p.gates = c.cost();
    for (int i = 0; i < c.inputs; ++i)
        p.node_depth.emplace(i, 0);
    if (c.has_zero)
        require(p.node_depth.emplace(c.zero_node, 0).second, "zero node aliases input");
    std::unordered_map<int, int> fanout;
    for (const auto &op : c.ops) {
        auto a = p.node_depth.find(op[1]), b = p.node_depth.find(op[2]);
        require(!p.node_depth.count(op[0]) && a != p.node_depth.end() && b != p.node_depth.end(),
                "duplicate or unavailable depth node");
        p.node_depth.emplace(op[0], 1 + std::max(a->second, b->second));
        ++fanout[op[1]];
        ++fanout[op[2]];
    }
    for (const auto &[node, count] : fanout)
        p.max_fanout = std::max(p.max_fanout, count);
    p.output_depth.reserve(c.outputs.size());
    for (auto node : c.outputs) {
        auto it = p.node_depth.find(node);
        require(it != p.node_depth.end(), "unavailable depth output");
        p.output_depth.push_back(it->second);
        p.depth = std::max(p.depth, it->second);
    }
    return p;
}

int depth_lower_bound(const Rows &targets) {
    int bound = 0;
    for (const auto &row : targets) {
        std::size_t weight = std::size_t(row.count());
        int levels = 0;
        for (std::size_t covered = 1; covered < weight; covered *= 2)
            ++levels;
        bound = std::max(bound, levels);
    }
    return bound;
}

Circuit synthesize_schedule(const Rows &targets, int inputs, const ScheduleOptions &options,
                            std::uint64_t seed) {
    require(options.share_threshold >= 0, "share threshold must be nonnegative");
    GreedyCircuit c(inputs);
    std::unordered_map<Mask, int, MaskHash> by_mask;
    std::size_t seen = 0;
    index_new(c, by_mask, seen);
    std::unordered_map<int, int> depth;
    for (int i = 0; i <= inputs; ++i)
        depth.emplace(i, 0);
    std::vector<std::vector<int>> expressions;
    for (const auto &mask : targets) {
        require(mask.bits() <= std::size_t(inputs), "mask exceeds schedule input count");
        std::vector<int> row;
        for (int bit = 0; bit < inputs; ++bit)
            if (mask.bit(bit))
                row.push_back(bit);
        expressions.push_back(std::move(row));
    }
    std::unordered_map<std::uint64_t, int> freq;
    for (const auto &row : expressions)
        update_pairs(freq, row, 1);
    ringlab::SeedRandom random(seed);
    const int threshold = std::max(options.share_threshold, 2);
    while (options.share_threshold > 1) {
        int best = 0;
        for (const auto &[pair, count] : freq)
            best = std::max(best, count);
        if (best < threshold)
            break;
        std::vector<std::uint64_t> choices;
        for (const auto &[pair, count] : freq)
            if (count == best)
                choices.push_back(pair);
        std::sort(choices.begin(), choices.end());
        std::uint64_t pair = choices[0];
        if (options.depth_aware) {
            int arrival = INT_MAX;
            for (auto candidate : choices) {
                int a = int(candidate >> 32), b = int(std::uint32_t(candidate));
                int ready = 1 + std::max(depth.at(a), depth.at(b));
                if (ready < arrival) {
                    arrival = ready;
                    pair = candidate;
                }
            }
        } else if (seed) {
            pair = choices[random.below(choices.size())];
        }
        int a = int(pair >> 32), b = int(std::uint32_t(pair));
        int node = xor_node(c, by_mask, seen, a, b);
        record_depth(depth, node, 1 + std::max(depth.at(a), depth.at(b)));
        for (auto &row : expressions)
            if (std::binary_search(row.begin(), row.end(), a) &&
                std::binary_search(row.begin(), row.end(), b)) {
                update_pairs(freq, row, -1);
                row.erase(std::lower_bound(row.begin(), row.end(), a));
                row.erase(std::lower_bound(row.begin(), row.end(), b));
                auto it = std::lower_bound(row.begin(), row.end(), node);
                if (it != row.end() && *it == node)
                    row.erase(it);
                else
                    row.insert(it, node);
                update_pairs(freq, row, 1);
            }
    }
    std::vector<int> outputs;
    outputs.reserve(expressions.size());
    for (const auto &row : expressions) {
        Mask mask;
        for (auto node : row)
            mask ^= c.masks.at(node);
        auto it = by_mask.find(mask);
        if (it != by_mask.end()) {
            outputs.push_back(it->second);
            continue;
        }
        outputs.push_back(accumulate(c, by_mask, seen, depth, row));
    }
    auto result = c.export_circuit(outputs);
    result.verify(targets);
    return result;
}

Circuit rebalanced(const Circuit &c) {
    std::unordered_map<int, Mask> value;
    for (int i = 0; i < c.inputs; ++i)
        value.emplace(i, Mask::unit(i));
    if (c.has_zero)
        value.emplace(c.zero_node, Mask());
    std::unordered_map<int, std::array<int, 2>> operands;
    std::unordered_map<int, int> consumers;
    for (const auto &op : c.ops) {
        require(value.emplace(op[0], value.at(op[1]) ^ value.at(op[2])).second,
                "duplicate or unavailable rebalance node");
        operands.emplace(op[0], std::array<int, 2>{op[1], op[2]});
        ++consumers[op[1]];
        ++consumers[op[2]];
    }
    std::set<int> outputs(c.outputs.begin(), c.outputs.end());
    // A gate whose value only one consumer needs, and that is not itself an
    // output, can be folded into that consumer's accumulation tree.
    auto foldable = [&](int node) {
        return operands.count(node) && consumers[node] == 1 && !outputs.count(node);
    };
    GreedyCircuit c2(c.inputs);
    std::unordered_map<Mask, int, MaskHash> by_mask;
    std::size_t seen = 0;
    index_new(c2, by_mask, seen);
    std::unordered_map<int, int> depth;
    for (int i = 0; i <= c.inputs; ++i)
        depth.emplace(i, 0);
    std::unordered_map<int, int> mapped;
    for (int i = 0; i < c.inputs; ++i)
        mapped.emplace(i, i);
    if (c.has_zero)
        mapped.emplace(c.zero_node, c.inputs);
    std::vector<int> pool, work;
    for (const auto &op : c.ops) {
        int node = op[0];
        if (foldable(node) || mapped.count(node))
            continue;
        pool.clear();
        work.assign({operands.at(node)[0], operands.at(node)[1]});
        while (!work.empty()) {
            int current = work.back();
            work.pop_back();
            if (foldable(current)) {
                work.push_back(operands.at(current)[0]);
                work.push_back(operands.at(current)[1]);
                continue;
            }
            int fresh = mapped.at(current);
            auto duplicate = std::find(pool.begin(), pool.end(), fresh);
            if (duplicate == pool.end())
                pool.push_back(fresh);
            else
                pool.erase(duplicate);
        }
        mapped.emplace(node, accumulate(c2, by_mask, seen, depth, pool));
    }
    std::vector<int> result;
    result.reserve(c.outputs.size());
    for (auto node : c.outputs) {
        auto it = mapped.find(node);
        require(it != mapped.end(), "rebalance output is not available");
        require(c2.masks.at(it->second) == value.at(node), "rebalance changed an output mask");
        result.push_back(it->second);
    }
    auto out = c2.export_circuit(result);
    require(out.cost() <= c.cost(), "rebalance increased the gate count");
    return out;
}
} // namespace xd
