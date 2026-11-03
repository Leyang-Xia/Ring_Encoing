#include "xd/depth.hpp"
#include <algorithm>
#include <array>
#include <set>

namespace xd {
namespace {
using Deadline = std::chrono::steady_clock::time_point;

void check_deadline(Deadline deadline) {
    if (deadline != Deadline::max() && std::chrono::steady_clock::now() >= deadline)
        throw SearchDeadlineExceeded();
}

void index_new(const GreedyNetwork &c, std::unordered_map<BitMask, int, BitMaskHash> &by_mask,
               std::size_t &seen) {
    for (; seen < c.masks.size(); ++seen)
        by_mask.emplace(c.masks[seen], int(seen));
}

int xor_node(GreedyNetwork &c, std::unordered_map<BitMask, int, BitMaskHash> &by_mask,
             std::size_t &seen, int a, int b) {
    int node = c.xor_nodes(a, b);
    index_new(c, by_mask, seen);
    return node;
}

void record_depth(std::unordered_map<int, int> &depth, int node, int value) {
    auto it = depth.find(node);
    if (it == depth.end())
        depth.emplace(node, value);
    else
        it->second = std::max(it->second, value);
}

int accumulate(GreedyNetwork &c, std::unordered_map<BitMask, int, BitMaskHash> &by_mask,
               std::size_t &seen, std::unordered_map<int, int> &depth, std::vector<int> pool,
               Deadline deadline) {
    while (pool.size() > 1) {
        check_deadline(deadline);
        std::sort(pool.begin(), pool.end(), [&](int a, int b) {
            return depth.at(a) != depth.at(b) ? depth.at(a) < depth.at(b) : a < b;
        });
        int a = pool[0], b = pool[1];
        pool.erase(pool.begin(), pool.begin() + 2);
        int node = xor_node(c, by_mask, seen, a, b);
        record_depth(depth, node, 1 + std::max(depth.at(a), depth.at(b)));
        auto duplicate = std::find(pool.begin(), pool.end(), node);
        if (duplicate == pool.end())
            pool.push_back(node);
        else
            pool.erase(duplicate);
    }
    return pool.empty() ? c.inputs : pool.front();
}
} // namespace

Circuit rebalanced(const Circuit &c, std::chrono::steady_clock::time_point deadline) {
    check_deadline(deadline);
    std::unordered_map<int, BitMask> value;
    for (int i = 0; i < c.inputs; ++i)
        value.emplace(i, BitMask::unit(i));
    if (c.has_zero)
        value.emplace(c.zero_node, BitMask());

    std::unordered_map<int, std::array<int, 2>> operands;
    std::unordered_map<int, int> consumers;
    for (const auto &op : c.ops) {
        check_deadline(deadline);
        require(value.emplace(op[0], value.at(op[1]) ^ value.at(op[2])).second,
                "duplicate or unavailable rebalance node");
        operands.emplace(op[0], std::array<int, 2>{op[1], op[2]});
        ++consumers[op[1]];
        ++consumers[op[2]];
    }

    std::set<int> outputs(c.outputs.begin(), c.outputs.end());
    auto foldable = [&](int node) {
        return operands.count(node) && consumers[node] == 1 && !outputs.count(node);
    };

    GreedyNetwork c2(c.inputs);
    std::unordered_map<BitMask, int, BitMaskHash> by_mask;
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
        check_deadline(deadline);
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
        mapped.emplace(node, accumulate(c2, by_mask, seen, depth, pool, deadline));
    }

    std::vector<int> result;
    result.reserve(c.outputs.size());
    for (auto node : c.outputs) {
        check_deadline(deadline);
        auto it = mapped.find(node);
        require(it != mapped.end(), "rebalance output is not available");
        require(c2.masks.at(it->second) == value.at(node), "rebalance changed an output mask");
        result.push_back(it->second);
    }
    auto out = c2.export_circuit(result);
    check_deadline(deadline);
    require(out.cost() <= c.cost(), "rebalance increased the gate count");
    return out;
}
} // namespace xd
