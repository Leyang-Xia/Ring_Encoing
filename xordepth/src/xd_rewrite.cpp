#include "xd/depth.hpp"
#include <algorithm>
#include <set>

namespace xd {
namespace {
using Deadline = std::chrono::steady_clock::time_point;

void check_deadline(Deadline deadline) {
    if (deadline != Deadline::max() && std::chrono::steady_clock::now() >= deadline)
        throw SearchDeadlineExceeded();
}

struct RewriteState {
    Circuit circuit;
    BitValues values;
    std::unordered_map<int, int> latest;
    std::unordered_map<int, int> uses;
    std::set<int> outputs;
};

RewriteState inspect(const Circuit &input, const BitRows &targets, int max_depth,
                     Deadline deadline) {
    check_deadline(deadline);
    RewriteState state;
    state.circuit = input.pruned();
    state.values = state.circuit.verify(targets);
    auto profile = depth_profile(state.circuit);
    require(profile.depth <= max_depth, "rewrite input exceeds depth limit");
    for (auto node : state.circuit.outputs) {
        state.outputs.insert(node);
        state.latest[node] = max_depth;
    }
    for (auto it = state.circuit.ops.rbegin(); it != state.circuit.ops.rend(); ++it) {
        check_deadline(deadline);
        int deadline_at = state.latest.at((*it)[0]) - 1;
        for (int parent : {(*it)[1], (*it)[2]}) {
            auto found = state.latest.find(parent);
            if (found == state.latest.end() || deadline_at < found->second)
                state.latest[parent] = deadline_at;
        }
        ++state.uses[(*it)[1]];
        ++state.uses[(*it)[2]];
    }
    return state;
}

void add_signal(const BitValues &values, std::unordered_map<int, int> &alias,
                std::unordered_map<int, int> &ready,
                std::unordered_map<BitMask, int, BitMaskHash> &by_mask,
                std::vector<int> &available, int node, int node_depth) {
    alias[node] = node;
    ready[node] = node_depth;
    available.push_back(node);
    auto found = by_mask.find(values.at(node));
    if (found == by_mask.end() || node_depth < ready.at(found->second))
        by_mask[values.at(node)] = node;
}

bool rebuild_without(const RewriteState &state, int victim, int max_depth, Circuit &result,
                     Deadline deadline) {
    Circuit next;
    next.inputs = state.circuit.inputs;
    next.has_zero = state.circuit.has_zero;
    next.zero_node = state.circuit.zero_node;
    std::unordered_map<int, int> alias, ready;
    std::unordered_map<BitMask, int, BitMaskHash> by_mask;
    std::vector<int> available;
    for (int input = 0; input < next.inputs; ++input)
        add_signal(state.values, alias, ready, by_mask, available, input, 0);
    if (next.has_zero)
        add_signal(state.values, alias, ready, by_mask, available, next.zero_node, 0);

    for (const auto &op : state.circuit.ops) {
        check_deadline(deadline);
        int node = op[0];
        if (node == victim)
            continue;
        const BitMask &target = state.values.at(node);
        auto equal = by_mask.find(target);
        if (equal != by_mask.end() && ready.at(equal->second) <= state.latest.at(node)) {
            alias[node] = equal->second;
            continue;
        }

        int left = -1, right = -1;
        int best_depth = state.latest.at(node) + 1;
        if (alias.count(op[1]) && alias.count(op[2])) {
            left = alias.at(op[1]);
            right = alias.at(op[2]);
            best_depth = 1 + std::max(ready.at(left), ready.at(right));
        } else {
            for (int candidate : available) {
                check_deadline(deadline);
                auto complement = by_mask.find(state.values.at(candidate) ^ target);
                if (complement == by_mask.end())
                    continue;
                int other = complement->second;
                if (candidate == other)
                    continue;
                int candidate_depth = 1 + std::max(ready.at(candidate), ready.at(other));
                if (candidate_depth < best_depth) {
                    left = candidate;
                    right = other;
                    best_depth = candidate_depth;
                }
            }
        }
        if (left < 0 || right < 0 || best_depth > state.latest.at(node))
            return false;
        next.ops.push_back({node, left, right});
        add_signal(state.values, alias, ready, by_mask, available, node, best_depth);
    }
    for (int output : state.circuit.outputs) {
        auto found = alias.find(output);
        if (found == alias.end())
            return false;
        next.outputs.push_back(found->second);
    }
    result = next.pruned();
    BitRows expected;
    expected.reserve(state.circuit.outputs.size());
    for (int output : state.circuit.outputs)
        expected.push_back(state.values.at(output));
    result.verify(expected);
    require(depth_profile(result).depth <= max_depth, "rewrite exceeds depth limit");
    return true;
}
} // namespace

Circuit rewrite_depth_gates(const Circuit &circuit, const BitRows &targets, int max_depth,
                            int max_trials, Deadline deadline, int *used_trials) {
    if (used_trials) *used_trials = 0;
    require(max_depth >= 0, "negative rewrite depth limit");
    require(max_trials >= 0, "negative rewrite trial budget");
    RewriteState state = inspect(circuit, targets, max_depth, deadline);
    Circuit current = state.circuit;
    if (max_trials == 0 || current.ops.empty())
        return current;

    int trials = 0;
    for (int pass = 0; pass < 4 && trials < max_trials; ++pass) {
        check_deadline(deadline);
        bool changed = false;
        std::vector<int> victims;
        for (const auto &op : state.circuit.ops)
            if (!state.outputs.count(op[0]))
                victims.push_back(op[0]);
        std::stable_sort(victims.begin(), victims.end(), [&](int a, int b) {
            return state.uses[a] < state.uses[b];
        });
        for (int victim : victims) {
            if (trials >= max_trials)
                break;
            ++trials;
            // Track partial usage so a caller that catches a deadline exception
            // can still account for the trials actually spent.
            if (used_trials) *used_trials = trials;
            Circuit candidate;
            if (!rebuild_without(state, victim, max_depth, candidate, deadline))
                continue;
            candidate.verify(targets);
            if (candidate.cost() < current.cost()) {
                current = std::move(candidate);
                state = inspect(current, targets, max_depth, deadline);
                changed = true;
                break;
            }
        }
        if (!changed)
            break;
    }
    current.verify(targets);
    require(depth_profile(current).depth <= max_depth, "rewrite final depth limit");
    return current;
}
} // namespace xd
