#pragma once
#include "ringlab/core.hpp"
#include <cstdint>
#include <unordered_map>
#include <vector>
#include <string>

namespace xd {
using ringlab::Circuit;
using ringlab::Mask;
using ringlab::Rows;

// Critical path of a XOR network, counted in two-input XOR gates from a primary
// input (or the constant zero) to an output. Wiring and aliases cost no gate.
struct DepthProfile {
    int gates = 0;
    int depth = 0;
    int max_fanout = 0;
    std::unordered_map<int, int> node_depth;
    std::vector<int> output_depth;
};
DepthProfile depth_profile(const Circuit &circuit);

// One gate can at most double the number of distinct inputs it covers, so an
// output of weight w needs at least ceil(log2 w) levels. This bound ignores
// sharing and fanout, so it is weak: it is a floor, not a target.
int depth_lower_bound(const Rows &targets);

struct ScheduleOptions {
    // A pair of signals must appear in at least this many output rows before it
    // is materialized as a shared node. 0 or 1 disables sharing entirely.
    int share_threshold = 2;
    // Break equal-frequency ties towards the earliest-arriving pair.
    bool depth_aware = false;
};

struct BoundedScheduleOptions {
    int max_depth = 0;
    int beam_width = 64;
    int branch_width = 16;
    int max_representations = 8;
    int restarts = 8;
    std::uint64_t seed = 0;
    double timeout_seconds = 180.0;
};

struct BoundedScheduleResult {
    std::string status;
    Circuit circuit;
    int depth_lower_bound = 0;
    int requested_depth = 0;
    std::size_t explored_states = 0;
    std::size_t pruned_depth = 0;
    std::size_t pruned_bound = 0;
};

BoundedScheduleResult synthesize_bounded(const Rows &targets, int inputs,
                                         const BoundedScheduleOptions &options);

// Plank-style common-subexpression scheduling with a depth-aware selection pass
// and a balanced (Huffman by arrival time) final accumulation. Every output,
// dependency and live gate is re-verified before returning.
Circuit synthesize_schedule(const Rows &targets, int inputs, const ScheduleOptions &options,
                            std::uint64_t seed = 0);

// Re-associate the accumulation trees of a pruned circuit. Only gates with a
// single consumer that are not outputs are folded away, so each affected cone
// has |frontier| - 1 gates before and after: the gate count never grows.
// Every output mask is compared against the original before returning.
Circuit rebalanced(const Circuit &circuit);
} // namespace xd
