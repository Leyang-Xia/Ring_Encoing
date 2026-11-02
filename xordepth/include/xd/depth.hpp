#pragma once
#include "xd/core.hpp"
#include <chrono>
#include <cstdint>
#include <exception>
#include <unordered_map>
#include <vector>
#include <string>

namespace xd {

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
int depth_lower_bound(const BitRows &targets);

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
    // Optional selective lookahead. Zero disables the extra probing pass;
    // otherwise the first K candidates are simulated for L future steps.
    int lookahead_width = 0;
    int lookahead_depth = 0;
    std::uint64_t seed = 0;
    double timeout_seconds = 180.0;
    // Internal monotonic deadline used by SearchState; not serialized as a
    // result parameter. The public default leaves primitive tests unbounded.
    std::chrono::steady_clock::time_point deadline =
        std::chrono::steady_clock::time_point::max();
};

struct BoundedScheduleResult {
    std::string status;
    Circuit circuit;
    int depth_lower_bound = 0;
    int requested_depth = 0;
    std::size_t explored_states = 0;
    std::size_t pruned_depth = 0;
    std::size_t pruned_bound = 0;
    // Cost of the best externally supplied, verified incumbent that was
    // feasible at the requested depth, or -1 when none was supplied.
    int reference_cost = -1;
    Circuit reference_circuit;
    bool has_reference = false;
};

// Internal control flow for the bounded search budget. Public callers normally
// use the default, unbounded schedule helpers and never observe this type.
struct SearchDeadlineExceeded final : std::exception {
    const char *what() const noexcept override { return "xordepth search deadline exceeded"; }
};

BoundedScheduleResult synthesize_bounded(const BitRows &targets, int inputs,
                                         const BoundedScheduleOptions &options,
                                         const std::vector<Circuit> &incumbents = {});

// Plank-style common-subexpression scheduling with a depth-aware selection pass
// and a balanced (Huffman by arrival time) final accumulation. Every output,
// dependency and live gate is re-verified before returning.
Circuit synthesize_schedule(const BitRows &targets, int inputs, const ScheduleOptions &options,
                            std::uint64_t seed = 0,
                            std::chrono::steady_clock::time_point deadline =
                                std::chrono::steady_clock::time_point::max());

// Re-associate the accumulation trees of a pruned circuit. Only gates with a
// single consumer that are not outputs are folded away, so each affected cone
// has |frontier| - 1 gates before and after: the gate count never grows.
// Every output mask is compared against the original before returning.
Circuit rebalanced(
    const Circuit &circuit,
    std::chrono::steady_clock::time_point deadline =
        std::chrono::steady_clock::time_point::max());
} // namespace xd
