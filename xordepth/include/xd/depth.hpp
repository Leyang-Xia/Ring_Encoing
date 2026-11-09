#pragma once
#include "xd/core.hpp"
#include <chrono>
#include <cstdint>
#include <exception>
#include <functional>
#include <unordered_map>
#include <utility>
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

// MW_SQ-style support summary over node-id rows: counts every unordered pair of
// node ids within each row and returns (maximum support, number of pairs that
// attain it). Empty rows and rows with fewer than two entries contribute
// nothing, so an all-trivial input yields (0, 0). Cost is O(sum t(t-1)/2).
std::pair<std::size_t, std::size_t> mwsq_support_summary(
    const std::vector<std::vector<int>> &rows);

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
    // Total completable-greedy attempts, alternating plain and lookahead.
    int greedy_restarts = 2;
    // MW_SQ-style lookahead rescoring of greedy tie sets. Off by default; when
    // off every greedy attempt is bit-identical to the previous behavior.
    bool greedy_mwsq = false;
    // Optional selective lookahead. Zero disables the extra probing pass;
    // otherwise the first K candidates are simulated for L future steps.
    int lookahead_width = 0;
    int lookahead_depth = 0;
    int residual_pairs = 0;
    int residual_targets = 0;
    bool allow_empty_toggle = false;
    int rewrite_trials = 0;
    // Repeatable external incumbent files whose verified circuits join the
    // unified candidate pool with zero search time.
    std::vector<std::string> incumbent_paths;
    // Worker threads for the greedy portfolio and per-layer Beam expansion.
    // 1 keeps the historical sequential path bit-for-bit; larger values are
    // deterministic as well: results are collected and merged in the same
    // order the sequential path produces (wall-clock timeout runs excepted,
    // as they already are machine-dependent).
    int threads = 1;
    std::size_t max_states = 10000;
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
    std::size_t toggle_candidates = 0;
    std::size_t bridge_candidates = 0;
    std::size_t diversity_buckets = 0;
    std::size_t state_limit = 0;
    std::size_t round_limit = 0;
    std::size_t max_rounds = 0;
    std::size_t published_candidates = 0;
    std::size_t incumbents_used = 0;
    std::size_t incumbents_skipped = 0;
    std::size_t rewrite_trials_used = 0;
    std::string best_source;
};

// Internal control flow for the bounded search budget. Public callers normally
// use the default, unbounded schedule helpers and never observe this type.
struct SearchDeadlineExceeded final : std::exception {
    const char *what() const noexcept override { return "xordepth search deadline exceeded"; }
};

BoundedScheduleResult synthesize_bounded(const BitRows &targets, int inputs,
                                         const BoundedScheduleOptions &options);

// Called once for every successful greedy attempt so a caller can publish each
// result immediately instead of waiting for the whole portfolio to finish.
using CircuitObserver = std::function<void(const Circuit &)>;

// Strictly term-reducing common-subexpression elimination. Every accepted pair
// preserves a feasible earliest-ready completion within max_depth. Restarts
// alternate plain and bounded-lookahead attempts with the same seed; each pair
// of attempts advances the seed by one.
Circuit synthesize_completable(
    const BitRows &targets, int inputs, int max_depth, int restarts, std::uint64_t seed = 0,
    std::chrono::steady_clock::time_point deadline =
        std::chrono::steady_clock::time_point::max(),
    CircuitObserver observer = {},
    // MW_SQ-style rescoring of tied greedy candidates; false keeps the
    // historical random choice among ties.
    bool mwsq_tiebreak = false,
    // Worker threads for independent restarts; attempts are published in
    // restart order either way.
    int threads = 1);

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

// Try to remove non-output shared gates by rewriting all consumers in one
// transaction. The returned circuit preserves the target masks and hard depth
// limit, and is accepted only when its live XOR count is strictly smaller.
// used_trials reports the exact number of trials consumed so callers can share
// one total rewrite budget across several invocations.
Circuit rewrite_depth_gates(const Circuit &circuit, const BitRows &targets, int max_depth,
                            int max_trials = 400,
                            std::chrono::steady_clock::time_point deadline =
                                std::chrono::steady_clock::time_point::max(),
                            int *used_trials = nullptr);
} // namespace xd
