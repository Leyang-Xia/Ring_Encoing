#include "xd/search.hpp"
#include "xd/experiments.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
using namespace xd;

namespace fs = std::filesystem;
namespace {
BoundedScheduleOptions options(int depth) {
    BoundedScheduleOptions o;
    o.max_depth = depth;
    o.beam_width = 4;
    o.branch_width = 8;
    o.max_representations = 8;
    o.restarts = 1;
    o.greedy_restarts = 0;
    return o;
}
bool has_pair(const std::vector<Candidate> &candidates, int left, int right) {
    for (const auto &candidate : candidates)
        if (candidate.left == left && candidate.right == right)
            return true;
    return false;
}

void exercise_json_uint64() {
    require(Json::parse("0").uint64() == 0, "uint64 JSON accepts zero");
    require(Json::parse("18446744073709551615").uint64() ==
                std::numeric_limits<std::uint64_t>::max(),
            "uint64 JSON accepts UINT64_MAX exactly");
    require(Json(std::numeric_limits<std::uint64_t>::max()).uint64() ==
                std::numeric_limits<std::uint64_t>::max(),
            "uint64 JSON roundtrips UINT64_MAX");
    for (const std::string text : {"-1", "-0", "1.0", "1e0", "18446744073709551616"}) {
        bool rejected = false;
        try {
            (void)Json::parse(text).uint64();
        } catch (const std::exception &) {
            rejected = true;
        }
        require(rejected, "uint64 JSON rejects non-integral or overflowing value: " + text);
    }
}

void exercise_mwsq_support_summary() {
    auto empty = mwsq_support_summary({});
    require(empty.first == 0 && empty.second == 0,
            "mwsq summary of empty rows is (0,0)");
    auto single = mwsq_support_summary({{5}});
    require(single.first == 0 && single.second == 0,
            "mwsq summary of a single-element row is (0,0)");
    auto disjoint = mwsq_support_summary({{0, 1}, {2, 3}});
    require(disjoint.first == 1 && disjoint.second == 2,
            "mwsq summary counts every unordered pair once without repeats");
    auto tied = mwsq_support_summary({{0, 1}, {0, 1}, {2, 3}, {2, 3}, {0, 2}});
    require(tied.first == 2 && tied.second == 2,
            "mwsq summary reports the maximum support and the ties at the maximum");
}

// MW_SQ-style greedy tie rescore: verified, depth-feasible, seed-deterministic.
void exercise_mwsq_greedy() {
    const auto no_deadline = std::chrono::steady_clock::time_point::max();
    const BitRows small{BitMask(7), BitMask(11), BitMask(), BitMask(7), BitMask(1)};
    auto first = synthesize_completable(small, 4, 2, 2, 1, no_deadline, {}, true);
    auto second = synthesize_completable(small, 4, 2, 2, 1, no_deadline, {}, true);
    first.verify(small);
    second.verify(small);
    require(depth_profile(first).depth <= 2, "mwsq greedy small-matrix depth limit");
    require(first.cost() == second.cost(),
            "mwsq greedy is deterministic for a fixed seed");

    int inputs = 0;
    auto rows = read_matrix("data/matrices/vandermonde_8to17_direct.txt", inputs);
    require(inputs == 64 && rows.size() == 24, "Vandermonde mwsq dimensions");
    auto mwsq = synthesize_completable(rows, inputs, 5, 2, 1, no_deadline, {}, true);
    mwsq.verify(rows);
    auto profile = depth_profile(mwsq);
    require(profile.depth <= 5, "mwsq greedy Vandermonde respects D=5");
    require(mwsq.cost() <= 223, "mwsq greedy Vandermonde cost is at most 223 (actual " +
                                     std::to_string(mwsq.cost()) + ")");
    std::cout << "mwsq greedy Vandermonde: cost=" << mwsq.cost()
              << " depth=" << profile.depth << '\n';

    auto mwsq_options = options(5);
    mwsq_options.seed = 1;
    mwsq_options.greedy_restarts = 2;
    mwsq_options.max_states = 1;
    mwsq_options.timeout_seconds = 0;
    mwsq_options.greedy_mwsq = true;
    auto bounded_first = synthesize_bounded(rows, inputs, mwsq_options);
    auto bounded_second = synthesize_bounded(rows, inputs, mwsq_options);
    require(bounded_first.status == "verified", "bounded mwsq Vandermonde is verified");
    bounded_first.circuit.verify(rows);
    require(depth_profile(bounded_first.circuit).depth <= 5,
            "bounded mwsq Vandermonde respects D=5");
    require(bounded_first.circuit.cost() == bounded_second.circuit.cost(),
            "bounded mwsq is deterministic for a fixed seed");
    std::cout << "bounded mwsq Vandermonde: cost=" << bounded_first.circuit.cost()
              << " depth=" << depth_profile(bounded_first.circuit).depth << '\n';
}

void exercise_completable_greedy() {
    const BitRows small{BitMask(7), BitMask(11), BitMask(), BitMask(7), BitMask(1)};
    auto circuit = synthesize_completable(small, 4, 2, 2, 1);
    circuit.verify(small);
    auto profile = depth_profile(circuit);
    require(profile.depth <= 2, "greedy small-matrix depth limit");
    require(profile.output_depth.size() == small.size() &&
                profile.output_depth[2] == 0 && profile.output_depth[4] == 0,
            "greedy handles zero rows and input aliases");
    require(circuit.outputs.size() == small.size(), "greedy preserves duplicate outputs");

    BitMask wide;
    for (int bit = 0; bit < 80; ++bit) wide.toggle(bit);
    const BitRows wide_rows{wide, BitMask(), wide};
    auto wide_circuit = synthesize_completable(wide_rows, 80, 63, 2, 7);
    wide_circuit.verify(wide_rows);
    require(depth_profile(wide_circuit).depth <= 63,
            "greedy completion budget is safe at D>=63");

    bool deadline_seen = false;
    try {
        (void)synthesize_completable(small, 4, 2, 2, 1,
                                     std::chrono::steady_clock::now() - std::chrono::seconds(1));
    } catch (const SearchDeadlineExceeded &) {
        deadline_seen = true;
    }
    require(deadline_seen, "greedy observes an expired deadline");

    int inputs = 0;
    auto rows = read_matrix("data/matrices/vandermonde_8to17_direct.txt", inputs);
    require(inputs == 64 && rows.size() == 24, "Vandermonde regression dimensions");
    auto direct = synthesize_completable(rows, inputs, 5, 2, 1);
    direct.verify(rows);
    auto direct_profile = depth_profile(direct);
    require(direct_profile.depth <= 5, "direct greedy Vandermonde respects D=5");
    require(direct.cost() <= 223, "direct greedy Vandermonde cost is at most 223 (actual " +
                                      std::to_string(direct.cost()) + ")");
    std::cout << "direct greedy Vandermonde: cost=" << direct.cost()
              << " depth=" << direct_profile.depth << '\n';

    auto greedy_options = options(5);
    greedy_options.seed = 1;
    greedy_options.greedy_restarts = 2;
    greedy_options.max_states = 1;
    greedy_options.timeout_seconds = 0;
    auto greedy = synthesize_bounded(rows, inputs, greedy_options);
    require(greedy.status == "verified", "greedy Vandermonde D=5 is verified");
    greedy.circuit.verify(rows);
    auto greedy_profile = depth_profile(greedy.circuit);
    require(greedy_profile.depth <= 5, "greedy Vandermonde respects D=5");
    require(greedy.circuit.cost() <= 223, "bounded greedy Vandermonde cost is at most 223");
    require(greedy.explored_states <= 1, "greedy incumbent respects the Beam state budget");
    std::cout << "greedy Vandermonde: cost=" << greedy.circuit.cost()
              << " depth=" << greedy_profile.depth << '\n';

    auto control_options = greedy_options;
    control_options.greedy_restarts = 0;
    auto control = synthesize_bounded(rows, inputs, control_options);
    require(control.status == "timeout" && control.explored_states == 1,
            "zero greedy restarts preserves the pure-Beam control");
}

Circuit make_circuit8(std::vector<std::array<int, 3>> ops, std::vector<int> outputs) {
    Circuit c;
    c.inputs = 8;
    c.has_zero = false;
    c.zero_node = -1;
    c.ops = std::move(ops);
    c.outputs = std::move(outputs);
    return c;
}

std::string write_incumbent_file(const std::string &path,
                                 const std::vector<Circuit> &circuits) {
    auto arr = Json::array();
    for (const auto &c : circuits) arr.arr.push_back(c.json());
    Json j = Json::object({{"candidates", arr}});
    j.write(path);
    return path;
}

// Unified candidate pool: every offered circuit must be verified and within the
// hard depth bound, and the best is chosen by (cost, depth) lexicographic order.
void exercise_candidate_pool() {
    // All compute mask 255 = 0^1^...^7 over 8 inputs.
    auto balanced = make_circuit8(
        {{8, 0, 1}, {9, 2, 3}, {10, 4, 5}, {11, 6, 7}, {12, 8, 9}, {13, 10, 11},
         {14, 12, 13}},
        {14});
    auto medium = make_circuit8(
        {{8, 0, 1}, {9, 2, 3}, {10, 8, 9}, {11, 4, 5}, {12, 11, 10}, {13, 12, 6},
         {14, 13, 7}},
        {14});
    auto chain = make_circuit8(
        {{8, 0, 1}, {9, 8, 2}, {10, 9, 3}, {11, 10, 4}, {12, 11, 5}, {13, 12, 6},
         {14, 13, 7}},
        {14});
    auto wrong = make_circuit8(
        {{8, 0, 1}, {9, 8, 2}, {10, 9, 3}, {11, 10, 4}, {12, 11, 5}, {13, 12, 6}},
        {13});
    require(depth_profile(balanced).depth == 3, "balanced circuit depth 3");
    require(depth_profile(medium).depth == 5, "medium circuit depth 5");
    require(depth_profile(chain).depth == 7, "chain circuit depth 7");

    const auto path = (fs::temp_directory_path() / "xd_unified_pool.json").string();
    write_incumbent_file(path, {balanced, medium, chain, balanced, wrong});

    const BitRows targets{BitMask(std::uint64_t{255})};
    BoundedScheduleOptions o;
    o.max_depth = 5;
    o.greedy_restarts = 0;
    o.timeout_seconds = 1e-9;  // effectively pre-expired: only incumbents load
    o.incumbent_paths = {path};
    auto result = synthesize_bounded(targets, 8, o);

    require(result.status == "verified", "pool returns a verified incumbent");
    result.circuit.verify(targets);
    require(result.circuit.cost() == 7, "pool best cost is 7");
    require(depth_profile(result.circuit).depth == 3,
            "pool picks lower depth at equal cost");
    require(result.published_candidates == 2, "only valid incumbents published");
    require(result.incumbents_used == 2, "two distinct valid incumbents used");
    require(result.incumbents_skipped == 3, "duplicate / wrong / over-depth skipped");
    require(result.best_source == "incumbent", "best source recorded");
    fs::remove(path);
}

// A pre-expired deadline must not discard a verified external incumbent.
void exercise_preexpired_incumbent() {
    int inputs = 0;
    auto rows = read_matrix("data/matrices/generated/dense_100x100_seed123.txt", inputs);
    require(inputs == 100 && rows.size() == 100, "dense100 dimensions");
    BoundedScheduleOptions o;
    o.max_depth = 7;
    o.greedy_restarts = 0;
    o.timeout_seconds = 1e-9;  // pre-expired deadline
    o.incumbent_paths = {
        "data/results/greedy_portfolio_full_20261104/dense_100x100_seed123.v1.json"};
    auto result = synthesize_bounded(rows, inputs, o);
    require(result.status == "verified", "pre-expired deadline keeps V1 incumbent");
    result.circuit.verify(rows);
    require(result.circuit.cost() == 1809, "V1 incumbent cost is 1809");
    require(depth_profile(result.circuit).depth <= 7, "V1 incumbent respects D=7");
    require(result.incumbents_used >= 1, "V1 incumbent counted as used");
    require(result.best_source == "incumbent", "V1 incumbent is the best source");
}

// Regression: the original 5s full-budget configuration must now publish the
// cheap CSE baseline before the greedy portfolio eats the wall-clock budget.
void exercise_dense100_5s() {
    int inputs = 0;
    auto rows = read_matrix("data/matrices/generated/dense_100x100_seed123.txt", inputs);
    BoundedScheduleOptions o;
    o.max_depth = 7;
    o.beam_width = 4;
    o.branch_width = 4;
    o.max_representations = 4;
    o.restarts = 2;
    o.greedy_restarts = 2;
    o.rewrite_trials = 400;
    o.seed = 1;
    o.timeout_seconds = 5.0;
    o.max_states = 256;
    auto result = synthesize_bounded(rows, inputs, o);
    require(result.status == "verified", "dense100 5s is verified");
    result.circuit.verify(rows);
    auto profile = depth_profile(result.circuit);
    require(profile.depth <= 7, "dense100 5s respects D=7");
    require(result.circuit.cost() <= 1809,
            "dense100 5s cost is at most 1809 (actual " +
                std::to_string(result.circuit.cost()) + ")");
    std::cout << "dense100 5s: cost=" << result.circuit.cost()
              << " depth=" << profile.depth << " source=" << result.best_source << '\n';
}
}
int main() {
    exercise_json_uint64();
    require(combine_depth({0, 0, 0, 0}) == 2, "balanced four-leaf depth");
    require(!dominates(std::vector<int>{3, 2}, std::vector<int>{4, 1}), "term ordering");
    require(dominates(std::vector<int>{3, 2}, std::vector<int>{4, 3}), "term-and-depth dominance");
    Representation a{{3, 2}, {}, 2}, b{{4, 3}, {}, 3};
    require(dominates(a, b), "term-and-depth dominance");
    require(normalize_terms({4, 2, 4}) == std::vector<int>({2}), "duplicate terms cancel");
    require(candidate_depth(2, 2) == 3, "candidate depth");
    require(!candidate_allowed(2, 2, 2), "hard candidate depth filter");
    require(completion_round_limit(0, 0) == 0, "zero outputs need no completion rounds");
    require(completion_round_limit(1, 0) == 1, "one output at depth zero needs one round");
    require(completion_round_limit(2, 2) == 14, "completion rounds include both output trees");
    require(completion_round_limit(std::numeric_limits<std::size_t>::max(), 1) ==
                std::numeric_limits<std::size_t>::max(),
            "completion round multiplication saturates");
    require(completion_round_limit(1, std::numeric_limits<std::size_t>::digits) ==
                std::numeric_limits<std::size_t>::max(),
            "completion round exponentiation saturates");
    ParetoFrontier p;
    pareto_insert(p, Representation{{1, 2}, {}, 2});
    pareto_insert(p, Representation{{1}, {}, 3});
    require(p.entries.size() == 2, "pareto retains depth tradeoff");
    std::unordered_map<BitMask, int, BitMaskHash> nodes;
    nodes.emplace(BitMask(3), 7);
    require(lookup_node(nodes, BitMask(3)) == 7, "duplicate mask reuses node");
    require(lookup_node(nodes, BitMask(5)) == -1, "missing mask lookup");
    validate_target_width({BitMask(3)}, 2);
    bool width_rejected = false;
    try { validate_target_width({BitMask(4)}, 2); } catch (...) { width_rejected = true; }
    require(width_rejected, "target width validated");

    auto state = SearchState::initial({BitMask(15)}, 4, options(2));
    auto candidates = state.candidates();
    require(has_pair(candidates, 0, 1), "initial pair candidate");
    auto pair = std::find_if(candidates.begin(), candidates.end(), [](const Candidate &candidate) {
        return candidate.left == 0 && candidate.right == 1;
    });
    require(pair != candidates.end(), "pair lookup");
    auto parent_signature = state.signature();
    auto child = state.apply(*pair);
    require(child.valid(), "transactional pair application");
    require(!child.frontier(0).entries.empty(), "pair replacement representation");
    require(state.signature() == parent_signature, "apply does not mutate parent");
    auto wrong_depth = state.apply(Candidate{0, 1, 0, 1});
    require(wrong_depth.valid(), "candidate depth is recomputed");
    require(wrong_depth.signature().find("3@1;") != std::string::npos,
            "candidate depth is not trusted");
    auto depth_two = child.apply(Candidate{2, 4, 2, 1});
    require(depth_two.valid(), "second pair application");
    auto depth_two_signature = depth_two.signature();
    require(!depth_two.apply(Candidate{3, 5, 3, 1}).valid(), "over-depth transition rejected");
    require(depth_two.signature() == depth_two_signature, "rejected transition preserves parent");
    auto absent_pair_state = SearchState::initial({BitMask(3)}, 4, options(2));
    require(!absent_pair_state.apply(Candidate{2, 3, 1, 1}).valid(),
            "absent pair transition rejected");

    auto shared = SearchState::initial({BitMask(7), BitMask(11)}, 4, options(2));
    auto shared_candidates = shared.candidates();
    auto shared_pair = std::find_if(shared_candidates.begin(), shared_candidates.end(),
                                    [](const Candidate &candidate) {
                                        return candidate.left == 0 && candidate.right == 1;
                                    });
    require(shared_pair != shared_candidates.end(), "shared pair indexed");
    auto shared_child = shared.apply(*shared_pair);
    require(shared_child.valid(), "shared pair remains feasible");
    require(shared_child.pair_support(*shared_pair) == 0, "pair support is updated");

    auto toggle_state = SearchState::initial({BitMask(15)}, 4, options(2));
    auto toggled = toggle_state.apply_toggle(0, 1, 0, 0);
    require(toggled.valid(), "residual toggle transition is valid");
    require(toggled.frontier(0).entries.front().mask == BitMask(15),
            "residual toggle preserves target mask");
    require(toggle_state.signature() != toggled.signature(),
            "residual toggle does not mutate parent state");
    auto scored_toggle = toggle_state.apply_toggle(0, 1, 0, 0, 3);
    require(scored_toggle.valid(), "scored residual toggle is valid");
    require(scored_toggle.frontier(0).entries.front().future_support == 3,
            "residual future support survives frontier insertion");
    auto toggle_too_shallow = SearchState::initial({BitMask(15)}, 4, options(1));
    require(!toggle_too_shallow.apply_toggle(0, 1, 0, 0).valid(),
            "residual toggle respects hard depth limit");
    auto residual_options = options(2);
    residual_options.residual_pairs = 2;
    residual_options.residual_targets = 1;
    auto residual_state = SearchState::initial({BitMask(15), BitMask(7)}, 4,
                                               residual_options);
    auto residual_candidates = residual_state.toggle_candidates();
    require(residual_candidates.size() <= 2, "residual target expansion is bounded");
    for (const auto &candidate : residual_candidates)
        require(candidate.toggle, "residual candidate is marked as toggle");
    for (const auto &candidate : residual_candidates) {
        require(candidate.affected_targets > 0, "residual candidate reports affected targets");
        require(candidate.depth_slack >= 0, "residual candidate reports depth slack");
        require(candidate.toggle_penalty >= 0, "residual candidate reports toggle penalty");
    }
    auto bridge_options = options(2);
    bridge_options.residual_pairs = 8;
    bridge_options.residual_targets = 1;
    auto bridge_state = SearchState::initial({BitMask(7), BitMask(11)}, 4,
                                             bridge_options);
    auto bridges = bridge_state.bridge_candidates();
    require(!bridges.empty(), "cross-representation bridge candidate generated");
    require(bridges.front().bridge && !bridges.front().toggle,
            "bridge candidate has distinct transition kind");
    auto bridged = bridge_state.apply_bridge(
        bridges.front().left, bridges.front().right,
        static_cast<std::size_t>(bridges.front().target),
        static_cast<std::size_t>(bridges.front().representation),
        static_cast<std::size_t>(bridges.front().target2),
        static_cast<std::size_t>(bridges.front().representation2));
    require(bridged.valid(), "cross-representation bridge transition is valid");
    auto empty_toggle_state = toggle_state.apply(Candidate{0, 1, 1, 1});
    require(empty_toggle_state.valid(), "state exposes a residual representation");
    auto empty_default = empty_toggle_state.apply_toggle(0, 1, 0, 0);
    require(!empty_default.valid(), "empty toggle is disabled by default");
    auto empty_options = options(2);
    empty_options.allow_empty_toggle = true;
    auto empty_enabled = SearchState::initial({BitMask(15)}, 4, empty_options)
                             .apply(Candidate{0, 1, 1, 1});
    require(empty_enabled.valid(), "empty toggle setup remains valid");
    auto empty_transition = empty_enabled.apply_toggle(0, 1, 0, 0);
    require(empty_transition.valid(), "empty toggle option reaches absent operands");
    require(residual_state.fingerprint() == residual_state.fingerprint(),
            "structure fingerprint is deterministic");
    require(residual_state.fingerprint() != residual_state.apply(*shared_pair).fingerprint(),
            "different search structures receive different fingerprints");

    // The second output has alternative producers. The unused producer of
    // mask 7 must not become a mandatory cost of completing mask 15.
    auto alternatives = SearchState::initial({BitMask(3), BitMask(15)},
                                              4, options(4));
    int ancestor_bound = alternatives.best_cost_lower_bound();
    for (auto candidate : std::vector<Candidate>{{0, 1}, {2, 4}, {3, 4}, {2, 6}}) {
        alternatives = alternatives.apply(candidate);
        require(alternatives.valid(), "alternative producer transition remains valid");
        ancestor_bound = std::max(ancestor_bound, alternatives.best_cost_lower_bound());
    }
    require(alternatives.complete(), "alternative producer path completes");
    auto alternative_circuit = alternatives.export_circuit();
    alternative_circuit.verify({BitMask(3), BitMask(15)});
    require(alternative_circuit.cost() == 3, "unused alternative producer is pruned");
    require(ancestor_bound <= alternative_circuit.cost(),
            "every ancestor lower bound is admissible for its completion");

    auto balanced = synthesize_bounded({BitMask(15)}, 4, options(2));
    require(balanced.status == "verified", "D=2 feasible");
    require(balanced.circuit.cost() == 3, "balanced XOR count");
    require(depth_profile(balanced.circuit).depth <= 2, "hard depth bound");
    balanced.circuit.verify({BitMask(15)});
    auto repeat = synthesize_bounded({BitMask(15)}, 4, options(2));
    require(repeat.status == "verified", "repeat feasible");
    require(repeat.circuit.json().dump() == balanced.circuit.json().dump(),
            "same seed is deterministic");

    auto impossible = synthesize_bounded({BitMask(15)}, 4, options(1));
    require(impossible.status == "infeasible_by_lower_bound", "lower-bound rejection");

    auto state_limited_options = options(2);
    state_limited_options.beam_width = 1;
    state_limited_options.branch_width = 1;
    state_limited_options.restarts = 1;
    state_limited_options.timeout_seconds = 0;
    state_limited_options.max_states = 1;
    auto state_limited = synthesize_bounded({BitMask(7), BitMask(15)}, 4,
                                            state_limited_options);
    require(state_limited.explored_states == 1,
            "state budget stops before the next state (actual " +
                std::to_string(state_limited.explored_states) + ")");
    require(state_limited.status == "timeout", "state budget exhaustion is a timeout");
    require(state_limited.round_limit == 14, "result records the completion round limit");
    require(state_limited.max_rounds == 1, "result records the actual maximum restart rounds");

    auto incumbent = synthesize_bounded({BitMask(15)}, 4, state_limited_options);
    require(incumbent.status == "verified",
            "state budget exhaustion preserves a verified baseline incumbent");
    require(incumbent.explored_states <= 1, "incumbent search respects the state budget");
    require(incumbent.circuit.outputs.size() == 1,
            "state budget exhaustion retains the incumbent circuit");
    incumbent.circuit.verify({BitMask(15)});

    auto restart_limited_options = options(1);
    restart_limited_options.beam_width = 1;
    restart_limited_options.branch_width = 1;
    restart_limited_options.restarts = 2;
    restart_limited_options.timeout_seconds = 0;
    restart_limited_options.max_states = 1;
    auto restart_limited = synthesize_bounded({BitMask(3)}, 2, restart_limited_options);
    require(restart_limited.explored_states == 1,
            "state budget remains cumulative across restarts (actual " +
                std::to_string(restart_limited.explored_states) + ")");
    require(restart_limited.status == "verified",
            "cross-restart state exhaustion preserves a verified incumbent");
    restart_limited.circuit.verify({BitMask(3)});

    auto sufficient_options = state_limited_options;
    sufficient_options.max_states = 32;
    auto sufficient = synthesize_bounded({BitMask(7), BitMask(15)}, 4, sufficient_options);
    require(sufficient.status == "verified", "sufficient state budget can complete");
    require(sufficient.explored_states > 1 && sufficient.explored_states <= 32,
            "completed search stays within its state budget");
    require(sufficient.max_rounds <= sufficient.round_limit,
            "completed search stays within its round budget");
    sufficient.circuit.verify({BitMask(7), BitMask(15)});

    auto final_round_options = options(0);
    final_round_options.beam_width = 1;
    final_round_options.branch_width = 1;
    final_round_options.restarts = 1;
    final_round_options.timeout_seconds = 0;
    final_round_options.max_states = 1;
    auto final_round = synthesize_bounded({BitMask(1)}, 1, final_round_options);
    require(final_round.status == "verified", "complete state is exported in the final round");
    require(final_round.round_limit == 1 && final_round.max_rounds == 1 &&
                final_round.explored_states == 1,
            "final-round completion updates round statistics");

    auto exhausted_options = options(2);
    exhausted_options.beam_width = 1;
    exhausted_options.branch_width = 1;
    auto exhausted = synthesize_bounded({BitMask(7), BitMask(11)}, 4,
                                        exhausted_options);
    require(exhausted.status == "search_failed" || exhausted.status == "verified",
            "beam exhaustion is not infeasibility");
    exercise_mwsq_support_summary();
    exercise_mwsq_greedy();
    exercise_completable_greedy();
    exercise_candidate_pool();
    exercise_preexpired_incumbent();
    exercise_dense100_5s();
    std::cout << "xd search primitives passed\n";
}
