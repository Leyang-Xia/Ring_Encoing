#include "xd/search.hpp"
#include <algorithm>
#include <iostream>
using namespace xd;

namespace {
BoundedScheduleOptions options(int depth) {
    BoundedScheduleOptions o;
    o.max_depth = depth;
    o.beam_width = 4;
    o.branch_width = 8;
    o.max_representations = 8;
    o.restarts = 1;
    return o;
}
bool has_pair(const std::vector<Candidate> &candidates, int left, int right) {
    for (const auto &candidate : candidates)
        if (candidate.left == left && candidate.right == right)
            return true;
    return false;
}
}
int main() {
    require(combine_depth({0, 0, 0, 0}) == 2, "balanced four-leaf depth");
    require(!dominates(std::vector<int>{3, 2}, std::vector<int>{4, 1}), "term ordering");
    require(dominates(std::vector<int>{3, 2}, std::vector<int>{4, 3}), "term-and-depth dominance");
    Representation a{{3, 2}, {}, 2}, b{{4, 3}, {}, 3};
    require(dominates(a, b), "term-and-depth dominance");
    require(normalize_terms({4, 2, 4}) == std::vector<int>({2}), "duplicate terms cancel");
    require(candidate_depth(2, 2) == 3, "candidate depth");
    require(!candidate_allowed(2, 2, 2), "hard candidate depth filter");
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

    Circuit incumbent;
    incumbent.inputs = 4;
    incumbent.has_zero = false;
    incumbent.ops = {{4, 0, 1}, {5, 2, 3}, {6, 4, 5}};
    incumbent.outputs = {6};
    incumbent.verify({BitMask(15)});
    auto incumbent_options = options(2);
    incumbent_options.lookahead_width = 2;
    incumbent_options.lookahead_depth = 1;
    auto referenced = synthesize_bounded({BitMask(15)}, 4, incumbent_options,
                                         {incumbent});
    require(referenced.status == "verified", "verified incumbent keeps bounded search feasible");
    require(referenced.has_reference && referenced.reference_cost == incumbent.cost(),
            "bounded result retains the verified reference circuit");
    require(referenced.circuit.cost() <= incumbent.cost(),
            "bounded search never regresses below verified incumbent");
    require(depth_profile(referenced.circuit).depth <= incumbent_options.max_depth,
            "lookahead result respects hard depth bound");

    Circuit too_deep;
    too_deep.inputs = 4;
    too_deep.has_zero = false;
    too_deep.ops = {{4, 0, 1}, {5, 4, 2}, {6, 5, 3}};
    too_deep.outputs = {6};
    too_deep.verify({BitMask(15)});
    auto filtered = synthesize_bounded({BitMask(15)}, 4, options(2), {too_deep});
    require(filtered.reference_cost < 0, "incumbent deeper than D is ignored");
    require(filtered.status == "verified" && depth_profile(filtered.circuit).depth <= 2,
            "depth filter does not weaken bounded feasibility");

    auto impossible = synthesize_bounded({BitMask(15)}, 4, options(1));
    require(impossible.status == "infeasible_by_lower_bound", "lower-bound rejection");

    auto exhausted_options = options(2);
    exhausted_options.beam_width = 1;
    exhausted_options.branch_width = 1;
    auto exhausted = synthesize_bounded({BitMask(7), BitMask(11)}, 4,
                                        exhausted_options);
    require(exhausted.status == "search_failed" || exhausted.status == "verified",
            "beam exhaustion is not infeasibility");
    std::cout << "xd search primitives passed\n";
}
