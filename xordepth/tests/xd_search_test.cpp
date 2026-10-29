#include "xd/search.hpp"
#include <algorithm>
#include <iostream>
using namespace xd;
using ringlab::require;
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
    std::unordered_map<ringlab::Mask, int, ringlab::MaskHash> nodes;
    nodes.emplace(ringlab::Mask(3), 7);
    require(lookup_node(nodes, ringlab::Mask(3)) == 7, "duplicate mask reuses node");
    require(lookup_node(nodes, ringlab::Mask(5)) == -1, "missing mask lookup");
    validate_target_width({ringlab::Mask(3)}, 2);
    bool width_rejected = false;
    try { validate_target_width({ringlab::Mask(4)}, 2); } catch (...) { width_rejected = true; }
    require(width_rejected, "target width validated");

    auto state = SearchState::initial({ringlab::Mask(15)}, 4, options(2));
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
    auto depth_two = child.apply(Candidate{4, 2, 2, 1});
    require(depth_two.valid(), "second pair application");
    auto depth_two_signature = depth_two.signature();
    require(!depth_two.apply(Candidate{5, 3, 3, 1}).valid(), "over-depth transition rejected");
    require(depth_two.signature() == depth_two_signature, "rejected transition preserves parent");

    auto shared = SearchState::initial({ringlab::Mask(7), ringlab::Mask(11)}, 4, options(2));
    auto shared_candidates = shared.candidates();
    auto shared_pair = std::find_if(shared_candidates.begin(), shared_candidates.end(),
                                    [](const Candidate &candidate) {
                                        return candidate.left == 0 && candidate.right == 1;
                                    });
    require(shared_pair != shared_candidates.end(), "shared pair indexed");
    auto shared_child = shared.apply(*shared_pair);
    require(shared_child.valid(), "shared pair remains feasible");
    require(shared_child.pair_support(*shared_pair) == 0, "pair support is updated");

    auto balanced = synthesize_bounded({ringlab::Mask(15)}, 4, options(2));
    require(balanced.status == "verified", "D=2 feasible");
    require(balanced.circuit.cost() == 3, "balanced XOR count");
    require(depth_profile(balanced.circuit).depth <= 2, "hard depth bound");
    balanced.circuit.verify({ringlab::Mask(15)});
    auto repeat = synthesize_bounded({ringlab::Mask(15)}, 4, options(2));
    require(repeat.status == "verified", "repeat feasible");
    require(repeat.circuit.json().dump() == balanced.circuit.json().dump(),
            "same seed is deterministic");

    auto impossible = synthesize_bounded({ringlab::Mask(15)}, 4, options(1));
    require(impossible.status == "infeasible_by_lower_bound", "lower-bound rejection");

    auto exhausted_options = options(2);
    exhausted_options.beam_width = 1;
    exhausted_options.branch_width = 1;
    auto exhausted = synthesize_bounded({ringlab::Mask(7), ringlab::Mask(11)}, 4,
                                        exhausted_options);
    require(exhausted.status == "search_failed" || exhausted.status == "verified",
            "beam exhaustion is not infeasibility");
    std::cout << "xd search primitives passed\n";
}
