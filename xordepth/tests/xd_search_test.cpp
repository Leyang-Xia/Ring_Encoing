#include "xd/search.hpp"
#include <iostream>
using namespace xd;
using ringlab::require;
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
    std::cout << "xd search primitives passed\n";
}
