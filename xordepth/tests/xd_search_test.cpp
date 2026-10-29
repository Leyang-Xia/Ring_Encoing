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
    ParetoFrontier p;
    pareto_insert(p, Representation{{1, 2}, {}, 2});
    pareto_insert(p, Representation{{1}, {}, 3});
    require(p.entries.size() == 2, "pareto retains depth tradeoff");
    BoundedScheduleOptions bad; bad.max_depth = -1;
    bool rejected = false; try { synthesize_bounded({ringlab::Mask(1)}, 1, bad); } catch (...) { rejected = true; }
    require(rejected, "negative options rejected");
    std::cout << "xd search primitives passed\n";
}
