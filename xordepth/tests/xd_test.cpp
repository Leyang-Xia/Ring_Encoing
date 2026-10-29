#include "xd/experiments.hpp"
#include <iostream>

using namespace xd;
using ringlab::require;

template <class F> void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, "invalid input accepted");
}

int main() {
    try {
        // A chain of three gates and a balanced tree over the same four terms.
        auto chained = ringlab::synthesize({Mask(15)}, 4);
        require(chained.cost() == 3 && depth_profile(chained).depth == 3, "chained depth");
        require(depth_lower_bound({Mask(15), Mask(3)}) == 2, "lower bound");
        require(depth_lower_bound({Mask(1)}) == 0, "single input bound");
        require(depth_lower_bound({Mask(0)}) == 0, "zero row bound");
        ScheduleOptions options;
        auto balanced = synthesize_schedule({Mask(15)}, 4, options);
        require(balanced.cost() == 3 && depth_profile(balanced).depth == 2, "balanced depth");
        balanced.verify({Mask(15)});

        // Without sharing, every row is its own balanced tree.
        options.share_threshold = 0;
        Rows disjoint{Mask(15), Mask(0xf0)};
        auto lone = synthesize_schedule(disjoint, 8, options);
        require(lone.cost() == (4 - 1) + (4 - 1), "no-sharing cost");
        lone.verify(disjoint);

        // Re-association keeps the gate count and shortens the chain.
        auto repaired = rebalanced(chained);
        require(repaired.cost() == chained.cost(), "rebalance kept the gate count");
        require(depth_profile(repaired).depth == 2, "rebalance shortened a chain");
        repaired.verify({Mask(15)});
        require(rebalanced(repaired).cost() == repaired.cost(), "rebalance is stable");

        Rows rows{Mask(15), Mask(9), Mask(6), Mask(0)};
        auto shared = ringlab::synthesize(rows, 4);
        auto fixed = rebalanced(shared);
        fixed.verify(rows);
        require(fixed.cost() <= shared.cost(), "rebalance cost is non-increasing");
        require(depth_profile(fixed).depth <= depth_profile(shared).depth,
                "rebalance depth is non-increasing");

        // Malformed circuits are rejected, not silently measured.
        rejects([] { depth_profile(Circuit()); });
        auto broken = chained;
        broken.ops[0][1] = 9999;
        rejects([&] { depth_profile(broken); });
        auto duplicated = chained;
        duplicated.ops.push_back(duplicated.ops[0]);
        rejects([&] { depth_profile(duplicated); });

        std::cout << "xd: depth metrics, balanced scheduling and re-association passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
