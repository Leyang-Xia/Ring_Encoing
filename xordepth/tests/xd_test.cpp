#include "xd/experiments.hpp"
#include <iostream>

using namespace xd;


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
        auto chained = synthesize({BitMask(15)}, 4);
        require(chained.cost() == 3 && depth_profile(chained).depth == 3, "chained depth");
        require(depth_lower_bound({BitMask(15), BitMask(3)}) == 2, "lower bound");
        require(depth_lower_bound({BitMask(1)}) == 0, "single input bound");
        require(depth_lower_bound({BitMask(0)}) == 0, "zero row bound");
        ScheduleOptions options;
        auto balanced = synthesize_schedule({BitMask(15)}, 4, options);
        require(balanced.cost() == 3 && depth_profile(balanced).depth == 2, "balanced depth");
        balanced.verify({BitMask(15)});

        // Without sharing, every row is its own balanced tree.
        options.share_threshold = 0;
        BitRows disjoint{BitMask(15), BitMask(0xf0)};
        auto lone = synthesize_schedule(disjoint, 8, options);
        require(lone.cost() == (4 - 1) + (4 - 1), "no-sharing cost");
        lone.verify(disjoint);

        // Re-association keeps the gate count and shortens the chain.
        auto repaired = rebalanced(chained);
        require(repaired.cost() == chained.cost(), "rebalance kept the gate count");
        require(depth_profile(repaired).depth == 2, "rebalance shortened a chain");
        repaired.verify({BitMask(15)});
        require(rebalanced(repaired).cost() == repaired.cost(), "rebalance is stable");

        BitRows rows{BitMask(15), BitMask(9), BitMask(6), BitMask(0)};
        auto shared = synthesize(rows, 4);
        auto fixed = rebalanced(shared);
        fixed.verify(rows);
        require(fixed.cost() <= shared.cost(), "rebalance cost is non-increasing");
        require(depth_profile(fixed).depth <= depth_profile(shared).depth,
                "rebalance depth is non-increasing");

        // A shared gate can be removed only by rewriting all of its consumers
        // as one transaction while preserving the hard depth limit.
        Circuit multi;
        multi.inputs = 4;
        multi.zero_node = 4;
        multi.ops = {{{5, 0, 1}}, {{6, 0, 2}}, {{7, 0, 3}}, {{8, 5, 2}}, {{9, 5, 3}}};
        multi.outputs = {6, 7, 8, 9};
        BitRows multi_rows{BitMask(5), BitMask(9), BitMask(7), BitMask(11)};
        multi.verify(multi_rows);
        auto rewritten = rewrite_depth_gates(multi, multi_rows, 2, 40);
        rewritten.verify(multi_rows);
        require(rewritten.cost() == 4, "shared-gate rewrite removes one gate");
        require(depth_profile(rewritten).depth <= 2, "shared-gate rewrite keeps depth");
        rejects([&] { rewrite_depth_gates(multi, multi_rows, 1, 40); });

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
