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

        // MW_SQ support summary: (max pair support, ties at max) over rows of
        // node ids. Golden values lock the exact counting semantics (the fast
        // flat-key implementation must reproduce them bit for bit).
        require(mwsq_support_summary({}) == std::pair<std::size_t, std::size_t>{0, 0},
                "empty summary");
        require(mwsq_support_summary({{1}, {}}) ==
                    std::pair<std::size_t, std::size_t>{0, 0},
                "no pairable rows");
        // A row with a repeated id contributes that (id, id) pair once.
        require(mwsq_support_summary({{1}, {}, {7, 7}}) ==
                    std::pair<std::size_t, std::size_t>{1, 1},
                "repeated id counts as one pair");
        // Rows {1,2},{2,3},{1,2}: pair (1,2) twice, (2,3) once -> max 2, ties 1.
        require(mwsq_support_summary({{1, 2}, {2, 3}, {1, 2}}) ==
                    std::pair<std::size_t, std::size_t>{2, 1},
                "max support with one winner");
        // Rows {0,1},{0,1},{1,2},{1,2}: (0,1)=2 and (1,2)=2 -> max 2, ties 2.
        require(mwsq_support_summary({{0, 1}, {0, 1}, {1, 2}, {1, 2}}) ==
                    std::pair<std::size_t, std::size_t>{2, 2},
                "ties at max");
        // One wide row: every pair appears exactly once.
        require(mwsq_support_summary({{3, 1, 2}}) ==
                    std::pair<std::size_t, std::size_t>{1, 3},
                "unsorted row with all singleton pairs");
        // Mixed: {1,2},{1,3},{2,3},{1,2}: (1,2)=2, (1,3)=1, (2,3)=1 -> (2,1).
        require(mwsq_support_summary({{1, 2}, {1, 3}, {2, 3}, {1, 2}}) ==
                    std::pair<std::size_t, std::size_t>{2, 1},
                "triangle plus duplicate");

        std::cout << "xd: depth metrics, balanced scheduling and re-association passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
