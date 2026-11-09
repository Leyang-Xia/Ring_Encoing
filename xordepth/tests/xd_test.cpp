#include "xd/experiments.hpp"
#include "xd/xsets.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>

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

        // X-Sets schedule import: the "Schedule:" trailer of a Plank X-Sets
        // run parses into a Circuit, verifies against the matrix, and feeds
        // rebalanced() directly.
        {
            const char *schedule_path = "/tmp/xd_xsets_test_schedule.txt";
            const char *matrix_path = "/tmp/xd_xsets_test_matrix.txt";
            {
                std::ofstream out(schedule_path);
                out << "Targets Left 300.  New Element 0111 Edge Weight: 2\n";
                out << "1011\n1101\n";
                out << "Schedule:\n";
                out << "s0 = x0\n";
                out << "s1 = x1\n";
                out << "s2 = x2\n";
                out << "s3 = x3\n";
                out << "s4 = s0 + s1\n";
                out << "s5 = s4 + s2\n";
                out << "s6 = s5 + s3\n";
                out << "b0 = s6\n";
                out << "Total XOR operations: 3\n";
            }
            {
                std::ofstream out(matrix_path);
                out << "1111\n";
            }
            auto parsed = parse_xsets_schedule(schedule_path, 4);
            require(parsed.circuit.inputs == 4, "xsets input count");
            require(parsed.circuit.cost() == 3, "xsets gate count");
            require(parsed.declared_xor == 3, "xsets declared xor trailer");
            // b0 = s6 (node 6).
            require(parsed.circuit.outputs == std::vector<int>{6},
                    "xsets bindings in b order");
            BitRows xsets_rows{BitMask(0xF)};
            parsed.circuit.verify(xsets_rows);
            require(depth_profile(parsed.circuit).depth == 3,
                    "xsets chain keeps its original depth");
            auto repaired = rebalanced(parsed.circuit);
            repaired.verify(xsets_rows);
            require(depth_profile(repaired).depth == 2,
                    "rebalance shortens the imported chain");

            // The one-shot command verifies before and after and reports both.
            Json report = rebalance_xsets(schedule_path, matrix_path);
            require(report["method"].string() == "XSETS_REBALANCE",
                    "rebalance result method");
            require(report["original"]["cost"].int32() == 3 &&
                        report["original"]["depth"].int32() == 3,
                    "rebalance report of the original");
            require(report["rebalanced"]["cost"].int32() == 3 &&
                        report["rebalanced"]["depth"].int32() == 2,
                    "rebalance report after re-association");
            require(report["original"]["depth_histogram"]["3"].int32() == 1,
                    "original depth histogram");
            require(report["rebalanced"]["depth_histogram"]["2"].int32() == 1,
                    "rebalanced depth histogram");
            auto from_report = Circuit::from_json(report["circuit"]);
            from_report.verify(xsets_rows);

            // The parser also consumes a raw stream, e.g. a piped X-Sets run.
            {
                std::istringstream piped(
                    "Schedule:\ns0 = x0\ns1 = x1\ns2 = x2\ns3 = x3\n"
                    "s4 = s0 + s1\ns5 = s4 + s2\ns6 = s5 + s3\nb0 = s6\n"
                    "Total XOR operations: 3\nSchedule Depth: 3\n"
                    "# technique=UBER_XSET thresh=0 nstart=0 giveup=0 rows=1 cols=4 xor=3 depth=3\n");
                auto streamed = parse_xsets_schedule(piped, 4);
                require(streamed.circuit.cost() == 3 && streamed.declared_xor == 3,
                        "streamed xsets parse");
                require(streamed.declared_depth == 3, "streamed depth trailer");
                streamed.circuit.verify(xsets_rows);
                // A disagreeing depth trailer is caught by the command.
                {
                    std::istringstream bad_stream(
                        "Schedule:\nb0 = s1\nSchedule Depth: 9\n");
                    rejects([&] {
                        Json report = rebalance_xsets(bad_stream, "bad", matrix_path);
                        (void)report;
                    });
                }
            }

            // An output gate is never folded away: a chain of two outputs
            // keeps its depth even after re-association.
            {
                std::ofstream out(schedule_path);
                out << "Schedule:\n";
                out << "s0 = x0\ns1 = x1\ns2 = x2\ns3 = x3\n";
                out << "s4 = s0 + s1\ns5 = s4 + s2\ns6 = s5 + s3\n";
                out << "b0 = s6\nb1 = s5\n";
            }
            {
                std::ofstream out(matrix_path);
                out << "1111\n1110\n";
            }
            auto chained_outputs = parse_xsets_schedule(schedule_path, 4);
            BitRows chain_rows{BitMask(0xF), BitMask(0x7)};
            chained_outputs.circuit.verify(chain_rows);
            auto kept = rebalanced(chained_outputs.circuit);
            kept.verify(chain_rows);
            require(depth_profile(kept).depth == 3,
                    "rebalance keeps output-gate chains");

            // Malformed schedules are rejected, not silently accepted.
            auto bad = [&](const char *text) {
                std::ofstream out(schedule_path);
                out << text;
                rejects([&] { parse_xsets_schedule(schedule_path, 4); });
            };
            bad("s0 = x0\n"); // no Schedule: header
            bad("Schedule:\ns0 = s1 + s2\n"); // use before definition
            bad("Schedule:\ns0 = x9\n"); // input index out of range
            bad("Schedule:\ns0 = s0 + x0\n"); // self reference before definition
            bad("Schedule:\ns0 = x0\ns0 = x1\n"); // redefinition
            bad("Schedule:\nb0 = x0\nb0 = x1\n"); // duplicate binding
            bad("Schedule:\ns0 = x0\ns0 = x0 ~ x1\n"); // unrecognized line
            bad("Schedule:\ns0 = x0 + x1\nTotal XOR operations: 9\n"); // trailer mismatch
            bad("Schedule:\nb1 = x0\n"); // binding gap: b0 missing

            // Target count must match the matrix in the command wrapper.
            {
                std::ofstream out(matrix_path);
                out << "1111\n1110\n0001\n";
            }
            {
                std::ofstream out(schedule_path);
                out << "Schedule:\nb0 = x0\nb1 = x1\n";
            }
            rejects([&] { rebalance_xsets(schedule_path, matrix_path); });

            std::remove(schedule_path);
            std::remove(matrix_path);
        }

        std::cout << "xd: depth metrics, balanced scheduling and re-association passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
