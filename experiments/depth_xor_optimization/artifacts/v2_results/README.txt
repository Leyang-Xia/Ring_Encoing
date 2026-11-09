V2 depth-constrained XOR gate minimization, 2026-11-03.

Apply full_integration.patch to a clean Leyang-Xia/Ring_Encoing snapshot at4b6a4e83b00c2cd66cfd9798e41959ca0361a25a. It includes earlier native OUR_DCCSE integration plus this V2, and passed git apply --check against the original snapshot. Do not apply the full patch over an already patched checkout; code/ contains the final edited core, CLI, tests and V2 documentation for review.
Build: make -j3; make test; make -C xordepth test.
C++17 only for all synthesis and optimization; report.py only audits existing circuits and produces comparison tables, and is not needed by the optimizer.

New usage:
./build/ringlab analyze MATRIX.txt --method OUR_DCCSE_V2 --depth D --seeds 8 --seed 1 --output result.json
Add --incumbent FRONTIER.json for optional verified historical candidates.
./build/ringlab verify --input result.json --samples 0

Core algorithm: combine depth-feasible greedy and bounded-lookahead CSE, then transactionally remove shared intermediate gates by rewriting all consumers together. Reverse-propagated output deadlines allow consuming noncritical timing slack. Strict gate-count improvement and hard global depth checks are required. At most400 deletion trials across4 passes. This is heuristic, not SAT or an optimality proof.
Default V2 --seeds8 means8 per strategy,16 total. Old original comparison used8 total. An additional v1r16 control uses16 greedy restarts to separate portfolio gains from a higher restart count. All have an external180s ceiling; no attempt is made to force early-finished algorithms to consume the full budget. Wall times include concurrent execution, not isolated speedup measurements.

Four exact native matrices and old candidate/beam evidence are bundled under four_matrices. The unchanged bounded algorithm is not rerun: its previous same-matrix/D, beam64,branch16,representations8,restarts8,seed1,timeout180 measurements are reused. Its timeout records remain inconclusive, not infeasibility proofs. Original historical candidate optimization budgets differ.

The standalone runs never load historical circuits. Warm runs select the minimum-gate feasible circuit from standalone and native frontier candidates, then run the same local optimization. Report these separately: selecting203/6 for Vandermonde is reuse, whereas RS140/6->139/6 is an actual new rewrite improvement. Original input/output transform, conversion costs and live gates are verified, including historical ring candidates.

Logs: root-tests.log and xordepth-tests.log passed; ASAN+UBSAN core test passed with ASAN_OPTIONS=detect_leaks=0 (no leak-sanitizer claim). CMake configuration already includes the edited source file, no new source-list changes needed; CMake/CTest not installed here. Patch whitespace and application checks passed.

run.sh and run_equal_restarts.sh reproduce the new CLI experiments. Set REPO_ROOT to the patched checkout directory; keep four_matrices and v2_results as sibling directories. All input matrices are unchanged upstream files. report.py independently validates GF2 outputs, dependencies, all live gates and depth ceilings; it also checks native verification records and monotonicity vs initial candidates.

Final findings:
Independent V2 vs old8: lower gate count in all8 tested matrix/depth pairs.
Independent V2 vs old16 (equal total restarts):6 wins,2 losses. RS old16 gives145/5 for both D5 andD6, versus V2 147/5 and146/6. Therefore portfolio allocation is not universally better than spending all16 seeds on greedy. Users should retain best verified circuits from both when available; pass the preferred circuit as --incumbent to optimize it.
Warm mode in this table imports only the historical native frontier, not the newly computed old16 control. Thus RS D5 warm147/5 is not the best circuit among every newly measured candidate (old16 gives145/5). Preserve this distinction when reporting global best-known results.
Warm RS D6 improves a loaded140/6 circuit to139/6. GF8 standalone D7 yields1210/7,9 fewer than historical1219/7. Vandermonde warm203/6 is preserved historical evidence, not a newly discovered improvement.
Vandermonde D6 independent run produced217/5, which also satisfies D5; the direct D5 run gave218/5. Heuristic search outcomes need not be monotone with D. A portfolio can reuse any verified circuit whose actual depth fits the requested bound. These are single-seed benchmark instances, not statistical/global-optimality claims.
All24 new final/control circuits passed native and independent exact verification. The16 initial_elimination_only runs are exploratory checkpoints predating the final portfolio and deadline-based rewrite; do not mix them into final comparison.csv.
