Four exact xordepth fixed matrices, no synthetic substitutions.
Vandermonde: vandermonde_8to17_direct.txt (24 outputs x64 inputs), D0=5.
RS: systematic_rs_7_4_direct.txt (24x32), D0=5.
Cauchy: cauchy_7_4_direct.txt (24x32), D0=5.
GF8: gf8_10x10_expanded.txt (80x80), D0=6.
D0=max ceil(log2 row weight). Each matrix tested at D0 and D0+1.

Native beam configuration matches xordepth/docs/depth-bounded-search.md:
--depth D --beam 64 --branch 16 --representations 8 --restarts 8 --seed 1 --timeout 180
External watchdog190s. OUR_DCCSE: --depth D --seeds8 --seed1, default refinement=true, lookahead=false, external watchdog180s.
This equalizes target, depth ceiling, configured restarts, seed and nominal maximum180s budget; no beam width analogue exists for greedy CSE, and actual time can differ. Neither method is run to optimality. Source main4b6a4e83b00c2cd66cfd9798e41959ca0361a25a; same local native integration patch as matched_depth/integration.patch. No source algorithm changes in this experiment.

Vandermonde and RS bounded/OURS results reused without alteration from the immediately preceding experiment. Cauchy and GF8 are new measurements. See run.sh for new commands; times include concurrent execution and should not be treated as isolated speedup ratios.

*_native_frontier.json contains the original repository candidates reproduced earlier with frontier --share2,3,0, importing structured_ring_example.json for Vandermonde and practical_mds_examples.json for RS (no imported evidence for Cauchy/GF8). Comparing best stored/native candidate is separate from comparing a freshly run bounded algorithm; archival optimization budgets differ and are not an equal-time experiment.

Every successful circuit checked natively and by independent GF2/dependency/live-gate/depth audit in report.py. No new synthesis implementation was added. Timeout means inconclusive, never proof of infeasibility. Circuit depths assume XOR2 unit delay, simultaneous inputs, free wires and fanout.

Final results in combined_comparison.csv. All12 successful newly compared circuits and60 original frontier circuits passed independent audit. Native verification also passed for successful circuits. No additional OOM event occurred in this batch. GF8 D6 reached external190s watchdog without exporting a circuit; GF8 D7 returned1219/7 after189.301s. Cauchy D5 returned native timeout; D6 returned150/6. OURS Cauchy147/5 at both ceilings; GF8 1247/6 and1224/7.
The no-share CSE baseline still permits equivalent-signal reuse; it is not strictly independent per-output trees.
The archive retains preceding matched_depth experiments (including generated128x96); four_matrices is the authoritative four-official-matrix comparison requested in this extension.
