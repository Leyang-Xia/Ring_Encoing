# Repository-aligned XOR benchmark results

Upstream: https://github.com/Leyang-Xia/Ring_Encoing
Pinned source: 4b6a4e83b00c2cd66cfd9798e41959ca0361a25a (remote main when checked).

Open report.html for the comparison. comparison.csv contains both depth limits; all_candidates.csv contains every old/new candidate, preserving labels, provenance and the Pareto result.

## What was actually reproduced

- Four official xordepth matrices and the exact documented `frontier --share 2,3,0` commands, including the original `--also` certificate files for Vandermonde and systematic RS.
- All 60 original candidate labels/costs/depths match archived measurements. Their original target masks match too.
- The 206-circuit audit reproduces 178 depth improvements and 59 removed gates.
- Two documented depth-6 Beam searches were actually run with beam=64, branch=16, representations=8, restarts=8, seed=1, timeout=180s. Cost/depth reproduce 221/6 and 149/6; explored-state counts differ because the native search is wall-clock bounded and the historical document also identifies an older source revision.
- Root README's systematic-rs-7-4/direct UBER_XSET profile with params=(0,3,0), timeout=60s completed at cost 126 and passed ringlab verify.
- Original root Make and xordepth Make tests, modified versions of both suites, new core/CLI tests and ASan/UBSan core tests passed. LeakSanitizer unsupported under ptrace, so its check is explicitly excluded (detect_leaks=0). CMake source lists/tests are updated, but cmake/ctest are not installed in this environment, so that build route was not run.

This scope is the official four-matrix count/depth benchmark and two bounded cases, not a rerun of all quotient-ring/family research searches.

## Native integration

`integration.patch` adds OUR_DCCSE to shared C++ core, exposes it in ringlab analyze, and optionally includes its candidates in the existing xordepth frontier. The existing candidate serialization, Pareto calculation and native verify command are reused. No Python experiment driver is required. Historical certificates, matrices and third-party code were not changed.

Main comparison: one restart, seed0, hard depth=D0 or D0+1, lookahead=false, local refinement=true. This aligns the fresh CSE frontier's single-run convention. Two separately labelled supplementary ringlab runs use depth6/restarts8/seed1. Equal restart labels are not equal computational budgets; archived best candidates also had different historical budgets.

All 68 integrated frontier candidates pass xordepth verify. Without the opt-in flag, all four result files remain BYTE IDENTICAL to outputs from the unmodified executable. Original entries in the expanded frontier preserve every field except their Pareto marker, which legitimately changes when a new candidate dominates them.

## Results at the minimum depth

| Matrix (outputs x inputs) | Depth | Original frontier best | OUR_DCCSE | Gate reduction |
|---|---:|---:|---:|---:|
| Vandermonde (24 x 64) | 5 | 336 | 224 | 33.33% |
| Systematic RS (24 x 32) | 5 | 168 | 152 | 9.52% |
| Cauchy (24 x 32) | 5 | 224 | 151 | 32.59% |
| GF(2^8) expanded (80 x 80) | 6 | 1830 | 1246 | 31.91% |

These update finite verified candidate frontiers, NOT global optima or all-algorithm records. At larger allowed depths, old candidates may still win: Vandermonde at depth6 is 203 gates versus ours222; GF8 at depth7 is 1219 versus ours1226. Full tradeoffs are retained in report.html and comparison.csv.

## Reproduce from the real upstream checkout

```bash
git clone https://github.com/Leyang-Xia/Ring_Encoing.git
cd Ring_Encoing
git checkout 4b6a4e83b00c2cd66cfd9798e41959ca0361a25a
git apply /absolute/path/to/integration.patch
bash docs/run-aligned-benchmark.sh /tmp/aligned-results
```

The patch passed `git apply --check` against the original source snapshot. The script runs Make tests, native frontiers and verifies, historical audit, the actual 180s Beam searches, the 8-restart OUR_DCCSE supplements, and the 60s README profile. RUN_BOUNDED=0 / RUN_PROFILE=0 optionally skip the longer stages for a quick rerun; neither stage was skipped in this delivered experiment.

Source implementation and limits are documented in METHOD.md and, after patching, docs/aligned-depth-cse.md. No remote repository write was performed.
