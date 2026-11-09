# Larger XOR matrices: native CLI extension benchmark

This is a synthetic SIZE-EXTENSION experiment using the same integrated C++ native CLI and verification paths. It is not a rerun of the four archived official matrices, and no old ring circuits were imported for these new transforms.

## Inputs

Nine matrices, outputs x inputs: 128x96, 192x128, 300x200, each with sparse/dense/blocks patterns. One fixed seed=123 per case: these are scale probes, not statistical performance estimates. sparse/dense use independent Bernoulli probabilities .15/.5. blocks: the largest multiple of four not exceeding half the columns forms disjoint four-input blocks, each selected with p=.5 per output; remaining columns are independent Bernoulli(.2). The original C++ SeedRandom generator and dynamic Mask type are used; generation is exposed through ringlab generate. Complete text matrices and exact mask manifests are included.

## Actual methods and parameters

- Native xordepth frontier --share 2,3,0: original CSE-chained, its rebalanced circuit, share2, share3, and no-share candidates. No OUR candidates are injected into this call, so a slow OUR solve cannot hide completed native baselines.
- Native ringlab analyze --method OUR_DCCSE, seeds=1, seed=0, refine=true, lookahead=false, depth=D0 or D0+1.
- Native MW_SQ and UBER_XSET: unchanged parameters (0,3,0), internal timeout=10s.
- Native bounded Beam: depth=D0+1, beam=8, branch=4, representations=4, restarts=1, seed=1, internal timeout=10s. This is explicitly a reduced-budget scale probe, not the previous official 64/16/8/8/180s configuration.
- Main command limit: 10 seconds. Native internally timed backends and Beam additionally have a 12-second external guard for initialization/cleanup. Their reported wall time can therefore exceed 10 seconds. The whole frontier task generates multiple candidates, whereas one OUR task generates one circuit: this is not an equal-CPU ranking.
- Four concurrent processes; XSETS_THREADS=1 for every task. The first calibration batch revealed that the backend otherwise defaults to hardware-concurrency threads; those oversubscribed timings/results were excluded from the primary report. The controlled batch was rerun in full.
- Supplement: ONLY previously unfinished OUR cases were rerun at 30 seconds (same seed, restart count and algorithm). These are labelled separately and do not overwrite the 10-second results.
- 300x200 dense full threshold scan did not complete in 10 seconds. A separate native frontier --share 0 task completed its chained, rebalanced and no-share baselines. share2/share3 remain unmeasured for this case under the recorded scan budget.

CSE-noshare disables explicit common-pair extraction but the repository's equal-mask node reuse can still share coincident intermediate functions. It is not the naive sum(row_weight-1) independent-tree count. Both the naive count and actual native baseline are preserved in metadata.

## Verification and results

63 controlled/supplement/fallback tasks: 28 completed verified tasks and 35 timeout tasks. Tasks can emit multiple circuits. All 62 completed circuit records passed the original native verify commands AND independent exact Python GF(2) audits of every gate, original matrix masks, all outputs, topology, live gate count and depth. Python is only used for auditing/reporting, never circuit synthesis. Incomplete tasks have no reported gate count or depth.

The main 10s table is in report.html, all_results.csv and summary.json. Thirty-second supplements are separate. Main qualitative findings:

- 128x96 dense: OUR 2301/6 versus native no-share 3285/6; OUR with one extra allowed layer gives 2214/7 versus CSE+rebalanced 2223/7.
- 192x128 blocks: OUR 1871/7 versus native no-share 2945/7; CSE+rebalanced is smaller but deeper at 1851/8.
- 300x200 blocks (30s supplement): OUR 4431/7 and 4308/8; native CSE+rebalanced 4285/9, share3 4745/8, no-share 6726/7.
- 300x200 sparse (30s supplement): OUR 4873/6 versus native CSE+rebalanced 4866/6; no gain and slower solver.
- 192x128 and 300x200 dense OUR attempts remained unfinished at 30s. No inference about feasibility or eventual gate count follows.
- Both X-Sets methods timed out in all nine main cases at the original (0,3,0) parameters. This does not establish that their eventual circuits are worse.

Timing is wall-clock command time in this environment, one matrix/seed per cell. The CSE frontier time covers all candidates in that command, not an isolated cost per algorithm. Do not treat these numbers as robust microbenchmarks or physical ASIC timing. The implementation currently rebuilds all depth-feasible pairs at each OUR step, unlike the native incremental CSE; optimizing this maintenance is a clear scaling priority.

## Reproduce

The included complete patch applies directly to upstream commit 4b6a4e83b00c2cd66cfd9798e41959ca0361a25a. It includes the earlier OUR_DCCSE native integration plus the new C++ input generator. Do NOT apply both the earlier patch and this full patch.

    git clone https://github.com/Leyang-Xia/Ring_Encoing.git
    cd Ring_Encoing
    git checkout 4b6a4e83b00c2cd66cfd9798e41959ca0361a25a
    git apply /path/to/integration_and_large_inputs.patch
    bash /path/to/run_large_benchmark.sh "$PWD" /tmp/large-xor-results

RUN_SUPPLEMENT=0 skips the 30s extensions. The helper run_task.sh must be next to run_large_benchmark.sh. Root make test (including new generation repeatability/dimension tests) passed. No third-party solver or historical certificate was changed. No remote repository write was performed.
