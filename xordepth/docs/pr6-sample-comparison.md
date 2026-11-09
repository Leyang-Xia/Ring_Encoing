# PR #6、xordepth 与 MW_SQ 的 sample 对比

## 口径

2026-11-04，使用 `data/matrices/sample.manifest` 的全部 12 个矩阵。每行格式为 `XOR 门数/实际深度`；`D` 取行权重深度下界。纯 Beam、PR V1/V2 和最终 xordepth 使用 seed=1、无墙钟截止；纯 Beam 与最终版本均为 `beam=4, branch=4, representations=4, restarts=2, max_states=256`，最终版本另加 `greedy_restarts=2, rewrite_trials=400`。PR V1 为 2 次普通 greedy，V2 为普通+lookahead 共 2 次。

MW_SQ 使用 `params=0 0 0`，后端墙钟上限 180 秒，再用 xordepth `rebalanced()` 重结合。MW_SQ 本身不使用硬深度 D，因此另列原始和重结合后的 `(cost/depth)`；重结合深度大于 D 时不能当作深度约束下的可比线路。强参数 `0 3 0` 在 dense 40×64 上 180 秒超时，故完整 sample 使用已完成的最强统一配置 `0 0 0`，超时不作为成本结论。

## 结果

| 矩阵 | D | 纯 Beam | PR V1 | PR V2 | 最终 xordepth | MW_SQ 原始 → 重结合 |
|---|---:|---:|---:|---:|---:|---:|
| dense 10×10 | 3 | 22/3 | 21/3 | 21/3 | **21/3** | 20/4 → 20/4 |
| dense 40×64 | 6 | timeout | 524/6 | 515/6 | **515/6** | 501/14 → 501/7 |
| sparse 12×20 | 2 | 24/2 | 24/2 | 24/2 | **24/2** | 24/2 → 24/2 |
| sparse 100×100 | 2 | 194/2 | 194/2 | 194/2 | **194/2** | 194/2 → 194/2 |
| blocks 24×40 | 3 | timeout | 56/3 | 56/3 | **56/3** | 55/5 → 55/4 |
| blocks 100×100 | 5 | timeout | 548/5 | 551/5 | **551/5** | 537/8 → 537/6 |
| modular 10×10 | 3 | 24/3 | 24/3 | 24/3 | **24/3** | 23/4 → 23/4 |
| modular 64×96 | 6 | timeout | 335/6 | 335/6 | **335/6** | 312/10 → 312/8 |
| Cauchy 7→4 | 5 | timeout | 147/5 | 150/5 | **150/5** | 147/9 → 147/7 |
| systematic RS 7→4 | 5 | timeout | 150/5 | 149/5 | **149/5** | 148/9 → 148/6 |
| Vandermonde 8→17 | 5 | timeout | 226/5 | 223/5 | **223/5** | 220/12 → 220/6 |
| GF(2⁸) 10×10 | 6 | timeout | 1247/6 | 1235/6 | **1235/6** | 1192/16 → 1192/7 |

汇总结论：

- 纯 Beam 只有 **4/12 verified**，PR V1、PR V2、最终 xordepth 和 MW_SQ 都是 **12/12 verified**。
- 最终 xordepth 在 **12/12** 上与 PR V2 的 `(cost, depth)` 完全一致；V1 在 blocks 100×100 和 Cauchy 上优于 V2，因此保留普通/lookahead portfolio 是合理的。
- MW_SQ+重结合的门数在 10 个样本上更少、2 个上持平，但只有 sparse 两例满足本轮 D；其余 10 例深度多 1–2 层。它提供了更强的门数参照，也显示最终方法换来了更浅线路，不能据此声称门数最优。
- Vandermonde D=5 的关键回归从纯 Beam 无结果变为 **223/5**；MW_SQ 为 **220/6**，两者是不同深度点上的有效上界。

## 耗时与限制

12 个样本的外部 `real` 时间求和：纯 Beam 12.61 秒、PR V1 2.89 秒、PR V2 3.16 秒、最终 xordepth 24.84 秒、MW_SQ 37.38 秒，重结合合计约 0.05 秒。最终 xordepth 的成本/深度与 PR V2 一致，但运行更慢，主要因为多消费者重写和 Beam 仍继续执行；这不是运行时性能追平的结论。GF(2⁸) 上 MW_SQ 占 32.56 秒，最终 xordepth 为 6.44 秒，PR V2 为 2.10 秒。

所有最终 xordepth、PR 和 MW_SQ 重结合线路均通过各自的符号功能验证；最终有 12 条 bounded 证据、12 条 MW_SQ 原始证据和 12 个包含重结合线路的 frontier 文件。逐矩阵参数、耗时、状态和完整电路见：

`data/results/pr6_sample_baseline_20261104/comparison.json`

这些数字是固定矩阵、固定 seed 和固定预算下的启发式上界，不是全局最优性证明。单一样本集不能证明其他尺寸或其他 seed 上仍与 PR 一致。
