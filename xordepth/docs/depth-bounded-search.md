# 深度硬约束调度搜索

`DEPTH_BOUNDED_XSET_BEAM` 是一个参考 X-Sets/Uber-XSet 思路的启发式调度器。它在固定 GF(2) 输入矩阵上搜索二输入 XOR 网络，并把最大输出深度 `D` 作为硬约束：任何中间节点和最终输出都不能超过 `D`。

## 搜索状态

每个状态保存已物化节点、节点深度、每个输出的候选表示，以及表示中的 pair 支持索引。每个输出保留局部 Pareto 表，表示的两个坐标是：

```text
(remaining_terms, arrival_depth)
```

如果一个表示在这两个坐标上都不差且至少一项更好，较差表示会被删除。这样可以保留“项数更多但深度更低”的表示。

候选 pair 的优先级依次参考：支持数、替换后的剩余项数、表示完成深度、新节点深度和稳定节点顺序。候选生成、表示更新和最终电路导出都会检查深度。Beam 被剪空或搜索预算耗尽不会被解释成不可行证明。

## 状态含义

| 状态 | 含义 |
|---|---|
| `verified` | 找到完整电路，逐门功能、门数和深度都通过验证 |
| `infeasible_by_lower_bound` | `D` 小于行权重给出的安全深度下界 |
| `search_failed` | 搜索空间被启发式剪枝耗尽，没有找到线路 |
| `timeout` | 达到搜索时间预算，没有找到线路 |

只有 `verified` 结果包含 `cost`、`depth` 和完整 `circuit`。其余状态的 `cost` 为 null，不能用于 XOR 数比较。

## 使用

```bash
cd xordepth
make -j4
./build/xordepth bounded data/matrices/vandermonde_8to17_direct.txt \
  --depth 6 --beam 64 --branch 16 --representations 8 \
  --restarts 8 --lookahead-width 8 --lookahead-depth 2 \
  --seed 1 --timeout 180 --output /tmp/vandermonde-d6.json
./build/xordepth verify /tmp/vandermonde-d6.json
```

`--depth D` 是最大输出深度；`--beam`、`--branch`、`--representations`、`--restarts`、`--seed` 和 `--timeout` 分别控制 Beam 宽度、每状态候选数、每目标 Pareto 表上限、重启次数、候选重排种子和墙钟预算。
`--lookahead-width K --lookahead-depth L` 会对每个状态中最多 `K` 个候选做 `L` 步的选择性试走，并使用试走后的下界重新排序；默认两个值都是 `0`，以保持基线搜索行为。它只改变候选优先级，不放宽深度检查。

## Completable greedy 与有限终止

`--greedy-restarts N` 在 Beam 前运行严格减少表达式项数的深度可行 greedy，普通和有限 lookahead 尝试按相同 seed 成对执行；默认 `N=2`，`N=0` 可恢复纯 Beam 对照。每个输出维护一个 GF(2) 表达式，候选只有在替换后仍满足固定叶子集合的完成深度预算时才接受。新门至少服务两个输出，复用门至少服务一个输出，因此每轮表达式项数严格下降，主循环有限。最终使用 earliest-ready 合并，并逐门验证功能和实际深度。

`--no-timeout` 只关闭墙钟截止，不关闭结构性预算。`--max-states N` 限制跨 restart 累计的 Beam 取状态数；每个 restart 还受 `outputs * (2^(D+1)-1)` 轮次上界约束，算术使用饱和无符号运算。预算或轮次耗尽不是不可行证明。若已经得到 verified incumbent，后续状态预算、墙钟或重写截止不会丢弃它。

## 多线程（--threads N）

`--threads N`（默认 1）为 greedy portfolio 的 restart 与 Beam 的层内扩展分配工作线程。确定性通过三段式层结构保持：顺序发布 complete 状态并快照 incumbent 剪枝阈值 → 各父状态独立并行扩展 → 按原顺序回放转移表查重。因此任意线程数的结果与串行逐门一致（`--timeout` 运行除外，其结果本就依赖机器速度）。sample-12 无时限 4 线程 11.35s → 6.83s（1.66×），全量 57 矩阵 64.73s → 42.59s（1.52×）；三方位相容校验（改动前 vs 新串行 vs 新 4 线程，含完整电路 ops）零失配，ThreadSanitizer 无竞态报告。单元测试 `exercise_threaded_determinism` 锁定该性质。

注意：bounded 管线内部的 greedy 阶段最初漏传 `options.threads`（恒为 1 线程），修复后大矩阵加速明显提升。MW_SQ 并列重评分循环本身也已并行化：每个 tied 候选的评分只依赖 key 多重集，因此按行分块并行生成 key（各块由工作线程排序）后按固定块序归并、串行游程扫描，再按 tied 顺序归约出与串行一致的 winner 列表；`threads=1` 走原串行路径。为避免超额订阅，重启并行与评分内部并行共享线程预算（每个 attempt 分得 `threads/重启数` 个工作线程）。实测 6 个最大矩阵合计 47.95s（1 线程）→ 24.16s（4 线程，**1.98×**），全量 57 矩阵 64.73s → 38.63s（**1.68×**）。任意线程数输出与改动前串行逐门一致（6 大矩阵 t1–t10 相互一致；全量 57 对改动前基线 57/57 逐门相同），TSan 复验无警告（含 restarts×内部双层并行）。加速上限仍受 Amdahl 制约：Beam 的 Phase 1/3 与 CSE 基线为串行。

2026-11-04 的 12 矩阵对照中，纯 Beam 仅 4/12 verified；加入 completable greedy 和 400 次重写预算后为 12/12，并在全部样本上达到 PR #6 V2 相同的 `(cost, depth)`。完整数据及 MW_SQ+重结合参照见 [PR sample 对比](pr6-sample-comparison.md)。

## 固定矩阵结果

实验使用源码修订 `04a5114f8541e5472c3201ecbe876f6c16e697c3`，参数均为 `beam=64`、`branch=16`、`representations=8`、`restarts=8`、`seed=1`、`timeout=180s`。

| 矩阵 | `D` | 下界 | 状态 | XOR | 实际深度 | 搜索状态数 |
|---|---:|---:|---|---:|---:|---:|
| Vandermonde 8→17 direct | 6 | 5 | verified | 221 | 6 | 10,897 |
| systematic RS [7,4,4] direct | 6 | 5 | verified | 149 | 6 | 17,761 |

这两个结果是启发式上界，不是全局 XOR 最优性证明。完整目标掩码、搜索参数和电路保存在：

- `data/results/xd_bounded_vandermonde_d6.json`
- `data/results/xd_bounded_rs_d6.json`

可用 `xordepth verify` 重新执行功能、门数和深度检查。

## MW_SQ 参照

为了给深度搜索提供同一矩阵上的计数参照，又在 `64×24` 的
`vandermonde_8to17_direct.txt` 上运行了当前工作区的 MW_SQ。`0 3 0` 是更高的
前瞻参数；`0 4 0` 在 180 秒预算内超时。

| 方法 | 参数 | 状态 | XOR | 原始深度 | 重结合后深度 | 墙钟时间 |
|---|---|---|---:|---:|---:|---:|
| MW_SQ | `0 0 0` | verified | 220 | 12 | 6 | 0.126 s |
| MW_SQ | `0 3 0` | verified | 219 | 12 | 6 | 5.684 s |
| MW_SQ | `0 4 0` | timeout | — | — | — | 180 s |

重结合只重新安排单消费者 XOR 树，保留 MW_SQ 已经选择的共享节点，因此门数不变；
`0 3 0` 的 219 门线路在深度 6 预算下仍是一个有效参照。当前 bounded Beam
在同一 `D=6` 配置下得到 221 门，所以这组结果显示了搜索还没有达到该参照上界。
超时配置不提供 XOR 数结论。

完整 verified 电路和目标掩码保存在：

- `data/results/mw_sq_vandermonde_8to17_p0_0_0.json`
- `data/results/mw_sq_vandermonde_8to17_p0_3_0.json`
- `data/results/mw_sq_vandermonde_8to17_p0_4_0_timeout.json`

## 运行结果差异与 canonical 输入

MW_SQ 的候选评分包含平局选择，目标行的输入顺序会影响候选访问顺序。因此，下面两条
命令虽然计算的是同一个目标集合，可能得到不同的启发式上界：

```bash
./build/xsets 0 3 0 MW_SQ --quiet < data/matrices/vandermonde_8to17_direct.txt
# 220 XOR

./build/xsets 0 3 0 MW_SQ --quiet < data/matrices/vandermonde_8to17_direct_canonical.txt
# 219 XOR
```

`ringlab analyze` 的 backend 适配层会用 `std::set<Mask>` 对目标行去重并排序，因而使用
canonical 顺序。为保证直接 X-Sets 调用与 `ringlab`、xordepth 的参照一致，仓库保留了
排序后的 [`vandermonde_8to17_direct_canonical.txt`](../data/matrices/vandermonde_8to17_direct_canonical.txt)。
在 canonical 输入上，直接 X-Sets 和 `ringlab analyze` 都得到 219 XOR、深度 12；这是
行顺序导致的启发式路径差异，不是电路验证或 XOR 计数差异。

## 正确性边界

- 深度按二输入 XOR 门的拓扑深度计算，输入和常量零深度为 0；
- 补零、接线和输出别名不计 XOR；
- Beam 宽度、表示上限和超时都会改变搜索覆盖范围；
- `search_failed` 和 `timeout` 不说明不存在深度不超过 `D` 的线路；
- 当前版本不复制同一状态中的重复 mask 节点，也不建模 LUT、布线、扇出时延或多输入 XOR 门。

## 安全下界修复与多矩阵比较

早期版本把一个状态中所有 Pareto 表示的依赖节点取并集，并将其作为已支付门数；但这些表示是互相替代的候选，未被最终选择的依赖可以被裁掉。当前版本将分支界限改为：对每个目标取一个合法表示的最小依赖闭包，再对目标取最大值。这个下界较弱，但不会因为互斥表示而高估成本；Beam 排序仍可使用包含替代结构的启发式优先级。

固定 `beam=64`、`branch=16`、`representations=8`、`restarts=4`、`seed=1`，得到如下直接矩阵结果：

| 矩阵 | `D` | xordepth | MW_SQ + 重结合 | UBER_XSET + 重结合 |
|---|---:|---:|---:|---:|
| Cauchy 7→4 | 8 | 150 / 6 | 132 / 8 | 141 / 11 |
| systematic RS 7→4 | 7 | 149 / 6 | 136 / 7 | 126 / 11 |
| Vandermonde 8→17 | 6 | 221 / 6 | 219 / 6 | 212 / 8 |

这些结果说明安全下界修复消除了错误剪枝，但当前候选生成仍主要复现 CSE 型结构，尚未达到 MW_SQ 或 UBER_XSET 的直接结果。完整目标、电路和参数保存在 `data/results/xd_series_compare_safe_bound_20261102.json`。GF(2⁸) 80×80 矩阵上的两种 X-Sets 配置在 60 秒内均超时，因此没有把它作为 XOR 质量结论。

## Residual candidate 与独立 core

当前 xordepth 已经使用自己的 `xd::BitMask`、`xd::Circuit`、`xd::Json`、`xd::GreedyNetwork` 和 `xd::SeedRandom`，Makefile 不再链接 ringlab 的 `../src`。`tests/independence_test.sh` 会检查源码、头文件和构建入口没有 ringlab 依赖。

Residual 模式通过 `R' = normalize(R XOR {a,b,a XOR b})` 暴露跨表示共享节点，并用结构 fingerprint 分桶保留多种 Beam 状态。默认 `residual_pairs=0`，所以基础结果保持兼容；实验配置为 `residual_pairs=8`、`residual_targets=2`、`lookahead_width=8`、`lookahead_depth=2`、`beam=64`、`branch=16`、`restarts=4`、`seed=1`。

| 矩阵 | `D` | xordepth residual/diversity | MW_SQ + 重结合 | UBER_XSET + 重结合 |
|---|---:|---:|---:|---:|
| Cauchy 7→4 | 8 | 150 / 6 | 132 / 8 | 141 / 11 |
| systematic RS 7→4 | 7 | 149 / 6 | 136 / 7 | 126 / 11 |
| Vandermonde 8→17 | 6 | 221 / 6 | 219 / 6 | 212 / 8 |

完整 bounded circuit、参数、toggle 计数、Beam bucket 计数和外部重结合参照保存在 `data/results/xd_series_residual_bound_20261102.json`。Residual 候选扩大了可达状态空间，但在这组预算下尚未降低三类矩阵的 XOR 数；后续应继续优化候选评分和增量状态更新。

## 219 门 MW_SQ trace 回放

为了区分“搜索没有选到”与“状态模型无法表达”，对
`data/results/mw_sq_vandermonde_8to17_p0_3_0.json` 的原始 `ops` 做了逐门回放。
按当前普通 pair 规则，第一处不可达门固定为：

```text
step 56: node 120 = input 3 XOR node 118
```

在该时刻，`input 3` 出现在 3 个 representation 中，`node 118` 出现在 1 个
representation 中，两者没有共存于同一个 representation。因此普通 `apply()` 无法
生成该门。新增的跨表示 bridge transition 可以在两个 representation 上同步执行
`R' = normalize(R XOR {3,118,120})`，测试已经确认该门在状态模型中可表达。

当前 Vandermonde 搜索仍得到 `221 XOR / depth 6`。bridge 候选已经生成，但还需要记录
`bridge_scored`、`bridge_expanded` 和 `bridge_pruned_by_beam`，以定位它在 lookahead、
branch width、transposition 或 fingerprint Beam 阶段的剩余损失。
