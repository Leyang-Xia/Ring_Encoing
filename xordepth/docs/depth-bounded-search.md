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

`--threads N`（默认 1）为 greedy portfolio 的 restart 与 Beam 的层内扩展分配工作线程。确定性通过三段式层结构保持：顺序发布 complete 状态并快照 incumbent 剪枝阈值 → 各父状态独立并行扩展 → 按原顺序回放转移表查重。因此任意线程数的结果与串行逐门一致（`--timeout` 运行除外，其结果本就依赖机器速度）。全量 57 矩阵无时限运行 64.73s（串行）→ 38.63s（4 线程，**1.68×**）；任意线程数输出与串行逐门一致，ThreadSanitizer 无竞态报告。单元测试 `exercise_threaded_determinism` 锁定该性质。

MW_SQ 并列重评分循环本身也已并行化：每个 tied 候选的评分只依赖 key 多重集，因此按行分块并行生成 key（各块由工作线程排序）后按固定块序归并、串行游程扫描，再按 tied 顺序归约出与串行一致的 winner 列表；`threads=1` 走原串行路径。为避免超额订阅，重启并行与评分内部并行共享线程预算（每个 attempt 分得 `threads/重启数` 个工作线程）。加速上限仍受 Amdahl 制约：Beam 的顺序发布/回放阶段与 CSE 基线为串行。

## 深度下界与安全剪枝

分支界限使用的成本下界为：对每个目标取一个合法表示的最小依赖闭包，再对目标取最大值。这个下界较弱，但不会因为互斥表示而高估成本（早期版本对所有 Pareto 表示的依赖取并集，会把互相替代的候选依赖重复计入，导致错误剪枝，已修复）；Beam 排序仍可使用包含替代结构的启发式优先级。

`infeasible_by_lower_bound` 是唯一的数学不可行判据：当 `D` 小于行权重深度下界 `⌈log₂(max row weight)⌉` 时，任何二输入 XOR 网络都不可能在深度 `D` 内完成。

## canonical 输入

MW_SQ 的候选评分包含平局选择，目标行的输入顺序会影响候选访问顺序。因此下面两条命令虽然计算同一个目标集合，可能得到不同的启发式上界：

```bash
./build/xsets 0 3 0 MW_SQ --quiet < data/matrices/vandermonde_8to17_direct.txt
# 220 XOR

./build/xsets 0 3 0 MW_SQ --quiet < data/matrices/vandermonde_8to17_direct_canonical.txt
# 219 XOR
```

`ringlab analyze` 的 backend 适配层会用 `std::set<Mask>` 对目标行去重并排序，因而使用 canonical 顺序。为保证直接 X-Sets 调用与 `ringlab`、xordepth 的参照一致，仓库保留了排序后的 [`vandermonde_8to17_direct_canonical.txt`](../data/matrices/vandermonde_8to17_direct_canonical.txt)。这是行顺序导致的启发式路径差异，不是电路验证或 XOR 计数差异。

## 正确性边界

- 深度按二输入 XOR 门的拓扑深度计算，输入和常量零深度为 0；
- 补零、接线和输出别名不计 XOR；
- Beam 宽度、表示上限和超时都会改变搜索覆盖范围；
- `search_failed` 和 `timeout` 不说明不存在深度不超过 `D` 的线路；
- 当前版本不复制同一状态中的重复 mask 节点，也不建模 LUT、布线、扇出时延或多输入 XOR 门。

## 独立 core

xordepth 使用自己的 `xd::BitMask`、`xd::Circuit`、`xd::Json`、`xd::GreedyNetwork` 和 `xd::SeedRandom`，Makefile 不链接 ringlab 的 `../src`。`tests/independence_test.sh` 会检查源码、头文件和构建入口没有 ringlab 依赖。

---

历史实验记录（PR sample 对比、逐版本全量基准、小矩阵 frontier/audit/bounded 测量）已从 `data/results/` 移除，完整数据保留在 git 历史中；当前结果以 [`data/results/xordepth_full_threads4_mwsqpar_rerun_20261109/`](../data/results/xordepth_full_threads4_mwsqpar_rerun_20261109/bench.tsv) 全量 57 矩阵运行为准。
