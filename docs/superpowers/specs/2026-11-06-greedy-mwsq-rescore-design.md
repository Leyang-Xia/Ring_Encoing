# MW_SQ 风格 greedy 候选重评分设计

日期：2026-11-06

## 背景与依据

### Plank DSN-2012 论文（Heuristics for Optimizing Matrix-Based Erasure Codes, Plank/Schuman/Robison, DSN'12）

- **X-Set 框架**：贪心迭代加元素；每个未完成 target 维护若干 X-Set（能异或出该 target 的已调度元素集合）。Theorem 4.1：大小为 2 的 X-Set 应立即生成 target。Threshold/L 两个参数扩展可发现带抵消的公共子式。
- **MW2（代码库 `MW_SQ` 的来源，X-Sets.cpp `EvaluateHeavyCandidate:1345`）**：对每个最大权候选边，克隆调度试放并跑廉价归约到不动点，度量归约后的最高边权 `newtop` 与顶权边数 `newn`，按 `(newtop↓, newn↓)` 选。本质是"试放后剩余共享结构越集中越好"的一步 lookahead 重评分。论文实测 MW2 在 CRS 编码/解码上产生最佳调度但最慢（每步对所有 top 边各克隆一次调度）。
- **Uber-XSet**：Uber-CSHR 与 X-Set 的混合——先只考虑 min X-Set 尺寸最小（最接近完成）的 targets，再在这些 targets 可缩小的候选中选总权重（能缩小的 target 数）最大者。论文综合评测中最优的通用启发式。xordepth beam 的 UBER 轻量排序即源于此思想，但无模拟。
- **对本设计的启示**：xordepth greedy 的 `tied` 随机选择正对应 MW 只按权选、不区分并列的情形；补上 MW2 式模拟重评分（对每个 tied 候选试替换、度量剩余 pair 支持度分布）即为论文中被验证有效的方向。`best_source` 数据显示最终结果几乎全部来自 greedy/cse_baseline/rewrite，beam 极少胜出，因此重评分挂在 greedy 选择处，而非 beam 排序。

## 方案

在 greedy 每轮选出 `tied`（rank 相同）候选后，若启用 `mwsq_tiebreak` 且 `tied.size() > 1`：

1. 对每个 tied 候选试替换：复制全部 expressions，对该候选 `uses` 的行应用 `replaced()`（node 不存在时用虚拟 id `nodes.masks.size()`，与 lookahead 的 `candidate_node` 做法一致，不物化节点）。
2. 对替换后的全部行统计无序 pair 支持度 `count[{a,b}]`，得到 MW_SQ 语义的 `(max_support, ties_at_max)`。
3. 按 `(max_support↓, ties_at_max↓)` 取最优；仍并列则维持原随机选择（种子确定性不变）。
4. 全程 `check_deadline`；不改变任何 gain>0/可行性判定，只重排 tie，因此严格有限性与正确性论证不变。

## 接口与 CLI

- `BoundedScheduleOptions` 增加 `bool greedy_mwsq = false;`。
- `synthesize_completable` 增加带默认值的尾部参数（或轻量 struct），默认关闭保持位相容基线。
- CLI：`--greedy-mwsq` 布尔开关；bounded JSON 增加 `"greedy_mwsq": bool`，旧 JSON 缺字段仍可验证。

## 可测试性

在 `xd` 命名空间暴露纯函数便于单元测试：

```cpp
// rows 为每输出表达式的节点 id 集合；返回 (max_support, pairs_at_max_support)
std::pair<std::size_t, std::size_t> mwsq_support_summary(
    const std::vector<std::vector<int>> &rows);
```

## 非目标

- 不改 beam `candidates()` 排序；不改 rewrite；不改 incumbent 池；不改默认行为（flag 关闭时所有输出与当前位相容）。
- 不承诺门数必然下降——由对比实验判定；重评分只改变 tie 走向。

## 验收

- `mwsq_support_summary` 单元测试覆盖空行/单元素/全并列。
- flag 关闭时：sample/full 结果与当前完全一致（回归）。
- flag 打开时：Vandermonde D=5 仍 <=223/5 且 verified；12 矩阵 sample 全 verified、确定性（同配置两次运行 cost 相同）；ASan/UBSan 通过。
- 对比实验：sample 12 矩阵、15s 预算、同 seed，`--greedy-mwsq` 开 vs 关，逐矩阵报告 cost/depth/耗时；有效果再跑 full 57 并与 `xordepth_full_merged_20261106` 汇总比较。
