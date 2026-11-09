# xordepth 统一 incumbent 候选池设计

## 背景与已复现事实

2026-11-06 的 UBER 风格轻量排序全量运行结果：`xordepth/data/results/xordepth_full_uber_20261106/summary.tsv` 共 57 行，52 个 verified、5 个 infeasible_by_lower_bound，无错误或深度违规；与 11-04 全量结果成本完全一致。

四个高参数无墙钟截止任务均确定性结束于 `max_states=10000`：Cauchy 147/5、RS 146/5、Vandermonde 223/5、GF(2⁸) 1235/6；耗时分别约 31.18、32.42、59.03、705.32 秒。高参数没有带来新的门数下界。

`dense_100x100_seed123` 的 5 秒全量结果为 1819/7，而 10 秒 recheck 为 1809/7。线路哈希显示 5 秒结果等于 V2 证据，10 秒结果等于 V1 证据。直接调用各候选生成器得到：

```text
plain greedy                 1819 / 7
greedy pair, no deadline     1818 / 7
CSE baseline                 1809 / 7
rebalanced CSE baseline      1809 / 7
greedy pair, 5s deadline     1819 / 7
```

当前 `synthesize_bounded()` 先运行 greedy portfolio，再运行 CSE baseline。5 秒截止在 greedy 内发生时，函数提前返回，1809 的 CSE baseline 从未进入候选池；10 秒时它来得及发布。因此这是候选发布顺序和统一池缺失问题，不是 UBER 排序失效，也不是需要更大搜索预算。

## 目标

1. 所有已验证候选在任何后续阶段超时前都立即进入统一池：外部 incumbent、CSE baseline、rebalanced baseline、每个 greedy attempt、Beam complete state、rewrite 输出。
2. 用单一 `(live XOR cost, actual depth)` 字典序选优，并在入池前逐门验证和检查硬深度。
3. 廉价候选（外部 incumbent 与 CSE baseline）先于耗时 greedy/Beam 生成，避免墙钟截止吞掉候选。
4. 保持原有 beam/branch/restart/greedy/rewrite/max_states/timeout 预算，不增加搜索时间；rewrite 使用共享总 trial 预算。
5. 提供 `--incumbent FILE` 与 `merge` 命令，使 V1/V2/历史 verified 线路可以零搜索时间加入最终候选池。

## 接口

- `BoundedScheduleOptions` 增加 `std::vector<std::string> incumbent_paths`。
- `BoundedScheduleResult` 增加 `published_candidates`、`incumbents_used`、`incumbents_skipped`、`rewrite_trials_used`、`best_source`。
- `synthesize_completable()` 增加可选 `CircuitObserver` 回调；每完成一个 attempt 就调用，不等到整个 portfolio 结束。
- `rewrite_depth_gates()` 增加可选 `used_trials` 输出，供共享预算记账；deadline 异常由调用方捕获并保留池内 incumbent。
- CLI 增加可重复 `--incumbent FILE` 和 `xordepth merge BASE.json CANDIDATE.json... --output FILE`。
- bounded JSON 序列化新计数；旧 JSON 缺字段仍可验证。

## 有限性与正确性边界

- 候选池只保存已 `verify(targets)` 且 `depth_profile().depth <= D` 的线路。
- 不同变换或重复线路跳过并计数；错误 JSON 直接失败，不静默吞掉。
- 超时、状态预算耗尽、搜索失败都不表示不可行；若池内已有 verified incumbent，仍返回 verified 和完整线路。
- rewrite 只在池内出现新的 best circuit 时尝试，总 trials 不超过 `rewrite_trials`。
- `merge` 不运行任何搜索，只读取精确 JSON、验证线路、选择 best 并写出完整 bounded 结果。
- 不修改历史证书；新结果写入独立目录。

## 测试

1. 单元：候选池对同 cost 更低 depth、错误变换、超深度、重复线路的选择/拒绝。
2. 回归：使用已存 V1 incumbent 和预过期 deadline，断言仍返回 V1 的 verified circuit。
3. 集成：`dense_100x100_seed123`、D=7、5 秒、原全量参数应返回 <=1809/7。
4. CLI：`--incumbent` 与 `merge` 的输出可由 `verify` 复核；tampered depth/cost 被拒绝。
5. 全量：同 5 秒原预算重跑 57 矩阵，统计改善/回退/错误；sample、full、根测试与 sanitizer 通过。

## 非目标

- 本阶段不实现 MW_SQ 风格候选重评分。
- 不增加 beam/branch/restart/greedy 次数或墙钟预算。
- 不改写 UBER 排序逻辑。
- 不自动提交、合并或关闭任何 PR。
