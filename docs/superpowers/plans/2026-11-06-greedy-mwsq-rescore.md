# MW_SQ 风格 greedy 候选重评分实施计划

Spec: `docs/superpowers/specs/2026-11-06-greedy-mwsq-rescore-design.md`

## Task 1: 支持度摘要纯函数

- 文件：`xordepth/include/xd/depth.hpp`、`xordepth/src/xd_greedy.cpp`、`tests/xd_search_test.cpp`。
- TDD RED：`mwsq_support_summary` 不存在（编译失败）。测试：空 rows、单元素行、无重复对、多行并列 max。
- GREEN：实现纯函数；O(每行 t(t-1)/2)。

## Task 2: greedy tie 重评分

- 文件：`xordepth/src/xd_greedy.cpp`、`include/xd/depth.hpp`。
- TDD：先加小矩阵测试断言 `--greedy-mwsq`（API 参数）开启时结果 verified/depth 合规且与关闭时可能不同但确定性（两次运行同 cost）；Vandermonde D=5 开启时 cost<=223。
- GREEN：在 `tied.size()>1` 且开关开启时，对每个 tied 候选试替换并打分 `(max_support, ties_at_max)`，取最优、仍并列维持随机；deadline 贯穿。

## Task 3: CLI 与 JSON

- 文件：`xordepth/apps/xordepth.cpp`、`src/xd_frontier.cpp`（序列化/验证）、`tests/cli_test.sh`。
- TDD RED：CLI 不识别 `--greedy-mwsq`；JSON 缺 `greedy_mwsq`。
- GREEN：`--greedy-mwsq` 开关、JSON 记录/验证（bool，旧 JSON 兼容）。
- 回归：flag 关闭时 sample 全部结果与当前一致。

## Task 4: 对比实验（主代理执行）

- sample 12：15s、seed1、`--greedy-mwsq` 开 vs 关，逐矩阵 cost/depth/秒；汇总胜/平/负。
- 若有净收益：full 57 同口径复测，并与 `xordepth_full_merged_20261106/summary.tsv` 合并比较（`xordepth merge`）。
- 报告：明确重评分有效/无效及量级；更新文档与日志；不 commit。
