# xordepth 统一 incumbent 候选池实施计划

Spec: `docs/superpowers/specs/2026-11-06-unified-incumbent-design.md`

## Task 1: 候选池与接口

- 文件：`xordepth/include/xd/depth.hpp`、`xordepth/include/xd/search.hpp`、`xordepth/src/xd_search.cpp`、`xordepth/src/xd_greedy.cpp`、`xordepth/src/xd_rewrite.cpp`。
- TDD RED：新增候选池测试和预过期 deadline + V1 incumbent 回归，看到当前返回 timeout/1819。
- GREEN：实现统一 `publish_candidate()`；所有来源逐门验证、深度检查、字典序选优；CSE baseline 与 rebalanced baseline 先于 greedy 发布；greedy 每完成一个 attempt 即回调发布；Beam complete state 和 rewrite 输出也入池；共享 rewrite trial 预算。
- 验证：`build/xd-search-test` 与 `dense_100x100_seed123` 5 秒原预算返回 <=1809/7。

## Task 2: 外部 incumbent 与离线 merge

- 文件：`xordepth/apps/xordepth.cpp`、`xordepth/src/xd_frontier.cpp`、`xordepth/include/xd/experiments.hpp`。
- TDD RED：CLI 使用 `--incumbent` 或 `merge` 时提示未知命令/选项。
- GREEN：支持可重复 `--incumbent FILE`；`merge BASE CANDIDATES... --output` 读取所有嵌套 circuit，按精确 target masks 验证，选择 best，保留完整 bounded JSON；不同变换/重复/超深度线路跳过或拒绝。
- 验证：CLI 回归、tampered JSON 拒绝、`xordepth verify` 复核。

## Task 3: 序列化与兼容性

- 文件：`xordepth/src/xd_frontier.cpp`、`xordepth/tests/cli_test.sh`。
- JSON 保存 `published_candidates`、`incumbents_used`、`incumbents_skipped`、`rewrite_trials_used`、`best_source`、`incumbent_files`。
- 旧结果缺字段仍可验证；unsigned 字段使用 `Json::uint64()`；不允许宣称超预算的 rewrite 或 explored states。
- 验证：旧 11-04/11-06 bounded JSON 复核通过。

## Task 4: 复测与报告

- 用原 5 秒全量配置重跑 57 矩阵，输出独立 `summary.tsv` 与比较统计。
- 只统计 verified 成本；infeasible 与 timeout 不计为门数结论。
- 运行 `make -C xordepth test`、`test-full`、`make test`、ASan/UBSan；报告 dense100 改善、全量改善/回退数量和耗时。
- 更新文档及工作区日志，不提交。
