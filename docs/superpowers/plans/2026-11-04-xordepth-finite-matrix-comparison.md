# xordepth 有限搜索、矩阵套件与 PR 对比实施计划

Spec: `docs/superpowers/specs/2026-11-04-xordepth-finite-matrix-comparison-design.md`

## Task 1: 静态矩阵与分层测试

- 文件：`xordepth/apps/xordepth_matrix_gen.cpp`、`xordepth/data/matrices/generated/*`、`sample.manifest`、`full.manifest`、`tests/matrix_suite_test.cpp`、`Makefile`、`docs/matrix-suite.md`。
- TDD：先让测试要求 sample=12、full=57、尺寸范围和四类结构；确认旧测试因 manifest/数据缺失失败；再用 C++17 生成器生成静态数据并实现 manifest 测试。
- 验证：`make -C xordepth test` 跑 sample，`make -C xordepth test-full` 跑 full，所有线路验证通过。

## Task 2: 有限性回归

- 文件：`include/xd/depth.hpp`、`include/xd/search.hpp`、`src/xd_search.cpp`、`src/xd_frontier.cpp`、`tests/xd_search_test.cpp`、`tests/cli_test.sh`。
- TDD RED：对 `{7,15}`、D=2、no-timeout、`max_states=1` 断言恰好探索 1 状态；当前实现探索 2。增加轮次上界算术、跨 restart 状态预算和最后一轮完成测试。
- GREEN：修正预算判断顺序；实现饱和 `completion_round_limit()` 和每 restart 轮次限制；序列化统计并保持旧 JSON 兼容。
- 验证：聚焦测试、CLI 和 `make -C xordepth test` 全绿；Vandermonde D=5 no-timeout 在外部 watchdog 内自行结束。

## Task 3: PR 基线对比

- 文件：独立结果 `xordepth/data/results/xd_pr6_sample_baseline_20261104.json` 或等价格式及对比文档。
- 使用 PR #6 隔离构建，不合入其代码。对 sample 在同一 D/seed/restart 口径下运行当前 Beam、OUR_DCCSE V1、V2；保存每条 verified circuit 或 PR 原始证据路径。
- 验证：所有 verified 结果由对应验证器复核；汇总成功数、cost/depth、耗时，不把失败当成本。

## Task 4: 按需移植 PR greedy

- 文件：`include/xd/depth.hpp`、`src/xd_depth.cpp` 或独立 `src/xd_greedy.cpp`、`src/xd_search.cpp`、`tests/xd_search_test.cpp`、`Makefile`。
- TDD RED：添加 Vandermonde D=5 可行回归、表达式规模严格下降、取消/零行/重复输出、实际 reused-node depth 和 deadline 测试。
- GREEN：实现有界 completable greedy，并作为额外 incumbent；所有候选维护完成深度，主循环以表达式项数总和为严格下降测度。
- 验证：关键回归达到 PR 同量级或明确记录差距；sample 复测，full correctness、根 `make test` 和 sanitizer 全绿。

## Task 5: 评审与交付

- 对照 spec 检查有限性、验证语义、数据可复现性和对比公平性，再查代码质量。
- 保存最终基准和方法说明；更新工作区日志。
- 不自动提交、合并或关闭 PR，向用户报告吸收内容、实测表格和剩余限制。
