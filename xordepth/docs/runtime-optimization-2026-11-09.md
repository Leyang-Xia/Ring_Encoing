# 无时限模式运行时长优化（结果位相容）

## 背景与依据

`xordepth bounded --no-timeout`（PR 对比协议：seed=1、beam=4、branch=4、
representations=4、restarts=2、greedy-restarts=2、--greedy-mwsq、
rewrite-trials=400、max-states=256）在 12 矩阵 sample 上耗时 25.30 秒，
约为 PR 纯 greedy 端口的 8 倍。采样（macOS `sample`）显示时间分布：

- greedy 组合（`synthesize_completable`）约一半，其中 `mwsq_support_summary`
  的 std::map 逐对插入约占全程 30%，`NodeStore::pair_info` 每次 BitMask 异或
  触发一次堆分配；
- Beam `SearchState::apply` 约一半，其中每次 apply 深拷贝父状态的
  `pair_index_`（随即被重建，纯浪费）并全量重建该索引，`BitMask` 的
  `std::vector` 表示让每节点/每掩码都是一次堆分配；
- `SearchState::candidates()` 的第二排序比较器每次比较做两次 std::map 查找。

## 改动清单（全部保持结果逐位一致）

1. **BitMask 小缓冲表示**（`include/xd/core.hpp`、`src/xd_core_mask.cpp`）：
   ≤2 个 64 位字（128 位输入以内）内联存储，超出才落堆。仍是精确多字掩码，
   值语义、`trim/toggle/^/==/</hash/decimal` 全部保持；移动构造清零源内联区
   以维持"size 之外恒为零"不变量。全局消灭热路径 malloc 风暴。
2. **mwsq 平铺计数**（`src/xd_greedy.cpp`）：`mwsq_support_summary` 由
   std::map 改为 uint64 键平铺 + 排序游程计数；tie-break 评估直接从表达式
   枚举键（不再物化 rows/复制 std::set）。`(max_support, ties_at_max)`
   按构造等价。
3. **greedy pair_info 按迭代缓存 + by_pair 平坦化**：NodeStore 在枚举期间
   不变，pair_info 结果按迭代缓存；可行 (pair,row) 出现记录平铺收集后一次
   排序分组，分组序 = 原 std::map 的 (left,right) 字典序，行序 = 升序，
   下游全序排序保证候选序列不变。
4. **SearchState::transition_base**（`include/xd/search.hpp`、
   `src/xd_search.cpp`）：apply/apply_toggle/apply_bridge 不再拷贝
   `pair_index_`/`residual_pair_support_`（每个转换都会全量重建；非法子
   状态从未被读取）。已证明 `pair_index_` 的哈希迭代顺序不可观测
   （candidates() 两轮排序为全序决定论、fingerprint() 对 supports 排序、
   pair_support/toggle/bridge 只计数）。
5. **apply 跳过未受影响 target**：不含 (left,right) 的 frontier 保留
   transition_base 的拷贝（把已 Pareto 的条目逐个重插回自身不可能增删条目，
   值恒等）；residual 调整仍按原不对称 -1/+1 顺序执行，residual 开启时
   语义逐位一致。
6. **candidates() 进度字段内联**：UBER 完成度累计存入 Candidate 自身字段，
   排序比较器不再做 map 查找；单次 lower_bound 同时完成查重与累计。
7. **闭包去重时间戳**：`best_cost_lower_bound`/`cost_priority` 的
   std::set 活集改为代计数戳（每个 entry 一个 generation / 全程一个），
   大小与成员关系精确复现；边界 require 先于 stamp 访问，非法节点仍以
   相同错误拒绝。

尝试并回退：`pair_index_` 改为排序平坦数组（每次 apply 全量排序 P 个元素
慢于哈希插入，sample 总时长 12.5→21.6 秒，已回退）。

## 验证（全部通过）

- **sample 12 矩阵（无时限协议）位相容**：`(cost, depth, status)` 与优化前
  基线逐位一致（基线 `notimeout_baseline_20261109` vs 最终
  `notimeout_opt7_20261109`）。
- **全量 57 矩阵位相容**：与 `xordepth_full_mwsq_20261106`（同配置、15s
  上限）**57/57 (cost, depth) 完全一致**（旧运行从未触到 15s 截止，故
  无时限等价）。结果目录 `xordepth_full_notimeout_20261109`。
- `make test`（xd/xd-search/trace-replay/matrix-suite/cli/independence）、
  xd_test 新增 `mwsq_support_summary` 金标准用例（含重复 id、空行、并列
  最大值等语义锁定）。
- ASan+UBSan（新代码重建 `build-sanitize/xordepth`）覆盖 bounded 主路径、
  residual toggle/bridge 路径、merge/verify，无报告。

## 提速效果

sample 12 矩阵无时限总时长：

| 版本 | 总时长 | 备注 |
|---|---:|---|
| 优化前（2026-11-09 基线） | 25.30s | gf8 11.89s、modular_64x96 7.71s |
| +mwsq 平铺/pair_info 缓存/transition_base | 20.87s | |
| +BitMask 小缓冲 | 15.34s | |
| +candidates() 进度内联 | 12.51s | |
| +greedy 平坦化 + 闭包时间戳 + 跳过未受影响 target | **11.60s** | gf8 4.92s、modular_64x96 3.65s |

净提速 **2.18×**（−54%）。全量 57 矩阵无时限合计 64.73s（含 5 个
infeasible_by_lower_bound 的 D 配置问题矩阵，秒级完成）。结果仍严格
支配 PR：sample 与 PR V2 12/12 一致，全量 26 优 26 平 0 劣。

基准脚本：`data/results/run_notimeout_bench.py`（sample 协议）、
`data/results/run_full_notimeout_bench.py`（全量、读参照 D 表）。
中间版本目录 `notimeout_opt1..opt6_20261109` 保留作实验轨迹。
