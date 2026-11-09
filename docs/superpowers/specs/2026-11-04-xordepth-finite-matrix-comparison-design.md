# xordepth 有限搜索、矩阵套件与 PR 对比设计

## 目标

1. 建立 10x10 至 100x100 的数十个确定性二进制矩阵，划分快速 sample 和完整 full set。
2. 保证 `bounded --no-timeout` 没有墙钟截止时仍由结构性预算有限终止；预算耗尽不表示不可行。
3. 在相同矩阵、深度、seed 和无墙钟截止口径下比较当前 Beam 与 PR #6 OUR_DCCSE V1/V2，记录结果质量与耗时。
4. 只有 PR 在 sample 上显示稳定优势时才移植其优化；最终线路必须逐门验证并满足硬深度。

## 已确认事实

- 当前未提交工作区已有 9 个代表矩阵、`--no-timeout`、`max_states` 和重写边界修复，必须保留。
- Vandermonde 8→17、64 输入、24 输出、D=5 是关键回归。PR V1/V2 分别报告 218/5、217/5；当前 Beam 长时间没有可行线路。
- 当前 `combine_depth()` 已对每个表示做固定叶子集合的完成深度检查，因此重复加入 Kraft/W 检查不是根本修复。
- 当前 Beam 的 residual/toggle、替代表示和新节点可产生非进展路径；`max_states` 目前多探索一个状态，且不覆盖 lookahead 的临时状态。
- PR greedy 每次接受的候选严格减少表达式项数总和，外层循环有限；固定表达式用 `W=sum 2^depth <= 2^D` 保持可完成性，最后用 earliest-ready 合并并验证。

## 矩阵套件

新增 C++17 确定性生成器，不引入 Python。生成 48 个静态矩阵：dense、sparse、blocks、modular 四类，尺寸 10x10、12x20、24x40、40x64、64x96、100x100，每类每尺寸两个 seed。名称统一为 outputs x inputs。

- `sample.manifest`：12 个代表，覆盖四类结构、小/中/大尺寸，并包含 Cauchy、RS、Vandermonde 和 GF(2^8)。
- `full.manifest`：48 个生成矩阵加现有 9 个代表矩阵，共 57 个。
- 默认 `make test` 只跑 sample；`make test-full` 跑 full。
- 每个矩阵检查解析、尺寸非空、行权重下界、无共享平衡调度、逐门功能、重写后成本和深度不回退。

## 有限性

- `max_states` 明确定义为实际从 Beam 取出并处理的累计状态数，跨 restart 共享；达到预算后不再处理下一个状态。lookahead 是有限临时探测，不伪装成全局工作预算。
- 每个 restart 另设 Beam 轮次上界 `outputs * (2^(D+1)-1)`，使用饱和 `size_t` 运算。这是每个输出独立深度 D 二叉树的保守节点数上界，共享只会减少所需节点。
- 到达最后一轮时仍导出已完成状态，但不扩展未完成状态。轮次耗尽返回 `search_failed`；状态预算耗尽沿用 `timeout`，表示计算预算结束而非不可行。
- 结果 JSON 记录轮次预算及实际最大轮次，旧 JSON 缺少新字段仍可验证。

## 公平对比与移植门槛

第一阶段不移植 PR，使用有限性修复后的当前实现作为基线。PR 在其隔离构建中运行。统一样本、矩阵行序、D、总 restart/seed 口径，不设置墙钟超时；分别记录 status、cost、depth、elapsed、状态和 restart 预算。超时/搜索失败不当作高成本线路参与门数平均。

若 PR 在关键回归和 sample 上有明确优势，第二阶段移植最小 greedy 路径：

1. 每个输出维护一个当前 GF(2) 表达式；
2. 候选只在替换后可完成深度不超过 D 时接受；
3. 只接受严格减少表达式项数总和的候选；
4. earliest-ready 完成各输出，prune、verify、depth_profile；
5. 有界 restart/lookahead，保留现有 deadline 检查；
6. 将它作为 bounded 的额外 baseline/incumbent，不删除现有 Beam 表达力。

重写仅在得到深度可行 incumbent 后执行。吸收后重新跑相同样本，区分“有限性修复带来的收益”和“PR greedy 带来的收益”。

## 非目标

- 不声称启发式结果全局最优。
- 不修改 third_party/xsets 或 third_party/uber。
- 不重写历史证书。
- 不用单一 64 位整数保存任意宽度 mask。
- 不把无结果、search_failed 或预算耗尽解释为不可行证明。
