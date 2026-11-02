# xordepth Safe Bound Implementation Plan

**Goal:** 修复互斥表示成本并集造成的错误剪枝，独立比较同 D 搜索质量。
**Architecture:** 将可证明下界与用于 Beam 排序的启发式估计分开；保持候选生成不变。所有电路用共享 C++17 core 验证。
**Spec:** docs/superpowers/specs/2026-10-31-depth-bounded-safe-lower-bound-design.md
**Execution:** Native，用户已要求同意并实现；保留工作区已有改动。

- [x] 增加 rows={3,15} 的四步回归：apply(0,1),(2,4),(3,4),(2,6)，历史下界不得超过完成电路的 3 门成本。旧代码第三步返回 4。
- [x] 小规模完整枚举，验证每个祖先 lower bound 不超过后继已验证完成线路成本。
- [x] best_cost_lower_bound 使用 max_target(min_representation(closure gates))；旧并集估计仅用于排序，命名明确为启发式分数，不用于安全剪枝。
- [x] 保持候选生成及 Beam 限额；复核 transposition 和 baseline，不输入外部 incumbent。
- [x] 复现 Cauchy D8/11、RS D7/11、Vandermonde D6/8；对照 MW_SQ/UBER_XSET(0,3,0)+rebalanced。参数、状态、完整电路和实际 wall time 放入新结果，不改历史记录。
- [x] 运行 xordepth tests、根 make test、diff checks；文档如实记录结果和不能声称精确的启发式删除。

第二阶段的候选生成、结构多样性 Beam 和增量 pair index 不在本计划内；当前多矩阵结果仍显示需要后续架构工作。

Review focus: 替代表示和死门；重复 mask 的复用；zero/alias；超时不构成最优证明；Pareto/Beam 删除仍是启发式。
