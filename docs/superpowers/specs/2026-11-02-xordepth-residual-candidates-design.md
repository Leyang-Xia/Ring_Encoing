# xordepth Residual Candidate 与结构多样性 Beam 设计

## 目标

扩展独立的 `xordepth` 深度受限搜索，使它能够生成当前 MW_SQ/UBER_XSET 路径中出现、但不会同时出现在同一 Pareto 表示中的共享节点。第一阶段安全 lower bound 已经修复；本阶段只改变候选生成、表示转移和 Beam 保留策略，不使用 MW_SQ/UBER_XSET 电路作为搜索输入。

目标仍是固定 GF(2) 矩阵和硬深度 `D` 下最小化完整固定输入到固定输出 XOR 门数。所有最终结果必须逐门验证，且不能把启发式搜索结果解释为全局最优。

## 已确认的现状问题

当前候选索引只扫描一个目标表示内部的两个 term：

```text
pair_index[edge] <- (target, representation, left, right)
```

`apply()` 也只在两个操作数同时存在时做替换。对已保存的 MW_SQ `219/6` 重结合线路按拓扑顺序诊断时，前 56 个门可以按当前规则生成；第一个不可达门为 `node 121 = input 3 XOR node 118`，两个操作数分散在不同目标表示中。贪心地尝试所有可用参考门最多可生成 194/219 个门，说明主要缺口在表示转移，而不是 lower bound。

X-Sets 的 MXS/XSet 状态允许通过全局 edge 和残余集合继续暴露后续候选；当前 xordepth 的 pair 共存限制没有表达这类转移。

## Residual/toggle transition

### 代数定义

对任意目标表示 `R` 和已物化节点 `a,b`，令：

```text
c = a XOR b
R' = normalize(R XOR {a,b,c})
```

其中集合 XOR 表示对 term 出现次数取模 2。因为 `mask(a) XOR mask(b) XOR mask(c) = 0`，所以 `R'` 与 `R` 计算同一个目标。

这允许以下情况：

- `a` 和 `b` 同时存在：退化为现有 pair replacement；
- 只有一个存在：通过 toggle 暴露另一个操作数和新节点；
- 两个都不存在：只在 residual 预计有明显收益时尝试，默认关闭以控制分支数。

### 合法性检查

每个 toggle 候选必须满足：

1. `a != b` 且两个节点存在；
2. `depth(c) = 1 + max(depth(a), depth(b)) <= D`；
3. 对受影响目标生成 `R'` 后，`arrival_depth(R') <= D`；
4. `R'` 归一化后仍然是合法表示；
5. 至少一个目标的 residual 项数、完成深度或可复用支持得到改善；
6. 同一状态中相同 `(a,b,c,target)` 只保留一次。

新节点按 mask 去重；已存在相同 mask 时复用节点，不增加 XOR 门。最终 circuit 仍由 `Circuit::verify()`、`depth_profile()` 和 `pruned()` 检查。

### 候选限制

不枚举所有节点对和所有目标的笛卡尔积。候选生成分两步：

1. 从当前全局节点 mask 中生成有限的高支持 pair；
2. 只对 pair 影响的目标表示生成 toggle residual。

新增选项：

```text
--residual-pairs N      每轮最多评估的全局 pair 数
--residual-targets N    每个 pair 最多扩展的目标表示数
--allow-empty-toggle    是否允许两个操作数都不在表示中
```

默认 `allow-empty-toggle=false`。第一版只允许至少一个操作数存在，先验证 MW 路径缺失门能否被重新暴露。

## Candidate 评分

每个普通 pair 或 residual candidate 计算整数评分：

```text
affected_targets       受影响目标数
term_reduction         residual 项数减少量
depth_slack            D - arrival_depth(R')
future_support         新节点在下一轮可形成的 pair 支持
toggle_penalty         新增 term 数惩罚
```

排序优先级：

1. `future_support` 降序；
2. `term_reduction` 降序；
3. `affected_targets` 降序；
4. `depth_slack` 降序；
5. `toggle_penalty` 升序；
6. mask 和节点编号稳定排序。

评分只用于候选排序，不能用于证明剪枝。安全下界仍使用第一阶段的 per-target minimum closure bound。

## 有限 MW_SQ 式 lookahead

对排序后的前 `K` 个候选，复制状态并执行最多 `L` 次后续闭合：

1. 应用候选或 toggle；
2. 只在受影响目标上更新 residual index；
3. 选取下一个 future_support 最大的合法候选；
4. 记录新 frontier 的最大支持、完成目标数和 residual 宽度。

选择候选时优先最大化后续 frontier 指标，再比较安全 lower bound。整个 lookahead 只产生临时状态，不改变父状态。超时检查必须贯穿复制、索引和评分过程。

## 结构多样性 Beam

当前 Beam 全局按一个 `cost_priority` 排序，会保留大量 CSE 型相似状态。新增结构指纹：

```text
fingerprint = hash(
    per-target minimum terms,
    per-target arrival-depth histogram,
    top pair-support histogram,
    number of toggle representations,
    live closure estimate
)
```

每轮先按 fingerprint 分桶，再按以下顺序在每桶内选代表：

1. 安全 lower bound 小；
2. future support 大；
3. 已完成目标数多；
4. 最大输出深度小；
5. signature 稳定排序。

Beam 总大小仍为 `beam_width`，每个非空桶至少保留一个状态；剩余名额按桶内评分循环分配。这样可以保留低门数、低深度和高共享潜力的不同结构。

## 接口和结果

在 `BoundedScheduleOptions` 中加入 residual 参数，并将其写入 bounded JSON；历史结果缺少这些字段时仍可验证。结果中新增：

```text
residual_pairs
residual_targets
allow_empty_toggle
diversity_buckets
toggle_candidates
```

`toggle_candidates` 只作为搜索计数，不代表新增门数。所有完整电路继续保存全部 ops、outputs、target masks 和验证元数据。

## 测试

### 单元测试

1. 构造一个只有一个操作数存在的表示，验证 `R XOR {a,b,c}` 保持目标 mask；
2. 验证 toggle 深度检查和 D 超限拒绝；
3. 验证相同 mask 复用不增加 gate；
4. 验证 empty toggle 默认关闭，开启后受 `residual_targets` 限制；
5. 验证结构指纹相同状态只保留一个代表，不同 fingerprint 在 Beam 中同时存在；
6. 验证临时 lookahead 不改变父状态；
7. 重新运行已有 lower-bound 反例和所有 `make test`。

### 多矩阵验证

不传 incumbent，固定：

```text
beam=64, branch=16, representations=8, restarts=4, seed=1
```

在 Cauchy、systematic RS、Vandermonde 三个矩阵上分别比较：

- 当前安全下界版本；
- residual candidate 版本；
- residual + diversity Beam 版本。

MW_SQ/UBER_XSET 只作为独立对照，先重结合再按相同 `D` 比较 XOR 和深度。GF(2⁸) 80×80 只作超时压力测试。

## 不做的事

- 不把 MW_SQ/UBER_XSET 电路作为候选或 incumbent；
- 不修改 `third_party/xsets/X-Sets.cpp` 或 `third_party/uber/Uber.cpp`；
- 不使用浮点评分、Python 或新的运行时依赖；
- 不在第一版加入任意高阶子集枚举；
- 不把 timeout、Beam 耗尽或单次实验结果解释为最优性证明。
