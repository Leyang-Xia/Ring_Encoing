# 深度受限 XOR 调度设计

## 目标

为 `xordepth` 增加一个受深度硬约束的启发式综合器。给定固定 GF(2) 目标矩阵和最大门级深度 `D`，输出必须满足所有输出深度 `<= D`；在搜索预算内尽量减少完整固定输入到固定输出变换的二输入 XOR 门数。

搜索失败只能报告为启发式搜索未找到解。只有 `D < depth_lower_bound` 时，才允许报告由下界证明的不可行。`D` 不应被转换成软惩罚项。

## 现状与边界

当前 `xordepth` 已有：

- `depth_profile()`：逐门计算深度、输出深度和 fanout；
- `rebalanced()`：不增加门数的累加树重结合；
- `synthesize_schedule()`：按共享 pair 频次的 Plank 式前向调度；
- `frontier()`：扫描共享阈值并保存完整、可复核的线路；
- 独立的 X-Sets/Uber 第三方源码，源码不能修改。

## 第一阶段源码审计结论

第一阶段对 `third_party/xsets/X-Sets.cpp` 和 `third_party/uber/Uber.cpp` 做了逐函数审计，并完成了保持候选选择语义的 X-Sets 性能补丁。X-Sets 的可复用核心是：每个目标维护多个 `MXS/XSet` 表示，全局 edge weight 统计一个候选 XOR 在多少目标中有用，`MW_SQ` 对候选做一步前瞻后比较后续 edge frontier，`UBER_XSET` 优先选择同时覆盖最小剩余集合的高权重 edge。Uber 的可复用核心是：对当前最容易完成的目标保存局部最小 Hamming 距离表示，再用小子集枚举寻找反向改进。

第一阶段补丁已验证 `MW_SQ` 和 `UBER_XSET` 的固定输入计数不变，同时显著减少运行时间。深度版本复用这些算法思想，但不直接修改第三方状态机：使用当前 C++17 `Mask` 作为 key，维护增量 pair 支持索引，并把 `MW_SQ` 的完整 Schedule 复制替换为受限 Beam 状态扩展。

审计还确认两个语义边界：X-Sets/Uber 的“最小项数”单一表示不足以表达深度约束；深度版本必须保留每个目标的一组局部非支配 `(remaining_terms, completion_depth)` 表示。第一版不允许重复 mask 节点，避免把节点复制的面积代价与深度约束混在一起；完成基础版本后，再以独立实验开关评估受控复制。

现有 `--depth-aware` 只影响相同频次 pair 的平局选择，不能形成深度硬约束。因此新算法必须使用独立的 `DEPTH_BOUNDED_BEAM` 路径，复用 `Mask`、`Circuit`、`GreedyCircuit` 和现有验证函数。

## 核心算法

### 1. 深度预算

对 `D` 做迭代搜索，从 `depth_lower_bound(targets)` 开始。对每个 `D`：

1. 先生成并重结合现有 CSE 基线，作为可行性上界和门数 incumbent；
2. 若基线深度 `<= D`，它保证该深度预算至少有一个可行解；
3. 受限搜索只接受能在 `D` 内完成的状态；
4. 在所有 verified 候选中选择门数最少者，门数相同再选择深度更小、fanout 更小者。

若 `D` 小于下界，直接写 `infeasible_by_lower_bound`。若 `D` 不小于下界但预算耗尽没有候选，写 `search_failed`，不能写不可行。

### 2. 搜索状态

每个 Beam 状态包含：

- 已物化节点：`Mask -> {node, depth}`；输入和零节点深度为 0；
- 每个输出的规范化剩余项列表；每项对应一个已物化节点；
- 当前 live XOR 门数；
- 每个节点的消费者数；
- 每个输出的局部非支配表示集合，每项包含剩余项数和预计完成深度；
- 用于增量更新的受影响输出集合。

状态比较和去重使用规范化 mask 表示，不依赖临时 node ID。相同规范化状态只保留门数、关键深度和 fanout 更优的状态。

### 3. 候选生成

候选来自两类来源：

1. **X-Sets/Uber 风格前向 pair**：统计当前剩余表达式中 pair 的跨目标支持数，保留支持数最高的前 `pair_limit` 个候选；相同支持数按预计到达深度、关键路径影响、mask 稳定顺序排序。
2. **局部反向拆分**：对剩余项数较大的目标，枚举 `target = A XOR B` 的有限候选拆分。优先选择已存在节点、低 `depth(A)`/`depth(B)`、以及能被多个目标复用的 `A`/`B`。拆分枚举受 `split_limit` 限制，不进行无界子集枚举。

已有相同 mask 的节点只复用，不新增门。新 pair 的深度为 `1 + max(depth(a), depth(b))`；若超过 `D`，立即丢弃。

### 4. 深度剪枝

每次候选应用后，只重新计算受影响输出：

- 使用当前节点到达深度做最早就绪优先的平衡累加估计；
- 若某个输出的估计完成深度超过 `D`，丢弃状态；
- 对每个剩余表达式使用 `ceil(log2(weight))` 作为安全的局部深度下界；
- 对剩余项数、当前门数和 incumbent 做廉价支配剪枝；
- 当所有输出均已完成时，再用真实 `Circuit` 导出和 `depth_profile()` 做最终确认。

估计器只能用于剪枝和排序，最终结果必须通过完整逐门验证。任何估计失败不得被解释为不可行证明。

### 5. Beam 与重启

每轮扩展状态后：

- 按整数词典序排序：预计新增门数、最大输出深度、关键输出深度余量、fanout、mask 顺序；
- 对完全相同的规范化状态去重；
- 只保留 `beam_width` 个状态；
- 每个 seed 使用确定性 `SeedRandom` 改变候选等价项顺序；
- `restarts` 次搜索共享同一深度预算和 incumbent，避免重复验证已知较差状态。

默认参数必须保证小矩阵可在常规测试时间内完成，并允许 CLI 覆盖。

## 接口与结果

在 `include/xd/depth.hpp` 中增加：

```cpp
struct BoundedScheduleOptions {
    int max_depth = 0;
    int beam_width = 64;
    int pair_limit = 32;
    int split_limit = 16;
    int restarts = 8;
    std::uint64_t seed = 0;
};

struct BoundedScheduleResult {
    std::string status; // verified, search_failed, infeasible_by_lower_bound
    Circuit circuit;
    int depth_lower_bound = 0;
    int requested_depth = 0;
    int explored_states = 0;
};

BoundedScheduleResult synthesize_bounded(
    const Rows& targets, int inputs, const BoundedScheduleOptions& options);
```

新增 `xordepth bounded MATRIX.txt` CLI 命令，支持：

```text
--depth D       必需，最大输出深度
--beam N        Beam 宽度
--pairs N       每轮保留的 pair 候选数
--splits N      每个目标保留的反向拆分数
--restarts N    确定性重启数
--seed N        随机种子
--output FILE   保存完整结果
```

结果保存：目标 masks、请求深度、下界、全部参数、状态、探索状态数、实际 cost/depth/max_fanout、完整 circuit 和验证说明。既有 JSON 结果不重写。

## 正确性与测试

必须覆盖：

1. 单行四输入：`D=2` 得到平衡树，`D=1` 由下界判定不可行；
2. 多输出共享：结果逐门计算与目标 mask 完全一致，所有输出深度不超过 `D`；
3. 已有节点复用：复用不重复增加 XOR 门；
4. 超深候选剪枝：任何中间节点深度超过 `D` 都不会进入最终线路；
5. Beam 去重和固定 seed：相同输入参数生成相同摘要；不同 seed 仍保持验证正确；
6. 搜索预算耗尽：返回 `search_failed`，不伪造 lower bound 或最佳性声明；
7. CLI JSON `verify`：核对 cost、depth、target masks 和状态；
8. `make test` 与现有 root 测试保持通过。

## 不做的事

- 不修改 `third_party/xsets/X-Sets.cpp` 或 `third_party/uber/Uber.cpp`；
- 不引入 Python、浮点评分或新的运行时依赖；
- 不声称启发式结果是全局 XOR 最优；
- 第一版不加入 LUT、布线、扇出时延或多输入 XOR 门模型；
- 第一版不复制高 fanout 节点，避免同时改变深度语义和门数模型。
