# xordepth 深度硬约束调度器重构设计

## 目标

将 `xordepth` 从当前的 CSE 调度、深度统计和重结合实验程序，重构为一个独立的深度受限 XOR 调度器。给定固定 GF(2) 输入矩阵、深度上限 `D` 和搜索预算，调度器寻找满足

```text
所有输出深度 <= D
```

的线路，并在搜索预算内尽量减少完整固定输入到固定输出变换的二输入 XOR 门数。

深度是硬约束，不转换为软惩罚项。搜索得到的线路必须同时通过：

1. 逐门 mask 语义验证；
2. XOR 门数核对；
3. 输出深度核对；
4. 输出节点可达性和 DAG 拓扑验证。

搜索耗尽只能报告 `search_failed` 或 `timeout`。只有 `D < depth_lower_bound` 时，才可以报告 `infeasible_by_lower_bound`。

## 现状和重构边界

当前 `xordepth` 包含：

- `depth_profile()`：从电路拓扑顺序计算节点深度和输出深度；
- `rebalanced()`：保持功能的累加树重结合；
- `synthesize_schedule()`：按 pair 频次的单路径 CSE 调度；
- `frontier()`、`audit()` 和历史结果验证；
- 只读的矩阵解析和完整 JSON 结果。

当前 `--depth-aware` 只影响等频 pair 的平局选择，最终累加树仍可能超过给定深度。它不再作为新搜索器的实现基础，但保留用于复现实验结果。

第三方 `third_party/xsets/X-Sets.cpp` 和 `third_party/uber/Uber.cpp` 不修改。新调度器只提取 X-Sets 的全局 edge support、MXS/XSet 表示和增量索引思想，不复用第三方可变状态机。

旧的 `frontier`、`audit`、`verify` 命令和历史 JSON 保持兼容。新算法通过独立的 `bounded` 命令和独立结果文件接入。

## 算法概览

算法名称为 `DEPTH_BOUNDED_XSET_BEAM`，分为四层：

1. 深度下界与可行初始解；
2. 带深度的 MXS/Pareto 表示；
3. X-Sets 风格的 pair 索引和候选评分；
4. 带转置表和 Beam 宽度限制的多状态搜索。

### 1. 深度预算和初始解

对每个目标行使用现有的 `depth_lower_bound()` 计算安全下界：权重为 `w` 的输出至少需要 `ceil(log2(w))` 层。

若 `D` 小于下界，直接返回 `infeasible_by_lower_bound`。

否则先尝试：

1. 完全不共享的平衡累加树；
2. 当前 CSE 调度器的结果；
3. 当前 CSE 结果经过 `rebalanced()` 的结果。

所有深度不超过 `D` 的结果都作为 incumbent。若没有初始可行解，搜索仍然可以继续；若搜索结束时没有完整线路，则返回 `search_failed`。

### 2. 搜索状态

每个状态包含：

```text
State:
    nodes                 已物化节点
    node_by_mask          mask 到节点候选的索引
    node_depth[v]         节点深度
    target_frontier[i]    第 i 个输出的局部 Pareto 表
    pair_index[e]         edge e 出现在哪些表示中
    gates                 已生成的 live XOR 门数
    max_depth             已生成节点的最大深度
    consumers[v]          节点消费者数
    history               生成顺序和操作数
```

输入节点深度为 0，常量零节点深度为 0。新节点的深度固定为：

```text
d(a XOR b) = 1 + max(d(a), d(b))
```

新节点只有在 `d(a XOR b) <= D` 时才允许进入状态。

为了保证状态可规范化，节点由 `(mask, depth, producer)` 标识，状态签名按 mask、深度和目标表示排序后生成，不使用临时 node id 作为等价判断依据。

第一版对同一个状态禁止重复 mask 节点。若同一个 mask 已有节点，则优先复用已有节点；后续版本可以加入同一 mask 的深度 Pareto 节点，但这不属于第一版的正确性范围。

### 3. 目标的 MXS/Pareto 表示

每个输出不再只保留一个剩余表达式，而保留若干表示：

```text
Representation:
    terms[]              有序、无重复的已物化节点集合
    mask                 terms 的 XOR
    remaining_terms      terms.size()
    arrival_depth        固定 terms 下完成该输出的最小二叉合并深度
```

`arrival_depth` 使用固定叶子到达深度的最早就绪合并计算。该值是当前表示的安全深度下界：

```text
while terms.size() > 1:
    取两个最小到达深度 a,b
    插入 1 + max(a,b)
```

如果 `arrival_depth > D`，表示直接删除。

对同一目标，若表示 `R1` 满足：

```text
R1.remaining_terms <= R2.remaining_terms
R1.arrival_depth  <= R2.arrival_depth
```

且至少一项严格，则 `R2` 被 `R1` 支配并删除。这样保留的是局部 `(remaining_terms, arrival_depth)` Pareto 集，而不是只保留项数最少的表示。

为避免表示数量失控，每个目标最多保留 `max_representations` 个表示。超过上限时按以下稳定顺序截断：

```text
remaining_terms 小者优先
arrival_depth 小者优先
全局 pair support 大者优先
mask 字典序
```

截断是启发式剪枝，结果不能声称全局最优。

### 4. X-Sets 风格候选索引

对每个表示中的节点对 `(a,b)` 生成 edge：

```text
e = mask[a] XOR mask[b]
```

维护增量索引：

```text
pair_index[e] = {
    (target_id, representation_id), ...
}
support[e] = pair_index[e].size()
```

候选生成只访问包含该 edge 的表示，不重新扫描整张目标矩阵。

候选 `e` 的指标为：

```text
support(e)          被多少目标表示使用
new_depth(e)        新节点深度
terms_after(e)      替换 pair 后的剩余项数
arrival_after(e)    替换后的预计完成深度
affected_targets    会被更新的目标数量
```

候选必须通过：

```text
new_depth(e) <= D
```

并且对每个受影响表示执行：

```text
R' = R - {a,b} + {e}
arrival_depth(R') <= D
```

如果某个受影响目标没有任何合法表示，候选被丢弃。

候选优先级采用 X-Sets/UBER-XSET 的结构：

```text
1. support(e) 最大
2. 所属表示的 remaining_terms 最小
3. arrival_after(e) 最小
4. new_depth(e) 最小
5. edge mask 字典序
```

第一版只使用前向 pair 候选，不做无界的目标反向子集枚举。反向拆分会作为后续独立扩展，因为它需要额外维护任意子集的 mask 索引，容易改变搜索空间和结果解释。

## Beam 搜索和成本目标

单路径贪心会过早固定共享结构，因此主搜索保留一组状态：

```text
beam = {initial_state}

while beam 非空且未超时:
    记录所有已完成状态
    对每个状态生成前 branch_width 个候选
    应用候选并做深度可行性检查
    用转置表按规范化状态去重
    保留 beam_width 个状态
```

完成状态是每个目标至少有一个表示满足：

```text
terms.size() == 0       输出常量零
terms.size() == 1       输出节点别名
```

状态排序使用：

```text
estimated_total_cost
max(arrival_depth)
总剩余项数
-total_support
-completed_targets
规范化状态签名
```

`estimated_total_cost` 只用于排序和启发式剪枝，不作为严格最优性证明。安全但较弱的门数下界为：

```text
gates_so_far + max_i(max(0, terms_i.size() - 1))
```

它考虑共享可能同时服务多个目标，因此不把各目标剩余门数直接相加。

每个 seed 使用确定性随机数改变等价候选顺序。多个 restart 共享当前最佳 incumbent，但不会改变同一 seed 下的候选排序。

## 正确性边界

深度检查必须出现于三个位置：

1. 新节点创建时检查 `1 + max(depth[a], depth[b]) <= D`；
2. 表示更新时检查 `arrival_depth <= D`；
3. 导出 Circuit 后调用 `depth_profile()`，检查每个输出深度和全局深度均不超过 `D`。

最终结果还必须调用现有 `Circuit::verify(targets)`，并核对实际 `cost()` 与 JSON 中的 cost。

搜索状态被剪空不等于不可行。以下状态含义必须区分：

```text
infeasible_by_lower_bound  深度下界直接证明不可行
search_failed              搜索空间被启发式剪枝耗尽
timeout                    搜索预算耗尽
verified                   找到并完成验证
```

## 接口和结果文件

在 `xordepth/include/xd/depth.hpp` 中增加：

```cpp
struct BoundedScheduleOptions {
    int max_depth = 0;
    int beam_width = 64;
    int branch_width = 16;
    int max_representations = 8;
    int restarts = 8;
    std::uint64_t seed = 0;
    double timeout_seconds = 180.0;
};

struct BoundedScheduleResult {
    std::string status;
    Circuit circuit;
    int depth_lower_bound = 0;
    int requested_depth = 0;
    std::size_t explored_states = 0;
    std::size_t pruned_depth = 0;
    std::size_t pruned_bound = 0;
};

BoundedScheduleResult synthesize_bounded(
    const Rows& targets, int inputs, const BoundedScheduleOptions& options);
```

新增命令：

```text
xordepth bounded MATRIX.txt --depth D [options] --output FILE
```

选项：

```text
--depth D       必需，最大输出深度
--beam N        Beam 宽度
--branch N      每个状态扩展的候选数
--representations N
--restarts N
--seed N
--timeout SEC
--output FILE
```

新结果文件必须保存：

- 输入矩阵和 target masks；
- `depth_limit` 和 `depth_lower_bound`；
- 搜索参数、seed 和超时预算；
- status、cost、depth、max_fanout；
- explored/pruned 统计；
- verified circuit 和完整操作列表。

既有历史 JSON 不重写，新实验写入 `xordepth/data/results/` 的独立文件。

## 重构后的文件边界

```text
xordepth/include/xd/depth.hpp          深度指标和公开搜索接口
xordepth/include/xd/search.hpp         搜索状态、表示和候选接口
xordepth/src/xd_depth.cpp              深度指标、下界和已有重结合
xordepth/src/xd_search.cpp             深度受限 X-Sets Beam 搜索
xordepth/src/xd_frontier.cpp           旧 frontier/audit/verify 兼容层
xordepth/apps/xordepth.cpp             bounded 命令和结果序列化
xordepth/tests/xd_search_test.cpp      状态转移、剪枝和确定性测试
```

搜索核心不直接依赖 CLI 或 JSON。所有 mask、线路和验证使用现有 C++17 `ringlab::Mask`、`Circuit` 和 `GreedyCircuit`。

## 测试计划

必须覆盖：

1. 单输出四输入，`D=2` 找到三门平衡树；`D=1` 返回 `infeasible_by_lower_bound`；
2. 两个共享输出，验证共享节点只计一次；
3. 中间候选深度超过 `D` 时不会进入状态；
4. 表示 Pareto 支配删除不删除深度更优表示；
5. 相同 seed、输入和参数得到相同 circuit 摘要；
6. 不同 seed 的结果仍通过逐门 mask 和深度验证；
7. Beam 被剪空时返回 `search_failed`，不报告不可行证明；
8. 超时结果没有 circuit 和 cost；
9. `xordepth bounded ...` 的 JSON 可以被 `verify` 重新核验；
10. 64×24 Vandermonde、systematic RS、Cauchy 和 GF(2^8) 展开矩阵的固定预算实验；
11. `make test` 以及根目录 `make test` 均通过。

## 不做的事

- 不修改 X-Sets/Uber 第三方源代码；
- 不引入 Python 或浮点评分；
- 不把启发式结果称为全局 XOR 最优；
- 不把搜索超时解释为不可行；
- 第一版不做 LUT、布线、寄存器、扇出时序和多输入 XOR 建模；
- 第一版不复制高扇出节点，也不允许同一状态重复 mask 节点；
- 不重写历史 JSON 或改变现有 `frontier`/`audit` 结果的含义。
