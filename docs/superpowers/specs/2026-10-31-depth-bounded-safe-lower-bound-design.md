# xordepth 安全下界与多矩阵搜索设计

## 目标

改进 `xordepth` 的独立深度受限搜索，使其在固定深度预算 `D` 下能够公平比较自身与 MW_SQ、UBER_XSET 重结合结果。第一阶段只修复搜索状态的成本下界和剪枝语义，不引入 incumbent，也不修改第三方后端。

成功标准：

- 所有最终线路继续满足逐门 GF(2) 验证和 `depth <= D`；
- 分支界限使用的 lower bound 不得高估任何可完成线路的 XOR 成本；
- 在 Cauchy、系统化 RS、Vandermonde 三个直接矩阵上，固定参数的结果可复现；
- 记录修复前后的 explored states、pruned states、cost、depth 和预算；
- 搜索结果仍然只是启发式上界，不声称全局最优。

## 现状问题

当前 `SearchState` 为每个目标保存多个 Pareto 表示。`apply()` 保留旧表示并加入替换后的新表示，因此同一状态中的表示是互相替代的方案，而不是需要同时出现在最终线路中的节点集合。

当前 `best_cost_lower_bound()` 却把所有表示的依赖闭包合并为 `live` 集合，再用这个并集的门数作为已支付成本。由于最终导出只会为每个目标选择一个完成表示，其他表示的生产门可以全部被裁掉；因此该并集成本可能大于最终线路成本。将其用于：

```text
child.lower_bound >= incumbent_cost
```

会错误剪枝仍可能达到更低 XOR 数的状态。转置表也使用同一个不安全下界，因此问题不只影响一个剪枝点。

## 第一阶段设计：安全但可计算的下界

### 表示依赖闭包

对一个表示 `R`，递归收集 `R.terms` 所依赖的已物化 XOR 节点，去重后得到 `closure(R)`。输入和常量零不计门数。闭包计算只读当前状态，不修改节点表。

对每个表示定义一个保守的表示下界：

```text
representation_lb(R) = |closure(R)|
```

第一阶段不把 `terms.size()-1` 直接加到闭包门数上，避免把已经存在、可复用的中间节点重复计数。深度合法性由表示的 `arrival_depth <= D` 单独保证。

### 目标下界

对目标 `i`：

```text
target_lb(i) = min representation_lb(R)
                   over all retained legal representations R
```

全局下界取所有目标下界的最大值：

```text
state_lb = max_i target_lb(i)
```

这是安全下界：最终线路必须为每个目标选择一个表示，其生产闭包至少包含该目标所选表示的闭包；全局线路成本至少不小于任一目标所需的闭包门数。跨目标共享只会降低实际并集成本，因此取最大值不会高估。

### 可选的安全项数下界

在第一阶段验证通过后，可以加入一个独立的安全项数下界，但必须通过小规模穷举确认不高估：

```text
term_lb(i) = min_R max(0, |R.terms| - 1)
```

它只能与 `state_lb` 取 `max`，不能直接相加，因为多个目标的最后合并门可能共享。若无法证明某个增强项的安全性，保留较弱的 `state_lb`。

## 状态和剪枝变化

`SearchState::best_cost_lower_bound()` 改为只计算上述安全表示闭包下界，不再把所有 Pareto 表示的节点并集当作已支付成本。

分支界限继续使用：

```text
if state_lb >= best_verified_cost:
    prune
```

因为目标是寻找严格低于当前最好结果的线路，等于 incumbent 的状态没有必要继续扩展。转置表的比较也统一使用新的 `state_lb`。

结果和历史 JSON 格式不变；只增加实验报告中的搜索参数和计数，不重写已有证书。

## 第二阶段候选生成（本阶段不实现）

如果安全下界修复后仍无法接近 X-Sets 结果，再进行候选生成升级：

1. 对目标表示按“实际闭包增量”排序，而不只按 pair support；
2. 为 Beam 保留不同共享结构的代表状态，避免全局排序后只剩 CSE 型路径；
3. 增加经过逐门验证的二层 XOR 宏候选；
4. 维护增量 pair index，避免每次转移完整重建。

这些改变必须在第一阶段基线稳定后单独评估，避免把剪枝正确性与候选启发式混在一起。

## 测试与实验

### 单元测试

新增测试覆盖：

1. 含两个互相替代 Pareto 表示的状态，确认 lower bound 不超过任一合法完成线路的真实成本；
2. 旧表示和新表示同时存在时，未被最终选择的闭包不会计入强制成本；
3. 转置表使用新下界后仍保持确定性；
4. 深度非法表示仍被硬剪枝；
5. 小型四输入和多输出矩阵的完整线路、门数和深度继续通过验证。

### 多矩阵基准

固定：

```text
MW_SQ/UBER_XSET: params = 0 3 0
xordepth: beam=64, branch=16, representations=8, restarts=4, seed=1
```

矩阵：

- Cauchy 7→4 direct，比较 `D=8` 和 `D=11`；
- systematic RS 7→4 direct，比较 `D=7` 和 `D=11`；
- Vandermonde 8→17 direct，比较 `D=6` 和 `D=8`；
- GF(2^8) 10×10 expanded 作为超时压力测试，不把 timeout 当作质量结论。

每个结果保存：矩阵、方法、参数、深度预算、status、cost、实际深度、墙钟预算、探索状态数、剪枝计数和完整已验证 circuit。MW_SQ/UBER_XSET 的深度统一用 `xordepth rebalanced()` 计算，比较时不传入 xordepth 搜索。

## 不做的事

- 不把 MW_SQ、UBER_XSET 或历史结果作为 xordepth 的 incumbent；
- 不修改 `third_party/xsets/X-Sets.cpp` 或 `third_party/uber/Uber.cpp`；
- 不把 timeout、search exhaustion 或 bounded result 当作最优性证明；
- 不引入 Python、浮点评分或额外运行时依赖；
- 不改变历史 JSON 证书的字节内容。
