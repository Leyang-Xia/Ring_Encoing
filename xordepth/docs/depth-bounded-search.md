# 深度硬约束调度搜索

`DEPTH_BOUNDED_XSET_BEAM` 是一个参考 X-Sets/Uber-XSet 思路的启发式调度器。它在固定 GF(2) 输入矩阵上搜索二输入 XOR 网络，并把最大输出深度 `D` 作为硬约束：任何中间节点和最终输出都不能超过 `D`。

## 搜索状态

每个状态保存已物化节点、节点深度、每个输出的候选表示，以及表示中的 pair 支持索引。每个输出保留局部 Pareto 表，表示的两个坐标是：

```text
(remaining_terms, arrival_depth)
```

如果一个表示在这两个坐标上都不差且至少一项更好，较差表示会被删除。这样可以保留“项数更多但深度更低”的表示。

候选 pair 的优先级依次参考：支持数、替换后的剩余项数、表示完成深度、新节点深度和稳定节点顺序。候选生成、表示更新和最终电路导出都会检查深度。Beam 被剪空或搜索预算耗尽不会被解释成不可行证明。

## 状态含义

| 状态 | 含义 |
|---|---|
| `verified` | 找到完整电路，逐门功能、门数和深度都通过验证 |
| `infeasible_by_lower_bound` | `D` 小于行权重给出的安全深度下界 |
| `search_failed` | 搜索空间被启发式剪枝耗尽，没有找到线路 |
| `timeout` | 达到搜索时间预算，没有找到线路 |

只有 `verified` 结果包含 `cost`、`depth` 和完整 `circuit`。其余状态的 `cost` 为 null，不能用于 XOR 数比较。

## 使用

```bash
cd xordepth
make -j4
./build/xordepth bounded data/matrices/vandermonde_8to17_direct.txt \
  --depth 6 --beam 64 --branch 16 --representations 8 \
  --restarts 8 --seed 1 --timeout 180 --output /tmp/vandermonde-d6.json
./build/xordepth verify /tmp/vandermonde-d6.json
```

`--depth D` 是最大输出深度；`--beam`、`--branch`、`--representations`、`--restarts`、`--seed` 和 `--timeout` 分别控制 Beam 宽度、每状态候选数、每目标 Pareto 表上限、重启次数、候选重排种子和墙钟预算。

## 固定矩阵结果

实验使用源码修订 `04a5114f8541e5472c3201ecbe876f6c16e697c3`，参数均为 `beam=64`、`branch=16`、`representations=8`、`restarts=8`、`seed=1`、`timeout=180s`。

| 矩阵 | `D` | 下界 | 状态 | XOR | 实际深度 | 搜索状态数 |
|---|---:|---:|---|---:|---:|---:|
| Vandermonde 8→17 direct | 6 | 5 | verified | 221 | 6 | 10,897 |
| systematic RS [7,4,4] direct | 6 | 5 | verified | 149 | 6 | 17,761 |

这两个结果是启发式上界，不是全局 XOR 最优性证明。完整目标掩码、搜索参数和电路保存在：

- `data/results/xd_bounded_vandermonde_d6.json`
- `data/results/xd_bounded_rs_d6.json`

可用 `xordepth verify` 重新执行功能、门数和深度检查。

## 正确性边界

- 深度按二输入 XOR 门的拓扑深度计算，输入和常量零深度为 0；
- 补零、接线和输出别名不计 XOR；
- Beam 宽度、表示上限和超时都会改变搜索覆盖范围；
- `search_failed` 和 `timeout` 不说明不存在深度不超过 `D` 的线路；
- 当前版本不复制同一状态中的重复 mask 节点，也不建模 LUT、布线、扇出时延或多输入 XOR 门。
