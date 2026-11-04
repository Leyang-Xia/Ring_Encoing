# 深度约束后的共享门重写

`rewrite_depth_gates()` 把 PR 中最适合吸收到 xordepth 的局部优化引入到独立模块
`src/xd_rewrite.cpp`。它接收一条已经满足最大深度 `D` 的完整线路，尝试删除一个非输出共享门，再一次性重建它的所有消费者。

重建时先从输出反向传播每个节点的最晚允许到达层数 `latest`。对每个候选删除点，算法按原拓扑顺序重建线路：优先复用已有等价 mask；如果原两个操作数已经有别名，则直接重建；否则从当前可用信号中寻找 GF(2) 补集。任何无法在 `latest[node]` 内完成的候选都会被放弃。

候选线路必须同时满足：

- 输出 mask 与目标矩阵完全一致；
- 所有节点依赖和拓扑顺序有效；
- 最大门级深度不超过 `D`；
- live XOR 门数严格少于当前线路。

每轮只接受严格降门的候选，最多执行四轮和 `max_trials` 次删除尝试。默认 bounded 搜索的 `rewrite_trials` 为 0，因此既有行为保持不变；CLI 可用 `--rewrite-trials N` 开启。

## 固定矩阵试验

以下试验使用相同输入、相同深度、相同 bounded 参数、seed=1 和 25 秒预算，只有第二次开启 `--rewrite-trials 400`。完整目标、参数和两条完整 verified circuit 保存在
`data/results/xd_rewrite_compare_20261104.json`。

| 矩阵 | 原始 bounded | 加重写 | 变化 |
|---|---:|---:|---:|
| Vandermonde 8→17 | 221 / 6 | 221 / 6 | 无变化 |
| systematic RS 7→4 | 149 / 6 | 149 / 6 | 无变化 |
| Cauchy 7→4 | 150 / 6 | 149 / 6 | 少 1 门 |
| GF(2⁸) 10×10 | 1219 / 7 | 1218 / 7 | 少 1 门 |

格式为 `(XOR 门数, 实际深度)`。四条重写结果都通过了 `xordepth verify` 的功能、门数和深度检查。这个结果说明 PR 的多消费者事务式重写与 latest-arrival slack 能在当前 xordepth 上带来局部收益，但它没有改变 Vandermonde 和 RS 在这组搜索预算下的结果，也不是全局优化证明。

## 与 `rebalanced()` 的区别

`rebalanced()` 只展开单消费者、非输出节点并重排累加树，门数不会增加但共享决策不变。`rewrite_depth_gates()` 可以删除多消费者共享节点并同时重建所有消费者，因此可能减少门数；它需要目标 mask 和深度上限作为输入，并且只接受严格改进。
