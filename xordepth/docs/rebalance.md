# xordepth rebalance 模块

`rebalanced()` 是 xordepth 的独立后处理模块，实现在
`src/xd_rebalance.cpp`，公开接口仍位于 `include/xd/depth.hpp`。它接收一条已经完成共享选择的 XOR circuit，重新安排可以安全展开的累加树，降低关键路径深度。

除了库调用，CLI 也直接暴露该变换：

```bash
xordepth rebalance (XSETS_OUTPUT.txt|-) --matrix MATRIX.txt --output OUT.json
# 或直接接管道（"-" 从 stdin 读，跳过中间文件）：
./X-Sets 0 0 0 UBER_XSET < MATRIX.txt | xordepth rebalance - --matrix MATRIX.txt --output OUT.json
```

它解析 Plank X-Sets（UBER_XSET/MW_SQ）输出的 `Schedule:` 段
（`src/xd_xsets.cpp` 的 `parse_xsets_schedule()`：`sK = xJ` 别名、
`sK = sA + sB` 门、`bT = sK` 目标绑定，绑定行允许任意顺序，输出按 b 序），
先 `verify()` 逐门校验其确实计算 `--matrix` 矩阵，再 `rebalanced()` 并在
结果 JSON 中给出前后 `(cost, depth, max_output_depth, depth_histogram)`。
直方图按层级统计每个输出的深度分布（KP4 300×200 的 UBER_XSET 000 档
rebalance 后为 5 层 9 个 / 6 层 10 个 / 7 层 143 个 / 8 层 138 个）。
先引用后定义、重复绑定、trailer 门数不符、无法识别的行、目标数与矩阵
不符都会报错。输出线路与 `target_masks` 存在 JSON 中，供 `audit` 等
后续工具读取。实测示例（KP4 300×200，UBER_XSET 000 档）：8713 门 /
深度 34 → rebalance 后 8713 门 / 深度 8。

对导入后的线路再跑 `rewrite_depth_gates()`（100/400 trials）实测无增益：
X-Sets 的共享结构在 CSE 意义上已饱和，所有高频对都已物化为共享门，
重写阶段找不到可删除的共享门。rebalance 不改共享决策，因此导入线路的
门数已经接近该共享方案的 下限；要再降门数需换共享方案本身。

## 算法

模块先逐门重建原线路的 GF(2) mask，并统计每个节点的消费者数量。一个内部节点满足以下条件时可以折叠：

```text
consumers[node] == 1 && node 不是输出
```

对每个输出以及不可折叠的共享节点，算法沿输入锥展开这些单消费者节点，直到遇到 primary input、零节点、多消费者节点或输出。展开后的边界节点按到达深度排序，每次取最早就绪的两个节点做 XOR，形成 Huffman 式累加树。

重结合不会改变共享决策。一个包含 `k` 个边界信号的累加锥在展开前后都需要至多 `k - 1` 个 XOR；相同 mask 的节点复用，重复信号则在 GF(2) 中抵消，因此结果不会增加门数。

## 验证不变量

模块返回结果前会检查：

- 每个输出的 mask 与输入线路完全相同；
- 所有节点引用和拓扑顺序有效；
- 输出线路的 XOR 门数不超过输入线路；
- deadline 在扫描、展开和累加过程中都生效。

因此 rebalance 是可独立测试的功能保持变换，而不是新的全局共享搜索器。它不能拆分多消费者共享节点，也不会主动在门数和深度之间搜索 Pareto 权衡。

## 优化空间

当前实现的主要成本在 `accumulate()` 和输入锥展开：每次合并都会重新排序 `pool`，重复节点通过线性查找完成 GF(2) 抵消。后续可以按以下顺序优化：

1. 用按深度排序的 priority queue 替代每轮全量排序；
2. 用 parity set 或节点计数表替代 `pool` 上的线性查找和删除；
3. 缓存每个节点展开后的不可折叠 frontier，避免重复遍历输入锥；
4. 对高深度、多消费者节点增加有限的复制/拆分候选，探索共享结构本身的深度优化。

前 3 项主要改善运行时间，保持当前门数和功能不变量。第 4 项会改变共享决策，需要额外的门数、深度和完整线路验证，不能直接视为现有 rebalance 的等价优化。

## 结果边界

rebalance 的深度是门级 XOR 深度，未计入 LUT 映射、扇出驱动、布线和寄存器。它适合作为计数线路的低风险深度后处理；如果深度仍受多消费者节点限制，应进入 bounded search 或专门的共享拆分实验。
