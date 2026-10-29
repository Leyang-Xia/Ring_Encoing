# xordepth：XOR 门数与关键路径深度

固定 GF(2) 线性变换的两种成本可以同时度量：

- **门数**：线路上实际的二输入 XOR 门个数，决定软件实现的工作量与硬件实现的面积；
- **深度**：从输入到输出的最长 XOR 门链，决定组合逻辑的关键路径。

现有 XOR 最小化工作（Plank 式 CSE、X-Sets/Uber 一类搜索）只优化前者。这个实验把两者放进同一个读数里，回答三个具体问题：

1. 现有线路在 `(门数, 深度)` 平面上落在哪里；
2. 在**不增加任何门**的前提下，深度还能压低多少；
3. 要再往下压，需要付出多少门。

这是一个独立实验：它不参与 `ringlab` 的命令行、证据文件或验证流程。它只复用 `../src` 下的四个核心文件（精确 JSON、比特掩码、线路构造、Plank 式 CSE），不引用 ringlab 的实验层。

## 构建

```bash
cd xordepth
make -j4          # 生成 build/xordepth
make test         # 单元检查 + CLI 检查
```

只依赖 C++17 编译器和 GNU Make。若要独立分发，把 `../include` 与 `../src` 的四个核心文件一起带走即可（`CORE=` / `INCLUDE=` 可覆盖路径）。

## 用法

```bash
# 1. 一个矩阵的计数–深度前沿（--also 把已保存的线路也作为候选）
./build/xordepth frontier data/matrices/vandermonde_8to17_direct.txt \
  --share 2,3,0 --also ../data/certificates/structured_ring_example.json \
  --output /tmp/frontier.json
./build/xordepth verify /tmp/frontier.json

# 2. 任何已保存线路的深度审计（只读，不改文件）
./build/xordepth audit ../data/certificates/structured_ring_example.json \
  --output /tmp/audit.json

# 3. 深度硬约束搜索
./build/xordepth bounded data/matrices/vandermonde_8to17_direct.txt \
  --depth 6 --beam 64 --branch 16 --representations 8 \
  --restarts 8 --seed 1 --timeout 180 --output /tmp/bounded.json
./build/xordepth verify /tmp/bounded.json
```

`MATRIX.txt`：每行一个输出，每列一个输入，只允许 `0`/`1` 与空白。

- `--share LIST`：一对信号至少出现在多少个输出行里才物化成共享门；`0` 表示完全不共享。
- `--depth-aware`：频次相同的候选里优先选最早就绪的一对。
- `--also FILE`：从任意结果/证据 JSON 中挑出**能算出同一个变换**的线路，作为额外候选（不能算出同一变换的会被跳过并计数）。
- `bounded` 的 `--depth D` 是所有输出的硬深度上限；`--beam`、`--branch`、`--representations`、`--restarts`、`--seed` 和 `--timeout` 控制启发式搜索预算。

## 主要结果

`(门数, 深度)`，只列 Pareto 前沿；完整数据在 `data/results/`。

| 矩阵 | 输入×输出 | 前沿 |
|---|---|---|
| Vandermonde 8→17 直接 | 64×24 | (203, 6)、(336, 5) |
| 系统化 RS `[7,4,4]` 直接 | 32×24 | (124, 14)、(126, 11)、(127, 9)、(133, 7)、(140, 6)、(168, 5) |
| Cauchy 7→4 直接 | 32×24 | (150, 6)、(224, 5) |
| GF(2⁸) 10×10 展开 | 80×80 | (1219, 7)、(1830, 6) |

两个与直觉相反的观察：

- **门数最优不等于深度最优，但差距常常可以免费消掉。** 8→17 上 `UBER_XSET` 是 212 门/15 层，`MW_SQ` 是 219 门/12 层；对两者各做一次不改共享决策的重结合，变成 212/8 与 219/6。
- **重结合之后的瓶颈换成共享结构。** 系统化 RS 上最强的 124 门线路只能从 15 层压到 14 层——它的深度卡在多消费者共享节点上，局部重排动不了，必须改共享决策。

对 206 条已保存线路做审计：**178 条深度可降低，总门数还减少 59**。最强的 8→17 完整线路 `203` 门从 10 层降到 6 层，并在此过程中支配了同问题中所有计数不超过 220 门的候选（只有「完全不共享」的 336 门/5 层仍留在前沿上）。

## 目录

```text
include/xd/          深度度量与调度接口
src/xd_depth.cpp     深度、下界、平衡调度、重结合
src/xd_search.cpp    深度硬约束 X-Sets Beam Search
src/xd_frontier.cpp  前沿、审计、证据复核、矩阵解析
apps/xordepth.cpp    命令行
tests/               单元与 CLI 检查
data/matrices/       实验输入（01 矩阵）
data/results/        测量结果（每条候选都带完整线路与目标掩码）
docs/                方法、证明与完整数据
```

方法与完整结果见 [计数–深度前沿](docs/count-depth-frontier.md) 和 [深度硬约束搜索](docs/depth-bounded-search.md)。

## 边界

- 深度是**门级逻辑深度**：二输入 XOR、单位延迟、扇出无代价。真实现件还要 LUT 映射、扇出驱动、布线与流水，本实验的结论只是第一层近似。
- 所有数字都是启发式上界。`depth_lower_bound` 很弱（只看单行权重），`深度 = 下界 + 1` 不代表接近最优。
- 重结合不改变任何输出的功能、不增加门数，因此不会改变原有效果的任何计数结论。
- `cost` 不计输入/输出接线、补零与输出别名。
