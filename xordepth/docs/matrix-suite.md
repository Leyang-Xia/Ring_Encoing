# xordepth 矩阵测试集合

单个 64×24 Vandermonde 结果不能代表深度约束 XOR 搜索的性能。当前测试集合覆盖手工反例、确定性合成矩阵、四类编码结构和多个尺寸。

## 确定性生成矩阵

`apps/xordepth_matrix_gen.cpp` 是独立的 C++17 生成器，不依赖 Python。运行 `make generate` 会在 `data/matrices/generated/` 中重新生成 48 个静态矩阵，并更新两个 manifest。矩阵采用四类结构：

- `dense`：每个位置独立等概率取 0/1，并保证每行非零；
- `sparse`：每行恰好选择 3 个不同输入；
- `blocks`：输出和输入划分为 4 个对应块，每行只在自己的块内取随机位；
- `modular`：以固定模数和由 seed 初始化的行系数生成算术模式，并保证每行非零。

每类都覆盖 `10x10`、`12x20`、`24x40`、`40x64`、`64x96`、`100x100`（名称均为输出数 × 输入数），每个组合使用固定 seed `123` 和 `456`，合计 `4 × 6 × 2 = 48` 个文件。生成器只写入 `generated/` 和 manifest，不改写现有代表矩阵。

## 分层 manifest

- `data/matrices/sample.manifest`：12 个矩阵，包括 8 个覆盖四类结构及小/中/大尺寸的生成矩阵，以及 Cauchy、systematic RS、Vandermonde 和 GF(2⁸) 四个现有代表。默认 `make test` 只运行这一层。
- `data/matrices/full.manifest`：48 个生成矩阵加下表全部 9 个现有矩阵，共 57 个。`make test-full` 先运行默认测试，再运行完整矩阵集合。

manifest 中的路径相对于 manifest 所在目录，允许空行和 `#` 注释。测试拒绝重复项、缺失文件、错误数量、不完整结构/尺寸/seed 组合，以及文件名与实际尺寸不一致的生成矩阵。

现有 9 个代表矩阵保持不变：

| 文件 | 结构 | 输入 | 输出 | 用途 |
|---|---|---:|---:|---|
| `tiny_hand_2x5.txt` | 手工共享反例 | 5 | 2 | 检查多消费者重写 |
| `dense_8x16_seed123.txt` | 确定性稠密 | 16 | 8 | 小规模快速回归 |
| `dense_16x32_seed123.txt` | 确定性稠密 | 32 | 16 | 中小规模扩展 |
| `blocks_32x64_seed123.txt` | 确定性块结构 | 64 | 32 | 结构化共享压力 |
| `cauchy_7_4_direct.txt` | Cauchy 直接变换 | 32 | 24 | 稀疏/中等密度 MDS |
| `systematic_rs_7_4_direct.txt` | 系统化 RS 直接变换 | 32 | 24 | 编码矩阵与共享结构 |
| `vandermonde_8to17_direct.txt` | Vandermonde 直接变换 | 64 | 24 | 当前主深度基准 |
| `vandermonde_8to17_direct_canonical.txt` | 同一矩阵的 canonical 行序 | 64 | 24 | 检查行序对启发式的影响 |
| `gf8_10x10_expanded.txt` | GF(2⁸) 展开矩阵 | 80 | 80 | 更宽、更密的压力样本 |

`tests/matrix_suite_test.cpp` 对每个矩阵执行：解析和非空尺寸检查、生成矩阵尺寸范围及数量检查、行权重深度下界、完全不共享的平衡调度、普通/lookahead completable greedy、逐门功能验证、硬深度检查，以及共享门重写后的门数和深度不回退检查。这个测试验证跨结构和尺寸的正确性；它不把一次运行的门数当作最优性证明。

性能实验应另外固定矩阵、深度、Beam/branch/representation、重启、seed 和是否无限时间，并保存完整线路。小矩阵可以使用 `--no-timeout`，但结果必须记录实际运行时间；超时、手动中断和无结果不能当成较差或较优的门数结论。
