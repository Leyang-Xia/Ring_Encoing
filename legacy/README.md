# 历史演示源码

这里保留迁移前的早期程序，便于追溯研究想法；后续实验使用根目录 README 中的 `ringlab`。

| 源码 | 原职责 | 依赖或局限 | 当前入口 |
|---|---|---|---|
| `ring_based_XOR.cpp` | 固定八位域到十维环的系数提升、矩阵输出和 Uber 调用 | Jerasure、GF-Complete；固定参数和相对文件路径 | 共享 `build_case` / `expand` 及 `profile`、`analyze` |
| `matrix_xor_analyzer.cpp` | 二进制矩阵的随机生成与后端比较 | 引用仓库未提供的 `cseopt_run`，不能独立链接 | `ringlab analyze` |
| `rs_xor_analysis.c` | 使用 Jerasure 展开小型校验矩阵 | Jerasure；部分分支使用自定义模式，不保证为真正系统化 RS | 系统化 RS 的 `profile` 和 `verify --group practical` |

这些文件不进入默认构建，也不代表当前完整成本与 MDS 核验标准。原来的运行位置、外部依赖和硬编码后端路径均有历史假设，直接编译它们需要自行适配；保留原始内容用于阅读，不在这里继续开发。

新实现仅依赖 C++17 标准库和 POSIX 接口，矩阵、线路和实验不再依赖上述外部库。用户自定义商环与矩阵可通过 `profile --input` 的定义格式运行，见 [使用说明](../docs/usage.md)。
