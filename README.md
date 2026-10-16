# Ring Encoding：有限域编码的环提升与 XOR 优化

这个仓库研究一个问题：将有限域上的固定编码变换提升到更高维的商环，在环内复用中间结果，再投影回原域，能否减少**完整编码所需的 XOR 次数**？

项目提供 C++17 实验工具、可逐门检查的 XOR 线路证据，以及数学说明。当前重点是未缩放 Vandermonde、系统化 Reed–Solomon 和 Cauchy 编码矩阵。它是研究工具，尚未实现面向实际数据块的高性能纠删码编码器。

数学概念、比特表示、二进制矩阵展开及环提升原理由 [数学基础](docs/mathematical-foundations.md) 介绍。

## 构建与快速开始

需要 C++17 编译器、GNU Make 和 `sha256sum` 或 `shasum`；也可用下文的 CMake 构建。后端进程管理支持 Linux / macOS。构建和运行不需要 Python、Boost、Jerasure 或 GF-Complete。

```bash
make -j3
./build/ringlab --help
./build/ringlab verify --exhaustive
make test
```

默认 `-O3` 构建。完整验证会核对保存的全部线路证据、矩阵恒等式、MDS 子式、CSE 历史种子及每例 1,000 组随机数据；`--exhaustive` 额外穷举四位域小例的 65,536 个输入。验证在 Release 模式中同样有效。

也支持 CMake 3.16 或更新版本：

```bash
cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build-cmake -j3
(cd build-cmake && ctest --output-on-failure)
```

## 常用操作

只验证某组实验：

```bash
./build/ringlab verify --group practical
./build/ringlab verify --group family --samples 1000
```

对已有矩阵重新运行 `MW_SQ` 或 `UBER_XSET`，保存参数、状态和完整线路：

```bash
./build/ringlab profile systematic-rs-7-4 \
  --path direct --method UBER_XSET --params 0 3 0 \
  --timeout 60 --output /tmp/rs-direct.json
./build/ringlab verify --input /tmp/rs-direct.json
```

`--path` 可选 `direct`（原域完整变换）、`main`（环主计算）、`projection`（回到原域的投影及必要后处理）。单阶段 profile 的计数不能直接当作完整环路线的总成本。

继续细化已保存的候选配置，默认只搜索未缩放 Vandermonde：

```bash
./build/ringlab search --limit 3 --timeout 60 \
  --output /tmp/vandermonde-search.json
./build/ringlab verify --input /tmp/vandermonde-search.json
```

分析自己的 01 矩阵，每行是一个输出，每列是一个输入，列从第 0 位开始：

```bash
./build/ringlab analyze matrix.txt --method CSE --seeds 100
./build/ringlab analyze matrix.txt --method MW_SQ --timeout 60
```

[命令与证据格式](docs/usage.md) 说明全部入口、例名、输出及继续搜索的方法。

## 已有结果

以下数字比较固定原始输入到固定原域输出的完整线路，包含进入环和返回原域的实际代价；当前所用输入进入环仅补零，XOR 成本为 0。

| 矩阵实例 | 已完成直接基线 | 完整环线路 | 结论 |
|---|---:|---:|---|
| 原始 8→17、8 数据 / 3 校验 | 212 | 203 | 相对所列基线减少 9 次 |
| 未缩放 Vandermonde、`0x11D`、12 数据 | 338 | 312 | 相对所列基线减少 26 次 |
| 系统化 RS `[7,4,4]` | 126 | 124 | 相对所列基线减少 2 次 |
| 四位域 Cauchy 小例 | 32 | 35 | 当前环线路没有收益 |

这些是有限启发式搜索得到的可实现上界，不能当作最优性证明或吞吐提升。其他 Cauchy 实例保留负结果；大型系统化 RS 的强直接搜索尚未完成，仍标为未决。[文档导航](docs/README.md) 给出全部正例、负例、参数和适用范围。

## 目录

```text
apps/                  统一命令行 ringlab
include/ringlab/        C++ 公共接口
src/                   精确 JSON、比特向量、数学、线路、CSE 与实验实现
tests/                 核心、证据、CLI 和内存错误回归
data/certificates/     原样保留的实验记录与完整 XOR 线路
docs/                  数学基础、验证方案、实验解释与使用说明
third_party/           保留版权声明的 X-Sets、Uber 后端
legacy/                早期演示源码与其依赖、局限说明
```

新功能使用 C++，优先复用共享核心；开发和计数约定见 [AGENTS.md](AGENTS.md)。第三方来源及 Uber 的内存访问修复见 [third_party/README.md](third_party/README.md)。历史演示仍保留，推荐使用统一 CLI 开展后续实验。
