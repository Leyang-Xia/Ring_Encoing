# X-Sets 优化后基准：100×50 01 矩阵

本记录用于验证第一阶段 X-Sets 性能补丁在较大输入上的运行时间。矩阵固定为 100 行、50 列，每个元素由一个确定性的 LCG 生成：

```text
state_0 = 73129
state_{n+1} = (1009 * state_n + 9176) mod 1000003
A[i,j] = 1 iff state_n < 500002
```

矩阵文件为 [`matrix_100x50_lcg.txt`](../data/benchmarks/matrix_100x50_lcg.txt)，100 行均不重复，平均行权重为 25.23。运行环境为当前工作区，X-Sets 源码 SHA-256 为 `8fe109bc7e02b65a756c5ab1495f28a0cbc9e8fd3d0063120c3a9c772f45ef45`。

## 结果

每次运行均使用 `ringlab analyze`，后端为仓库内 `xsets`，超时上限为 30 秒；完成线路再用 `ringlab verify --samples 0` 做符号核验。

| 方法 | 参数 `thresh nstart give-up-when` | 状态 | 墙钟时间 | XOR |
|---|---|---|---:|---:|
| MW_SQ | `0 0 0` | verified | 23.297 s | 905 |
| UBER_XSET | `0 0 0` | verified | 0.160 s | 979 |
| MW_SQ | `0 3 0` | timeout | 30.002 s | — |
| UBER_XSET | `0 3 0` | timeout | 30.002 s | — |

完整结果分别保存在：

- [`benchmark_100x50_mw_sq_baseline.json`](../data/results/benchmark_100x50_mw_sq_baseline.json)
- [`benchmark_100x50_uber_xset_baseline.json`](../data/results/benchmark_100x50_uber_xset_baseline.json)
- [`benchmark_100x50_mw_sq_strong.json`](../data/results/benchmark_100x50_mw_sq_strong.json)
- [`benchmark_100x50_uber_xset_strong.json`](../data/results/benchmark_100x50_uber_xset_strong.json)

`baseline` 表示 `0 0 0`，不是 XOR 质量的理论基线；`0 3 0` 是更强但更昂贵的前瞻配置。timeout 只说明在本次预算内没有完成，不说明不可行或一定劣于另一方法。

## 与补丁前源码的时间对照

在同一矩阵、同一 `0 0 0` 参数下，用补丁前源码临时编译的可执行文件复测：`MW_SQ` 约 23.93 s、905 XOR；`UBER_XSET` 约 0.17 s、979 XOR。当前补丁版本约为 23.30 s 和 0.16 s。这个规模上，优化主要保持了结果并减少了常数开销；此前 8→17 矩阵的 `0 3 0` 基准中，`MW_SQ` 从约 27.9 s 降到 8.1 s，`UBER_XSET` 从约 4.5 s 降到 1.2 s。不同输入规模的收益不能直接外推。
## 候选评估并行化复测

在同一 100×50 矩阵、`MW_SQ 0 0 0` 和当前源码上，使用 `XSETS_THREADS=2` 的运行约 81.35 秒完成，得到与既有结果一致的 905 次 XOR。当前环境中 `XSETS_THREADS=1` 和 `4` 在 90 秒预算内未完成，因此不把超时视为 XOR 质量结论。小矩阵的串行、2 线程和 64 线程输出逐字节一致；大矩阵的并行线路保留完整输出，后续可在更快机器上补充同一输入的串行完整核对。

线程数由环境变量控制，例如：

```bash
XSETS_THREADS=1 ./build/xsets 0 0 0 MW_SQ --quiet < data/benchmarks/matrix_100x50_lcg.txt
XSETS_THREADS=2 ./build/xsets 0 0 0 MW_SQ --quiet < data/benchmarks/matrix_100x50_lcg.txt
```
