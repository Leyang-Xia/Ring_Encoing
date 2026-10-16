# 文档导航

建议先阅读数学基础，再查看验证框架和具体实验。运行命令统一使用 C++ 工具，见 [使用说明](usage.md)。

## 理论与方法

1. [数学基础](mathematical-foundations.md)：从有限域元素的比特表示讲起，解释乘法矩阵、商环、提升、投影、共享中间结果及完整成本。
2. [可行性验证方案](feasibility-validation-plan.md)：正确性、基线、公平比较、搜索和工程实现的验证层次。

## 实验与证据

| 阅读顺序 | 文档 | 主要内容 | 完整证据 |
|---|---|---|---|
| 1 | [四位域小例](existence-example.md) | 可穷举的构造；更强直接基线后来超过环路线 | [JSON](../data/certificates/ring_lift_example.json)、[X-Sets 比较](../data/certificates/xsets_comparison.json) |
| 2 | [8→17 原始实例](structured-ring-example.md) | 152 次主计算 + 51 次投影，对比 212 次直接基线 | [JSON](../data/certificates/structured_ring_example.json) |
| 3 | [Vandermonde 主实验](vandermonde-validation.md) | 固定未缩放幂矩阵，包含 `0x11D` 域表示 | [JSON](../data/certificates/structured_ring_family.json) |
| 4 | [系统化 RS 与 Cauchy](systematic-rs-cauchy-examples.md) | 实际系统化 RS 正例、Cauchy 负例及大例未决状态 | [JSON](../data/certificates/practical_mds_examples.json) |
| 补充 | [构造家族](structured-ring-family.md) | 缩放扩展、更多参数和更强的域分解基线 | [JSON](../data/certificates/structured_ring_family.json) |

## 实现与维护

- [命令与证据格式](usage.md)：构建、核验、profile、搜索、二进制矩阵分析和消毒器检查。
- [C++ 整理与迁移约定](development/cpp-migration.md)：共享模块与迁移核验原则。
- [迁移计划与完成记录](superpowers/plans/2026-10-16-cpp-migration.md)。
- [第三方后端](../third_party/README.md)、[历史演示](../legacy/README.md)、[开发约定](../AGENTS.md)。
