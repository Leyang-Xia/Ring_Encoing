# C++ Repository Migration Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan inline. Steps use checkbox syntax for tracking.

**Goal:** 将所有 Python 实现迁为共享 C++17 核心，保持研究证据，统一运行和文档入口。

**Architecture:** 精确 JSON 与任意长度比特向量支撑共享数学和线路模块。实验模块读取原始证据；CLI 调用同一套重构、校验、CSE 与后端执行逻辑。

**Tech Stack:** C++17、GNU Make、POSIX 进程接口，无额外运行依赖。

**Spec:** [迁移约定](../../development/cpp-migration.md)。

## Global Constraints

新代码使用 C++；保留所有 JSON 原始字节及数学文档；X-Sets 源码保持原样；所有转换和存活 XOR 成本计入；保留负例及未决状态。

## Review Focus

- 超过 64 位的 JSON 整数、重复/零/单位输出必须保持精确。
- 非法节点、死门、错误输出必须被校验器拒绝，Release 模式也如此。
- 后端失败与超时须区分，超时后结束子进程，保留状态。
- 复合线路按真实存活门计费；保留更强的分解域基线。
- CSE 种子历史可复现；默认搜索保持未缩放 Vandermonde。

## Task 1: Exact core and CSE

Files: `include/ringlab/{json,core}.hpp`, `src/{json,bitvector,field,circuit,cse}.cpp`, `tests/core_test.cpp`, `Makefile`。

Interfaces: `Json::parse/dump`, `Mask::decimal/unit`, `multiply/remainder/expand`，`Circuit::verify`，`GreedyCircuit::synthesize`。

- [x] 先写并运行核心回归：112 位整数往返、投影、错误线路拒绝及历史 CSE 线路/直方图；观察缺失实现导致失败。
- [x] 实现上述接口；运行核心检查和旧程序基准，核对具体线路与历史成本。

## Task 2: Reconstruct, verify and run experiments

Files: `include/ringlab/experiments.hpp`, `src/{matrices,verification,search,backend}.cpp`, `apps/ringlab.cpp`, `tests/cli_test.sh`, `tests/experiment_test.cpp`。

Interfaces: `build_case` 输出域/环/投影矩阵；`verify_all` 输出 JSON 报告；`run_backend` 返回 verified/timeout/cancelled/error；CLI 使用同一接口。

- [x] 先运行 CLI 测试，确认未实现入口失败。
- [x] 重构全部五个旧程序的验证、profile、有限搜索和比较功能；提供任意二进制矩阵分析。
- [x] 逐一比较历史输出，保留大例未决状态；测试超时、后端错误、计数重用与默认矩阵族。

## Task 3: Organization and documentation

Files: `README.md`, `AGENTS.md`, `.gitignore`, `docs/README.md`, `docs/*.md`, `third_party/README.md`, `legacy/README.md`。

- [x] 整理目录，更新所有可执行说明和证据链接；删除五个 Python 文件。
- [x] README 提供快速开始、实际能力、结果边界及文档导航；记录后续 C++ 约定。
- [x] 检查内部链接、无 Python 依赖、证据字节不变；从干净构建运行全套测试。

## Task 4: Review and delivery

- [x] 检查 Release/调试构建及消毒器；运行代表性后端并记录迁移性能。
- [x] 完成整分支审查并修复重要问题，再运行全套验证。
- [x] 提交和推送，核对远端提交及工作区状态。
