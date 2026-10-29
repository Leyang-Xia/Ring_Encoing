# X-Sets 并行候选评估 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (recommended) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 并行执行 X-Sets heavyweight 候选的独立前瞻模拟，降低高参数运行时间，同时保持单线程实现的候选选择和最终线路完全一致。

**Architecture:** 在 `Schedule::GetBest` 的 `MW_SQ`/`MW_SMALLEST_SUM` heavyweight 分支中，先按现有有序集合快照候选，再由固定数量的 worker 为每个候选创建独立 `Schedule` 副本并计算纯评分结果。主线程等待任务完成后按候选原始顺序执行原有评分和 tie-break；其他算法继续走现有路径。线程数通过 `XSETS_THREADS` 环境变量控制，默认使用硬件并发数，`1` 表示严格串行。

**Tech Stack:** C++17 standard library (`std::thread`, `std::mutex`, `std::condition_variable`, `std::exception_ptr`), existing X-Sets data structures, POSIX shell regression tests, Makefile.

**Spec:** `docs/superpowers/specs/2026-10-28-xsets-parallel-candidate-evaluation-design.md`

## Global Constraints

- 保持候选集合、候选遍历顺序、评分函数、平局处理、`maxelts` 截断和最终调度线路不变。
- 不修改 `Schedule::Copy`、`Add_Element`、`Find_Two` 的状态转移语义。
- `UBER_XSET`、`MW`、`MW_SS`、`MW_MATCHING`、`SUBEX` 继续使用原有路径。
- 不引入 Python 或运行时依赖；保持 C++17 和第三方版权声明。
- 异常或任务失败不得返回部分候选；线程必须在所有路径 join 后再释放共享状态。
- 只有经过完整验证的线路才可用于性能和结果结论；超时不作为 XOR 质量结论。

## Review Focus

- `XSETS_THREADS=1` 与未设置环境变量的行为：必须复现原有线路和 XOR 数；由 Task 1/Task 3 的单线程黄金输出测试覆盖。
- 多线程 worker 完成顺序不同：归并必须按候选位置而非完成顺序；由 Task 2 的固定候选顺序测试覆盖。
- 只有一个候选、候选数少于线程数、空评分集合：不得死锁或访问越界；由 Task 2 的边界测试覆盖。
- worker 中 `Schedule` 深拷贝和临时 `Bitset` 的所有权：无泄漏、double free 或共享指针写入；由 Task 2 的 ASan/重复运行验证覆盖。
- 非 heavyweight 方法与 `UBER_XSET`：不得被并行路径改变；由 Task 3 的全套回归测试覆盖。

### Task 1: 定义线程配置与候选评估结果边界

**Files:**
- Modify: `third_party/xsets/X-Sets.cpp:36-80, 1290-1460`
- Test: `tests/xsets_parallel.sh`

**Interfaces:**
- Produces an internal `XSET_THREADS` parser/helper and a per-candidate evaluation record containing the original candidate position and the heavyweight score (`newtop/newn` or `sum`).

- [ ] **Step 1: Write the deterministic regression harness**

  Create `tests/xsets_parallel.sh` that runs the same deterministic matrix through `build/xsets` with `XSETS_THREADS=1` and `XSETS_THREADS=2`, captures stdout/stderr separately, and compares the complete quiet output and XOR count. Include a small matrix with one candidate and a matrix that exercises `MW_SQ` scoring. The serial run is the golden reference; because the requested change is intentionally semantics-preserving, this acceptance test can pass before parallel workers are wired.

- [ ] **Step 2: Record the pre-change golden output**

  Run `make build/xsets` followed by `bash tests/xsets_parallel.sh build/xsets`.

  Expected: the current serial implementation produces the same output for both environment values. Save this as the semantic baseline; later changes must remain byte-identical.

- [ ] **Step 3: Add internal configuration and result types**

  In `X-Sets.cpp`, add a helper that reads `XSETS_THREADS` with `std::getenv`, accepts only a positive integer, clamps it to the candidate count later, and falls back to `std::thread::hardware_concurrency()` (or `1` if unavailable). Add an internal result structure that stores candidate index, candidate key, success/exception state, and the two score fields needed by the existing branches.

- [ ] **Step 4: Run the focused harness**

  Rebuild and run `bash tests/xsets_parallel.sh build/xsets`.

  Expected: the configuration and output comparison portions pass; the script may still report the heavyweight path as serial until Task 2 wires workers.

### Task 2: Implement deterministic worker evaluation and ordered reduction

**Files:**
- Modify: `third_party/xsets/X-Sets.cpp:1290-1460`
- Modify: `Makefile:25-31`
- Test: `tests/xsets_parallel.sh`

**Interfaces:**
- Consumes the Task 1 thread configuration and result record.
- Produces a helper used only by `GetBest`, conceptually `evaluate_candidate(const Schedule&, const std::string&, std::size_t)`; it owns its `Schedule` clone and returns only value-type score data.

- [ ] **Step 1: Implement the serial reference through the same helper**

  Move the existing `candidate -> Copy -> Add_Element -> Find_Two` sequence into the helper. Preserve the existing score calculation exactly, including `Edges_By_Weight.size()-1`, `newn`, and the `sets` traversal for `MW_SMALLEST_SUM`. Use RAII (`std::unique_ptr` or an equivalent local cleanup guard) for `Schedule` and candidate objects.

- [ ] **Step 2: Add a bounded worker pool**

  For more than one candidate and more than one configured thread, create at most `min(thread_count, candidate_count)` workers. Use an atomic next-index counter or mutex-protected index queue; each worker writes only to its own indexed result slot. Capture the first `std::exception_ptr`, let all workers finish, join every thread, then rethrow on the caller thread. Do not update `pots`, `maxtop`, `maxn`, or `minxsz` from workers.

- [ ] **Step 3: Preserve the original reduction order**

  After all results are available, iterate result slots from index `0` upward and apply the exact original comparison and `pots` insertion logic. Keep the existing early return for `top == 1`; leave all non-heavyweight branches untouched. Ensure candidate keys are reconstructed into `Bitset` objects only during the ordered reduction.

- [ ] **Step 4: Add build support for threads**

  Add the platform-appropriate thread linker/compiler flag to the X-Sets and project targets without changing optimization level or C++ standard. Keep `build/xsets` and `build/ringlab` using the same source hash mechanism.

- [ ] **Step 5: Extend the regression harness**

  Compare `XSETS_THREADS=1`, `XSETS_THREADS=2`, and an intentionally oversized value such as `XSETS_THREADS=64` for both `MW_SQ` and `UBER_XSET`. Assert byte-identical quiet output for each method. Add a repeated parallel run to check determinism and a timeout-bound invocation to ensure worker cleanup on interruption/failure.

- [ ] **Step 6: Run focused tests and sanitizers**

  Run:

  ```bash
  make clean
  make build/xsets
  bash tests/xsets_parallel.sh build/xsets
  ```

  If the compiler supports it, rebuild with AddressSanitizer and rerun the script. Expected: all serial/parallel outputs match byte-for-byte and no sanitizer report appears.

### Task 3: Full regression and high-parameter benchmark

**Files:**
- Modify: `docs/xsets-optimization-benchmark.md`
- Modify: `third_party/README.md` (document the new runtime switch and semantic guarantee)

**Interfaces:**
- Consumes the parallel X-Sets binary from Task 2.
- Produces a reproducible benchmark table and regression evidence; no API or result format changes.

- [ ] **Step 1: Run the complete project suite**

  Run `make test` and confirm core tests, CLI tests, quiet-mode tests, exhaustive verification, and existing certificate checks pass.

- [ ] **Step 2: Compare fixed-input semantics**

  On the existing deterministic `data/benchmarks/matrix_100x50_lcg.txt` and the saved 8-to-17 ring example, run serial and parallel `MW_SQ` with identical parameters. Compare verified status, complete circuit, XOR count, and output mapping. Repeat once to check deterministic output.

- [ ] **Step 3: Measure wall-clock and memory behavior**

  Record `XSETS_THREADS=1,2,4` and the default worker count for a completed high-parameter configuration. Record timeout status separately when the configuration does not finish within the budget. Do not compare XOR counts across different methods or unfinished runs as quality claims.

- [ ] **Step 4: Update documentation**

  Add the `XSETS_THREADS` description, default behavior, strict semantic-equivalence guarantee, and benchmark command to `third_party/README.md`. Append the measured serial/parallel timings and unchanged XOR counts to `docs/xsets-optimization-benchmark.md`, preserving prior measurements.

- [ ] **Step 5: Final verification**

  Run `git diff --check`, `make test`, and the parallel harness one final time. Inspect `git diff` to ensure no historical JSON or unrelated backend logic changed.

