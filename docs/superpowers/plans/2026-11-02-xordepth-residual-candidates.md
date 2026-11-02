# xordepth Residual Candidates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 让独立 xordepth 搜索生成当前 pair 共存规则无法表达的跨表示共享节点，并用结构多样性 Beam 提高在硬深度 D 下接近 MW_SQ/UBER_XSET 重结合结果的概率。

**Architecture:** 在现有 `SearchState` 上增加受限 residual/toggle 转移，不改变安全 lower bound；候选索引从单表示 pair 扩展为有限的全局 pair 和受影响目标表示；Beam 按结构指纹分桶保留多种共享结构。MW_SQ、UBER_XSET 只作为独立外部对照，不进入搜索状态。

**Tech Stack:** C++17、现有 `ringlab::Mask`/`Circuit`、GNU Make、Shell CLI tests；不引入 Python、浮点评分或第三方源码修改。

**Spec:** `docs/superpowers/specs/2026-11-02-xordepth-residual-candidates-design.md`

## Global Constraints

- 所有新逻辑和测试使用 C++17 或更高版本。
- 所有 mask 使用 `ringlab::Mask`，不得使用单个 64 位整数承载任意二进制线性形式。
- 最终 circuit 必须通过 `Circuit::verify()`、`depth_profile()` 和 live gate 检查。
- `D` 是所有中间节点和输出的硬深度限制。
- lower bound 只能用于安全剪枝；启发式 priority 不能用于正确性证明。
- 不修改 `third_party/xsets/X-Sets.cpp` 或 `third_party/uber/Uber.cpp`。
- 不把 MW_SQ、UBER_XSET 或历史线路作为 xordepth candidate 或 incumbent。
- 保留已有 JSON 证书字节内容；新实验写入独立结果文件。
- 交付前运行 `make -C xordepth test` 和根目录 `make test`。

## Review Focus

- **代数等价性：** toggle 必须保持目标 mask；由 Task 1 的直接 mask 测试固定。
- **深度约束：** 新节点和 residual 表示都不能超过 D；由 Task 1 的 D 超限测试固定。
- **重复 mask：** 相同 mask 必须复用已有节点并保持 gate count；由 Task 1 固定。
- **候选爆炸：** 全局 pair 和 target 扩展必须受参数限制；由 Task 2 的计数和 CLI metadata 测试固定。
- **Beam 多样性：** 相同 fingerprint 去重、不同 fingerprint 共存；由 Task 3 的结构测试固定。

---

### Task 1: Residual/toggle 状态转移

**Files:**
- Modify: `xordepth/include/xd/search.hpp`
- Modify: `xordepth/src/xd_search.cpp`
- Test: `xordepth/tests/xd_search_test.cpp`

**Interfaces:**
- Consumes: `SearchState`, `Representation`, `Candidate`, `BoundedScheduleOptions`。
- Produces: 一个公开测试可调用的 `SearchState::apply_toggle(int left, int right, std::size_t target, std::size_t representation)`，以及内部 residual transition helper。

- [ ] **Step 1: Write failing tests for algebraic toggle.**

在 `xd_search_test.cpp` 增加一个四输入目标，例如 `Mask(15)` 的表示 `{0,2,3}`，调用 toggle `a=0,b=1`，检查新表示的 XOR mask 仍为 `Mask(15)`，且 `arrival_depth <= max_depth`。增加 D=1 的失败测试、重复 mask 复用测试和父状态不变测试。

- [ ] **Step 2: Run the focused test and verify it fails for the missing API.**

运行：

```bash
make -C xordepth build/xd-search-test
```

预期：编译失败，提示 `apply_toggle` 未定义；不得修改测试来绕过失败。

- [ ] **Step 3: Implement the minimal toggle transition.**

在 `SearchState` 中实现对一个目标表示的 symmetric difference：

```text
terms' = normalize_terms(terms XOR {left, right, new_node})
```

新节点使用现有 `node_by_mask_` 复用逻辑；至少一个操作数必须出现在表示中；新节点深度和表示 arrival depth 都必须通过 `max_depth`。父状态保持不变，失败返回 `valid_=false` 的 child。

- [ ] **Step 4: Run focused tests and then the xordepth test suite.**

运行：

```bash
make -C xordepth test
```

预期：toggle、深度、mask 复用、父状态不变测试通过，已有搜索测试保持通过。

- [ ] **Step 5: Commit Task 1.**

```bash
git add xordepth/include/xd/search.hpp xordepth/src/xd_search.cpp xordepth/tests/xd_search_test.cpp
git commit -m "feat: add depth-safe residual toggle transitions"
```

### Task 2: 全局 residual 候选索引与评分

**Files:**
- Modify: `xordepth/include/xd/depth.hpp`
- Modify: `xordepth/include/xd/search.hpp`
- Modify: `xordepth/src/xd_search.cpp`
- Modify: `xordepth/apps/xordepth.cpp`
- Modify: `xordepth/src/xd_frontier.cpp`
- Modify: `xordepth/include/xd/experiments.hpp`
- Test: `xordepth/tests/xd_search_test.cpp`
- Test: `xordepth/tests/cli_test.sh`

**Interfaces:**
- Consumes: Task 1 toggle transition and existing pair support index。
- Produces: `BoundedScheduleOptions::residual_pairs`, `residual_targets`, `allow_empty_toggle`；候选统计字段 `toggle_candidates`；有限 residual candidate generator。

- [ ] **Step 1: Add failing option and metadata tests.**

在单元测试中设置 `residual_pairs=2`、`residual_targets=1`，检查候选数量不超过限制；设置 `allow_empty_toggle=false`，检查两个操作数都不在表示中的候选不出现。CLI 测试检查 bounded JSON 包含三个参数和 `toggle_candidates`。

- [ ] **Step 2: Run tests to observe missing fields and behavior.**

运行：

```bash
make -C xordepth test
```

预期：新字段不存在或候选计数断言失败。

- [ ] **Step 3: Implement bounded global pair enumeration.**

从当前节点表生成有限 pair，优先使用现有 `pair_index_` 支持和节点 mask 稳定顺序；保留前 `residual_pairs` 个。对每个 pair 只检查前 `residual_targets` 个受影响表示，按目标、表示、节点编号稳定排序，避免无界笛卡尔积。

- [ ] **Step 4: Implement residual candidate scoring and application.**

为普通 pair 和 toggle candidate 计算整数指标：affected targets、term reduction、depth slack、future support、toggle penalty。只用这些指标排序，不修改安全 lower bound。至少一个操作数存在时调用 Task 1 的 toggle；`allow_empty_toggle` 为真时才允许两个操作数都不存在。

- [ ] **Step 5: Wire CLI and bounded JSON metadata.**

在 `apps/xordepth.cpp` 增加：

```text
--residual-pairs N
--residual-targets N
--allow-empty-toggle
```

在 `bounded_result()` 写入参数和 `toggle_candidates`。`verify_bounded_result()` 检查非负参数和非负计数；历史结果缺失字段时保持兼容。

- [ ] **Step 6: Run focused and complete tests.**

运行：

```bash
make -C xordepth test
```

预期：候选限制、mask、深度、CLI metadata 和旧测试全部通过。

- [ ] **Step 7: Commit Task 2.**

```bash
git add xordepth/include/xd/depth.hpp xordepth/include/xd/search.hpp xordepth/src/xd_search.cpp xordepth/apps/xordepth.cpp xordepth/src/xd_frontier.cpp xordepth/include/xd/experiments.hpp xordepth/tests/xd_search_test.cpp xordepth/tests/cli_test.sh
git commit -m "feat: generate bounded residual candidates"
```

### Task 3: 结构多样性 Beam 与有限 lookahead

**Files:**
- Modify: `xordepth/include/xd/search.hpp`
- Modify: `xordepth/src/xd_search.cpp`
- Modify: `xordepth/include/xd/depth.hpp`
- Modify: `xordepth/src/xd_frontier.cpp`
- Modify: `xordepth/apps/xordepth.cpp`
- Test: `xordepth/tests/xd_search_test.cpp`
- Test: `xordepth/tests/cli_test.sh`
- Modify: `xordepth/docs/depth-bounded-search.md`
- Create: `xordepth/data/results/xd_series_residual_bound_20261102.json`

**Interfaces:**
- Consumes: Task 2 candidate scores and residual counters。
- Produces: structure fingerprint、分桶 Beam、有限 MW_SQ 式 lookahead、三矩阵独立对比结果。

- [ ] **Step 1: Write failing Beam diversity tests.**

构造两个 lower bound 相同但 pair-support histogram 不同的状态；检查 fingerprint 不同且 Beam 保留两个状态。构造相同 fingerprint 的两个状态，检查只保留 priority 更好的状态。增加 lookahead 父状态 signature 不变测试。

- [ ] **Step 2: Run tests to verify the current global Beam fails diversity assertions.**

运行：

```bash
make -C xordepth test
```

预期：当前全局排序会丢掉结构不同状态，或缺少 fingerprint API。

- [ ] **Step 3: Implement deterministic structure fingerprints and bucket quotas.**

按每目标最小 terms、arrival-depth histogram、top support histogram、toggle representation 数和 live closure estimate 生成稳定 fingerprint。每轮 Beam 先按 fingerprint 分桶，每个非空桶保留一个，再按桶内 priority 循环填满 `beam_width`。

- [ ] **Step 4: Implement finite future-frontier lookahead.**

对前 `lookahead_width` 个候选复制状态，最多执行 `lookahead_depth` 次 future_support 最大的合法转移；记录完成目标数、最大 future support 和 residual 宽度。lookahead 只改变候选排序，不改变父状态和安全剪枝。

- [ ] **Step 5: Add metadata, documentation and independent comparison.**

在结果 JSON 中写入 `diversity_buckets`、`toggle_candidates` 和 lookahead 参数；在 `depth-bounded-search.md` 说明候选生成限制。运行三组矩阵：Cauchy、systematic RS、Vandermonde；不传 incumbent，MW_SQ/UBER_XSET 只作为重结合后的外部对照。保存完整 verified circuits 和参数到 `xd_series_residual_bound_20261102.json`。

- [ ] **Step 6: Run all verification commands.**

运行：

```bash
make -C xordepth test
make test
git diff --check
```

预期：所有 C++、CLI、历史证书和 X-Sets 并行测试通过；结果中明确标注 timeout 和 heuristic upper bound。

- [ ] **Step 7: Commit Task 3.**

```bash
git add xordepth/include/xd/search.hpp xordepth/src/xd_search.cpp xordepth/include/xd/depth.hpp xordepth/src/xd_frontier.cpp xordepth/apps/xordepth.cpp xordepth/tests/xd_search_test.cpp xordepth/tests/cli_test.sh xordepth/docs/depth-bounded-search.md xordepth/data/results/xd_series_residual_bound_20261102.json
git commit -m "feat: diversify depth-bounded search beam"
```

## Final Review

Before claiming completion, inspect the final diff for accidental incumbent use, verify that all new candidate paths enforce `depth <= D`, rerun the small residual algebra tests, and confirm the three-matrix JSON contains complete verified circuits rather than only metrics.
