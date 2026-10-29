# Depth-Bounded xordepth Rebuild Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (\`- [ ]\`) syntax for tracking.

**Goal:** Rebuild \`xordepth\` around a deterministic X-Sets-inspired Beam Search that minimizes XOR gates subject to the hard constraint that every output depth is at most \`D\`.

**Architecture:** Keep depth metrics, rebalancing, frontier, audit, and historical JSON verification as compatibility code. Add a separate search core with immutable-style \`SearchState\` transitions, per-target Pareto representations, an incremental pair-support index, deterministic candidate ordering, and Beam/transpose pruning. Expose the core through a new \`bounded\` CLI command and a result verifier that treats timeout/search exhaustion separately from lower-bound infeasibility.

**Tech Stack:** C++17, existing \`ringlab::Mask\`, \`ringlab::Circuit\`, \`ringlab::GreedyCircuit\`, exact JSON, GNU Make, shell CLI tests. No Python or floating-point scoring.

**Spec:** \`docs/superpowers/specs/2026-10-29-depth-bounded-xordepth-rebuild-design.md\`

## Global Constraints

- All new executable logic, experiments, search tools and tests use C++17 or later.
- All binary masks use \`ringlab::Mask\`; no floating point and no single 64-bit mask representation.
- Every returned circuit must pass \`Circuit::verify()\` and \`depth_profile()\` before it is marked \`verified\`.
- \`D\` is a hard output-depth limit; search timeout and Beam exhaustion are not infeasibility proofs.
- Existing \`frontier\`, \`audit\`, \`verify\`, historical JSON, and third-party X-Sets/Uber sources remain compatible and unchanged.
- New measurements are separate JSON files containing inputs, parameters, status, metrics, and the complete verified circuit.
- Before delivery run \`make test\` in \`xordepth\` and the root repository test target.

## Review Focus

- A target can have a larger-term representation with lower arrival depth; the Pareto test must retain it. Test in Task 1 with \`(3,2)\` versus \`(4,1)\`.
- A candidate may improve some representations but make another target infeasible under \`D\`; candidate application must be transactional. Test in Task 2 with a shared pair and \`D=2\`.
- A duplicate mask may already exist at a different depth; first-version state normalization must reuse the existing node and never silently violate \`D\`. Test in Task 1.
- Beam exhaustion is not a lower-bound proof. Test in Task 3 with an intentionally tiny beam and a feasible matrix.
- Timeout/error JSON must not contain a circuit or cost, while verified JSON must contain both and pass re-verification. Test in Task 4.

### Task 1: Add the search domain model and Pareto/depth primitives

**Files:**
- Create: \`xordepth/include/xd/search.hpp\`
- Create: \`xordepth/src/xd_search.cpp\`
- Create: \`xordepth/tests/xd_search_test.cpp\`
- Modify: \`xordepth/include/xd/depth.hpp\` to declare \`BoundedScheduleOptions\`, \`BoundedScheduleResult\`, and \`synthesize_bounded()\`.
- Modify: \`xordepth/Makefile\` to compile \`src/xd_search.cpp\` and \`tests/xd_search_test.cpp\`.

**Interfaces:**
- \`search.hpp\` exposes \`Representation\`, \`ParetoFrontier\`, \`Candidate\`, \`BoundedScheduleOptions\`, \`BoundedScheduleResult\`, and \`synthesize_bounded()\`.
- \`xd_search.cpp\` owns representation normalization, fixed-leaf arrival depth, Pareto dominance, canonical state signatures, and state transition helpers.
- Later tasks consume \`BoundedScheduleResult\` without accessing private state containers.

- [ ] **Step 1: Write failing domain tests**

Add tests with direct constructors/helpers for:

~~~cpp
require(combine_depth({0, 0, 0, 0}) == 2, "balanced four-leaf depth");
require(!dominates({3, 2}, {4, 1}), "depth-faster representation is retained");
require(dominates({3, 2}, {4, 3}), "term-and-depth dominance");
require(normalize_terms({4, 2, 4}) == std::vector<int>({2}), "duplicate terms cancel");
~~~

Also assert that a candidate node with \`1 + max(2, 2)\` is rejected for \`D=2\` and that a duplicate mask lookup returns the existing node index.

- [ ] **Step 2: Run the focused test to verify failure**

Run:

~~~bash
cd xordepth && make build/xd-search-test
~~~

Expected: compilation fails because the new search interfaces and helpers do not exist.

- [ ] **Step 3: Implement the domain model**

Define in \`search.hpp\`:

~~~cpp
struct Representation {
    std::vector<int> terms;
    ringlab::Mask mask;
    int arrival_depth = 0;
};

struct ParetoFrontier {
    std::vector<Representation> entries;
};

struct BoundedScheduleOptions {
    int max_depth = 0;
    int beam_width = 64;
    int branch_width = 16;
    int max_representations = 8;
    int restarts = 8;
    std::uint64_t seed = 0;
    double timeout_seconds = 180.0;
};

struct BoundedScheduleResult {
    std::string status;
    Circuit circuit;
    int depth_lower_bound = 0;
    int requested_depth = 0;
    std::size_t explored_states = 0;
    std::size_t pruned_depth = 0;
    std::size_t pruned_bound = 0;
};

BoundedScheduleResult synthesize_bounded(
    const Rows& targets, int inputs, const BoundedScheduleOptions& options);
~~~

Implement integer-only \`combine_depth\`, sorted XOR term normalization, Pareto dominance, and a canonical representation comparator. Reject negative options and invalid depth limits with \`ringlab::require\`.

- [ ] **Step 4: Run the focused tests to verify they pass**

Run:

~~~bash
cd xordepth && make build/xd-search-test && build/xd-search-test
~~~

Expected: PASS for depth calculation, duplicate cancellation, dominance, and depth-filter tests.

- [ ] **Step 5: Commit the domain model**

~~~bash
git add xordepth/include/xd/search.hpp xordepth/include/xd/depth.hpp \
  xordepth/src/xd_search.cpp xordepth/tests/xd_search_test.cpp xordepth/Makefile
git commit -m "feat: add depth-bounded search domain model"
~~~

### Task 2: Implement transactional X-Sets-style state transitions

**Files:**
- Modify: \`xordepth/include/xd/search.hpp\` with private/internal \`SearchState\`, \`NodeRecord\`, \`PairRef\`, and \`Candidate\` declarations.
- Modify: \`xordepth/src/xd_search.cpp\` with state construction, pair indexing, candidate generation, candidate application, and feasibility checks.
- Modify: \`xordepth/tests/xd_search_test.cpp\` with state transition tests.

**Interfaces:**
- \`SearchState::initial(targets, inputs, options)\` creates input nodes and one target representation per output.
- \`SearchState::candidates()\` returns deterministically ordered candidates after hard depth filtering.
- \`SearchState::apply(const Candidate&)\` returns a new state or an invalid result without mutating the parent.
- \`SearchState::complete()\` and \`SearchState::best_cost_lower_bound()\` are consumed by Task 3.

- [ ] **Step 1: Write failing transition tests**

Cover:

~~~cpp
auto state = SearchState::initial({Mask(15)}, 4, options(2));
auto candidates = state.candidates();
require(contains_pair(candidates, 0, 1), "initial pair candidate");
require(!contains_depth_three_candidate(candidates), "hard depth filter");

auto child = state.apply(find_pair(candidates, 0, 1));
require(child.valid(), "transactional pair application");
require(child.frontier(0).has_representation({node_for(0,1)}), "pair replacement");
~~~

Add a two-output case where applying a shared pair makes one target's only representation exceed \`D\`; assert that \`apply()\` returns invalid and leaves the parent signature unchanged.

- [ ] **Step 2: Run the transition tests to verify failure**

Run:

~~~bash
cd xordepth && make build/xd-search-test && build/xd-search-test
~~~

Expected: FAIL in state construction or candidate application before implementation.

- [ ] **Step 3: Implement state transitions**

Implement:

1. Initial target representations as sorted input node IDs.
2. \`pair_index\` mapping each edge mask to \`(target, representation)\` references.
3. Candidate generation from pair support, with \`new_depth <= max_depth\`.
4. Transactional replacement of \`{a,b}\` by the new node, including duplicate cancellation.
5. Recompute only affected representations, Pareto-prune them, and reject the child if any affected target has no legal representation.
6. A first-version \`node_by_mask\` rule that reuses an existing node and does not create duplicate mask nodes.
7. Canonical signatures that sort node masks/depths and all target representations.

Do not mutate the parent state during \`apply()\`. Use ordinary value containers first; optimize allocations only after the transition tests are stable.

- [ ] **Step 4: Run transition tests and existing xordepth tests**

Run:

~~~bash
cd xordepth && make test
~~~

Expected: all old depth/frontier tests pass, and the new tests verify pair support, Pareto retention, transactional rejection, and canonical state equality.

- [ ] **Step 5: Commit state transitions**

~~~bash
git add xordepth/include/xd/search.hpp xordepth/src/xd_search.cpp \
  xordepth/tests/xd_search_test.cpp
git commit -m "feat: add depth-bounded X-Sets state transitions"
~~~

### Task 3: Add Beam Search, incumbent selection, and circuit export

**Files:**
- Modify: \`xordepth/src/xd_search.cpp\` with Beam expansion, restarts, timeout handling, transposition pruning, and Circuit export.
- Modify: \`xordepth/tests/xd_search_test.cpp\` with end-to-end bounded synthesis tests.
- Modify: \`xordepth/Makefile\` only if a separate benchmark/test binary is needed.

**Interfaces:**
- \`synthesize_bounded(targets, inputs, options)\` returns only the public \`BoundedScheduleResult\`.
- A verified result has \`status == "verified"\`, a valid circuit, and exact cost/depth metrics.
- A timeout or search exhaustion has no circuit payload and no cost claim.

- [ ] **Step 1: Write failing end-to-end tests**

Add tests for:

~~~cpp
auto balanced = synthesize_bounded({Mask(15)}, 4, options_with_depth(2));
require(balanced.status == "verified", "D=2 feasible");
require(balanced.circuit.cost() == 3, "balanced XOR count");
require(depth_profile(balanced.circuit).depth <= 2, "hard depth bound");
balanced.circuit.verify({Mask(15)});

auto impossible = synthesize_bounded({Mask(15)}, 4, options_with_depth(1));
require(impossible.status == "infeasible_by_lower_bound", "lower-bound rejection");

auto exhausted = synthesize_bounded(shared_case, 4, tiny_beam_options);
require(exhausted.status == "search_failed" || exhausted.status == "verified",
        "beam exhaustion is not infeasibility");
~~~

Run the same feasible case twice with the same seed and compare the circuit JSON or its stable hash.

- [ ] **Step 2: Run the end-to-end tests to verify failure**

Run:

~~~bash
cd xordepth && make build/xd-search-test && build/xd-search-test
~~~

Expected: FAIL because \`synthesize_bounded()\` does not yet expand Beam states or export a circuit.

- [ ] **Step 3: Implement the search loop**

Implement:

1. Lower-bound rejection before creating a Beam.
2. Feasible CSE/rebalanced incumbents as optional initial best circuits.
3. \`beam_width\` states per round and \`branch_width\` candidates per state.
4. Transposition table keyed by canonical state signature; retain the state with lower gates, then lower maximum arrival depth.
5. Stable priority tuple: estimated cost, maximum arrival depth, total remaining terms, negative support, completed targets, signature.
6. Restart seeds using \`ringlab::SeedRandom\`, while preserving deterministic ordering for equal candidates.
7. Monotonic timeout checks using \`steady_clock\`.
8. Completion only when every target has a representation of zero or one term.
9. Circuit export from node history, followed by \`Circuit::verify()\` and \`depth_profile()\`.
10. Status separation: \`verified\`, \`search_failed\`, \`timeout\`, \`infeasible_by_lower_bound\`.

The weak cost bound may prune only when it is strictly greater than the incumbent cost. It must never be used to claim optimality.

- [ ] **Step 4: Run focused and existing tests**

Run:

~~~bash
cd xordepth && make test
~~~

Expected: the four-input depth tests pass, repeated seeds produce the same circuit, the \`D=1\` case is rejected by the lower bound, and no verified result exceeds \`D\`.

- [ ] **Step 5: Commit Beam Search**

~~~bash
git add xordepth/src/xd_search.cpp xordepth/tests/xd_search_test.cpp
git commit -m "feat: implement depth-bounded X-Sets Beam Search"
~~~

### Task 4: Add the bounded CLI, JSON result schema, and verification

**Files:**
- Modify: \`xordepth/apps/xordepth.cpp\` with the \`bounded\` command and options.
- Modify: \`xordepth/src/xd_frontier.cpp\` with bounded-result serialization/verification helpers.
- Modify: \`xordepth/include/xd/experiments.hpp\` with public result helper declarations.
- Modify: \`xordepth/tests/cli_test.sh\` with success, timeout, lower-bound, and malformed-result cases.
- Modify: \`xordepth/README.md\` with bounded command usage and status semantics.

**Interfaces:**
- CLI syntax: \`xordepth bounded MATRIX.txt --depth D [--beam N] [--branch N] [--representations N] [--restarts N] [--seed N] [--timeout SEC] --output FILE\`.
- \`xordepth verify FILE.json\` recognizes both old frontier JSON and new bounded JSON.
- \`verify\` requires target masks, declared cost/depth consistency, circuit verification, and \`depth <= depth_limit\` for bounded results.

- [ ] **Step 1: Write failing CLI tests**

Add shell cases:

~~~sh
"$tool" bounded "$tmp/four.txt" --depth 2 --beam 8 --branch 8 \
  --restarts 1 --timeout 5 --output "$tmp/bounded.json"
grep -q '"status":"verified"' "$tmp/bounded.json"
grep -q '"depth_limit":2' "$tmp/bounded.json"
"$tool" verify "$tmp/bounded.json" | grep -q '"status":"verified"'

expect_fail "$tool" bounded "$tmp/four.txt" --depth 1 --output "$tmp/no.json"
~~~

Also test that a timeout result has \`"cost": null\` and no \`"circuit"\` field, and that a tampered depth or cost is rejected by \`verify\`.

- [ ] **Step 2: Run CLI tests to verify failure**

Run:

~~~bash
cd xordepth && make build/xordepth && bash tests/cli_test.sh build/xordepth
~~~

Expected: FAIL because the \`bounded\` command and result schema are not wired in.

- [ ] **Step 3: Implement CLI and JSON**

Parse and validate positive \`--depth\`, \`--beam\`, \`--branch\`, \`--representations\`, \`--restarts\`, and finite positive \`--timeout\`. Serialize complete verified circuits only. For non-verified status, serialize status, parameters, counters, and null cost without a circuit.

Add a bounded-result verifier that reconstructs the matrix from \`target_masks\`, checks \`Circuit::verify()\`, checks declared \`cost\`, checks \`depth_profile()\`, and enforces \`profile.depth <= depth_limit\`.

- [ ] **Step 4: Run CLI tests and full xordepth tests**

Run:

~~~bash
cd xordepth && make test
~~~

Expected: old frontier/audit behavior remains unchanged, bounded verified results round-trip through \`verify\`, and invalid status/cost/depth combinations fail.

- [ ] **Step 5: Commit CLI and result verification**

~~~bash
git add xordepth/apps/xordepth.cpp xordepth/src/xd_frontier.cpp \
  xordepth/include/xd/experiments.hpp xordepth/tests/cli_test.sh xordepth/README.md
git commit -m "feat: expose bounded xordepth search and verification"
~~~

### Task 5: Add fixed-matrix experiments, documentation, and release validation

**Files:**
- Create: \`xordepth/data/results/xd_bounded_vandermonde_8to17_d*.json\` for completed experiments only.
- Create: \`xordepth/docs/depth-bounded-search.md\` with algorithm, invariants, and result limits.
- Modify: \`xordepth/README.md\` with reproduction commands and status interpretation.
- Modify: \`xordepth/Makefile\` only if a fixed benchmark target is useful.

**Interfaces:**
- Experiments invoke the compiled \`bounded\` command and preserve all parameters and complete verified circuits.
- Documentation distinguishes verified upper bounds, lower-bound infeasibility, search exhaustion, and timeout.

- [ ] **Step 1: Add documentation tests/checks**

Run shell checks that every saved verified result contains \`depth_limit\`, \`depth_lower_bound\`, \`cost\`, \`depth\`, \`target_masks\`, and \`circuit\`, then invoke \`verify\` on each saved result.

- [ ] **Step 2: Run fixed matrices**

Run at least:

~~~bash
cd xordepth
./build/xordepth bounded data/matrices/vandermonde_8to17_direct.txt \
  --depth 6 --beam 64 --branch 16 --representations 8 --restarts 8 \
  --seed 1 --timeout 180 --output data/results/xd_bounded_vandermonde_d6.json
./build/xordepth bounded data/matrices/systematic_rs_7_4_direct.txt \
  --depth 6 --beam 64 --branch 16 --representations 8 --restarts 8 \
  --seed 1 --timeout 180 --output data/results/xd_bounded_rs_d6.json
~~~

Do not replace an existing result file; use a new filename for each parameter set.

- [ ] **Step 3: Document the measured limits**

Record matrix shape, \`D\`, Beam parameters, seed, status, cost, depth, wall time, source revision, and whether the result is a verified upper bound or only search exhaustion. Do not compare timeout with a completed XOR count as if it were a result.

- [ ] **Step 4: Run repository validation**

Run:

~~~bash
cd xordepth && make test
cd .. && make test
~~~

If the search core changes memory ownership or process handling, also run the sanitizer commands documented in \`docs/usage.md\` and record their result separately.

- [ ] **Step 5: Commit documentation and measured results**

~~~bash
git add xordepth/docs/depth-bounded-search.md xordepth/README.md \
  xordepth/data/results
git commit -m "docs: record depth-bounded xordepth experiments"
~~~

## Plan Self-Review

- Spec coverage: depth hard constraint is covered by Tasks 1–3; Pareto representations by Tasks 1–2; pair indexing by Task 2; Beam/restarts/transposition by Task 3; CLI/result semantics by Task 4; fixed experiments and reproducibility by Task 5.
- Placeholder scan: the plan contains no TODO, TBD, or unspecified implementation step.
- Type consistency: \`BoundedScheduleOptions\` and \`BoundedScheduleResult\` are introduced in Task 1 and consumed unchanged by Tasks 3–5; \`synthesize_bounded()\` has one signature throughout.
- Review focus coverage: all five review-focus risks have explicit tests in Tasks 1–4.
- Compatibility check: no task edits third-party sources, rewrites historical JSON, or removes existing frontier/audit behavior.
