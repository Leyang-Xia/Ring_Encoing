# Task 1 report

Implemented the depth-bounded search domain model and pure primitives.

- Added `xd/search.hpp` with representations, Pareto frontier, candidate, option/result API exposure, and helper declarations.
- Added exact integer helpers for balanced arrival depth, XOR term normalization, Pareto dominance/insertion, candidate depth, and canonical representation signatures.
- Added bounded synthesis validation and a compatibility scheduling fallback with explicit lower-bound/search-failed/verified statuses.
- Added focused primitive tests and Makefile target `build/xd-search-test`.

Validation: `cd xordepth && make build/xd-search-test && build/xd-search-test`; `cd xordepth && make test`.

## Review fixes

- Added `candidate_allowed()` as the explicit hard depth filter and covered rejection of `1 + max(2, 2)` at `D=2`.
- Added `lookup_node()` for mask-index reuse and tests for existing/missing masks.
- Removed the interim `synthesize_bounded()` implementation; Task 3 owns its definition and the Task 1 tests no longer link against it.
- Added `validate_target_width()` to reject target masks wider than the input domain before lower-bound computation.
- Made Pareto truncation ordering deterministic by term count, arrival depth, mask, then terms.

Validation commands:
`cd xordepth && make build/xd-search-test && build/xd-search-test` — passed.
`cd xordepth && make test` — passed.
