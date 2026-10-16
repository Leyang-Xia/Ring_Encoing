# Repository conventions

- All new executable logic, experiments, search tools and tests must use C++17 or later. Do not add Python implementations or a Python runtime dependency. Shell may coordinate builds and CLI tests.
- Put reusable mathematics, bit vectors, circuit construction and verification in `include/ringlab/` and `src/`. Expose experiments through `apps/ringlab.cpp`; avoid duplicate standalone implementations.
- Keep README focused on project purpose, build, examples, result limits and navigation. Put mathematical explanations and worked examples under `docs/`.
- Preserve existing JSON certificates byte-for-byte unless deliberately updating an experiment; do not reserialize historical records for formatting. New measurements belong in separate result files and must retain parameters, method, path, status and full verified circuits.
- Never use floating point or a single 64-bit integer to hold a binary linear-form mask. JSON numeric literals must retain exact values.
- Verify every XOR, output, dependency and live-gate count. Release builds must reject invalid evidence; do not rely on `assert` for correctness checks.
- Count the full fixed-input-to-fixed-output transform, including both mappings and all necessary field postprocessing. Zero padding, wiring and output aliases cost no XOR. Compose and prune circuits before reporting a total.
- Compare against the strongest completed direct baseline, including stored field factorizations. Timeout is inconclusive, execution error is an error, and neither is evidence of a lower bound or positive gain.
- Default candidate search to actual unscaled Vandermonde matrices; identify scaled extensions explicitly. Reconstruct systematic RS as `V * inverse(V_data)` and check its identity data block.
- Keep `third_party/xsets/X-Sets.cpp` unchanged unless explicitly changing the backend; preserve third-party copyright notices. Record patches to other vendored programs.
- `legacy/` is historical reference. New capabilities must go into the supported C++ core and CLI.
- Before delivery run `make test`; for memory/process changes run the sanitizer checks described in `docs/usage.md`. Keep `CMakeLists.txt` and `Makefile` consistent.
