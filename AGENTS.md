## Cursor Cloud specific instructions

This is a C++17 / C research project for ring-based XOR optimization in erasure coding. There is no web UI, no package manager, and no service to keep running — it produces CLI binaries.

### System dependencies

`libgf-complete-dev` and `libjerasure-dev` (installed via `apt-get`). These are required before any compilation.

### Build commands

The Jerasure headers ship in `/usr/include/jerasure/`, so every compile that uses Jerasure needs `-I/usr/include/jerasure`. The shared library is `libJerasure` (capital J), so link with `-lJerasure` not `-ljerasure`.

```bash
g++ -std=c++17 -O3 Uber.cpp -o Uber
g++ -std=c++17 -O3 X-Sets.cpp -o X-Sets
g++ -std=c++17 -O3 -I/usr/include/jerasure ring_based_XOR.cpp -lgf_complete -lJerasure -o ring_xor
gcc -O3 -I/usr/include/jerasure rs_xor_analysis.c -lJerasure -lgf_complete -o rs_xor_analysis
```

### Running

- `ring_xor` is the main application. It prompts for an optimization level (1-5); pipe input to avoid blocking: `echo "2" | ./ring_xor`
- `ring_xor` calls `./Uber` at runtime via `popen`, so the `Uber` binary must be in the working directory.
- `rs_xor_analysis` runs a longer sweep across multiple RS configurations and calls both `./Uber` and `./X-Sets`.
- `matrix_xor_analyzer.cpp` references an extern `cseopt_run()` not included in the repo, so it cannot be linked standalone.

### Lint / tests

There is no linter or automated test suite configured. Correctness is validated by compiling without errors and running the binaries to verify output.
