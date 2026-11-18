#!/usr/bin/env bash
set -euo pipefail
helper_dir="$(cd "$(dirname "$0")" && pwd)"
LARGE_REPO="$(cd "${1:?Usage: bash run_large_benchmark.sh REPO OUTPUT}" && pwd)"
mkdir -p "${2:?Output directory required}"
out="$(cd "$2" && pwd)"
make -C "$LARGE_REPO" -j3 test > "$out/tests.log" 2>&1
make -C "$LARGE_REPO/xordepth" -j3 test > "$out/xordepth-tests.log" 2>&1
LARGE_INPUT_DIR="$out/inputs"; LARGE_RESULTS_DIR="$out/results_controlled"
mkdir -p "$LARGE_INPUT_DIR" "$LARGE_RESULTS_DIR"
export LARGE_REPO LARGE_INPUT_DIR LARGE_RESULTS_DIR XSETS_THREADS=1
: > "$out/tasks.txt"
for size in 128x96 192x128 300x200; do
 rows=${size%x*}; cols=${size#*x}
 for pattern in sparse dense blocks; do
  stem="${size}_${pattern}"
  "$LARGE_REPO/build/ringlab" generate "$LARGE_INPUT_DIR/$stem.txt" --rows "$rows" --columns "$cols" --pattern "$pattern" --seed 123 --output "$LARGE_INPUT_DIR/$stem.json" > /dev/null
  for method in frontier ours_min ours_plus1 MW_SQ UBER_XSET beam; do printf '%s %s\n' "$stem" "$method" >> "$out/tasks.txt"; done
 done
done
xargs -P 4 -n 2 bash "$helper_dir/run_task.sh" < "$out/tasks.txt" > "$out/run_controlled.log" 2>&1
bash "$helper_dir/run_task.sh" 300x200_dense frontier_basic > "$out/fallback.log" 2>&1
if [[ "${RUN_SUPPLEMENT:-1}" == 1 ]]; then
 : > "$out/supplement_tasks.txt"
 for stem in 192x128_dense 300x200_sparse 300x200_dense 300x200_blocks; do
  printf '%s ours_min30\n%s ours_plus130\n' "$stem" "$stem" >> "$out/supplement_tasks.txt"
 done
 xargs -P 4 -n 2 bash "$helper_dir/run_task.sh" < "$out/supplement_tasks.txt" > "$out/run_supplement.log" 2>&1
fi
printf 'Large-matrix native CLI comparison complete: %s\n' "$out"
