#!/usr/bin/env bash
set -u
base=/workspace/scratch/14e94f4a5e35
repo=$base/strict_alignment/repo
out=$base/four_matrices
for name in cauchy gf8; do
 if [[ $name == cauchy ]]; then matrix=$repo/xordepth/data/matrices/cauchy_7_4_direct.txt; depths='5 6'; else matrix=$repo/xordepth/data/matrices/gf8_10x10_expanded.txt; depths='6 7'; fi
 cp "$matrix" "$out/$name.txt"
 for depth in $depths; do
  { TIMEFORMAT='%R'; time timeout 180 "$repo/build/ringlab" analyze "$matrix" --method OUR_DCCSE --depth "$depth" --seeds 8 --seed 1 --output "$out/ours_${name}_d${depth}.json" > "$out/ours_${name}_d${depth}.stdout" 2> "$out/ours_${name}_d${depth}.stderr"; } 2> "$out/ours_${name}_d${depth}.seconds"
  printf '%s\n' "$?" > "$out/ours_${name}_d${depth}.exit"
  if [[ -f "$out/ours_${name}_d${depth}.json" ]]; then "$repo/build/ringlab" verify --input "$out/ours_${name}_d${depth}.json" --samples 0 > "$out/ours_${name}_d${depth}.verify.json"; fi
 done
done
cat /sys/fs/cgroup/memory.events > "$out/memory_before.txt"
for name in cauchy gf8; do
 if [[ $name == cauchy ]]; then depths='5 6'; else depths='6 7'; fi
 for depth in $depths; do
  (
  { TIMEFORMAT='%R'; time timeout 190 "$repo/xordepth/build/xordepth" bounded "$out/$name.txt" --depth "$depth" --beam 64 --branch 16 --representations 8 --restarts 8 --seed 1 --timeout 180 --output "$out/beam_${name}_d${depth}.json" > "$out/beam_${name}_d${depth}.stdout" 2> "$out/beam_${name}_d${depth}.stderr"; } 2> "$out/beam_${name}_d${depth}.seconds"
  printf '%s\n' "$?" > "$out/beam_${name}_d${depth}.exit"
  if [[ -f "$out/beam_${name}_d${depth}.json" ]]; then "$repo/xordepth/build/xordepth" verify "$out/beam_${name}_d${depth}.json" > "$out/beam_${name}_d${depth}.verify.json"; fi
  ) &
 done
done
wait
cat /sys/fs/cgroup/memory.events > "$out/memory_after.txt"
printf 'DONE\n'
