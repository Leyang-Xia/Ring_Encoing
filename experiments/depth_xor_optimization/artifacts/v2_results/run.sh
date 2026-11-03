#!/usr/bin/env bash
set -euo pipefail
base="$(cd "$(dirname "$0")/.." && pwd)"
repo="${REPO_ROOT:-$base/strict_alignment/repo}"
out=$base/v2_results
for mode in standalone warm; do
 (
 for name in vandermonde rs cauchy gf8; do
  if [[ $name == gf8 ]]; then depths='6 7'; else depths='5 6'; fi
  for depth in $depths; do
   extra=(); if [[ $mode == warm ]]; then extra=(--incumbent "$base/four_matrices/${name}_native_frontier.json"); fi
   stem=$out/${mode}_${name}_d${depth}
   { TIMEFORMAT='%R'; time timeout 180 "$repo/build/ringlab" analyze "$base/four_matrices/$name.txt" --method OUR_DCCSE_V2 --depth "$depth" --seeds 8 --seed 1 "${extra[@]}" --output "$stem.json" > "$stem.stdout" 2> "$stem.stderr"; } 2> "$stem.seconds"
   "$repo/build/ringlab" verify --input "$stem.json" --samples 0 > "$stem.verify.json"
  done
 done
 ) &
done
wait
printf 'DONE\n'
