#!/usr/bin/env bash
set -u
base="$(cd "$(dirname "$0")/.." && pwd)"
repo="${REPO_ROOT:-$base/strict_alignment/repo}"
out=$base/pareto_extension
for method in v1 v2; do
 (
 for name in vandermonde rs cauchy gf8; do
  for offset in 2 3; do
   lower=5; if [[ $name == gf8 ]]; then lower=6; fi
   depth=$((lower+offset));stem=$out/${method}_${name}_d${depth}
   if [[ $method == v1 ]]; then alg=OUR_DCCSE; seeds=16; else alg=OUR_DCCSE_V2;seeds=8;fi
   { TIMEFORMAT='%R'; time timeout 180 "$repo/build/ringlab" analyze "$base/four_matrices/$name.txt" --method "$alg" --depth "$depth" --seeds "$seeds" --seed 1 --output "$stem.json" > "$stem.stdout" 2> "$stem.stderr"; } 2> "$stem.seconds"
   code=$?;printf '%s\n' "$code" > "$stem.exit"
   if [[ $code == 0 ]]; then "$repo/build/ringlab" verify --input "$stem.json" --samples 0 > "$stem.verify.json";printf '%s\n' "$?" > "$stem.verify.exit";fi
  done
 done
 ) &
done
wait
printf 'DONE\n'
