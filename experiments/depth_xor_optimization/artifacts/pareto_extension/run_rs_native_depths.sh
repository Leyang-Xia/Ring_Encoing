#!/usr/bin/env bash
set -euo pipefail
base="$(cd "$(dirname "$0")/.." && pwd)"
repo="${REPO_ROOT:-$base/strict_alignment/repo}"
for method in v1 v2; do
 for depth in 9 11 14; do
  if [[ $method == v1 ]];then alg=OUR_DCCSE;seeds=16;else alg=OUR_DCCSE_V2;seeds=8;fi
  stem=$base/pareto_extension/${method}_rs_d${depth}
  { TIMEFORMAT='%R';time timeout 180 "$repo/build/ringlab" analyze "$base/four_matrices/rs.txt" --method "$alg" --depth "$depth" --seeds "$seeds" --seed 1 --output "$stem.json" > "$stem.stdout" 2> "$stem.stderr";} 2> "$stem.seconds"
  printf '0\n' > "$stem.exit"
  "$repo/build/ringlab" verify --input "$stem.json" --samples 0 > "$stem.verify.json"
 done
done
