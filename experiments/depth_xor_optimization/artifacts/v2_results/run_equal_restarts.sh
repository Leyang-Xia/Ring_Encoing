#!/usr/bin/env bash
set -euo pipefail
base="$(cd "$(dirname "$0")/.." && pwd)"
repo="${REPO_ROOT:-$base/strict_alignment/repo}"
for offset in 0 1; do
 (
 for name in vandermonde rs cauchy gf8; do
  depth=$((5+offset)); if [[ $name == gf8 ]]; then depth=$((6+offset)); fi
  stem=$base/v2_results/v1r16_${name}_d${depth}
  { TIMEFORMAT='%R'; time timeout 180 "$repo/build/ringlab" analyze "$base/four_matrices/$name.txt" --method OUR_DCCSE --depth "$depth" --seeds 16 --seed 1 --output "$stem.json" > "$stem.stdout" 2> "$stem.stderr"; } 2> "$stem.seconds"
  "$repo/build/ringlab" verify --input "$stem.json" --samples 0 > "$stem.verify.json"
 done
 ) &
done
wait
