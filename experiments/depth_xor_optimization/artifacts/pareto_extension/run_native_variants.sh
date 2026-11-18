#!/usr/bin/env bash
set -euo pipefail
base="$(cd "$(dirname "$0")/.." && pwd)"
repo="${REPO_ROOT:-$base/strict_alignment/repo}"
for name in vandermonde rs cauchy gf8; do
 "$repo/xordepth/build/xordepth" frontier "$base/four_matrices/$name.txt" --share 2,3,0 --depth-aware --output "$base/pareto_extension/native_depth_aware_$name.json" > "$base/pareto_extension/native_depth_aware_$name.stdout"
 "$repo/xordepth/build/xordepth" verify "$base/pareto_extension/native_depth_aware_$name.json" > "$base/pareto_extension/native_depth_aware_$name.verify.json"
done
