#!/usr/bin/env bash
set -euo pipefail

xsets=${1:?usage: xsets_parallel.sh path-to-xsets}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/small.txt" <<'MATRIX'
1111
0011
0101
1001
MATRIX

for matrix in "$tmp/small.txt"; do
  label=$(basename "$matrix" .txt)
  for method in MW_SQ UBER_XSET; do
    for threads in 1 2 64; do
      timeout 45 env XSETS_THREADS="$threads" "$xsets" 0 0 0 "$method" --quiet \
        < "$matrix" > "$tmp/${label}_${method}_${threads}.out"
      grep -q 'Total XOR operations:' "$tmp/${label}_${method}_${threads}.out"
    done
    cmp "$tmp/${label}_${method}_1.out" "$tmp/${label}_${method}_2.out"
    cmp "$tmp/${label}_${method}_1.out" "$tmp/${label}_${method}_64.out"
  done
done

printf 'X-Sets: serial and parallel candidate evaluation outputs match\n'
