#!/usr/bin/env bash
set -euo pipefail

xsets=${1:?usage: xsets_quiet.sh path-to-xsets}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
cat > "$tmp/matrix.txt" <<'MATRIX'
1111
0011
0101
1001
MATRIX

for method in MW_SQ UBER_XSET; do
  "$xsets" 0 0 0 "$method" --quiet < "$tmp/matrix.txt" > "$tmp/$method.out"
  grep -q 'Total XOR operations:' "$tmp/$method.out"
  if grep -q '^Targets Left' "$tmp/$method.out"; then
    echo "$method emitted progress output in quiet mode" >&2
    exit 1
  fi
done

grep -q 'Total XOR operations: 5' "$tmp/MW_SQ.out"
grep -q 'Total XOR operations: 5' "$tmp/UBER_XSET.out"
printf 'X-Sets: quiet mode preserves verified schedule output\n'
