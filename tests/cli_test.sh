#!/usr/bin/env bash
set -euo pipefail
lab="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
expect_fail() { if "$@" >"$tmp/out" 2>"$tmp/err"; then printf 'unexpected success\n' >&2; exit 1; fi; }
"$lab" verify --group practical --samples 4 --output "$tmp/practical.json" > /dev/null
# These statuses must not turn into spurious positive results during migration.
grep -q '"net_saved": null' "$tmp/practical.json"
grep -q '"ring": 124' "$tmp/practical.json"
expect_fail "$lab" verify --samples -1
expect_fail "$lab" profile systematic-rs-7-4 --timeout 0 --output "$tmp/bad.json"
expect_fail "$lab" profile systematic-rs-7-4
expect_fail "$lab" search --family unusual --output "$tmp/bad.json"
"$lab" profile systematic-rs-7-4 --path main --method UBER_XSET --output "$tmp/profile.json" > /dev/null
grep -q '"cost": 112' "$tmp/profile.json"
"$lab" verify --input "$tmp/profile.json" --samples 0 > /dev/null
"$lab" profile original --method MW_SQ --params 10 2 0 --timeout 0.001 --output "$tmp/timeout.json" > /dev/null
grep -q '"status": "timeout"' "$tmp/timeout.json"
printf '#!/bin/sh\nexit 7\n' > "$tmp/fail-backend"
chmod +x "$tmp/fail-backend"
expect_fail "$lab" profile original --backend "$tmp/fail-backend" --output "$tmp/failed.json"
grep -q '"status": "error"' "$tmp/failed.json"
printf '000\n100\n110\n110\n' > "$tmp/matrix.txt"
"$lab" analyze "$tmp/matrix.txt" --method UBER_XSET --output "$tmp/matrix.json" > /dev/null
grep -q '"cost": 1' "$tmp/matrix.json"
printf '01q\n' > "$tmp/invalid.txt"
expect_fail "$lab" analyze "$tmp/invalid.txt"
"$lab" search --limit 1 --timeout 0.001 --output "$tmp/search.json" > /dev/null
grep -q '"matrix_family": "vandermonde"' "$tmp/search.json"
"$lab" verify --input "$tmp/search.json" --samples 4 > /dev/null
printf 'CLI: verification, profiles, aliases, timeouts, backend errors and finite search passed\n'
"$lab" profile --input data/examples/original-8-to-10.json --path main \
  --method UBER_XSET --output "$tmp/custom.json" > /dev/null
"$lab" verify --input "$tmp/custom.json" --samples 0 > /dev/null
"$lab" verify --input "$tmp/matrix.json" --samples 0 > /dev/null
# Cancelling the parent must reap the isolated optimizer, even without timeout.
printf '#!/bin/sh\necho $$ > "$0.pid"\nexec sleep 30\n' > "$tmp/slow-backend"
chmod +x "$tmp/slow-backend"
"$lab" profile original --backend "$tmp/slow-backend" --output "$tmp/cancelled.json" > /dev/null &
parent=$!
for ((i=0;i<200;++i)); do
  if test -s "$tmp/slow-backend.pid"; then break; fi
  sleep 0.01
done
test -s "$tmp/slow-backend.pid"
child=$(cat "$tmp/slow-backend.pid")
kill -TERM "$parent"
set +e
wait "$parent"
status=$?
set -e
test "$status" -eq 143
if kill -0 "$child" 2>/dev/null; then kill -KILL "$child"; exit 1; fi
grep -q '"status": "cancelled"' "$tmp/cancelled.json"
printf 'CLI: custom definitions and cancellation cleanup passed\n'
