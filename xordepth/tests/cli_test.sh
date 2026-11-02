#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
tool="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
expect_fail() { if "$@" >"$tmp/out" 2>"$tmp/err"; then printf 'unexpected success\n' >&2; exit 1; fi; }

printf '0111\n1011\n1100\n' > "$tmp/matrix.txt"
printf '1111\n' > "$tmp/four.txt"
"$tool" bounded "$tmp/four.txt" --depth 2 --beam 8 --branch 8 \
  --representations 8 --restarts 1 --seed 0 --timeout 5 --output "$tmp/bounded.json" > /dev/null
grep -q '"status": "verified"' "$tmp/bounded.json"
grep -q '"depth_limit": 2' "$tmp/bounded.json"
"$tool" verify "$tmp/bounded.json" | grep -q '"status": "verified"'

cat > "$tmp/incumbent.json" <<'EOF'
{"inputs":4,"zero_node":null,"ops":[[4,0,1],[5,2,3],[6,4,5]],"outputs":[6],"cost":3}
EOF
"$tool" bounded "$tmp/four.txt" --depth 2 --beam 1 --branch 1 --restarts 1 \
  --lookahead-width 2 --lookahead-depth 1 --incumbent "$tmp/incumbent.json" \
  --timeout 0.000001 --output "$tmp/referenced.json" > /dev/null
grep -q '"status": "verified"' "$tmp/referenced.json"
grep -q '"reference_cost": 3' "$tmp/referenced.json"
"$tool" verify "$tmp/referenced.json" | grep -q '"status": "verified"'

"$tool" bounded "$tmp/four.txt" --depth 2 --beam 4 --branch 4 --restarts 1 \
  --residual-pairs 2 --residual-targets 1 --timeout 5 --output "$tmp/residual.json" > /dev/null
grep -q '"residual_pairs": 2' "$tmp/residual.json"
grep -q '"residual_targets": 1' "$tmp/residual.json"
grep -q '"toggle_candidates":' "$tmp/residual.json"
cp "$tmp/bounded.json" "$tmp/tampered.json"
perl -0pi -e 's/"depth":\s*2/"depth": 99/' "$tmp/tampered.json"
expect_fail "$tool" verify "$tmp/tampered.json"
cp "$tmp/bounded.json" "$tmp/tampered-bound.json"
perl -0pi -e 's/"depth_lower_bound":\s*2/"depth_lower_bound": 99/' "$tmp/tampered-bound.json"
expect_fail "$tool" verify "$tmp/tampered-bound.json"

"$tool" bounded "$root/data/matrices/gf8_10x10_expanded.txt" --depth 6 --beam 8 --branch 8 \
  --restarts 1 --timeout 0.000001 --output "$tmp/timeout.json" > /dev/null || true
grep -q '"status": "timeout"' "$tmp/timeout.json"
grep -q '"cost": null' "$tmp/timeout.json"
! grep -q '"circuit"' "$tmp/timeout.json"

"$tool" bounded "$tmp/four.txt" --depth 1 --output "$tmp/no.json" > /dev/null
grep -q '"status": "infeasible_by_lower_bound"' "$tmp/no.json"

"$tool" frontier "$tmp/matrix.txt" --share 2,0 --output "$tmp/frontier.json" > /dev/null
grep -q '"pareto_front"' "$tmp/frontier.json"
grep -q '"target_masks"' "$tmp/frontier.json"
grep -q '"depth_lower_bound"' "$tmp/frontier.json"
"$tool" verify "$tmp/frontier.json" > /dev/null
"$tool" verify "$tmp/frontier.json" | grep -q '"candidates_verified"'

# The audit pass reads back its own result without rewriting it.
"$tool" audit "$tmp/frontier.json" --output "$tmp/audit.json" > /dev/null
grep -q '"rebalanced_depth"' "$tmp/audit.json"
before=$(shasum "$tmp/frontier.json" | cut -d' ' -f1)
"$tool" audit "$tmp/frontier.json" > /dev/null
test "$before" = "$(shasum "$tmp/frontier.json" | cut -d' ' -f1)"

# A stored-circuit file may contribute candidates only when it computes the
# same transform; unrelated circuits are skipped instead of failing.
"$tool" frontier "$tmp/matrix.txt" --share 2 --also "$root/data/results/xd_frontier_cauchy_7_4.json" \
  --output "$tmp/also.json" > /dev/null
grep -q '"skipped"' "$tmp/also.json"

expect_fail "$tool" frontier "$tmp/matrix.txt"
expect_fail "$tool" frontier "$tmp/matrix.txt" --share abc --output "$tmp/bad.json"
expect_fail "$tool" frontier "$tmp/matrix.txt" --unknown-flag --output "$tmp/bad.json"
expect_fail "$tool" frontier "$tmp/absent.txt" --output "$tmp/bad.json"
expect_fail "$tool" frontier "$tmp/matrix.txt" --also "$tmp/absent.json" --output "$tmp/bad.json"
printf '01x\n' > "$tmp/bad.txt"
expect_fail "$tool" frontier "$tmp/bad.txt" --output "$tmp/bad.json"
printf '010\n0110\n' > "$tmp/ragged.txt"
expect_fail "$tool" frontier "$tmp/ragged.txt" --output "$tmp/bad.json"
expect_fail "$tool" audit
expect_fail "$tool" verify "$tmp/matrix.txt"
expect_fail "$tool" bogus-command

printf 'CLI: frontier, verify, audit and error paths passed\n'
