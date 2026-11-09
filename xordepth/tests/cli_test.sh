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
  --representations 8 --restarts 1 --greedy-restarts 2 --seed 0 --timeout 5 \
  --output "$tmp/bounded.json" > /dev/null
grep -q '"status": "verified"' "$tmp/bounded.json"
grep -q '"depth_limit": 2' "$tmp/bounded.json"
grep -q '"greedy_restarts": 2' "$tmp/bounded.json"
grep -q '"greedy_mwsq": false' "$tmp/bounded.json"
"$tool" verify "$tmp/bounded.json" | grep -q '"status": "verified"'

# The MW_SQ greedy tie rescore is recorded in the result and verifies; legacy
# JSON without the field still verifies.
"$tool" bounded "$tmp/four.txt" --depth 2 --beam 8 --branch 8 \
  --representations 8 --restarts 1 --greedy-restarts 2 --seed 0 --timeout 5 \
  --greedy-mwsq --output "$tmp/mwsq.json" > /dev/null
grep -q '"status": "verified"' "$tmp/mwsq.json"
grep -q '"greedy_mwsq": true' "$tmp/mwsq.json"
"$tool" verify "$tmp/mwsq.json" | grep -q '"status": "verified"'
cp "$tmp/mwsq.json" "$tmp/legacy-mwsq.json"
perl -0pi -e 's/\s*"greedy_mwsq":\s*(true|false),//' "$tmp/legacy-mwsq.json"
! grep -q '"greedy_mwsq"' "$tmp/legacy-mwsq.json"
"$tool" verify "$tmp/legacy-mwsq.json" | grep -q '"status": "verified"'
cp "$tmp/mwsq.json" "$tmp/tampered-mwsq.json"
perl -0pi -e 's/"greedy_mwsq":\s*true/"greedy_mwsq": 1/' "$tmp/tampered-mwsq.json"
expect_fail "$tool" verify "$tmp/tampered-mwsq.json"

"$tool" bounded "$tmp/four.txt" --depth 2 --beam 4 --branch 4 --restarts 1 \
  --residual-pairs 2 --residual-targets 1 --timeout 5 --output "$tmp/residual.json" > /dev/null
grep -q '"residual_pairs": 2' "$tmp/residual.json"
grep -q '"residual_targets": 1' "$tmp/residual.json"
grep -q '"toggle_candidates":' "$tmp/residual.json"
grep -q '"bridge_candidates":' "$tmp/residual.json"
grep -q '"diversity_buckets":' "$tmp/residual.json"

"$tool" bounded "$tmp/four.txt" --depth 2 --beam 4 --branch 4 --restarts 1 \
  --timeout 5 --rewrite-trials 40 --output "$tmp/rewrite.json" > /dev/null
grep -q '"rewrite_trials": 40' "$tmp/rewrite.json"
"$tool" verify "$tmp/rewrite.json" | grep -q '"status": "verified"'

"$tool" bounded "$tmp/four.txt" --depth 2 --beam 4 --branch 4 --restarts 1 \
  --no-timeout --output "$tmp/no-timeout.json" > /dev/null
grep -q '"timeout_seconds": 0' "$tmp/no-timeout.json"

"$tool" bounded "$tmp/four.txt" --depth 63 --beam 1 --branch 1 --restarts 1 \
  --greedy-restarts 0 --no-timeout --max-states 18446744073709551615 \
  --output "$tmp/saturation.json" > /dev/null
grep -q '"max_states": 18446744073709551615' "$tmp/saturation.json"
grep -q '"round_limit": 18446744073709551615' "$tmp/saturation.json"
"$tool" verify "$tmp/saturation.json" | grep -q '"status": "verified"'

printf '0111\n1111\n' > "$tmp/state-limited.txt"
"$tool" bounded "$tmp/state-limited.txt" --depth 2 --beam 1 --branch 1 --restarts 1 \
  --greedy-restarts 0 --no-timeout --max-states 1 --output "$tmp/state-limited.json" > /dev/null
grep -q '"status": "timeout"' "$tmp/state-limited.json"
grep -q '"greedy_restarts": 0' "$tmp/state-limited.json"
grep -q '"explored_states": 1' "$tmp/state-limited.json"
grep -q '"round_limit": 14' "$tmp/state-limited.json"
grep -q '"max_rounds": 1' "$tmp/state-limited.json"
"$tool" verify "$tmp/state-limited.json" | grep -q '"result_status": "timeout"'

cp "$tmp/state-limited.json" "$tmp/legacy-bounded.json"
perl -0pi -e 's/\s*"max_states":\s*\d+,//; s/\s*"round_limit":\s*\d+,//; s/\s*"max_rounds":\s*\d+,//' "$tmp/legacy-bounded.json"
"$tool" verify "$tmp/legacy-bounded.json" | grep -q '"status": "verified"'
cp "$tmp/state-limited.json" "$tmp/tampered-explored-states.json"
perl -0pi -e 's/"explored_states":\s*1/"explored_states": 2/' "$tmp/tampered-explored-states.json"
expect_fail "$tool" verify "$tmp/tampered-explored-states.json"
cp "$tmp/state-limited.json" "$tmp/tampered-round-limit.json"
perl -0pi -e 's/"round_limit":\s*14/"round_limit": 15/' "$tmp/tampered-round-limit.json"
expect_fail "$tool" verify "$tmp/tampered-round-limit.json"

cp "$tmp/bounded.json" "$tmp/tampered-greedy-restarts.json"
perl -0pi -e 's/"greedy_restarts":\s*2/"greedy_restarts": -1/' "$tmp/tampered-greedy-restarts.json"
expect_fail "$tool" verify "$tmp/tampered-greedy-restarts.json"
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
expect_fail "$tool" bounded "$tmp/four.txt" --depth 2 --greedy-restarts nope --output "$tmp/bad.json"
expect_fail "$tool" bounded "$tmp/four.txt" --depth 2 --max-states -1 --output "$tmp/bad.json"
expect_fail "$tool" frontier "$tmp/absent.txt" --output "$tmp/bad.json"
expect_fail "$tool" frontier "$tmp/matrix.txt" --also "$tmp/absent.json" --output "$tmp/bad.json"
printf '01x\n' > "$tmp/bad.txt"
expect_fail "$tool" frontier "$tmp/bad.txt" --output "$tmp/bad.json"
printf '010\n0110\n' > "$tmp/ragged.txt"
expect_fail "$tool" frontier "$tmp/ragged.txt" --output "$tmp/bad.json"
expect_fail "$tool" audit
expect_fail "$tool" verify "$tmp/matrix.txt"
expect_fail "$tool" bogus-command

# --- unified incumbent candidate pool: --incumbent and merge ---------------
# A verified bounded result is offered as an external incumbent. The near-zero
# timeout means only the zero-search-time incumbent path can publish.
"$tool" bounded "$tmp/four.txt" --depth 2 --beam 4 --branch 4 --restarts 1 \
  --timeout 0.000000001 --incumbent "$tmp/bounded.json" \
  --output "$tmp/with-inc.json" > /dev/null
grep -q '"status": "verified"' "$tmp/with-inc.json"
grep -q '"incumbents_used": 1' "$tmp/with-inc.json"
grep -q '"best_source": "incumbent"' "$tmp/with-inc.json"
grep -q '"published_candidates":' "$tmp/with-inc.json"
grep -q '"rewrite_trials_used":' "$tmp/with-inc.json"
grep -q '"incumbent_files":' "$tmp/with-inc.json"
"$tool" verify "$tmp/with-inc.json" | grep -q '"status": "verified"'

# --incumbent is repeatable; a duplicate circuit is skipped, not counted twice.
"$tool" bounded "$tmp/four.txt" --depth 2 --beam 4 --branch 4 --restarts 1 \
  --timeout 0.000000001 --incumbent "$tmp/bounded.json" \
  --incumbent "$tmp/bounded.json" --output "$tmp/with-inc2.json" > /dev/null
grep -q '"incumbents_skipped": 1' "$tmp/with-inc2.json"

# merge: zero-search-time recombination of already-verified results.
"$tool" merge "$tmp/bounded.json" "$tmp/with-inc.json" \
  --output "$tmp/merged.json" > /dev/null
grep -q '"status": "verified"' "$tmp/merged.json"
grep -q '"published_candidates":' "$tmp/merged.json"
grep -q '"best_source":' "$tmp/merged.json"
grep -q '"incumbent_files":' "$tmp/merged.json"
"$tool" verify "$tmp/merged.json" | grep -q '"status": "verified"'

# BASE target_masks may be exact uint64 numbers as well as decimal strings.
cp "$tmp/bounded.json" "$tmp/base-num.json"
perl -0pi -e 's/"target_masks":\s*\[\s*"\d+"\s*\]/"target_masks": [15]/' "$tmp/base-num.json"
grep -q '"target_masks": \[15\]' "$tmp/base-num.json"
"$tool" merge "$tmp/base-num.json" "$tmp/bounded.json" \
  --output "$tmp/merged-num.json" > /dev/null
"$tool" verify "$tmp/merged-num.json" | grep -q '"status": "verified"'

# Tampered merge output and a malformed incumbent file are rejected.
cp "$tmp/merged.json" "$tmp/tampered-merge.json"
perl -0pi -e 's/"depth":\s*\d+/"depth": 99/' "$tmp/tampered-merge.json"
expect_fail "$tool" verify "$tmp/tampered-merge.json"
cp "$tmp/merged.json" "$tmp/tampered-rewrite.json"
perl -0pi -e 's/"rewrite_trials_used":\s*0/"rewrite_trials_used": 999999/' "$tmp/tampered-rewrite.json"
expect_fail "$tool" verify "$tmp/tampered-rewrite.json"
printf '{ malformed json' > "$tmp/bad-incumbent.json"
expect_fail "$tool" bounded "$tmp/four.txt" --depth 2 --timeout 0.000000001 \
  --incumbent "$tmp/bad-incumbent.json" --output "$tmp/bad-inc.json"

printf 'CLI: frontier, verify, audit, incumbent, merge and error paths passed\n'
