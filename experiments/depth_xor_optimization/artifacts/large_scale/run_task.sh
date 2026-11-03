#!/usr/bin/env bash
set -u
root="$(cd "$(dirname "$0")/.." && pwd)"
stem="$1"; method="$2"
input_dir="${LARGE_INPUT_DIR:-$root/large_scale/inputs}"
matrix="$input_dir/$stem.txt"
meta="$input_dir/$stem.json"
result_dir="${LARGE_RESULTS_DIR:-$root/large_scale/results}"
mkdir -p "$result_dir"
base="$result_dir/${stem}_${method}"
export XSETS_THREADS=1
repo_dir="${LARGE_REPO:-$root/strict_alignment/repo}"
lab="$repo_dir/build/ringlab"
xd="$repo_dir/xordepth/build/xordepth"
D=$(sed -n 's/.*"depth_lower_bound": \([0-9]*\).*/\1/p' "$meta")
limit=10
case "$method" in
 frontier) cmd=("$xd" frontier "$matrix" --share 2,3,0);;
 frontier_basic) cmd=("$xd" frontier "$matrix" --share 0);;
 ours_min30) limit=30; cmd=("$lab" analyze "$matrix" --method OUR_DCCSE --depth "$D" --seeds 1 --seed 0);;
 ours_plus130) limit=30; cmd=("$lab" analyze "$matrix" --method OUR_DCCSE --depth "$((D+1))" --seeds 1 --seed 0);;
 ours_min) cmd=("$lab" analyze "$matrix" --method OUR_DCCSE --depth "$D" --seeds 1 --seed 0);;
 ours_plus1) cmd=("$lab" analyze "$matrix" --method OUR_DCCSE --depth "$((D+1))" --seeds 1 --seed 0);;
 MW_SQ|UBER_XSET) cmd=("$lab" analyze "$matrix" --method "$method" --params 0 3 0 --timeout 10); limit=12;;
 beam) cmd=("$xd" bounded "$matrix" --depth "$((D+1))" --beam 8 --branch 4 --representations 4 --restarts 1 --seed 1 --timeout 10); limit=12;;
 *) exit 2;;
esac
TIMEFORMAT='%3R'
{ time timeout --kill-after=2s "${limit}s" "${cmd[@]}" --output "$base.json" > "$base.stdout.json"; } 2> "$base.time.txt"
status=$?
printf '%s\n' "$status" > "$base.exit"
if [[ "$status" == 0 ]] && grep -q '"status": "verified"' "$base.json"; then
 if [[ "$method" == frontier* || "$method" == beam ]]; then
  "$xd" verify "$base.json" > "$base.verify.json" 2> "$base.verify.stderr"
 else
  "$lab" verify --input "$base.json" --samples 0 > "$base.verify.json" 2> "$base.verify.stderr"
 fi
 printf '%s\n' "$?" > "$base.verify.exit"
fi
printf '%s %s exit=%s\n' "$stem" "$method" "$status"
exit 0
