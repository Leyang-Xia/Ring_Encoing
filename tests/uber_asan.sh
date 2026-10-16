#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
build=${1:-"$root/build-sanitize"}
mkdir -p "$build"
"${CXX:-c++}" -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  "$root/third_party/uber/Uber.cpp" -o "$build/uber-asan"
# The historical implementation owns process-lifetime allocations; this check
# covers invalid access and undefined behavior, rather than exit-time leaks.
ASAN_OPTIONS=detect_leaks=0 "$build/uber-asan" I 3 < "$root/tests/uber-smoke.txt" > "$build/uber-smoke.log"
printf 'Uber: completed-target use-after-free regression passed\n'
