#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
if rg -n 'ringlab::|ringlab/(core|json)\.hpp|\.\./src|\.\./include' \
  "$root/include" "$root/src" "$root/apps" "$root/Makefile"; then
  echo "xordepth unexpectedly depends on ringlab" >&2
  exit 1
fi
test -x "$root/build/xordepth"
"$root/build/xordepth" --help >/dev/null
printf 'xordepth core independence passed\n'
