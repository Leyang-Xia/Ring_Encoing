#!/usr/bin/env bash
set -euo pipefail
folder="$(cd "$(dirname "$0")/.." && pwd)"
repo="$(git -C "$folder" rev-parse --show-toplevel)"
ref="${BASE_REF:-4b6a4e83b00c2cd66cfd9798e41959ca0361a25a}"
git -C "$repo" cat-file -e "$ref^{commit}"
work="$(mktemp -d "${TMPDIR:-/tmp}/xor-depth-build.XXXXXX")"
git -C "$repo" archive "$ref" | tar -x -C "$work"
(cd "$work" && git apply --check "$folder/integration.patch" && git apply "$folder/integration.patch")
make -C "$work" -j3
make -C "$work" test > "$work/root-tests.log" 2>&1
make -C "$work/xordepth" test > "$work/xordepth-tests.log" 2>&1
printf 'Verified independent build: %s\n' "$work"
