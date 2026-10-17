#!/usr/bin/env bash
# Correctness integration on the exact large inputs; explicitly partial graphs.
set -euo pipefail
bin=${1:?path to isogeny-class executable}
out=${2:?fresh output directory}
study=$(cd -- "$(dirname -- "$0")" && pwd)
mkdir -- "$out"
for name in A B C D E; do
  code=0
  "$bin" build --curve "$study/inputs/$name.json" --out "$out/$name.sqlite" --steps 1024 > "$out/$name-build.json" 2> "$out/$name-build.stderr" || code=$?
  [ "$code" = 3 ] || { echo "expected partial exit 3, got $code for $name" >&2; exit 1; }
  printf '%s\n' "$code" > "$out/$name-build.exit"
  code=0
  "$bin" verify --input "$out/$name.sqlite" > "$out/$name-verify.json" 2> "$out/$name-verify.stderr" || code=$?
  [ "$code" = 3 ] || { echo "expected partial replay exit 3, got $code for $name" >&2; exit 1; }
  printf '%s\n' "$code" > "$out/$name-verify.exit"
done
