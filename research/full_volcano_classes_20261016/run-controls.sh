#!/usr/bin/env bash
# Thin orchestration only. The executable performs all mathematics and checks.
set -euo pipefail
bin=${1:?path to isogeny-class executable}
out=${2:?fresh output directory}
study=$(cd -- "$(dirname -- "$0")" && pwd)
mkdir -- "$out"
for name in control-2 control-3 control-mixed; do
  "$bin" plan --curve "$study/inputs/$name.json" > "$out/$name-plan.json"
  "$bin" build --curve "$study/inputs/$name.json" --out "$out/$name.sqlite" --steps all > "$out/$name-build.json"
  "$bin" verify --input "$out/$name.sqlite" > "$out/$name-verify.json"
  "$bin" vertex --input "$out/$name.sqlite" --id 0 > "$out/$name-root.json"
done
