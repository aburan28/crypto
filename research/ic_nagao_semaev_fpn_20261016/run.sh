#!/usr/bin/env bash
# Registered run of PROTOCOL.md. Thin orchestration only: one harness process
# per (cell, kind, target, arm) under a 6 GB address-space cap, so an
# out-of-memory or timeout loses one row and is itself recorded as a row.
#
#   research/ic_nagao_semaev_fpn_20261016/run.sh <out-dir>
set -u
OUT=${1:?out dir}
BIN=${BIN:-target/release/examples/nagao_semaev_fpn}
SEED=20261016
P=65521
BUDGET=1800
CAP=6000000000
mkdir -p "$OUT"
{
  echo "commit $(git rev-parse HEAD)"
  echo "binary_sha256 $(sha256sum "$BIN" | cut -d' ' -f1)"
  echo "host $(uname -n) $(nproc) cores $(free -g | awk '/Mem:/{print $2}') GB"
  echo "start $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$OUT/manifest.txt"

run_one() { # n m kind target arm
  local n=$1 m=$2 kind=$3 t=$4 arm=$5 tmp err rc
  tmp=$(mktemp) err=$(mktemp)
  prlimit --as=$CAP timeout $((BUDGET + 600)) "$BIN" --n "$n" --m "$m" --p $P --seed $SEED \
    --targets $((t + 1)) --only "$t" --kind "$kind" --budget-s $BUDGET --arms "$arm" \
    >"$tmp" 2>"$err"
  rc=$?
  if grep -q '"row":"target"' "$tmp"; then
    grep '"row":"target"' "$tmp" >>"$OUT/rows.jsonl"
  else
    local why=crash
    [ $rc -eq 124 ] && why=wall_timeout
    { [ $rc -eq 137 ] || grep -q "memory allocation" "$err"; } && why=oom
    jq -cn --argjson n "$n" --argjson m "$m" --arg kind "$kind" --argjson t "$t" --arg arm "$arm" \
      --argjson rc $rc --arg why "$why" --arg tail "$(tail -c 300 "$err")" \
      '{row:"failure",n:$n,m:$m,kind:$kind,target:$t,arm:$arm,exit:$rc,reason:$why,stderr_tail:$tail}' \
      >>"$OUT/rows.jsonl"
  fi
  grep '"row":"setup"' "$tmp" | head -1 >>"$OUT/setup.jsonl"
  rm -f "$tmp" "$err"
}

for cell in "2 2 8" "3 2 8" "3 3 8" "4 3 8" "5 4 4" "4 4 4"; do
  set -- $cell
  for kind in planted random; do
    for ((t = 0; t < $3; t++)); do
      for arm in sem rr; do
        run_one "$1" "$2" "$kind" "$t" "$arm"
      done
    done
  done
done
echo "end $(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"$OUT/manifest.txt"
