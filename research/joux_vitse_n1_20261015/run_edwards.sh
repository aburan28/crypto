#!/bin/sh
# Thin orchestration of PROTOCOL_EDWARDS.md: sh research/joux_vitse_n1_20261015/run_edwards.sh cost|dlp
set -eu
D=research/joux_vitse_n1_20261015/data/edwards
BIN=./target/release/examples/joux_vitse_edwards
mkdir -p "$D"
case "$1" in
cost)
  for p in 271 521; do for s in 1 2; do
    out="$D/cost_p${p}_s${s}.jsonl"; rm -f "$out"
    $BIN cost --p $p --seed $s --residuals 32 --json "$out" > /dev/null
  done; done ;;
dlp)
  for ps in 31:1 41:1 61:1 71:1 31:2 41:2; do
    p=${ps%:*}; s=${ps#*:}
    out="$D/dlp_p${p}_s${s}.jsonl"; rm -f "$out"
    timeout 36000 $BIN dlp --p $p --seed $s --json "$out" > /dev/null || echo "exit $? p=$p s=$s" >> "$D/exits.txt"
  done ;;
esac
