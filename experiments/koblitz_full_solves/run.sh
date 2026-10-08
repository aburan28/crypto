#!/bin/bash
# One measured full solve: n a2 a6 l V  (V = mono | g<seed> | <seed>)
# Appends one JSON line to results.jsonl; cells already present are skipped.
set -u
cd "$(dirname "$0")/../.."
read -r n a2 a6 l v <<<"$*"
key="$n $a2 $a6 $l $v"
out=experiments/koblitz_full_solves/results.jsonl
grep -qF "\"key\":\"$key\"" "$out" 2>/dev/null && exit 0
seed=$([ "$v" = mono ] && echo "" || echo "$v")
row=$(KOBLITZ_PILOT_MAX_TRIALS=200000 ./target/release/examples/koblitz_ic_pilot $n $a2 $a6 2 $l 50 full $seed | tail -1)
echo "{\"key\":\"$key\",\"row\":$row}" >> "$out"
