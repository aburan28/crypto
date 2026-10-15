#!/bin/bash
# Q1 of PROTOCOL.md: the exact Möbius-class comparison of tower sets, over
# p ∈ [100, 3000] in chunks (resumable: a chunk whose end was recorded with
# exit 0 is skipped); q1/structure.jsonl is their concatenation in order.
#   BIN=path/to/examples/pkm_ss_tower research/pkm_supersingular_aux_20261015/run_q1.sh
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:?set BIN to the pkm_ss_tower binary}
mkdir -p "$here/q1/chunks"
for lo in 100 400 700 1000 1300 1600 1900 2200 2500 2800; do
  hi=$((lo + 299)); [ $hi -gt 3000 ] && hi=3000
  name=q1-$lo-$hi
  if grep -q "end $name exit 0" "$here/progress.txt" 2>/dev/null; then continue; fi
  echo "$(date -u +%FT%TZ) start $name" >> "$here/progress.txt"
  nice -n 10 prlimit --as=4294967296 -- "$bin" structure --p-min $lo --p-max $hi --t 2,3,4,5 \
    --ordinary 3 --brute 131 --out "$here/q1/chunks/$name.jsonl" 2> "$here/q1/chunks/$name.log"
  echo "$(date -u +%FT%TZ) end $name exit $?" >> "$here/progress.txt"
done
cat "$here"/q1/chunks/q1-{100,400,700,1000,1300,1600,1900,2200,2500,2800}-*.jsonl > "$here/q1/structure.jsonl"
