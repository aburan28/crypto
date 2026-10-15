#!/bin/bash
# Q1 of PROTOCOL.md: the exact Möbius-class comparison of tower sets.
#   BIN=path/to/examples/pkm_ss_tower research/pkm_supersingular_aux_20261015/run_q1.sh
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:?set BIN to the pkm_ss_tower binary}
echo "$(date -u +%FT%TZ) q1 start" >> "$here/progress.txt"
nice -n 10 prlimit --as=4294967296 -- "$bin" structure --p-min 100 --p-max 3000 --t 2,3,4,5 \
  --ordinary 3 --brute 131 --out "$here/q1/structure.jsonl" 2> "$here/q1/structure.log"
echo "$(date -u +%FT%TZ) q1 end exit $?" >> "$here/progress.txt"
