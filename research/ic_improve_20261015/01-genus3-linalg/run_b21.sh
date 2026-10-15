#!/bin/sh
# Extension (amendment A4): b = 21, main stream, LA arms inc and lpsgewd. Usage: sh run_b21.sh results/run
OUT=$1
for c in 0 1 2; do echo "$OUT - --bits 21 --curve $c --cx 8 --las inc,lpsgewd --out $OUT/b21_c${c}_main.jsonl"; done | xargs -P ${PAR:-1} -L 1 sh job.sh
