#!/bin/sh
# kappa-tuned arms (test curves 0-2): kappa = KAPPA, LA arms inc and lpsgewd (cx = 8).
# Usage: KAPPA=<kappa*> PAR=2 sh run_tuned.sh results/run
OUT=$1; mkdir -p $OUT
{
  for b in 20 19 18 17 16 15 14 13 12 11 10; do for c in 0 1 2; do
    echo "$OUT - --bits $b --curve $c --cx 8 --kappa $KAPPA --streams tuned --las inc,lpsgewd --out $OUT/b${b}_c${c}_tuned.jsonl"
  done; done
} | xargs -P ${PAR:-2} -L 1 sh job.sh
