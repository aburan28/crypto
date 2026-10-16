#!/bin/sh
# Diagnosis and R1 ratio on TRAINING curves 5-7 (never in a test fit): relation phase only
# (--las none --max-checks 1), shadow primitives. Usage: PAR=2 sh run_diag.sh results/train
OUT=$1; mkdir -p $OUT
{
for b in 12 13 14 15 16 17 18 19 20 21; do for c in 5 6 7; do
  echo "$OUT --bits $b --curve $c --kappa 2.0 --cx 8 --las none --max-checks 1 --no-gp-check --out $OUT/rel_b${b}_c${c}_k2.0.jsonl"
done; done
for b in 16 18; do for k in 1.0 1.19 1.41 1.68 2.38; do for c in 5 6 7; do
  echo "$OUT --bits $b --curve $c --kappa $k --cx 8 --las none --max-checks 1 --no-gp-check --out $OUT/rel_b${b}_c${c}_k${k}.jsonl"
done; done; done
for c in 5 6 7; do
  echo "$OUT --bits 20 --curve $c --kappa 1.19 --cx 8 --las none --max-checks 1 --no-gp-check --out $OUT/rel_b20_c${c}_k1.19.jsonl"
done
} | xargs -P ${PAR:-2} -L 1 sh job.sh
