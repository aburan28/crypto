#!/bin/sh
# Tuning of kappa for ic1lp (q^{5/7}) and ic2lpt (q^{2/3}) on training curves 5,6,7.
# Usage: PAR=2 sh run_tune.sh results/tune
OUT=$1; mkdir -p $OUT
BIN=${BIN:-../../../target/release/examples/ic_improve_04_genus3_2lp}
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
for b in 12 13 14; do for c in 5 6 7; do for k in 0.5 1 2 4; do
  echo "--bits $b --curve $c --arms ic1lp,ic2lpt --kappa1 $k --kappa2 $k --no-gp-check --out $OUT/tune_b${b}_c${c}_k${k}.jsonl"
done; done; done | xargs -P ${PAR:-2} -L 1 sh -c 'nice prlimit --as=3000000000 -- '"$BIN"' "$@" 2>>'"$OUT"'/progress.log' _
