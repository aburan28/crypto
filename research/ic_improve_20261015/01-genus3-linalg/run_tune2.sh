#!/bin/sh
# Tuning stage 2 (training curves 5-7): kappa, with LA arms inc and lpsgewd (cx = 8).
# Usage: PAR=2 sh run_tune2.sh results/tune2
OUT=$1; mkdir -p $OUT
BIN=${BIN:-../../../target/release/examples/ic_improve2_01_genus3_linalg}
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
{
  for k in 1 1.41 2 2.83; do for c in 5 6 7; do
    echo "--bits 20 --curve $c --kappa $k --cx 8 --las inc,lpsgewd --no-gp-check --out $OUT/t2_b20_c${c}_k${k}.jsonl"; done; done
  for k in 0.71 1 1.41 2 2.83; do for c in 5 6 7; do
    echo "--bits 18 --curve $c --kappa $k --cx 8 --las inc,lpsgewd --no-gp-check --out $OUT/t2_b18_c${c}_k${k}.jsonl"; done; done
} | xargs -P ${PAR:-2} -L 1 sh -c 'nice timeout 3600 prlimit --as=4000000000 -- '"$BIN"' "$@" 2>>'"$OUT"'/progress.log || echo "FAILED rc=$? $*" >>'"$OUT"'/progress.log' _
