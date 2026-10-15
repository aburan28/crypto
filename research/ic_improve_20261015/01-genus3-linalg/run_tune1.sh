#!/bin/sh
# Tuning stage 1 (training curves 5-7): LA arm and SGE constant cx at kappa = 2.
# Usage: PAR=3 sh run_tune1.sh results/tune1
OUT=$1; mkdir -p $OUT
BIN=${BIN:-../../../target/release/examples/ic_improve2_01_genus3_linalg}
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
LAS=inc,wied,sgelz,sgewd,bfslz,bfswd,bfssgelz,bfssgewd,lpsgelz,lpsgewd
for b in 18 16 14; do for c in 5 6 7; do for cx in 0 2 8; do
  echo "--bits $b --curve $c --kappa 2 --cx $cx --las $LAS --no-gp-check --out $OUT/t1_b${b}_c${c}_cx${cx}.jsonl"
done; done; done | xargs -P ${PAR:-3} -L 1 sh -c 'nice timeout 3600 prlimit --as=4000000000 -- '"$BIN"' "$@" 2>>'"$OUT"'/progress.log || echo "FAILED rc=$? $*" >>'"$OUT"'/progress.log' _
