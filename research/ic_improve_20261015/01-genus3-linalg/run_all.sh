#!/bin/sh
# Recorded run (test curves 0-2), kappa = 2, cx = 8. Usage: PAR=2 sh run_all.sh results/run
OUT=$1; mkdir -p $OUT
BIN=${BIN:-../../../target/release/examples/ic_improve2_01_genus3_linalg}
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
ALL=inc,wied,sgelz,sgewd,bfslz,bfswd,bfssgelz,bfssgewd,lpsgelz,lpsgewd
{
  for c in 0 1 2; do echo "--bits 20 --curve $c --cx 8 --las $ALL --out $OUT/b20_c${c}_main.jsonl"; done
  for c in 0 1 2; do echo "--bits 20 --curve $c --cx 8 --streams aa --las inc,lpsgewd --out $OUT/b20_c${c}_aa.jsonl"; done
  for b in 19 18 17 16 15 14 13 12 11 10; do for c in 0 1 2; do
    echo "--bits $b --curve $c --cx 8 --las $ALL --out $OUT/b${b}_c${c}_main.jsonl"
    echo "--bits $b --curve $c --cx 8 --streams aa --las inc,lpsgewd --out $OUT/b${b}_c${c}_aa.jsonl"
  done; done
} | xargs -P ${PAR:-2} -L 1 sh -c 'nice timeout 3600 prlimit --as=4000000000 -- '"$BIN"' "$@" 2>>'"$OUT"'/progress.log || echo "FAILED rc=$? $*" >>'"$OUT"'/progress.log' _
