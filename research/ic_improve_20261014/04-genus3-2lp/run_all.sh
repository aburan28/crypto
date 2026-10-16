#!/bin/sh
# Recorded run for item 04 (paired mode). Usage: PAR=3 sh run_all.sh <out_dir>
# Each job is single-threaded, nice, capped at 3 GB and 60 minutes.
OUT=$1; mkdir -p $OUT
BIN=${BIN:-../../../target/release/examples/ic_improve_04_genus3_2lp}
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
K1="--kappa1 2"
{
  # longest first
  for c in 0 1 2; do echo "--bits 18 --curve $c --arms basic --out $OUT/b18_c${c}_basic.jsonl"; done
  for c in 0 1 2; do echo "--bits 20 --curve $c --arms ic1lp $K1 --out $OUT/b20_c${c}_ic1lp.jsonl"; done
  for b in 19 20; do for c in 0 1 2; do
    echo "--bits $b --curve $c --arms ic2lp,aa --out $OUT/b${b}_c${c}_ic2lp.jsonl"
    echo "--bits $b --curve $c --arms rho --rho-targets 1 --out $OUT/b${b}_c${c}_rho.jsonl"
  done; done
  for c in 0 1 2; do echo "--bits 19 --curve $c --arms ic1lp $K1 --out $OUT/b19_c${c}_ic1lp.jsonl"; done
  for b in 16 17 18; do for c in 0 1 2; do
    [ $b -lt 18 ] && echo "--bits $b --curve $c --arms basic --out $OUT/b${b}_c${c}_basic.jsonl"
    echo "--bits $b --curve $c --arms ic2lp,aa,ic1lp $K1 --out $OUT/b${b}_c${c}_lp.jsonl"
    echo "--bits $b --curve $c --arms rho --rho-targets 2 --out $OUT/b${b}_c${c}_rho.jsonl"
  done; done
  for b in 10 11 12 13 14 15; do for c in 0 1 2; do
    echo "--bits $b --curve $c --arms basic,ic2lp,aa,ic1lp,rho --rho-targets 2 $K1 --out $OUT/b${b}_c${c}_all.jsonl"
  done; done
} | xargs -P ${PAR:-3} -L 1 sh -c 'nice timeout 3600 prlimit --as=3000000000 -- '"$BIN"' "$@" 2>>'"$OUT"'/progress.log || echo "FAILED rc=$? $*" >>'"$OUT"'/progress.log' _
