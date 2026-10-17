#!/bin/sh
# Recorded run for item 7. Usage: sh run_all.sh <out_dir> <bits...>
OUT=$1; shift
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
for b in "$@"; do
  for c in 0 1 2; do
    if [ "$b" -le 14 ]; then RT=3; else RT=2; fi
    if [ "$b" -le 13 ]; then ICARMS=ic2lp,basic; else ICARMS=ic2lp; fi
    echo "--bits $b --curve $c --arms $ICARMS --out $OUT/ic_b${b}_c${c}.jsonl"
    echo "--bits $b --curve $c --arms rho --rho-targets $RT --out $OUT/rho_b${b}_c${c}.jsonl"
  done
done | xargs -P ${PAR:-3} -L 1 sh -c 'nice prlimit --as=1000000000 -- python3 -I ghs3.py "$@" 2>>'"$OUT"'/progress.log' _
