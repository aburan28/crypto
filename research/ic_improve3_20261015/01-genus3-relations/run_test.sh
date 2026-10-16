#!/bin/sh
# Test collections (held-out curves). Usage:
#   PAR=2 sh run_test.sh results/run sched "12 ... 21" "0 1 2"
#   PAR=2 sh run_test.sh results/run aa    "12 ... 21" "0 1 2"
#   PAR=2 sh run_test.sh results/run r2    "12 ... 21" "0 1 2"   (needs r2.env: amendment C1)
OUT=$1; STREAM=$2; SIZES=$3; CURVES=$4; mkdir -p $OUT
. ./sched.env
{
for b in $SIZES; do for c in $CURVES; do
  CX=$(cx_of $b)
  case $STREAM in
    sched|aa) K=$(awk -v b=$b -v a=$SCHED_ALPHA -v be=$SCHED_BETA 'BEGIN{printf "%.4f", exp(a + be*(b-20))}') ;;
    r2) . ./r2.env
        K=$(awk -v b=$b -v a=$R2_ALPHA -v be=$R2_BETA 'BEGIN{printf "%.4f", exp(a + be*(b-20))}') ;;
  esac
  echo "$OUT --bits $b --curve $c --kappa $K --cx $CX --streams $STREAM --las lpsgebw2@$CX --out $OUT/b${b}_c${c}_${STREAM}.jsonl"
done; done
} | xargs -P ${PAR:-2} -L 1 sh job.sh
