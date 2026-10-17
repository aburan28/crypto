#!/bin/sh
# Test collections. Usage:
#   PAR=1 sh run_test.sh results/run main  "12 13 ... 22" "0 1 2"
#   PAR=1 sh run_test.sh results/run aa    "12 ... 22" "0 1 2"
#   PAR=1 sh run_test.sh results/run sched "12 ... 22" "0 1 2"   (needs sched.env: kappa(b), cx(b), m*)
OUT=$1; STREAM=$2; SIZES=$3; CURVES=$4; mkdir -p $OUT
{
for b in $SIZES; do for c in $CURVES; do
  case $STREAM in
    main) echo "$OUT --bits $b --curve $c --kappa 2 --cx 8 --streams main --las lpsgewd,lpsgebw2@8,lpsgebw3@8 --out $OUT/b${b}_c${c}_main.jsonl" ;;
    aa)   echo "$OUT --bits $b --curve $c --kappa 2 --cx 8 --streams aa --las lpsgewd --out $OUT/b${b}_c${c}_aa.jsonl" ;;
    sched) . ./sched.env
          K=$(awk -v b=$b -v a=$SCHED_ALPHA -v be=$SCHED_BETA 'BEGIN{printf "%.4f", exp(a + be*(b-20))}')
          CX=$(cx_of $b)
          echo "$OUT --bits $b --curve $c --kappa $K --cx $CX --streams sched --las lpsgewd@$CX,lpsgebw${MSTAR}@$CX --out $OUT/b${b}_c${c}_sched.jsonl" ;;
  esac
done; done
} | xargs -P ${PAR:-1} -L 1 sh job.sh
