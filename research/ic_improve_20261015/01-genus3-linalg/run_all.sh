#!/bin/sh
# Recorded run (test curves 0-2), kappa = 2, cx = 8. Usage: PAR=2 sh run_all.sh results/run
# Arms of one collection are independent of which other arms share it (each arm's cost is booked
# at its own successful check), so at b >= 19 the LA arms are split over two jobs to stay
# inside the 60-minute cap on the loaded host. job.sh skips a job whose skip-file is non-empty.
OUT=$1; mkdir -p $OUT
ALL=inc,wied,sgelz,sgewd,bfslz,bfswd,bfssgelz,bfssgewd,lpsgelz,lpsgewd
# at b >= 19 the four arms dominated in tuning stage 1 (sgelz, bfslz, bfswd, bfssgelz) are not run
G1=inc,wied,sgewd; G2=bfssgewd,lpsgelz,lpsgewd
{
  for b in 20 19; do for c in 0 2 1; do
    SK=-; [ $b = 20 ] && [ $c = 1 ] && SK=$OUT/b20_c1_main.jsonl
    echo "$OUT $SK --bits $b --curve $c --cx 8 --las $G1 --out $OUT/b${b}_c${c}_main_g1.jsonl"
    echo "$OUT $SK --bits $b --curve $c --cx 8 --las $G2 --out $OUT/b${b}_c${c}_main_g2.jsonl"
    echo "$OUT - --bits $b --curve $c --cx 8 --streams aa --las inc,lpsgewd --out $OUT/b${b}_c${c}_aa.jsonl"
  done; done
  for b in 18 17 16 15 14 13 12 11 10; do for c in 0 1 2; do
    echo "$OUT - --bits $b --curve $c --cx 8 --las $ALL --out $OUT/b${b}_c${c}_main.jsonl"
    echo "$OUT - --bits $b --curve $c --cx 8 --streams aa --las inc,lpsgewd --out $OUT/b${b}_c${c}_aa.jsonl"
  done; done
} | xargs -P ${PAR:-2} -L 1 sh job.sh
