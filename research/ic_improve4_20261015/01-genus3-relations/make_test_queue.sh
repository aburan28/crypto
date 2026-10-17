#!/bin/sh
# Test queue (held-out curves 0-4) in the protocol's priority order. Needs d2.env (amendment D2):
#   C3_ARM (LA arm name with cx, e.g. lpsgebw2o70c@16), C4_ARM, C4_ALPHA, C4_BETA, W
. ./sched.env
. ./d2.env
k_c4() { awk -v b=$1 -v a=$C4_ALPHA -v be=$C4_BETA 'BEGIN{printf "%.4f", exp(a + be*(b-20))}'; }
line() { # stream bits curve
  b=$2; c=$3; CX=$(cx_of $b)
  case $1 in
    r2) echo "results/run --bits $b --curve $c --kappa $(k_r2 $b) --cx $CX --streams r2 --las lpsgebw2@$CX,$C3_ARM --out results/run/b${b}_c${c}_r2.jsonl" ;;
    c4) echo "results/run --bits $b --curve $c --kappa $(k_c4 $b) --cx $CX --streams c4 --walks $W --las $C4_ARM --out results/run/b${b}_c${c}_c4.jsonl" ;;
    w)  echo "results/run --bits $b --curve $c --kappa $(k_r2 $b) --cx $CX --streams w --walks $W --las lpsgebw2@$CX --out results/run/b${b}_c${c}_w.jsonl" ;;
    aa) echo "results/run --bits $b --curve $c --kappa $(k_r2 $b) --cx $CX --streams aa --las lpsgebw2@$CX --out results/run/b${b}_c${c}_aa.jsonl" ;;
  esac
}
# big jobs early inside each priority block keeps the 3 slots busy
for s in r2 c4; do
  for c in 0 1 2; do line $s 22 $c; done
  for b in 21 20 19 18 17 16 15 14 13 12; do for c in 0 1 2; do line $s $b $c; done; done
done
for b in 21 20 19 18 17 16 15 14 13 12; do for c in 0 1 2; do line w $b $c; done; done
for b in 21 20 19 18 17 16 15 14 13 12; do for c in 0 1 2; do line aa $b $c; done; done
for s in r2 c4; do for b in 21 20 19 18; do for c in 3 4; do line $s $b $c; done; done; done
for c in 0 1 2; do line w 22 $c; done
