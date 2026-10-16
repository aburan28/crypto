#!/bin/sh
# Training (curves 5-7): kappa grid x LA arms. Usage: PAR=3 sh run_train.sh results/train "16 18" | "20 21"
OUT=$1; SIZES=$2; mkdir -p $OUT
{
for b in $SIZES; do
  if [ $b -le 18 ]; then KS="1.0 1.19 1.41 1.68 2.0 2.38"; LAS="lpsgewd@8,lpsgewd@16,lpsgewd@32,lpsgebw2@8,lpsgebw3@8,lpsgebw4@8"
  else KS=${KS21:-"1.19 1.41 1.68 2.0"}; LAS="lpsgewd@8,lpsgewd@16,lpsgebw2@8,lpsgebw3@8"; fi
  [ $b -eq 20 ] && KS="1.19 1.41 1.68 2.0"
  for k in $KS; do for c in 5 6 7; do
    echo "$OUT --bits $b --curve $c --kappa $k --cx 8 --las $LAS --no-gp-check --out $OUT/tr_b${b}_c${c}_k${k}.jsonl"
  done; done
done
} | xargs -P ${PAR:-3} -L 1 sh job.sh
