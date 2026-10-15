#!/bin/bash
# Held-out test split (PROTOCOL.md): seeds 1-3, item 1's sizes plus the p = 1301 stretch.
# Needs raw/frozen.json (from `tune`), written before this runs.
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=/home/user/crypto/target/release/examples/ic_improve_01_d4_2lp
cd /home/user/crypto/research/ic_improve_20261014/01-d4-2lp
for p in 101 181 269 397 521 769 1033 1301; do
  R=0; [ $p = 1301 ] && R=4   # rho-negation reference only where item 1 has none
  for s in 1 2 3; do
    timeout 3600 prlimit --as=3000000000 -- nice $B test --p $p --seed $s --frozen raw/frozen.json --rho-runs $R --json raw/test.jsonl > /dev/null 2>> raw/test.stderr || echo "FAIL p=$p s=$s rc=$?" >> raw/test.stderr
    echo "done p=$p s=$s $(date +%T)" >> raw/test.progress
  done
done
