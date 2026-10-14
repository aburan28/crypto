#!/bin/bash
# Frozen ladder: PROTOCOL.md sizes and seeds.
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=/home/user/crypto/target/release/examples/ic_next10_01_fp4_ladder
cd /home/user/crypto/research/ic_next10_20261014/01-fp4-ladder
for p in 101 181 269 397 521 769 1033; do
  for s in 1 2 3; do
    timeout 2700 prlimit --as=2000000000 -- nice $B ladder --p $p --seed $s --rho-runs 16 --json raw/ladder.jsonl > /dev/null 2>> raw/ladder.stderr || echo "FAIL p=$p s=$s rc=$?" >> raw/ladder.stderr
    echo "done p=$p s=$s $(date +%T)" >> raw/ladder.progress
  done
done
