#!/bin/bash
# Training split (PROTOCOL.md): seeds 4, 5; all 60 grid configs + d4 + null on one stream.
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=/home/user/crypto/target/release/examples/ic_improve_01_d4_2lp
cd /home/user/crypto/research/ic_improve_20261014/01-d4-2lp
for p in 101 181 269 397 521 769 1033; do
  for s in 4 5; do
    timeout 3600 prlimit --as=3000000000 -- nice $B train --p $p --seed $s --json raw/train.jsonl > /dev/null 2>> raw/train.stderr || echo "FAIL p=$p s=$s rc=$?" >> raw/train.stderr
    echo "done p=$p s=$s $(date +%T)" >> raw/train.progress
  done
done
