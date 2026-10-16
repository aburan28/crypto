#!/bin/bash
# Training split (PROTOCOL.md): seeds 4, 5 at item 01's 7 sizes; two lanes at most.
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=/home/user/crypto/target/release/examples/ic_improve2_04_d4_2lp_fix
cd /home/user/crypto/research/ic_improve_20261015/04-d4-2lp-fix
for p in 1033 769 521 397 269 181 101; do for s in 4 5; do echo "$p $s"; done; done |
  xargs -P 2 -n 2 sh -c 'timeout 3600 prlimit --as=3000000000 -- nice '"$B"' train --p $0 --seed $1 --json raw/train.jsonl > /dev/null 2>> raw/train.stderr || echo "FAIL p=$0 s=$1 rc=$?" >> raw/train.stderr; echo "done p=$0 s=$1 $(date +%T)" >> raw/train.progress'
