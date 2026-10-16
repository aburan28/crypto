#!/bin/bash
# Held-out test split (PROTOCOL.md): seeds 1-3 at item 01's 8 sizes plus the p = 1601 stretch.
# Needs raw/frozen.json (from `fit`), written before this runs.  Two lanes at most.
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=/home/user/crypto/target/release/examples/ic_improve2_04_d4_2lp_fix
cd /home/user/crypto/research/ic_improve_20261015/04-d4-2lp-fix
for p in 1601 1301 1033 769 521 397 269 181 101; do for s in 1 2 3; do R=0; [ $p = 1601 ] && R=2; echo "$p $s $R"; done; done |
  xargs -P 2 -n 3 sh -c 'timeout 3600 prlimit --as=3000000000 -- nice '"$B"' test --p $0 --seed $1 --rho-runs $2 --frozen raw/frozen.json --json raw/test.jsonl > /dev/null 2>> raw/test.stderr || echo "FAIL p=$0 s=$1 rc=$?" >> raw/test.stderr; echo "done p=$0 s=$1 $(date +%T)" >> raw/test.progress'
