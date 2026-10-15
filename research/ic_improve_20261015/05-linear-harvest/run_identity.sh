#!/bin/bash
# Default identity: harvest off (override) and KIC_LINEAR_HARVEST=0 (env) on
# item 03's test targets t = 100.., compared with item 03's committed i_* counters.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 KIC_LINEAR_HARVEST=0
B=target/release/examples/ic_improve2_05_linear_harvest
D=research/ic_improve_20261015/05-linear-harvest/runs
P="nice prlimit --as=3000000000 --"
A="--first 100 --phase identity --arms i:off/e:env"
( for c in "13 4 3 12" "17 6 3 12" "19 7 3 12" "17 9 2 12" "23 12 2 12"; do set -- $c
    $P $B --n $1 --ell $2 --m $3 --targets $4 $A; done ) > $D/identity.jsonl
