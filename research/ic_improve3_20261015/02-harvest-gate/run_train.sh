#!/bin/bash
# Training (PROTOCOL "Training"): t = 0..11, 12 cells, bench mode, arms i, g (G0), r, i2.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=target/release/examples/ic_improve3_02_harvest_gate
D=research/ic_improve3_20261015/02-harvest-gate/runs
P="nice prlimit --as=6000000000 --"
A="--targets 12 --first 0 --phase train --arms i:off/g:root,gpv=0.9607/r:root/i2:off"
( for c in "13 4 3" "17 6 3" "19 7 3" "17 9 2" "19 10 2" "23 12 2" "29 15 2" "29 10 3"; do set -- $c
    $P $B --n $1 --ell $2 --m $3 $A; done ) > $D/train_a.jsonl &
( for c in "23 8 3" "31 16 2" "31 10 3" "37 19 2"; do set -- $c
    $P $B --n $1 --ell $2 --m $3 $A; done ) > $D/train_b.jsonl &
wait
