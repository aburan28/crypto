#!/bin/bash
# Test (PROTOCOL "Test" and "Production path"): held-out t = 2000.., frozen spec from runs/train_fit.jsonl.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
unset KIC_LINEAR_HARVEST
B=target/release/examples/ic_improve3_02_harvest_gate
D=research/ic_improve3_20261015/02-harvest-gate/runs
P="nice prlimit --as=6000000000 --"
C=$(grep -o '"spec":"[^"]*"' $D/train_fit.jsonl | head -1 | sed 's/"spec":"//;s/"$//')
echo "frozen candidate: $C" >&2
AB="--first 2000 --phase test --mode bench --arms i:off/g:root,gpv=0.9607/c:$C/r:root/i2:off"
AP="--first 2000 --phase test --mode prod --arms i:off/g:root,gpv=0.9607/c:$C/i2:off"
( for c in "13 4 3 12" "17 6 3 12" "19 7 3 12" "17 9 2 12" "19 10 2 12" "23 12 2 12" "29 15 2 8" "31 16 2 8" "23 8 3 12"; do set -- $c
    $P $B --n $1 --ell $2 --m $3 --targets $4 $AB; $P $B --n $1 --ell $2 --m $3 --targets $4 $AP; done ) > $D/test_a.jsonl &
( for c in "29 10 3 8" "31 10 3 8" "37 19 2 8"; do set -- $c
    $P $B --n $1 --ell $2 --m $3 --targets $4 $AB; $P $B --n $1 --ell $2 --m $3 --targets $4 $AP; done ) > $D/test_b.jsonl &
wait
