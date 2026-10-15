#!/bin/bash
# Test phase (PROTOCOL "Test"): held-out targets t = 300.., frozen θ from runs/train_fit.jsonl.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
unset KIC_LINEAR_HARVEST
B=target/release/examples/ic_improve2_05_linear_harvest
D=research/ic_improve_20261015/05-linear-harvest/runs
P="nice prlimit --as=3000000000 --"
TH=$(grep -o '"train":true,"theta":[0-9.]*' $D/train_fit.jsonl | sed 's/.*://')
# g: frozen candidate; r: root ungated; u: root, gpv<1 (exploratory); h: every node; i2: A/A.
# (The recorded run spelled h as "1", which meant every node before the default change.)
A="--first 300 --phase test --arms i:off/g:root,gpv=$TH/r:root/u:root,gpv=1.0/h:every/i2:off"
( for c in "13 4 3 12" "17 6 3 12" "19 7 3 12" "17 9 2 12" "23 12 2 12" "29 15 2 8"; do set -- $c
    $P $B --n $1 --ell $2 --m $3 --targets $4 $A; done ) > $D/test_a.jsonl &
( $P $B --n 23 --ell 8 --m 3 --targets 12 $A; $P $B --n 31 --ell 10 --m 3 --targets 8 $A ) > $D/test_b.jsonl &
( $P $B --n 29 --ell 10 --m 3 --targets 8 $A ) > $D/test_c.jsonl &
wait
