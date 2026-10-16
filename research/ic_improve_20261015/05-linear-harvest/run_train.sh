#!/bin/bash
# Training phase (PROTOCOL "Training"): arms i, r, i2 on t = 0..7.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=target/release/examples/ic_improve2_05_linear_harvest
D=research/ic_improve_20261015/05-linear-harvest/runs
P="nice prlimit --as=3000000000 --"
A="--targets 8 --first 0 --phase train --arms i:off/r:root/i2:off"
( $P $B --n 13 --ell 4 --m 3 $A; $P $B --n 17 --ell 6 --m 3 $A ) > $D/train_m3.jsonl &
( $P $B --n 17 --ell 9 --m 2 $A; $P $B --n 19 --ell 10 --m 2 $A ) > $D/train_m2.jsonl &
wait
