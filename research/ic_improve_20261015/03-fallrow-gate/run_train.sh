#!/bin/bash
# Training phase (PROTOCOL §5): train-mode rows on the training cells, t = 0..7.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=target/release/examples/ic_improve2_03_fallrow_gate
D=research/ic_improve_20261015/03-fallrow-gate/runs
P="nice prlimit --as=3000000000 --"
( $P $B --family pipe --n 13 --ell 4 --m 3 --targets 8 --mode train
  $P $B --family pipe --n 17 --ell 6 --m 3 --targets 8 --mode train ) > $D/train_pipe3.jsonl &
( $P $B --family pipe --n 17 --ell 9 --m 2 --targets 8 --mode train
  $P $B --family pipe --n 19 --ell 10 --m 2 --targets 8 --mode train ) > $D/train_pipe2.jsonl &
( for n in 5 7 9; do $P $B --family barrel --n $n --targets 8 --mode train; done ) > $D/train_barrel.jsonl &
wait
