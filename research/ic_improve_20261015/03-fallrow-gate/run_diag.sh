#!/bin/bash
# Diagnostic arms (PROTOCOL §8): the "why" data; no verdict.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=target/release/examples/ic_improve2_03_fallrow_gate
D=research/ic_improve_20261015/03-fallrow-gate/runs
P="nice prlimit --as=3000000000 --"
( for c in "13 4 3" "17 6 3" "19 7 3" "23 12 2"; do set -- $c
    $P $B --family pipe --n $1 --ell $2 --m $3 --targets 4 --mode diag --first 200; done ) > $D/diag_pipe.jsonl &
( for n in 7 9 11; do $P $B --family barrel --n $n --targets 4 --mode diag --first 200; done ) > $D/diag_barrel.jsonl &
( for c in "13 4" "17 6" "19 7"; do set -- $c
    $P $B --family pipe --layout plain --n $1 --ell $2 --m 3 --targets 6 --mode train --first 200; done ) > $D/diag_plain_layout.jsonl &
wait
