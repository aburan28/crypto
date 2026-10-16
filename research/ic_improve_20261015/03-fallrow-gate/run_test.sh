#!/bin/bash
# Test phase (PROTOCOL §3, §7): measure-mode rows, fresh targets t = 100.., frozen gate.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=target/release/examples/ic_improve2_03_fallrow_gate
D=research/ic_improve_20261015/03-fallrow-gate/runs
P="nice prlimit --as=3000000000 --"
TH=$(python3 -c "import json;print(json.load(open('$D/gate.json'))['theta'])")
PI=$(python3 -c "import json;print(json.load(open('$D/gate.json'))['pi0'])")
G="--mode measure --theta $TH --nonfire $PI --first 100"
( for c in "13 4 12" "17 6 12" "19 7 12" "23 8 12" "29 10 8" "31 10 8"; do set -- $c
    $P $B --family pipe --n $1 --ell $2 --m 3 --targets $3 $G; done ) > $D/test_pipe3.jsonl &
( for c in "17 9 12" "23 12 12" "29 15 8"; do set -- $c
    $P $B --family pipe --n $1 --ell $2 --m 2 --targets $3 $G; done
  for c in "5 12" "7 12" "9 12"; do set -- $c
    $P $B --family barrel --n $1 --targets $2 $G; done ) > $D/test_pipe2_barrel_small.jsonl &
( $P $B --family barrel --n 11 --targets 8 $G
  $P $B --family barrel --n 13 --targets 6 $G ) > $D/test_barrel_large.jsonl &
wait
