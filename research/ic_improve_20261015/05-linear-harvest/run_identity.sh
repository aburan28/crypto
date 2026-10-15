#!/bin/bash
# Default identity: with the opt-out KIC_LINEAR_HARVEST=0 in the environment,
# arm e (env) and arm i (thread override off) on item 03's test targets t = 100..
# must reproduce item 03's committed incumbent counters exactly.  Arm d runs
# with no environment variable and no override (the new default) and arm g the
# explicit frozen spec: d must equal g counter for counter.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=target/release/examples/ic_improve2_05_linear_harvest
D=research/ic_improve_20261015/05-linear-harvest/runs
P="nice prlimit --as=3000000000 --"
A="--first 100 --phase identity"
( for c in "13 4 3 12" "17 6 3 12" "19 7 3 12" "17 9 2 12" "23 12 2 12"; do set -- $c
    KIC_LINEAR_HARVEST=0 $P $B --n $1 --ell $2 --m $3 --targets $4 $A --arms i:off/e:env
  done ) > $D/identity.jsonl
( for c in "19 7 3 6" "17 9 2 6" "23 12 2 6" "29 15 2 4"; do set -- $c
    env -u KIC_LINEAR_HARVEST $P $B --n $1 --ell $2 --m $3 --targets $4 $A --arms g:root,gpv=0.9607/d:env
  done ) > $D/identity_default.jsonl
