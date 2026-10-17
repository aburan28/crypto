#!/bin/bash
# Default identity (PROTOCOL "If C wins"): targets t = 100.. (item 03's test targets).
#  A: KIC_LINEAR_HARVEST=0 in the environment: arm e (env) and i (thread override off)
#     must equal item 03's committed incumbent counters (pre-harvest engine).
#  B: KIC_LINEAR_HARVEST=gpv: arm e (env) must equal the explicit item 05 spec g, and
#     item 05's committed default arm d (runs/identity_default.jsonl there).
#  C: variable unset: arm d (env, the new default) must equal the explicit frozen spec c.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=target/release/examples/ic_improve3_02_harvest_gate
D=research/ic_improve3_20261015/02-harvest-gate/runs
P="nice prlimit --as=6000000000 --"
X="--first 100 --phase identity"
( for c in "13 4 3 12" "17 6 3 12" "19 7 3 12" "17 9 2 12" "23 12 2 12"; do set -- $c
    KIC_LINEAR_HARVEST=0 $P $B --n $1 --ell $2 --m $3 --targets $4 $X --arms i:off/e:env
  done ) > $D/identity_off.jsonl &
( for c in "19 7 3 6" "17 9 2 6" "23 12 2 6" "29 15 2 4" "17 6 3 6" "23 8 3 6" "31 16 2 4"; do set -- $c
    KIC_LINEAR_HARVEST=gpv $P $B --n $1 --ell $2 --m $3 --targets $4 $X --arms g:root,gpv=0.9607/e:env
  done
  for c in "19 7 3 6" "17 9 2 6" "23 12 2 6" "29 15 2 4" "17 6 3 6" "23 8 3 6" "31 16 2 4" "13 4 3 6"; do set -- $c
    env -u KIC_LINEAR_HARVEST $P $B --n $1 --ell $2 --m $3 --targets $4 $X --arms c:root,gpv=0.9607,rank=331/d:env
  done ) > $D/identity_gpv_default.jsonl &
wait
