#!/bin/bash
# End-to-end (PROTOCOL "End-to-end"): m = 2, (23,12), (29,15) and (added) (31,16);
# 3 replicate rounds, arm order rotated per round, one process per arm.
# Arms (KIC_LINEAR_HARVEST): 0 (pre-harvest engine), gpv (item 05 gate, the incumbent),
# unset (the new default, the frozen candidate).
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=target/release/examples/ic_improve2_05_e2e
D=research/ic_improve3_20261015/02-harvest-gate/runs
P="nice prlimit --as=6000000000 --"
run() { if [ "$1" = unset ]; then env -u KIC_LINEAR_HARVEST $P $B "${@:2}"; else KIC_LINEAR_HARVEST=$1 $P $B "${@:2}"; fi; }
ARMS=(0 gpv unset)
for c in "23 12 2" "29 15 2" "31 16 2"; do set -- $c
  for r in 0 1 2; do
    for j in 0 1 2; do run ${ARMS[$(( (j + r) % 3 ))]} --n $1 --ell $2 --m $3; done
  done
done > ${OUT:-$D/e2e.jsonl}
