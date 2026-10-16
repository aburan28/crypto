#!/bin/bash
# End-to-end (PROTOCOL "End-to-end"): m = 2, two sizes, 3 replicate pairs,
# alternating order, one process per arm.  Arms: KIC_LINEAR_HARVEST=0 (the
# incumbent, opt-out) and the variable unset (the new default = frozen candidate).
# runs/e2e_optin.jsonl is the earlier run, before the default change, with the
# candidate selected as KIC_LINEAR_HARVEST="root,gpv=0.9607" (the same configuration).
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=target/release/examples/ic_improve2_05_e2e
D=research/ic_improve_20261015/05-linear-harvest/runs
P="nice prlimit --as=3000000000 --"
for c in "19 10 2" "23 12 2"; do set -- $c
  for r in 1 2 3; do
    if [ $((r % 2)) = 1 ]; then
      KIC_LINEAR_HARVEST=0 $P $B --n $1 --ell $2 --m $3; env -u KIC_LINEAR_HARVEST $P $B --n $1 --ell $2 --m $3
    else
      env -u KIC_LINEAR_HARVEST $P $B --n $1 --ell $2 --m $3; KIC_LINEAR_HARVEST=0 $P $B --n $1 --ell $2 --m $3
    fi
  done
done > ${OUT:-$D/e2e.jsonl}
