#!/bin/bash
# End-to-end (PROTOCOL §7a): one process per arm, 3 replicate pairs, alternating order.
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
B=target/release/examples/ic_improve2_03_e2e
D=research/ic_improve_20261015/03-fallrow-gate/runs
P="nice prlimit --as=3000000000 --"
for c in "17 6 3" "19 7 3"; do set -- $c
  for r in 1 2 3; do
    if [ $((r % 2)) = 1 ]; then
      env -u KIC_FALLROW_GATE $P $B --n $1 --ell $2 --m $3; KIC_FALLROW_GATE=1 $P $B --n $1 --ell $2 --m $3
    else
      KIC_FALLROW_GATE=1 $P $B --n $1 --ell $2 --m $3; env -u KIC_FALLROW_GATE $P $B --n $1 --ell $2 --m $3
    fi
  done
done > $D/e2e.jsonl
