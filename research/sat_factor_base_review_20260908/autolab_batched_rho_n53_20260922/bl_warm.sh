#!/bin/bash
set -u
LAB=$(cd "$(dirname "$0")" && pwd); OUT=$LAB/runs/bl_warm; mkdir -p $OUT; cd $OUT
IC=$LAB/../autolab_n53_eta_sweep_20260912/frozen/current/koblitz_rank_fixture
KS=$LAB/frozen/v3/koblitz_rho_batch_ks
C=n53-bl-warm-1024-v1
uptime > ic.load; KIC_BATCH_CORPUS=$C KIC_SHARED_PUBLIC_FIXTURE_DOMAIN=1 KIC_SHARED_FACTOR_LOG_PRECOMPUTATION=1 KIC_REQUIRED_SURPLUS_RELATIONS=0 /usr/bin/time -l $IC 53 0 1 10 531310 signed_expanded independent pair_pair_guided_256 1024 > ic.jsonl 2> ic.time; echo "ic rc=$? $(grep ' real' ic.time)"
for P in 242000 969000 3875000; do
  uptime > bl_$P.load; KIC_RHO_DP_BITS=6 KIC_RHO_PRECOMPUTE_WALKS=$P KIC_RHO_FREEZE_TABLE=1 KIC_RHO_BATCH_CORPUS=$C /usr/bin/time -l $KS 53 0 signed_frobenius 1024 531310 > bl_$P.jsonl 2> bl_$P.time; echo "bl P=$P rc=$? $(grep ' real' bl_$P.time)"
done
echo DONE
