#!/bin/bash
set -u
LAB=$(cd "$(dirname "$0")" && pwd)
OUT=$LAB/runs/panel_v2; mkdir -p $OUT; cd $OUT
IC=$LAB/../autolab_n53_eta_sweep_20260912/frozen/current/koblitz_rank_fixture
KS=$LAB/frozen/v2/koblitz_rho_batch_ks
run_ic() { echo "$(uptime | sed 's/.*averages: //')" > ic$1_b$2.load; KIC_BATCH_CORPUS=n53-ks-scaling-$1-v1 KIC_SHARED_PUBLIC_FIXTURE_DOMAIN=1 KIC_SHARED_FACTOR_LOG_PRECOMPUTATION=1 KIC_REQUIRED_SURPLUS_RELATIONS=0 /usr/bin/time -l $IC 53 0 1 10 531310 signed_expanded independent pair_pair_guided_256 $1 > ic$1_b$2.jsonl 2> ic$1_b$2.time; echo "ic L=$1 b=$2 rc=$? $(grep ' real' ic$1_b$2.time) load=$(cat ic$1_b$2.load)"; }
run_ks() { echo "$(uptime | sed 's/.*averages: //')" > ks$1_b$2.load; KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=n53-ks-scaling-$1-v1 /usr/bin/time -l $KS 53 0 signed_frobenius $1 531310 > ks$1_b$2.jsonl 2> ks$1_b$2.time; echo "ks L=$1 b=$2 rc=$? $(grep ' real' ks$1_b$2.time) load=$(cat ks$1_b$2.load)"; }
for b in 0 1 2; do
  for L in 1024 2048 4096 8192; do
    if [ $((b % 2)) = 0 ]; then run_ic $L $b; run_ks $L $b; else run_ks $L $b; run_ic $L $b; fi
  done
done
echo DONE
