#!/bin/bash
set -u
cd /private/tmp/claude-501/-Volumes-SSD990-cairn--claude-worktrees-crypto-autoresearcher-e2e-macos-gui-70374d/437b5898-3f94-4bd6-b829-6e9a1e5d3668/scratchpad/runs/paired 
IC=/Volumes/SSD990/crypto/research/sat_factor_base_review_20260908/autolab_n53_eta_sweep_20260912/frozen/current/koblitz_rank_fixture
KS=/private/tmp/claude-501/-Volumes-SSD990-cairn--claude-worktrees-crypto-autoresearcher-e2e-macos-gui-70374d/437b5898-3f94-4bd6-b829-6e9a1e5d3668/scratchpad/target/release/examples/koblitz_rho_batch_ks
run_ic() { KIC_BATCH_CORPUS=$2 KIC_SHARED_PUBLIC_FIXTURE_DOMAIN=1 KIC_SHARED_FACTOR_LOG_PRECOMPUTATION=1 KIC_REQUIRED_SURPLUS_RELATIONS=0 /usr/bin/time -l $IC 53 0 1 10 531310 signed_expanded independent pair_pair_guided_256 $1 > ic$1_b$3.jsonl 2> ic$1_b$3.time; echo "ic L=$1 b=$3 rc=$? $(grep ' real' ic$1_b$3.time) load=$(uptime | sed 's/.*averages: //')"; }
run_ks() { KIC_RHO_BATCH_CORPUS=$2 /usr/bin/time -l $KS 53 0 signed_frobenius $1 531310 > ks$1_b$3.jsonl 2> ks$1_b$3.time; echo "ks L=$1 b=$3 rc=$? $(grep ' real' ks$1_b$3.time)"; }
for b in 0 1 2; do
  for spec in "32 n53-shared-factor-logs-crossover-v1" "1024 n53-shared-factor-logs-1024-v1"; do
    set -- $spec
    if [ $((b % 2)) = 0 ]; then run_ic $1 $2 $b; run_ks $1 $2 $b; else run_ks $1 $2 $b; run_ic $1 $2 $b; fi
  done
done
echo DONE
