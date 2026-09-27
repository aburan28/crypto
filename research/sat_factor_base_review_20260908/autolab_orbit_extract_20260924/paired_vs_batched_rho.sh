#!/bin/bash
# Paired n=53 panel: compact-orbit shared-log DLP vs frozen Kuhn-Struik batched rho
# on identical published targets. Batched rho generates the corpus; its published
# scalars are fed to the compact arm. Block parity alternates which arm runs first.
set -u
LAB=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$LAB/../../.." && pwd)
OUT=${PANEL_OUT:-$LAB/paired_vs_batched_rho}; mkdir -p "$OUT"; cd "$OUT"
KS=$ROOT/research/sat_factor_base_review_20260908/autolab_batched_rho_n53_20260922/frozen/v2/koblitz_rho_batch_ks
IC=$ROOT/target/release/examples/koblitz_orbit_dlp_fast
BASE=$ROOT/research/sat_factor_base_review_20260908/autolab_shared_log_scaling_20260912/runs/n53/block_00/direct/stdout.txt
head -n 1 "$BASE" > base_n53.jsonl
shasum -a 256 "$KS" "$IC" base_n53.jsonl > SHA256SUMS
corpus() {
  if [ ! -s ks_corpus_$1.jsonl ]; then
    KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=n53-ks-scaling-$1-v1 "$KS" 53 0 signed_frobenius $1 531310 > ks_corpus_$1.jsonl 2>/dev/null
    python3 -c "import json,sys; [print(d['published_fixture_scalar']) for d in map(json.loads, open('ks_corpus_$1.jsonl')) if d.get('kind')=='rho_ks_batch_fixture']" > scalars_$1.txt
  fi
}
run_ks() {
  if [ -s ks$1_b$2.jsonl ] && [ -n "${REUSE_KS:-}" ]; then echo "ks L=$1 b=$2 reused"; return; fi
  uptime | sed 's/.*averages: //' > ks$1_b$2.load
  KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=n53-ks-scaling-$1-v1 /usr/bin/time -l "$KS" 53 0 signed_frobenius $1 531310 > ks$1_b$2.jsonl 2> ks$1_b$2.time
  echo "ks L=$1 b=$2 rc=$? $(grep ' real' ks$1_b$2.time | tr -s ' ') load=$(cat ks$1_b$2.load)"
}
run_ic() {
  K=$(eval echo "\${K_$1:-}")
  if [ -n "$K" ]; then SRC=construct:53:0:$K; TAG=K${K}_; export KIC_DUMP_BASE=$OUT/base_K${K}_L$1_b$2.jsonl; else SRC=base_n53.jsonl; TAG=; unset KIC_DUMP_BASE; fi
  uptime | sed 's/.*averages: //' > ic${TAG}$1_b$2.load
  /usr/bin/time -l "$IC" $SRC scalars_$1.txt 7 ic${TAG}$1_b$2.jsonl > ic${TAG}$1_b$2.summary.json 2> ic${TAG}$1_b$2.time
  echo "ic ${TAG}L=$1 b=$2 rc=$? $(grep ' real' ic${TAG}$1_b$2.time | tr -s ' ') load=$(cat ic${TAG}$1_b$2.load)"
}
for L in ${PANEL_L:-32 1024 4096}; do corpus $L; done
for b in 0 1 2; do
  for L in ${PANEL_L:-32 1024 4096}; do
    if [ $((b % 2)) = 0 ]; then run_ic $L $b; run_ks $L $b; else run_ks $L $b; run_ic $L $b; fi
  done
done
echo DONE
