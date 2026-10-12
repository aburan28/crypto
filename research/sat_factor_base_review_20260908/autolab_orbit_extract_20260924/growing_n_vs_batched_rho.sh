#!/bin/bash
# Growing-n (n=37/41/53) paired panel at fixed L: compact-orbit shared-log DLP
# vs frozen Kuhn-Struik batched rho on identical published targets.
# Phase 1 picks K per n on a disjoint tuning corpus (lowest process wall).
# Phase 2 evaluates on a fresh corpus, 3 blocks, alternating arm order.
set -u
LAB=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$LAB/../../.." && pwd)
OUT=$LAB/growing_n_vs_batched_rho; mkdir -p "$OUT"; cd "$OUT"
KS=${KS_OVERRIDE:-$ROOT/research/sat_factor_base_review_20260908/autolab_batched_rho_n53_20260922/frozen/v2/koblitz_rho_batch_ks}
IC=$ROOT/target/release/examples/koblitz_orbit_dlp_fast
L=${PANEL_L_FIXED:-1024}
shasum -a 256 "$KS" "$IC" > SHA256SUMS
candidates() { case $1 in 37) echo "7 14 28 42";; 41) echo "85 170 255 340";; 53) echo "220 440 660 880";; 61) echo "150 300 450 600";; esac; }
corpus() {  # $1=n $2=name
  if [ ! -s scalars_$2.txt ]; then
    KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=$2 "$KS" $1 0 signed_frobenius $L 531310 > corpus_$2.jsonl 2>/dev/null
    python3 -c "import json; [print(d['published_fixture_scalar']) for d in map(json.loads, open('corpus_$2.jsonl')) if d.get('kind')=='rho_ks_batch_fixture']" > scalars_$2.txt
  fi
}
wall() { grep ' real' $1 | awk '{print $1}'; }
for n in ${PANEL_N:-37 41 53}; do
  tune=n$n-ks-growing-tune-$L-v1; corpus $n $tune
  if [ ! -s chosen_K_n$n.txt ]; then
    best=; best_wall=
    for K in $(candidates $n); do
      /usr/bin/time -l "$IC" construct:$n:0:$K scalars_$tune.txt 7 tune_n${n}_K$K.jsonl > tune_n${n}_K$K.summary.json 2> tune_n${n}_K$K.time
      w=$(wall tune_n${n}_K$K.time); echo "tune n=$n K=$K wall=$w"
      if [ -z "$best" ] || awk "BEGIN{exit !($w < $best_wall)}"; then best=$K; best_wall=$w; fi
    done
    echo $best > chosen_K_n$n.txt; echo "chosen n=$n K=$best"
  fi
done
run_ks() { uptime | sed 's/.*averages: //' > ks_n$1_b$2.load; KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=n$1-ks-growing-$L-v1 /usr/bin/time -l "$KS" $1 0 signed_frobenius $L 531310 > ks_n$1_b$2.jsonl 2> ks_n$1_b$2.time; echo "ks n=$1 b=$2 rc=$? wall=$(wall ks_n$1_b$2.time) load=$(cat ks_n$1_b$2.load)"; }
run_ic() { K=$(cat chosen_K_n$1.txt); uptime | sed 's/.*averages: //' > ic_n$1_b$2.load; KIC_DUMP_BASE=$OUT/base_n$1_K$K.jsonl /usr/bin/time -l "$IC" construct:$1:0:$K scalars_n$1-ks-growing-$L-v1.txt 7 ic_n$1_b$2.jsonl > ic_n$1_b$2.summary.json 2> ic_n$1_b$2.time; echo "ic n=$1 K=$K b=$2 rc=$? wall=$(wall ic_n$1_b$2.time) load=$(cat ic_n$1_b$2.load)"; }
for n in ${PANEL_N:-37 41 53}; do corpus $n n$n-ks-growing-$L-v1; done
for b in 0 1 2; do
  for n in ${PANEL_N:-37 41 53}; do
    if [ $((b % 2)) = 0 ]; then run_ic $n $b; run_ks $n $b; else run_ks $n $b; run_ic $n $b; fi
  done
done
echo DONE
