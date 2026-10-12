#!/bin/bash
# n=61 paired panel at L=4096. K tuned on a disjoint 4096-target corpus.
# No L=32 panels.
set -u
LAB=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$LAB/../../.." && pwd)
OUT=$LAB/growing_n_n61_L4096
mkdir -p "$OUT"
cd "$OUT"
KS=${KS_OVERRIDE:-$ROOT/target/release/examples/koblitz_rho_batch_ks_v2_n61}
IC=$ROOT/target/release/examples/koblitz_orbit_dlp_fast
L=4096
N=61
shasum -a 256 "$KS" "$IC" > SHA256SUMS
echo "KS=$KS IC=$IC L=$L N=$N $(date -u +%FT%TZ)" | tee panel.log

candidates="600 800 1000 1200"
wall() { grep ' real' "$1" | awk '{print $1}'; }

corpus() {
  local name=$1
  if [ ! -s "scalars_$name.txt" ]; then
    echo "corpus $name" | tee -a panel.log
    KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=$name \
      "$KS" $N 0 signed_frobenius $L 531310 > "corpus_$name.jsonl" 2>/dev/null
    python3 -c "import json; [print(d['published_fixture_scalar']) for d in map(json.loads, open('corpus_$name.jsonl')) if d.get('kind')=='rho_ks_batch_fixture']" \
      > "scalars_$name.txt"
    wc -l "scalars_$name.txt" | tee -a panel.log
  fi
}

tune=n61-ks-growing-tune-$L-v1
eval_name=n61-ks-growing-$L-v1
corpus "$tune"
corpus "$eval_name"

if [ ! -s chosen_K_n61.txt ]; then
  best=; best_wall=
  for K in $candidates; do
    echo "tune K=$K start $(date -u +%FT%TZ)" | tee -a panel.log
    /usr/bin/time -l "$IC" "construct:$N:0:$K" "scalars_$tune.txt" 7 "tune_n61_K$K.jsonl" \
      > "tune_n61_K$K.summary.json" 2> "tune_n61_K$K.time"
    w=$(wall "tune_n61_K$K.time")
    echo "tune n=61 K=$K wall=$w" | tee -a panel.log
    if [ -z "$best" ] || awk "BEGIN{exit !($w < $best_wall)}"; then
      best=$K
      best_wall=$w
    fi
  done
  echo "$best" > chosen_K_n61.txt
  echo "chosen n=61 K=$best wall=$best_wall" | tee -a panel.log
fi

K=$(cat chosen_K_n61.txt)
run_ks() {
  local b=$1
  uptime | sed 's/.*averages: //' > "ks_n61_b$b.load"
  echo "ks b=$b start $(date -u +%FT%TZ) load=$(cat ks_n61_b$b.load)" | tee -a panel.log
  KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=$eval_name \
    /usr/bin/time -l "$KS" $N 0 signed_frobenius $L 531310 \
    > "ks_n61_b$b.jsonl" 2> "ks_n61_b$b.time"
  echo "ks n=61 b=$b rc=$? wall=$(wall ks_n61_b$b.time) load=$(cat ks_n61_b$b.load)" | tee -a panel.log
}
run_ic() {
  local b=$1
  uptime | sed 's/.*averages: //' > "ic_n61_b$b.load"
  echo "ic b=$b K=$K start $(date -u +%FT%TZ) load=$(cat ic_n61_b$b.load)" | tee -a panel.log
  KIC_DUMP_BASE=$OUT/base_n61_K$K.jsonl \
    /usr/bin/time -l "$IC" "construct:$N:0:$K" "scalars_$eval_name.txt" 7 "ic_n61_b$b.jsonl" \
    > "ic_n61_b$b.summary.json" 2> "ic_n61_b$b.time"
  echo "ic n=61 K=$K b=$b rc=$? wall=$(wall ic_n61_b$b.time) load=$(cat ic_n61_b$b.load)" | tee -a panel.log
}

for b in 0 1 2; do
  if [ $((b % 2)) = 0 ]; then run_ic $b; run_ks $b; else run_ks $b; run_ic $b; fi
done

echo "replay start $(date -u +%FT%TZ)" | tee -a panel.log
REPLAY_BASE=$OUT/base_n61_K$K.jsonl python3 "$LAB/independent_replay.py" \
  --dlp "$OUT" "$OUT/independent_replay.json" 61 | tee -a panel.log
echo "DONE $(date -u +%FT%TZ)" | tee -a panel.log
