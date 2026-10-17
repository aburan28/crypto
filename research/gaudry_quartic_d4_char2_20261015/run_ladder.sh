#!/bin/sh
# Part A of PROTOCOL.md. From the repository root, after
#   cargo build --release --example gaudry_quartic_bin
set -e
B=target/release/examples/gaudry_quartic_bin
R=research/gaudry_quartic_d4_char2_20261015/raw
NEW="--frobenius --sparsest --replay"
run() { name=$1; shift; "$B" "$@" --json "$R/$name.jsonl" > "$R/$name.log" 2>&1; }
# Ablation, m = 8, 10 planted targets each.
run abl_plain_s       --m 8 --seed 1 --targets 10 --planted --sparsest
run abl_plain_sr      --m 8 --seed 1 --targets 10 --planted --sparsest --replay
run abl_weighted      --m 8 --seed 1 --targets 10 --planted --weighted
run abl_weighted_s    --m 8 --seed 1 --targets 10 --planted --weighted --sparsest
run abl_weighted_sr   --m 8 --seed 1 --targets 10 --planted --weighted --sparsest --replay
run abl_frob          --m 8 --seed 1 --targets 10 --planted --frobenius
run abl_frob_s        --m 8 --seed 1 --targets 10 --planted --frobenius --sparsest
run abl_frob_sr       --m 8 --seed 1 --targets 10 --planted $NEW
# Ladder.
for m in 5 6 7 8 10 12; do
  run m${m}_planted --m $m --seed 1 --targets 20 --planted $NEW
  run m${m}_uniform --m $m --seed 1 --targets 50 $NEW
done
run m15_uniform_cost_only --m 15 --seed 1 --targets 20 $NEW --no-mitm
