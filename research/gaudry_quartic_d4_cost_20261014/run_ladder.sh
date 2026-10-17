#!/bin/sh
# Part A of PROTOCOL.md: the decomposition ladder and the p = 269 ablation.
# Run from the repository root after
#   cargo build --release --example gaudry_quartic_d4
# Each line appends one JSON row per target to raw/<name>.jsonl and its log
# to raw/<name>.log.
set -e
B=target/release/examples/gaudry_quartic_d4
R=research/gaudry_quartic_d4_cost_20261014/raw
NEW="--weighted --chains --sparsest --replay"
run() { name=$1; shift; "$B" "$@" --json "$R/$name.jsonl" > "$R/$name.log" 2>&1; }
# Ablation, p = 269, seed 1, 10 planted targets.
run abl_committed            --p 269 --seed 1 --targets 10 --planted
run abl_s                    --p 269 --seed 1 --targets 10 --planted --sparsest
run abl_sc                   --p 269 --seed 1 --targets 10 --planted --sparsest --chains
run abl_scr                  --p 269 --seed 1 --targets 10 --planted --sparsest --chains --replay
run abl_w                    --p 269 --seed 1 --targets 10 --planted --weighted
run abl_ws                   --p 269 --seed 1 --targets 10 --planted --weighted --sparsest
run abl_wc                   --p 269 --seed 1 --targets 10 --planted --weighted --chains
run abl_wcs                  --p 269 --seed 1 --targets 10 --planted --weighted --chains --sparsest
run abl_wcsr                 --p 269 --seed 1 --targets 10 --planted $NEW
# Ladder.
run p29_new_planted          --p 29 --seed 1 --find 16 --targets 20 --planted $NEW
run p29_new_uniform          --p 29 --seed 1 --find 16 --targets 300 $NEW
run p29_committed_planted    --p 29 --seed 1 --find 16 --targets 20 --planted
run p29_committed_uniform    --p 29 --seed 1 --find 16 --targets 50
run p269_new_planted         --p 269 --seed 1 --targets 20 --planted $NEW
run p269_new_uniform         --p 269 --seed 1 --targets 50 $NEW
run p269_committed_planted   --p 269 --seed 1 --targets 10 --planted
run p1009_new_planted        --p 1009 --seed 1 --targets 20 --planted $NEW
run p1009_new_uniform        --p 1009 --seed 1 --targets 50 $NEW
run p1009_committed_planted  --p 1009 --seed 1 --targets 5 --planted
run p4049_new_planted        --p 4049 --seed 1 --targets 10 --planted $NEW
run p4049_new_uniform        --p 4049 --seed 1 --targets 20 $NEW
run p4049_committed_planted  --p 4049 --seed 1 --targets 3 --planted
