#!/bin/sh
# Thin orchestration of the frozen grid in PROTOCOL.md: one process per cell.
# Launch detached: setsid nohup ./research/pdp_split_dreg_20261011/run_grid.sh &
set -u
D=research/pdp_split_dreg_20261011
for CELL in 3:1:11:4 3:0:13:4 3:1:17:4 3:1:19:4 3:0:13:5 3:1:17:5 4:1:11:3 4:0:13:3; do
  OUT="$D/runs/$(echo "$CELL" | tr ':' '-').jsonl"
  start=$(date -u +%FT%TZ)
  timeout 7200 target/release/pdp_split_dreg --cells "$CELL" --guesses 0,2,4 \
    --targets 4 --seed 20261011 --budget-secs 300 --max-free 24 \
    --lambda-targets 256 --out "$OUT" 2>"$OUT.stderr"
  rc=$?
  echo "$CELL start=$start end=$(date -u +%FT%TZ) exit=$rc" >> "$D/runs/exits.txt"
done
echo "grid done $(date -u +%FT%TZ)" >> "$D/runs/exits.txt"
