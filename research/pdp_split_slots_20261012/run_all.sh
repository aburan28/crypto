#!/bin/sh
# Thin orchestration of PROTOCOL.md Parts A, B, C, in that order.
# Launch detached: setsid nohup ./research/pdp_split_slots_20261012/run_all.sh &
set -u
D=research/pdp_split_slots_20261012
mkdir -p "$D/sessions" "$D/runs-B" "$D/audits"
log() { echo "$1 $(date -u +%FT%TZ) $2" >> "$D/progress.txt"; }
for n in 11 13 17 19; do
  log "A-n$n" start
  target/release/ecbench run --wait --spec "$D/SPEC-A-n$n.json" --out "$D/sessions/A-n$n" --cpus none \
    > "$D/sessions/A-n$n.log" 2>&1
  log "A-n$n" "end exit=$?"
done
for CELL in 3:1:11:4 3:0:13:4 3:1:17:4 3:1:19:4 3:0:13:5 3:1:17:5; do
  OUT="$D/runs-B/$(echo "$CELL" | tr ':' '-').jsonl"
  log "B-$CELL" start
  timeout 7200 target/release/pdp_split_dreg --cells "$CELL" --base closure --slots-per-target 4 \
    --targets 4 --guesses 0,2,4 --seed 20261012 --budget-secs 300 --max-free 24 \
    --lambda-targets 256 --lambda-systems 512 --out "$OUT" 2>"$OUT.stderr"
  log "B-$CELL" "end exit=$?"
done
log "C-n17" start
target/release/ecbench run --wait --spec "$D/SPEC-C-n17.json" --out "$D/sessions/C-n17" --cpus none \
  > "$D/sessions/C-n17.log" 2>&1
log "C-n17" "end exit=$?"
for s in A-n11 A-n13 A-n17 A-n19 C-n17; do
  target/release/ecbench verify --dir "$D/sessions/$s" --replay-all --out "$D/audits/$s.audit.json" \
    > "$D/audits/$s.verify.log" 2>&1
  log "audit-$s" "end exit=$?"
done
log all done
