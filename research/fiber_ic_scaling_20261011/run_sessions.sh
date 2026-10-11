#!/bin/bash
# Runs every spec of the round into research/fiber_ic_scaling_20261011/sessions,
# one at a time under the benchmark lock, waiting for a quiet host before each
# (ecbench refuses a busy host with exit 2), then audits each with exact replays.
set -u
cd "$(dirname "$0")/../.."
B=./target/release/ecbench
D=research/fiber_ic_scaling_20261011
order="scaling_h4 scaling_h72"
while pgrep -f "ecbench run" > /dev/null; do sleep 15; done
for name in $order; do
  spec=$D/specs/$name.json
  out=$D/sessions/$name
  if [ -d "$out" ] && [ ! -f "$out/session.json" ]; then rm -rf "$out"; fi
  if [ ! -d "$out" ]; then
    echo "== $name: start $(date -u +%FT%TZ)"
    tries=0
    while :; do
      $B run --wait --spec "$spec" --out "$out" > "$D/sessions/$name.run.log" 2>&1
      rc=$?
      if [ $rc -eq 2 ] && grep -q "host is not quiet" "$D/sessions/$name.run.log" && [ $tries -lt 40 ]; then
        tries=$((tries+1)); rm -rf "$out"; sleep 20; continue
      fi
      break
    done
    echo "   exit $rc after $tries busy retries at $(date -u +%FT%TZ); $(tail -n 1 "$D/sessions/$name.run.log")"
  else
    echo "== $name: exists"
  fi
  if [ -f "$out/session.json" ] && [ ! -f "$D/sessions/$name.audit.json" ]; then
    $B verify --dir "$out" --replay 12 --out "$D/sessions/$name.audit.json" > "$D/sessions/$name.verify.log" 2>&1
    echo "   verify exit $?: $(tail -n 1 "$D/sessions/$name.verify.log")"
    $B table --dir "$out" > "$D/sessions/$name.table.md" 2>&1
  fi
done
echo "== all done $(date -u +%FT%TZ)"
