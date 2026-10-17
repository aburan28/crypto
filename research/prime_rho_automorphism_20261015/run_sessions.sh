#!/bin/bash
# Runs the cubic and quartic specs one at a time, then audits each with exact replays.
set -u
cd "$(dirname "$0")/../.."
B=./target/release/ecbench
D=research/prime_rho_automorphism_20261015
for name in cubic quartic; do
  out=$D/sessions/$name
  mkdir -p $D/sessions
  if [ ! -f "$out/session.json" ]; then
    rm -rf "$out"
    echo "== $name: start $(date -u +%FT%TZ)"
    $B run --wait --spec $D/specs/$name.json --out "$out" > "$D/sessions/$name.run.log" 2>&1
    echo "   exit $? at $(date -u +%FT%TZ); $(tail -n 1 "$D/sessions/$name.run.log")"
  fi
  if [ -f "$out/session.json" ]; then
    $B verify --dir "$out" --replay 12 --out "$D/sessions/$name.audit.json" > "$D/sessions/$name.verify.log" 2>&1
    echo "   verify exit $?: $(tail -n 1 "$D/sessions/$name.verify.log")"
    $B table --dir "$out" > "$D/sessions/$name.table.md" 2>&1
  fi
done
echo "== all done $(date -u +%FT%TZ)"
