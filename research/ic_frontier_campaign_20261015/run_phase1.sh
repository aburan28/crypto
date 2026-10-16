#!/bin/sh
# Phase 1 of PROTOCOL.md: every frozen spec into phase1/sessions/, counts only (L0),
# then a replay-all audit of each session.  Thin orchestration only.
set -u
cd "$(dirname "$0")"
E=${ECBENCH:-../../target/release/ecbench}
mkdir -p phase1/sessions phase1/audits
for s in P20 P24 P28 P32 P36 P40 PKM20 PKM22 PKM24 PKM26 PKM28 K13 K17 K19 K23 K37 K43 K47 K41; do
  [ -f phase1/sessions/$s/session.json ] && grep -q '"status": *"complete"' phase1/sessions/$s/session.json && continue
  rm -rf phase1/sessions/$s
  start=$(date -u +%FT%TZ)
  $E run --spec phase1/spec-$s.json --out phase1/sessions/$s --cpus none --allow-busy --quiet --lock /tmp/ic-frontier-phase1.lock
  rc=$?
  echo "$s started $start finished $(date -u +%FT%TZ) rc=$rc"
done
