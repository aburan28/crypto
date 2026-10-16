#!/bin/sh
# Phase 0 step 3: run every selection spec into sessions/, counts only (L0).
set -u
E=${ECBENCH:-../../../target/release/ecbench}
for n in 13 17 19 23 37 43 47 41; do
  [ -d sessions/K$n ] && continue
  start=$(date -u +%FT%TZ)
  $E run --spec spec-select-K$n.json --out sessions/K$n --cpus none --allow-busy --quiet --lock /tmp/ic-frontier-select.lock
  echo "K$n started $start finished $(date -u +%FT%TZ) rc=$?"
done
