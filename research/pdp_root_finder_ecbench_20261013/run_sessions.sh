#!/bin/bash
# The measured sessions of research/pdp_root_finder_ecbench_20261013
# (PROTOCOL.md), one per size, each audited with every run replayed, then
# the saved comparisons, the one-unit table and the resource vectors.
# From the repository root:  research/pdp_root_finder_ecbench_20261013/run_sessions.sh [n ...]
set -u
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
bin=$root/target/release/ecbench
sizes=${*:-"13 19 23 31 37"}
mkdir -p "$here/sessions" "$here/audits" "$here/tables"
for n in $sizes; do
  out=$here/sessions/n$n
  echo "$(date -u +%FT%TZ) n=$n run" >> "$here/progress.txt"
  "$bin" run --wait --spec "$here/specs/SPEC-n$n.json" --out "$out" --cpus auto > "$here/sessions/n$n.run.log" 2>&1
  echo "$(date -u +%FT%TZ) n=$n run exit $?" >> "$here/progress.txt"
  "$bin" verify --dir "$out" --replay-all --out "$here/audits/n$n.audit.json" > "$here/audits/n$n.verify.log" 2>&1
  echo "$(date -u +%FT%TZ) n=$n audit exit $?" >> "$here/progress.txt"
  for pair in "rho-strong ic-rootfind" "rho-strong ic-mitm" "rho-strong ic-f4" "rho-strong ic-rootfind-pairs" "ic-mitm ic-rootfind" "ic-f4 ic-rootfind" "rho-strong rho-strong-aa"; do
    set -- $pair
    "$bin" compare --dir "$out" --a "$1" --b "$2" --save > "$here/tables/n$n.compare.$1.vs.$2.txt" 2>&1 || rm -f "$here/tables/n$n.compare.$1.vs.$2.txt"
  done
  "$bin" table --dir "$out" > "$here/tables/n$n.table.md" 2>&1
  "$bin" resources --dir "$out" --out "$here/tables/n$n.resources.json" > /dev/null 2>&1
  echo "$(date -u +%FT%TZ) n=$n tables done" >> "$here/progress.txt"
done
echo "$(date -u +%FT%TZ) ALLDONE" >> "$here/progress.txt"
