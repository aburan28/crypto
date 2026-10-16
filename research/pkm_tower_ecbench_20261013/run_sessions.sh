#!/bin/bash
# The measured sessions of research/pkm_tower_ecbench_20261013 (PROTOCOL.md):
# one per size, each audited with every run replayed, then the saved
# comparisons, the one-unit table and the resource vectors.
# From the repository root:  ECBENCH=path/to/ecbench research/pkm_tower_ecbench_20261013/run_sessions.sh [bits ...]
set -u
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
bin=${ECBENCH:-$root/target/release/ecbench}
sizes=${*:-"20 22 24 26 28"}
mkdir -p "$here/sessions" "$here/audits" "$here/tables"
for b in $sizes; do
  out=$here/sessions/p$b
  echo "$(date -u +%FT%TZ) p$b run" >> "$here/progress.txt"
  "$bin" run --wait --spec "$here/specs/SPEC-p$b.json" --out "$out" --cpus auto > "$here/sessions/p$b.run.log" 2>&1
  echo "$(date -u +%FT%TZ) p$b run exit $?" >> "$here/progress.txt"
  "$bin" verify --dir "$out" --replay-all --out "$here/audits/p$b.audit.json" > "$here/audits/p$b.verify.log" 2>&1
  echo "$(date -u +%FT%TZ) p$b audit exit $?" >> "$here/progress.txt"
  for pair in "rho-neg ic-rootfind" "rho-neg ic-mitm" "rho-neg ic-rootfind-pairs" "rho-neg ic-abscissa" "ic-mitm ic-rootfind" "ic-abscissa ic-rootfind" "rho-neg rho-neg-aa"; do
    set -- $pair
    "$bin" compare --dir "$out" --a "$1" --b "$2" --save > "$here/tables/p$b.compare.$1.vs.$2.txt" 2>&1 || rm -f "$here/tables/p$b.compare.$1.vs.$2.txt"
  done
  "$bin" table --dir "$out" > "$here/tables/p$b.table.md" 2>&1
  "$bin" resources --dir "$out" --out "$here/tables/p$b.resources.json" > /dev/null 2>&1
  echo "$(date -u +%FT%TZ) p$b tables done" >> "$here/progress.txt"
done
echo "$(date -u +%FT%TZ) ALLDONE" >> "$here/progress.txt"
