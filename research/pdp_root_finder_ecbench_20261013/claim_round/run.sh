#!/bin/bash
# Claim round: sessions, replay-all audits, comparisons, a claim per workload.
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${ECBENCH:?set ECBENCH to the binary that measures the sessions}
mkdir -p "$here/sessions" "$here/audits" "$here/tables" "$here/claims"
for n in ${*:-23 31}; do
  out=$here/sessions/n$n
  echo "$(date -u +%FT%TZ) n=$n run" >> "$here/progress.txt"
  "$bin" run --wait --spec "$here/specs/SPEC-n$n.json" --out "$out" --cpus auto > "$here/sessions/n$n.run.log" 2>&1
  echo "$(date -u +%FT%TZ) n=$n run exit $?" >> "$here/progress.txt"
  "$bin" verify --dir "$out" --replay-all --out "$here/audits/n$n.audit.json" > "$here/audits/n$n.verify.log" 2>&1
  echo "$(date -u +%FT%TZ) n=$n audit exit $?" >> "$here/progress.txt"
  for pair in "rho-strong ic-rootfind" "rho-strong ic-mitm" "ic-mitm ic-rootfind" "rho-strong rho-strong-aa"; do
    set -- $pair
    "$bin" compare --dir "$out" --a "$1" --b "$2" --save > "$here/tables/n$n.compare.$1.vs.$2.txt" 2>&1
  done
  "$bin" table --dir "$out" > "$here/tables/n$n.table.md" 2>&1
  for w in $(python3 -c "import json,sys; [print(x['workload_id']) for x in json.load(open(sys.argv[1]))['workloads']]" "$out/plan.json"); do
    "$bin" claim build --dir "$out" --ic ic-rootfind --rho rho-strong --workload "$w" --out "$here/claims/n$n-$w.json" > "$here/claims/n$n-$w.log" 2>&1
    echo "n=$n $w build_exit=$?" >> "$here/claims/SUMMARY.txt"
    [ -f "$here/claims/n$n-$w.json" ] && "$bin" claim check --report "$here/claims/n$n-$w.json" > "$here/claims/n$n-$w.check.txt" 2>&1
  done
  echo "$(date -u +%FT%TZ) n=$n claims done" >> "$here/progress.txt"
done
echo "$(date -u +%FT%TZ) ALLDONE" >> "$here/progress.txt"
