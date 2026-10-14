#!/bin/bash
# After the sessions: a vs_rho claim report for every workload of every
# size (ic-rootfind against rho-strong, first measured round), the
# checker's verdict on each, and the database SQL of every session.
# From the repository root:  research/pdp_root_finder_ecbench_20261013/after_sessions.sh [n ...]
set -u
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
bin=$root/target/release/ecbench
sizes=${*:-"13 19 23 31 37"}
mkdir -p "$here/claims" "$here/db"
: > "$here/claims/SUMMARY.txt"
for n in $sizes; do
  dir=$here/sessions/n$n
  [ -f "$dir/plan.json" ] || continue
  for w in $(python3 -c "import json,sys; [print(x['workload_id']) for x in json.load(open(sys.argv[1]))['workloads']]" "$dir/plan.json"); do
    "$bin" claim build --dir "$dir" --ic ic-rootfind --rho rho-strong --workload "$w" --out "$here/claims/n$n-$w.json" --exit-code > "$here/claims/n$n-$w.log" 2>&1
    echo "n=$n $w build_exit=$? $(grep -iE 'fail|pass|verdict' "$here/claims/n$n-$w.log" | head -2 | tr '\n' ' ')" >> "$here/claims/SUMMARY.txt"
  done
  "$bin" db sql "$dir" "$here"/claims/n$n-*.json > "$here/db/n$n.sql" 2> "$here/db/n$n.sql.log"
  echo "n=$n db_sql_exit=$? lines=$(wc -l < "$here/db/n$n.sql")" >> "$here/claims/SUMMARY.txt"
done
cat "$here/claims/SUMMARY.txt"
