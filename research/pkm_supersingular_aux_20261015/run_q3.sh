#!/bin/bash
# Q3 of PROTOCOL.md: one ecbench session per curve, audited with every run
# replayed, then the saved comparisons, the table and the resource vectors.
# Operation counts only: --cpus none --allow-busy (L0), niced, 4 GB cap.
#   ECBENCH=path/to/ecbench research/pkm_supersingular_aux_20261015/run_q3.sh [s20 s22 s24]
set -u
here=$(cd "$(dirname "$0")" && pwd)/q3
bin=${ECBENCH:?set ECBENCH to the ecbench binary}
names=${*:-"s20 s22 s24"}
mkdir -p "$here/sessions" "$here/audits" "$here/tables"
run() { nice -n 10 prlimit --as=4294967296 -- "$bin" "$@"; }
pairs="rho-neg:ic-ss-rootfind rho-neg:ic-ord-rootfind ic-plain-ss:ic-ss-rootfind ic-plain-ord:ic-ord-rootfind ic-ord-rootfind:ic-ss-rootfind ic-plain-ord:ic-plain-ss ic-ss-mitm:ic-ss-rootfind rho-neg:rho-neg-aa"
for n in $names; do
  out=$here/sessions/$n
  echo "$(date -u +%FT%TZ) q3 $n run" >> "$here/../progress.txt"
  run run --wait --allow-busy --cpus none --spec "$here/specs/SPEC-$n.json" --out "$out" > "$here/sessions/$n.run.log" 2>&1
  echo "$(date -u +%FT%TZ) q3 $n run exit $?" >> "$here/../progress.txt"
  run verify --dir "$out" --replay-all --out "$here/audits/$n.audit.json" > "$here/audits/$n.verify.log" 2>&1
  echo "$(date -u +%FT%TZ) q3 $n audit exit $?" >> "$here/../progress.txt"
  for pair in $pairs; do
    a=${pair%%:*}; b=${pair##*:}
    run compare --dir "$out" --a "$a" --b "$b" --save > "$here/tables/$n.compare.$a.vs.$b.txt" 2>&1
  done
  run table --dir "$out" > "$here/tables/$n.table.md" 2>&1
  echo "$(date -u +%FT%TZ) q3 $n tables done" >> "$here/../progress.txt"
done
