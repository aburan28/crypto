#!/bin/bash
# The measured sessions of research/pkm_isogeny_tower_ecbench_20261015
# (PROTOCOL.md): one session per curve, each audited with every run
# replayed, then the saved comparisons, the one-unit table and the
# resource vectors.  Operation counts only: --cpus none --allow-busy (L0),
# one thread, niced, address space capped at 2.5 GB.
# From the repository root:
#   ECBENCH=path/to/ecbench research/pkm_isogeny_tower_ecbench_20261015/run_sessions.sh [name ...]
set -u
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
bin=${ECBENCH:-$root/target/release/ecbench}
names=${*:-"p20 p22 p24 p26 p28 q20 q22 q24"}
mkdir -p "$here/sessions" "$here/audits" "$here/tables"
run() { nice -n 10 prlimit --as=2684354560 -- "$bin" "$@"; }
for n in $names; do
  out=$here/sessions/$n
  echo "$(date -u +%FT%TZ) $n run" >> "$here/progress.txt"
  run run --wait --allow-busy --cpus none --spec "$here/specs/SPEC-$n.json" --out "$out" > "$here/sessions/$n.run.log" 2>&1
  echo "$(date -u +%FT%TZ) $n run exit $?" >> "$here/progress.txt"
  run verify --dir "$out" --replay-all --out "$here/audits/$n.audit.json" > "$here/audits/$n.verify.log" 2>&1
  echo "$(date -u +%FT%TZ) $n audit exit $?" >> "$here/progress.txt"
  case $n in
    p*) pairs="rho-neg:ic-iso-rootfind rho-neg:ic-iso-mitm rho-neg:ic-kummer-rootfind rho-neg:ic-plain ic-plain:ic-iso-rootfind ic-kummer-rootfind:ic-iso-rootfind ic-iso-mitm:ic-iso-rootfind rho-neg:rho-neg-aa" ;;
    *)  pairs="rho-neg:ic-iso-rootfind rho-neg:ic-iso-mitm rho-neg:ic-plain ic-plain:ic-iso-rootfind ic-iso-mitm:ic-iso-rootfind rho-neg:rho-neg-aa" ;;
  esac
  for pair in $pairs; do
    a=${pair%%:*}; b=${pair##*:}
    run compare --dir "$out" --a "$a" --b "$b" --save > "$here/tables/$n.compare.$a.vs.$b.txt" 2>&1 || rm -f "$here/tables/$n.compare.$a.vs.$b.txt"
  done
  run table --dir "$out" > "$here/tables/$n.table.md" 2>&1
  run resources --dir "$out" --out "$here/tables/$n.resources.json" > /dev/null 2>&1
  echo "$(date -u +%FT%TZ) $n tables done" >> "$here/progress.txt"
done
echo "$(date -u +%FT%TZ) ALLDONE" >> "$here/progress.txt"
