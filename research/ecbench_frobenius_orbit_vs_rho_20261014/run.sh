#!/usr/bin/env bash
# Frozen claim run (PROTOCOL.md §4). Usage: bash research/ecbench_frobenius_orbit_vs_rho_20261014/run.sh
set -euo pipefail
ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
HERE="$ROOT/research/ecbench_frobenius_orbit_vs_rho_20261014"
ECBENCH="${ECBENCH:-$ROOT/target/release/ecbench}"
mkdir -p "$HERE/sessions" "$HERE/audits" "$HERE/tables" "$HERE/logs"
for spec in "$ROOT"/docs/ecbench/specs/claim-koblitz-frob-orbit-full-*.json; do
  name=$(basename "$spec" .json); name=${name#claim-koblitz-frob-orbit-full-}
  out="$HERE/sessions/$name"
  [ -f "$out/session.json" ] && { echo "skip $name"; continue; }
  "$ECBENCH" run --wait --spec "$spec" --out "$out" --cpus none --allow-busy > "$HERE/logs/run-$name.log" 2>&1
  "$ECBENCH" verify --dir "$out" --replay-all --out "$HERE/audits/audit-$name.json" > "$HERE/logs/verify-$name.log" 2>&1
  "$ECBENCH" table --dir "$out" > "$HERE/tables/table-$name.md"
  renew=$(jq -r '.arms[] | select(.name|test("^ic-orbit-.*-renew$")) | .name' "$spec")
  walk=$(jq -r '.arms[] | select(.name|test("^ic-orbit-.*-walk$")) | .name' "$spec")
  negfold=$(jq -r '.arms[] | select(.name|test("^ic-negfold")) | .name' "$spec")
  pairs="rho-frob $renew  rho-strong $renew  $walk $renew  rho-frob rho-frob-aa  rho-frob rho-strong"
  [ -n "$negfold" ] && pairs="$pairs  $negfold $renew  rho-frob $negfold"
  set -- $pairs
  while [ $# -gt 0 ]; do
    "$ECBENCH" compare --dir "$out" --a "$1" --b "$2" --save > "$HERE/tables/compare-$name-$2-vs-$1.txt"
    shift 2
  done
  echo "done $name"
done
echo RUN-DONE
