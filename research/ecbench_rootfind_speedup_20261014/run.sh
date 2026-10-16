#!/usr/bin/env bash
# Usage: bash research/ecbench_rootfind_speedup_20261014/run.sh  (ECBENCH=path overrides)
set -euo pipefail
ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
HERE="$ROOT/research/ecbench_rootfind_speedup_20261014"
ECBENCH="${ECBENCH:-$ROOT/target/release/ecbench}"
one() {
  local name=$1 spec=$2; shift 2
  local out="$HERE/sessions/$name"
  "$ECBENCH" run --wait --spec "$spec" --out "$out" --cpus none --allow-busy
  "$ECBENCH" verify --dir "$out" --replay-all --out "$HERE/audits/audit-$name.json"
  "$ECBENCH" table --dir "$out" > "$HERE/tables/table-$name.md"
  while [ $# -gt 0 ]; do
    "$ECBENCH" compare --dir "$out" --a "$1" --b "$2" --save > "$HERE/tables/compare-$name-$2-vs-$1.txt"
    shift 2
  done
}
one claim-n13 "$ROOT/docs/ecbench/specs/claim-koblitz-rootfind-stride-n13.json" \
  ic-negfold-m2-d5-walk ic-negfold-m2-d5-stride  rho-strong ic-negfold-m2-d4-stride \
  rho-strong rho-strong-aa  rho-strong ic-negfold-m2-d5-stride
one check-n19 "$ROOT/docs/ecbench/specs/check-koblitz-rootfind-stride-n19.json" \
  ic-negfold-m2-d6-walk ic-negfold-m2-d6-stride  rho-strong ic-negfold-m2-d6-stride  rho-strong ic-negfold-m2-d6-walk
