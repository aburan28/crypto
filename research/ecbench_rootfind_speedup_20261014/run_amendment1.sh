#!/usr/bin/env bash
# Amendment 1. Usage: bash research/ecbench_rootfind_speedup_20261014/run_amendment1.sh
set -euo pipefail
ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
HERE="$ROOT/research/ecbench_rootfind_speedup_20261014"
ECBENCH="${ECBENCH:-$ROOT/target/release/ecbench}"
# Sessions A and B must replay unchanged under this binary.
"$ECBENCH" verify --dir "$HERE/sessions/claim-n13" --replay-all --out "$HERE/audits/reaudit-claim-n13.json"
"$ECBENCH" verify --dir "$HERE/sessions/check-n19" --replay-all --out "$HERE/audits/reaudit-check-n19.json"
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
one claim-renew-n13 "$ROOT/docs/ecbench/specs/claim-koblitz-rootfind-stride-renew-n13.json" \
  ic-negfold-m2-d5-walk ic-negfold-m2-d5-renew  rho-strong ic-negfold-m2-d4-renew \
  ic-negfold-m2-d5-stride ic-negfold-m2-d5-renew  rho-strong rho-strong-aa  rho-strong ic-negfold-m2-d5-renew
one claim-renew-n19 "$ROOT/docs/ecbench/specs/claim-koblitz-rootfind-stride-renew-n19.json" \
  ic-negfold-m2-d6-walk ic-negfold-m2-d6-renew  ic-negfold-m2-d6-stride ic-negfold-m2-d6-renew \
  rho-strong ic-negfold-m2-d6-renew  rho-strong ic-negfold-m2-d6-walk
echo AMENDMENT1-DONE
