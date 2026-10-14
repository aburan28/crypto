#!/usr/bin/env bash
# Parity-loop runner. Usage: ECBENCH=path/to/ecbench bash .../run.sh [cycle1|cycle2|claim|all]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
ECBENCH="${ECBENCH:-$ROOT/target-parity/release/ecbench}"
CPUS="${CPUS:---cpus none --allow-busy}"

run_one() {
  local name="$1" spec="$2" out="$3"
  echo "==> $name"
  rm -rf "$out"
  # shellcheck disable=SC2086
  nice -n 19 "$ECBENCH" run --wait --spec "$spec" --out "$out" $CPUS
  nice -n 19 "$ECBENCH" verify --dir "$out" --replay 8 --out "$HERE/audits/audit-$name.json"
  nice -n 19 "$ECBENCH" table --dir "$out" | tee "$HERE/tables/table-$name.md"
  nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-f4-m3-d5 --b ic-negfold-m2-d5 --save 2>/dev/null || true
  python3 "$HERE/analyze_resources.py" --records "$out/records.jsonl" | tee "$HERE/tables/amortized-$name.txt"
}

case "${1:-all}" in
  cycle1)
    run_one cycle1 "$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-negfold.json" \
      "$HERE/sessions/cycle1"
    ;;
  *)
    echo "usage: $0 cycle1|..."
    exit 1
    ;;
esac
