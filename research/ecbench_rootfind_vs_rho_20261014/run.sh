#!/usr/bin/env bash
# vs-rho chase runner. Usage: ECBENCH=path/to/ecbench bash .../run.sh <cycle>
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
ECBENCH="${ECBENCH:-$ROOT/target-faster/release/ecbench}"
CPUS="${CPUS:---cpus none --allow-busy --lock /tmp/crypto-bench-vsrho-1774.lock}"

run_one() {
  local name="$1" spec="$2" out="$3"
  echo "==> $name"
  rm -rf "$out"
  # shellcheck disable=SC2086
  nice -n 19 "$ECBENCH" run --wait --spec "$spec" --out "$out" $CPUS
  nice -n 19 "$ECBENCH" verify --dir "$out" --replay 8 --out "$HERE/audits/audit-$name.json"
  nice -n 19 "$ECBENCH" table --dir "$out" | tee "$HERE/tables/table-$name.md"
  if [[ -f "$ROOT/research/ecbench_rootfind_faster_20261014/analyze_resources.py" ]]; then
    python3 "$ROOT/research/ecbench_rootfind_faster_20261014/analyze_resources.py" \
      --records "$out/records.jsonl" | tee "$HERE/tables/amortized-$name.txt"
  fi
}

case "${1:-}" in
  cycle4)
    run_one cycle4 "$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-vsrho-cycle4.json" \
      "$HERE/sessions/cycle4"
    ;;
  cycle5)
    run_one cycle5 "$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-vsrho-cycle5.json" \
      "$HERE/sessions/cycle5"
    ;;
  *)
    echo "usage: $0 cycle4|cycle5"
    exit 1
    ;;
esac
