#!/usr/bin/env bash
# Improve-loop runner for research/ecbench_rootfind_improve_20261013.
# Usage: ECBENCH=path/to/ecbench bash research/ecbench_rootfind_improve_20261013/run.sh [select|claim|all]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ECBENCH="${ECBENCH:-$ROOT/target/release/ecbench}"
OUT="$(cd "$(dirname "$0")" && pwd)"
MODE="${1:-all}"
CPUS="${ECBENCH_CPUS:---cpus none --allow-busy}"

run_one() {
  local name="$1" spec="$2"
  local sess="$OUT/sessions/$name"
  if [[ -e "$sess" ]]; then
    echo "skip $name: $sess exists" >&2
    return 0
  fi
  echo "=== run $name ===" >&2
  # shellcheck disable=SC2086
  nice -n 10 "$ECBENCH" run --wait $CPUS --spec "$spec" --out "$sess" \
    2>&1 | tee "$OUT/logs/run-$name.log"
  nice -n 10 "$ECBENCH" verify --dir "$sess" --replay-all --exit-code \
    --out "$OUT/audits/audit-$name.json" \
    2>&1 | tee "$OUT/logs/audit-$name.log"
  nice -n 10 "$ECBENCH" table --dir "$sess" \
    >"$OUT/tables/table-$name.md"
  nice -n 10 "$ECBENCH" resources --dir "$sess" \
    >"$OUT/comparisons/$name.resources.json"
}

mkdir -p "$OUT"/{sessions,tables,audits,comparisons,logs,cycles}

case "$MODE" in
  select)
    run_one select "$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-select.json"
    ;;
  claim)
    run_one claim-n13 "$ROOT/docs/ecbench/specs/claim-koblitz-rootfind-improve-n13.json"
    ;;
  all)
    run_one select "$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-select.json"
    run_one claim-n13 "$ROOT/docs/ecbench/specs/claim-koblitz-rootfind-improve-n13.json"
    ;;
  *)
    echo "usage: $0 [select|claim|all]" >&2
    exit 2
    ;;
esac
