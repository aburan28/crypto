#!/usr/bin/env bash
# Faster-loop runner. Usage: ECBENCH=path/to/ecbench bash .../run.sh <cycle>
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
ECBENCH="${ECBENCH:-$ROOT/target-faster/release/ecbench}"
CPUS="${CPUS:---cpus none --allow-busy --lock /tmp/crypto-bench-faster-1774.lock}"

run_one() {
  local name="$1" spec="$2" out="$3"
  echo "==> $name"
  rm -rf "$out"
  # shellcheck disable=SC2086
  nice -n 19 "$ECBENCH" run --wait --spec "$spec" --out "$out" $CPUS
  nice -n 19 "$ECBENCH" verify --dir "$out" --replay 8 --out "$HERE/audits/audit-$name.json"
  nice -n 19 "$ECBENCH" table --dir "$out" | tee "$HERE/tables/table-$name.md"
  python3 "$HERE/analyze_resources.py" --records "$out/records.jsonl" | tee "$HERE/tables/amortized-$name.txt"
}

case "${1:-}" in
  cycle1|xor)
    run_one cycle1 "$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-faster-xor.json" \
      "$HERE/sessions/cycle1"
    ;;
  cycle3|negfold-d4)
    run_one cycle3 "$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-faster-negfold-d4.json" \
      "$HERE/sessions/cycle3"
    ;;
  claim|claim-v3)
    name=claim-n13-v3
    spec="$ROOT/docs/ecbench/specs/claim-koblitz-rootfind-faster-n13.json"
    out="$HERE/sessions/claim-n13-v3"
    echo "==> $name"
    rm -rf "$out"
    # shellcheck disable=SC2086
    nice -n 19 "$ECBENCH" run --wait --spec "$spec" --out "$out" $CPUS
    nice -n 19 "$ECBENCH" verify --dir "$out" --replay 12 --out "$HERE/audits/audit-$name.json"
    nice -n 19 "$ECBENCH" table --dir "$out" | tee "$HERE/tables/table-$name.md"
    nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-f4-m3-d5 --b ic-negfold-m2-d5 --save
    nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-f4-m3-d5 --b ic-c2-frob-counted-m2 --save
    nice -n 19 "$ECBENCH" compare --dir "$out" --a rho-strong --b ic-negfold-m2-d5 --save
    nice -n 19 "$ECBENCH" compare --dir "$out" --a rho-strong --b ic-c2-frob-counted-m2 --save
    python3 "$HERE/analyze_resources.py" --records "$out/records.jsonl" | tee "$HERE/tables/amortized-$name.txt"
    ;;
  *)
    echo "usage: $0 cycle1|cycle3|claim"
    exit 1
    ;;
esac
