#!/usr/bin/env bash
# Parity-loop runner. Usage: ECBENCH=path/to/ecbench bash .../run.sh [cycle1|cycle2|claim|all]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
ECBENCH="${ECBENCH:-$ROOT/target-parity/release/ecbench}"
# Private lock so L0 engineering smokes do not queue behind unrelated benches.
CPUS="${CPUS:---cpus none --allow-busy --lock /tmp/crypto-bench-parity-1774.lock}"

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

run_cycle2() {
  local name=cycle2
  local spec="$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-window-div.json"
  local out="$HERE/sessions/cycle2"
  echo "==> $name"
  rm -rf "$out"
  # shellcheck disable=SC2086
  nice -n 19 "$ECBENCH" run --wait --spec "$spec" --out "$out" $CPUS
  nice -n 19 "$ECBENCH" verify --dir "$out" --replay 8 --out "$HERE/audits/audit-$name.json"
  nice -n 19 "$ECBENCH" table --dir "$out" | tee "$HERE/tables/table-$name.md"
  nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-f4-m3-d5 --b ic-negfold-m2-d5 --save 2>/dev/null || true
  nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-negfold-m2-d5 --b ic-win-neg-d2 --save 2>/dev/null || true
  python3 "$HERE/analyze_resources.py" --records "$out/records.jsonl" | tee "$HERE/tables/amortized-$name.txt"
}

run_cycle4() {
  local name=cycle4
  local spec="$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-frob-orbit.json"
  local out="$HERE/sessions/cycle4"
  echo "==> $name"
  rm -rf "$out"
  # shellcheck disable=SC2086
  nice -n 19 "$ECBENCH" run --wait --spec "$spec" --out "$out" $CPUS
  nice -n 19 "$ECBENCH" verify --dir "$out" --replay 8 --out "$HERE/audits/audit-$name.json"
  nice -n 19 "$ECBENCH" table --dir "$out" | tee "$HERE/tables/table-$name.md"
  nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-f4-m3-d5 --b ic-frob-orbit-m2 --save 2>/dev/null || true
  nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-negfold-m2-d5 --b ic-frob-orbit-m2 --save 2>/dev/null || true
  nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-f4-m3-d5 --b ic-subtract-m2-d5 --save 2>/dev/null || true
  python3 "$HERE/analyze_resources.py" --records "$out/records.jsonl" | tee "$HERE/tables/amortized-$name.txt"
}

case "${1:-all}" in
  cycle1)
    run_one cycle1 "$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-negfold.json" \
      "$HERE/sessions/cycle1"
    ;;
  cycle2)
    run_cycle2
    ;;
  cycle4)
    run_cycle4
    ;;
  cycle4b)
    run_one cycle4b "$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-compact-frob.json" \
      "$HERE/sessions/cycle4b"
    ;;
  cycle5)
    run_one cycle5 "$ROOT/docs/ecbench/specs/smoke-koblitz-rootfind-compact-counted.json" \
      "$HERE/sessions/cycle5"
    ;;
  claim|cycle3)
    name=claim-n13
    spec="$ROOT/docs/ecbench/specs/claim-koblitz-rootfind-parity-n13.json"
    out="$HERE/sessions/claim-n13"
    echo "==> $name"
    rm -rf "$out"
    # shellcheck disable=SC2086
    nice -n 19 "$ECBENCH" run --wait --spec "$spec" --out "$out" $CPUS
    nice -n 19 "$ECBENCH" verify --dir "$out" --replay 12 --out "$HERE/audits/audit-$name.json"
    nice -n 19 "$ECBENCH" table --dir "$out" | tee "$HERE/tables/table-$name.md"
    nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-f4-m3-d5 --b ic-negfold-m2-d5 --save
    nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-rootfind-m2-d5 --b ic-negfold-m2-d5 --save
    nice -n 19 "$ECBENCH" compare --dir "$out" --a rho-strong --b ic-negfold-m2-d5 --save
    python3 "$HERE/analyze_resources.py" --records "$out/records.jsonl" | tee "$HERE/tables/amortized-$name.txt"
    ;;
  claim-v2|cycle6)
    name=claim-n13-v2
    spec="$ROOT/docs/ecbench/specs/claim-koblitz-rootfind-parity-n13-v2.json"
    out="$HERE/sessions/claim-n13-v2"
    echo "==> $name"
    rm -rf "$out"
    # shellcheck disable=SC2086
    nice -n 19 "$ECBENCH" run --wait --spec "$spec" --out "$out" $CPUS
    nice -n 19 "$ECBENCH" verify --dir "$out" --replay 12 --out "$HERE/audits/audit-$name.json"
    nice -n 19 "$ECBENCH" table --dir "$out" | tee "$HERE/tables/table-$name.md"
    nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-f4-m3-d5 --b ic-c2-frob-counted-m2 --save
    nice -n 19 "$ECBENCH" compare --dir "$out" --a ic-negfold-m2-d5 --b ic-c2-frob-counted-m2 --save
    nice -n 19 "$ECBENCH" compare --dir "$out" --a rho-strong --b ic-c2-frob-counted-m2 --save
    python3 "$HERE/analyze_resources.py" --records "$out/records.jsonl" | tee "$HERE/tables/amortized-$name.txt"
    python3 "$HERE/analyze_f4_sensitivity.py" --records "$out/records.jsonl" | tee "$HERE/tables/f4-sensitivity-$name.txt"
    ;;
  *)
    echo "usage: $0 cycle1|cycle2|cycle4|cycle4b|cycle5|claim|claim-v2|..."
    exit 1
    ;;
esac
