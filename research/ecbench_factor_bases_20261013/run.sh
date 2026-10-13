#!/usr/bin/env bash
# Build is assumed: cargo build --release --bin ecbench
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$ROOT/../.." && pwd)"
ECB="${ECB:-$REPO/target/release/ecbench}"
cd "$REPO"

if [[ ! -x "$ECB" ]]; then
  echo "missing $ECB; run: cargo build --release --bin ecbench" >&2
  exit 1
fi

for spec in koblitz17 koblitz23 koblitz23-cyc prime20; do
  out="$ROOT/sessions/$spec"
  if [[ -e "$out" ]]; then
    echo "skip existing $out"
    continue
  fi
  echo "=== run $spec ==="
  "$ECB" run --wait --allow-busy --lock /tmp/crypto-bench-fb-protocol.lock --spec "$ROOT/specs/$spec.json" --out "$out" \
    | tee "$ROOT/run-$spec.log"
  "$ECB" verify --dir "$out" --replay 8 --out "$ROOT/audit-$spec.json"
  "$ECB" table --dir "$out" | tee "$ROOT/table-$spec.md"
done

python3 "$ROOT/check_fb_ids.py"
echo "all sessions complete"
