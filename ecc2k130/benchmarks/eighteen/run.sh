#!/bin/bash
# Best remaining priced leftovers vs the 18 B/s target.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$(cd "$(dirname "$0")" && pwd)"
EXTRACT="$ROOT/benchmarks/top-clmad/extract.py"
QUIET="$ROOT/benchmarks/top-clmad/modal-quiet"
TMP="${EIGHTEEN_RECEIPT_DIR:-/tmp/eighteen-receipts}"
mkdir -p "$TMP"
export PATH="${HOME}/.local/bin:${PATH}"
export MODAL="${MODAL:-$QUIET}"
chmod +x "$EXTRACT" "$QUIET"
cd "$ROOT"

run_arm() {
  local name="$1"
  shift
  echo "$name: $*" >&2
  make "$@" 2>"$TMP/${name}.err" | tee "$TMP/${name}.log"
  python3 "$EXTRACT" "$TMP/${name}.log" "$TMP/${name}.json"
}

run_arm table bench-rtx-pro6000 RTX_PRO6000_WALK_TABLE=1 RTX_PRO6000_WORKERS=0
run_arm pivot bench-rtx-pro6000 RTX_PRO6000_WALK_TABLE=1 RTX_PRO6000_TABLE_PIVOT_BYTES=1 RTX_PRO6000_WORKERS=0

cp -f "$TMP/table.json" "$TMP/pivot.json" "$OUT/"
python3 - "$TMP" "$OUT" <<'PY'
from pathlib import Path
import json, sys
tmp, out = Path(sys.argv[1]), Path(sys.argv[2])
rows = []
for name in ("table", "pivot"):
    obj = json.loads((tmp / f"{name}.json").read_text())
    rates = []
    for sample in obj.get("samples") or []:
        rate = float(sample["rate"])
        rates.append(round(rate / 1000.0 if rate > 1000 else rate, 6))
    if not rates and obj.get("rate"):
        rate = float(obj["rate"])
        rates.append(round(rate / 1000.0 if rate > 1000 else rate, 6))
    med = sorted(rates)[len(rates)//2]
    rows.append(dict(name=name, ratesB=rates, medianB=med, over18=med > 18.0,
                     identity=dict(packedTopClmad=obj.get("packedTopClmad"),
                                   walk=name)))
    src = tmp / f"{name}.log"
    if src.exists():
        data = src.read_bytes()
        (out / f"{name}.log").write_bytes(data[-200_000:] if len(data) > 200_000 else data)
summary = dict(unit="B scalar updates / s", targetB=18.0,
               floorB=(23.0, 25.0), shippingReferenceB=15.115792,
               class="engineering", arms=rows,
               anyOver18=any(r["over18"] for r in rows))
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
PY
