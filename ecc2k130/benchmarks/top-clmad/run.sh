#!/bin/bash
# Paired RTX PRO 6000 benches: shipping product vs PACKED_TOP_CLMAD=1.
# Writes raw Modal logs and a recovered JSON receipt for each arm.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$(cd "$(dirname "$0")" && pwd)"
# Modal add_local_dir snapshots ecc2k130/. Writing receipts there during the
# image build aborts the run ("modified during build process").
TMP="${TOP_CLMAD_RECEIPT_DIR:-/tmp/top-clmad-receipts}"
mkdir -p "$TMP"
export PATH="${HOME}/.local/bin:${PATH}"
chmod +x "$OUT/modal-quiet" "$OUT/extract.py" "$OUT/summarize.py"
export MODAL="${MODAL:-$OUT/modal-quiet}"
cd "$ROOT"

run_arm() {
  local name="$1"
  shift
  echo "$name: $*" >&2
  make "$@" 2>"$TMP/${name}.err" | tee "$TMP/${name}.log"
  python3 "$OUT/extract.py" "$TMP/${name}.log" "$TMP/${name}.json"
}

if [ "${SKIP_CONTROL:-0}" != "1" ]; then
  run_arm control bench-rtx-pro6000
else
  echo "SKIP_CONTROL=1; expecting $TMP/control.json" >&2
  test -f "$TMP/control.json"
fi

if [ "${SKIP_CANDIDATE:-0}" != "1" ]; then
  run_arm candidate bench-rtx-pro6000 RTX_PRO6000_TOP_CLMAD=1
fi

cp -f "$TMP/control.json" "$TMP/candidate.json" "$OUT/"
python3 - "$TMP" "$OUT" <<'PY'
from pathlib import Path
import sys
tmp, out = Path(sys.argv[1]), Path(sys.argv[2])
for name in ("control.log", "candidate.log"):
    src = tmp / name
    if not src.exists():
        continue
    data = src.read_bytes()
    (out / name).write_bytes(data[-200_000:] if len(data) > 200_000 else data)
PY
python3 "$OUT/summarize.py"
