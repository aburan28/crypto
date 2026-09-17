#!/bin/bash
# Paired RTX PRO 6000 benches: shipping product vs PACKED_TOP_CLMAD=1.
# Writes raw Modal logs and the JSON object each entry point prints.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$(cd "$(dirname "$0")" && pwd)"
# Modal add_local_dir snapshots ecc2k130/. Writing receipts there during the
# image build aborts the run ("modified during build process").
TMP="${TOP_CLMAD_RECEIPT_DIR:-/tmp/top-clmad-receipts}"
mkdir -p "$TMP"
export PATH="${HOME}/.local/bin:${PATH}"
MODAL="${MODAL:-modal}"
cd "$ROOT"

extract_json() {
  python3 - "$1" "$2" <<'PY'
import json, sys
text = open(sys.argv[1], encoding='utf-8', errors='replace').read()
decoder = json.JSONDecoder()
found = None
i = 0
while i < len(text):
    start = text.find('{', i)
    if start < 0:
        break
    try:
        obj, end = decoder.raw_decode(text, start)
    except json.JSONDecodeError:
        i = start + 1
        continue
    if isinstance(obj, dict) and ('valid' in obj or 'rate' in obj or 'samples' in obj):
        found = obj
    i = start + end
if found is None:
    raise SystemExit('no benchmark JSON object in ' + sys.argv[1])
open(sys.argv[2], 'w', encoding='utf-8').write(json.dumps(found, indent=2) + '\n')
print('wrote', sys.argv[2], 'valid=', found.get('valid'), 'rate=', found.get('rate'))
PY
}

echo "control: make bench-rtx-pro6000" >&2
make bench-rtx-pro6000 2>&1 | tee "$TMP/control.log"
extract_json "$TMP/control.log" "$TMP/control.json"

echo "candidate: make bench-rtx-pro6000 RTX_PRO6000_TOP_CLMAD=1" >&2
make bench-rtx-pro6000 RTX_PRO6000_TOP_CLMAD=1 2>&1 | tee "$TMP/candidate.log"
extract_json "$TMP/candidate.log" "$TMP/candidate.json"

cp -f "$TMP/control.json" "$TMP/candidate.json" "$OUT/"
# Logs are large; keep the last 200 KiB in-tree as the frozen transcript.
python3 - "$TMP" "$OUT" <<'PY'
from pathlib import Path
import sys
tmp, out = Path(sys.argv[1]), Path(sys.argv[2])
for name in ("control.log", "candidate.log"):
    data = (tmp / name).read_bytes()
    (out / name).write_bytes(data[-200_000:] if len(data) > 200_000 else data)
PY
python3 "$OUT/summarize.py"
