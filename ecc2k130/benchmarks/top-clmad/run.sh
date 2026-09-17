#!/bin/bash
# Paired RTX PRO 6000 benches: shipping product vs PACKED_TOP_CLMAD=1.
# Writes raw Modal logs and the JSON object each entry point prints.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$(cd "$(dirname "$0")" && pwd)"
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
make bench-rtx-pro6000 2>&1 | tee "$OUT/control.log"
extract_json "$OUT/control.log" "$OUT/control.json"

echo "candidate: make bench-rtx-pro6000 RTX_PRO6000_TOP_CLMAD=1" >&2
make bench-rtx-pro6000 RTX_PRO6000_TOP_CLMAD=1 2>&1 | tee "$OUT/candidate.log"
extract_json "$OUT/candidate.log" "$OUT/candidate.json"

python3 "$OUT/summarize.py"
