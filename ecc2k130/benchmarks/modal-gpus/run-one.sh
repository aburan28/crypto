#!/usr/bin/env bash
# Rent one Modal GPU, time the shipping packed walk, freeze the receipt.
# See SURVEY.md.
#
# Logs go to /tmp until every Modal call has returned. The image mount is
# the ecc2k130 tree; writing into it while `modal run` is uploading raises
# ExecutionError.
set -euo pipefail
cd "$(dirname "$0")/../.."
GPU="${1:?usage: run-one.sh MODAL-TYPE}"
OUTDIR="${OUTDIR:-benchmarks/modal-gpus}"
TMP="${SURVEY_TMP:-/tmp/ecc2k130-survey}"
CATALOG="$OUTDIR/catalog.json"

slug=$(python3 -c "from pathlib import Path; import sys; sys.path.insert(0, 'benchmarks/modal-gpus'); from freeze import slug; print(slug(sys.argv[1]))" "$GPU")
clmad=$(python3 -c "
import json, sys
from pathlib import Path
cat = json.loads(Path(sys.argv[1]).read_text())
want = sys.argv[2]
for g in cat['gpus']:
    if g['modal'] == want:
        print(int(g['clmad']))
        break
else:
    raise SystemExit('unknown GPU ' + want)
" "$CATALOG" "$GPU")

mkdir -p "$TMP/$slug" "$OUTDIR"
echo "=== bench CLMAD=$clmad on $GPU ($(date -u +%Y-%m-%dT%H:%M:%SZ)) ==="
set +e
make bench-modal-gpu SURVEY_GPU="$GPU" SURVEY_CLMAD="$clmad" 2>&1 | tee "$TMP/$slug/bench.log"
rc=${PIPESTATUS[0]}
set -e
if [ "$rc" -ne 0 ]; then
  python3 - "$GPU" "$TMP/$slug/bench.log" "$OUTDIR/${slug}.json" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, "benchmarks/modal-gpus")
from freeze import write_unavailable
log = Path(sys.argv[2]).read_text(errors="replace")[-4000:]
write_unavailable(sys.argv[1], "modal run failed:\n" + log, Path(sys.argv[3]))
print("unavailable", sys.argv[3])
PY
  cp -f "$TMP/$slug/bench.log" "$OUTDIR/${slug}.log"
  exit "$rc"
fi
if ! python3 benchmarks/modal-gpus/freeze.py "$GPU" "$TMP/$slug/bench.log" "$TMP/$slug/receipt.json"; then
  python3 - "$GPU" "$TMP/$slug/bench.log" "$OUTDIR/${slug}.json" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, "benchmarks/modal-gpus")
from freeze import write_unavailable
write_unavailable(sys.argv[1], "freeze refused the log (identity or repeats)", Path(sys.argv[3]))
print("unavailable", sys.argv[3])
PY
  cp -f "$TMP/$slug/bench.log" "$OUTDIR/${slug}.log"
  exit 1
fi
cp -f "$TMP/$slug/bench.log" "$OUTDIR/${slug}.log"
cp -f "$TMP/$slug/receipt.json" "$OUTDIR/${slug}.json"
echo "receipt $OUTDIR/${slug}.json"
