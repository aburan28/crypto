#!/usr/bin/env bash
# Sequential packed-walk benches on every catalog GPU that is a pinned SKU.
# Failures are frozen as unavailable.json-style receipts; the loop continues.
set -u
cd "$(dirname "$0")/../.."
CATALOG="${CATALOG:-benchmarks/modal-gpus/catalog.json}"
mapfile -t GPUS < <(python3 -c "
import json
from pathlib import Path
cat = json.loads(Path('$CATALOG').read_text())
for g in cat['gpus']:
    if g.get('run'):
        print(g['modal'])
")
status=0
for gpu in "${GPUS[@]}"; do
  echo "######## $gpu ########"
  if ! bash benchmarks/modal-gpus/run-one.sh "$gpu"; then
    echo "FAILED $gpu (recorded unavailable, continuing)"
    status=1
  fi
done
python3 benchmarks/modal-gpus/summarize.py
exit $status
