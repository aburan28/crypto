#!/usr/bin/env bash
# Rent one Modal B200, prove the packed client on it, then time the shipping
# (or software) walk. See ../../B200.md.
#
# Logs go to /tmp until every Modal call has returned. The image mount is
# the ecc2k130 tree; writing into it while `modal run` is uploading raises
# ExecutionError (the first attempt died that way on validate.log).
set -euo pipefail
cd "$(dirname "$0")/../.."
OUTDIR="${OUTDIR:-benchmarks/b200}"
CLMAD="${B200_CLMAD:-1}"
SKIP_VALIDATE="${SKIP_VALIDATE:-0}"
TMP="${B200_TMP:-/tmp/ecc2k130-b200}"
mkdir -p "$TMP" "$OUTDIR"

if [ "$CLMAD" = 1 ]; then
  stem=clmad
else
  stem=software
fi

if [ "$SKIP_VALIDATE" != 1 ] && [ "$CLMAD" = 1 ]; then
  echo "=== validate on B200 ($(date -u +%Y-%m-%dT%H:%M:%SZ)) ==="
  make validate-b200-modal 2>&1 | tee "$TMP/validate.log"
fi

echo "=== bench CLMAD=$CLMAD on B200 ($(date -u +%Y-%m-%dT%H:%M:%SZ)) ==="
make bench-b200-modal B200_CLMAD="$CLMAD" 2>&1 | tee "$TMP/${stem}.log"
python3 benchmarks/b200/freeze.py "$stem" "$TMP/${stem}.log" "$TMP/${stem}.json"

cp -f "$TMP/${stem}.log" "$TMP/${stem}.json" "$OUTDIR/"
if [ -f "$TMP/validate.log" ]; then
  cp -f "$TMP/validate.log" "$OUTDIR/"
fi
echo "receipt $OUTDIR/${stem}.json"
