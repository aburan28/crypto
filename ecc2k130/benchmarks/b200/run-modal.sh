#!/usr/bin/env bash
# Rent one Modal B200, prove the packed client on it, then time the shipping
# (or software) walk. See ../../B200.md.
set -euo pipefail
cd "$(dirname "$0")/../.."
OUTDIR="${OUTDIR:-benchmarks/b200}"
CLMAD="${B200_CLMAD:-1}"
SKIP_VALIDATE="${SKIP_VALIDATE:-0}"
mkdir -p "$OUTDIR"

if [ "$CLMAD" = 1 ]; then
  stem=clmad
else
  stem=software
fi
log="$OUTDIR/${stem}.log"

if [ "$SKIP_VALIDATE" != 1 ] && [ "$CLMAD" = 1 ]; then
  echo "=== validate on B200 ($(date -u +%Y-%m-%dT%H:%M:%SZ)) ==="
  make validate-b200-modal 2>&1 | tee "$OUTDIR/validate.log"
fi

echo "=== bench CLMAD=$CLMAD on B200 ($(date -u +%Y-%m-%dT%H:%M:%SZ)) ==="
make bench-b200-modal B200_CLMAD="$CLMAD" 2>&1 | tee "$log"
python3 benchmarks/b200/freeze.py "$stem" "$log" "$OUTDIR/${stem}.json"
echo "receipt $OUTDIR/${stem}.json"
