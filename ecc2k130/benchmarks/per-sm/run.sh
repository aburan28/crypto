#!/usr/bin/env bash
# One Modal allocation, several oversubscribed grids. See WAVES.md.
# Logs go to /tmp until Modal returns (image mount is this tree).
set -euo pipefail
cd "$(dirname "$0")/../.."
GPU="${WAVES_GPU:-RTX-PRO-6000}"
LIST="${WAVES_LIST:-1,4,6,8}"
OUTDIR="${OUTDIR:-benchmarks/per-sm}"
TMP="${WAVES_TMP:-/tmp/ecc2k130-per-sm}"
slug=$(printf '%s' "$GPU" | tr 'A-Z' 'a-z' | tr -c 'a-z0-9' '-')
slug=${slug%-}
mkdir -p "$TMP/$slug" "$OUTDIR"
echo "=== waves $LIST on $GPU ($(date -u +%Y-%m-%dT%H:%M:%SZ)) ==="
make bench-waves-modal WAVES_GPU="$GPU" WAVES_LIST="$LIST" 2>&1 | tee "$TMP/$slug/waves.log"
python3 benchmarks/per-sm/freeze.py "$GPU" "$TMP/$slug/waves.log" "$TMP/$slug/waves.json"
cp -f "$TMP/$slug/waves.log" "$TMP/$slug/waves.json" "$OUTDIR/${slug}-waves.log" "$OUTDIR/${slug}-waves.json"
echo "receipt $OUTDIR/${slug}-waves.json"
