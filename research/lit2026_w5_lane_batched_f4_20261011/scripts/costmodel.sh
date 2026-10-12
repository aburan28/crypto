#!/bin/sh
# Re-record the cost-model diagnostic into OUT/raw/diag-costmodel.jsonl.
OUT=${OUT:?set OUT}
P=/Volumes/SSD990-2/crypto/targets/lit2026-w5/bin/f4_lane_probe
cd /Volumes/SSD990-2/crypto/targets/lit2026-w5/scratch || exit 1
f="$OUT/raw/diag-costmodel.jsonl"
echo "{\"uptime\":\"$(uptime | sed 's/.*load averages*: //')\",\"utc\":\"$(date -u +%Y-%m-%dT%H:%M:%SZ)\",\"note\":\"costmodel re-recorded after the JSON-format fix\"}" > "$f"
$P costmodel --a 1 --n 17 --dim 6 >> "$f"
$P costmodel --a 1 --n 23 --dim 7 >> "$f"
$P costmodel --a 0 --n 23 --dim 8 >> "$f"
$P costmodel --a 0 --n 31 --dim 9 >> "$f"
echo costmodel-done
