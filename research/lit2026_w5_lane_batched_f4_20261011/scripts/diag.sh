#!/bin/sh
# lit2026 W5 diagnostics with f4_lane_probe on the four production cells:
# costmodel (per-degree root Macaulay shapes), lockstep (node sharing across
# targets), profile (root builds against specialisations) and lanes (root
# phase, k = 2 … 64 targets, sequential echelon root against lane-batched).
#   OUT  output directory (raw/ is created under it)
OUT=${OUT:?set OUT}
P=/Volumes/SSD990-2/crypto/targets/lit2026-w5/bin/f4_lane_probe
mkdir -p "$OUT/raw"
cd /Volumes/SSD990-2/crypto/targets/lit2026-w5/scratch || exit 1
up() { echo "{\"uptime\":\"$(uptime | sed 's/.*load averages*: //')\",\"utc\":\"$(date -u +%Y-%m-%dT%H:%M:%SZ)\"}"; }
for mode in costmodel lockstep profile; do
  f="$OUT/raw/diag-$mode.jsonl"
  up > "$f"
  $P $mode --a 1 --n 17 --dim 6 --targets 8 >> "$f"
  $P $mode --a 1 --n 23 --dim 7 --targets 4 >> "$f"
  $P $mode --a 0 --n 23 --dim 8 --targets 4 >> "$f"
  $P $mode --a 0 --n 31 --dim 9 --targets 1 --budget 2000 >> "$f"
  up >> "$f"
done
# lockstep on K0/2^31 needs more than its one production target
f="$OUT/raw/diag-lockstep-k0-31-4targets.jsonl"
up > "$f"
$P lockstep --a 0 --n 31 --dim 9 --targets 4 --budget 2000 >> "$f"
for root in echelon rref; do
  f="$OUT/raw/diag-lanes-$root.jsonl"
  up > "$f"
  if [ "$root" = rref ]; then export KIC_F4_INHERIT_ROOT=rref; fi
  $P lanes --a 1 --n 17 --dim 6 --ks 2,4,8,16,32,64 --reps 5 >> "$f"
  $P lanes --a 1 --n 23 --dim 7 --ks 2,4,8,16,32,64 --reps 5 >> "$f"
  $P lanes --a 0 --n 23 --dim 8 --ks 2,4,8,16,32,64 --reps 5 >> "$f"
  $P lanes --a 0 --n 31 --dim 9 --ks 2,4,8,16,32,64 --reps 5 >> "$f"
  up >> "$f"
done
echo diag-done
