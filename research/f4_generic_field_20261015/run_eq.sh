#!/usr/bin/env bash
# Equality runs (PROTOCOL.md §3): every round-8 cell both engines run, N <= 18,
# once with the u32 engine and once with the generic engine at Fp32; and the
# Kummer p2 series once more at Fp64, Fp52 and FpWide<4> (same prime).
# Usage: run_eq.sh BIN OUTDIR
set -euo pipefail
BIN=$1; OUT=$2
mkdir -p "$OUT"
base=(--ladder-t none --planted 0 --random 2 --budget 3600 --t-min 1 --max-t 9 --max-t-m3 4 --max-t-m4 3)
series=(
  "kummer-p2|--p 3221225473 --kinds kummer --m 2,3 --controls tower,null"
  "dickson-p1|--p 2013265921 --kinds dickson --m 2,3,4 --controls tower"
  "torus-T31|--p 2046820351 --kinds dickson-torus --m 2,3,4 --controls tower,null --null-chain"
  "torus-T24|--p 14680063 --kinds dickson-torus --m 2,3,4 --controls tower,null --null-chain"
)
for s in "${series[@]}"; do
  IFS='|' read -r name args <<<"$s"
  for arm in "u32|--engine tower" "fp32|--engine generic --backend fp32"; do
    IFS='|' read -r a eng <<<"$arm"
    # shellcheck disable=SC2086
    RAYON_NUM_THREADS=2 nice -n 10 prlimit --as=4294967296 -- \
      "$BIN" $eng $args "${base[@]}" --out "$OUT/$name.$a.jsonl" 2>>"$OUT/$name.$a.log"
    echo "$(date -u +%FT%TZ) $name $a done" >>"$OUT/progress.txt"
  done
done
for b in fp64 fp52 wide4; do
  # shellcheck disable=SC2086
  RAYON_NUM_THREADS=2 nice -n 10 prlimit --as=4294967296 -- \
    "$BIN" --engine generic --backend $b --p 3221225473 --kinds kummer --m 2,3 --controls tower \
    "${base[@]}" --max-t 8 --out "$OUT/kummer-p2.$b.jsonl" 2>>"$OUT/kummer-p2.$b.log"
  echo "$(date -u +%FT%TZ) kummer-p2 $b done" >>"$OUT/progress.txt"
done
