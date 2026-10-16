#!/usr/bin/env bash
# Extension cells (PROTOCOL.md §8 B), one random target each.
# Usage: run_ext.sh BIN OUTDIR
set -euo pipefail
BIN=$1; OUT=$2
mkdir -p "$OUT"
run() {
  local name=$1; shift
  # resumable: a cell already finished is not run again
  if grep -q " $name ok$" "$OUT/progress.txt" 2>/dev/null; then return; fi
  echo "$(date -u +%FT%TZ) $name start" >>"$OUT/progress.txt"
  RAYON_NUM_THREADS=2 nice -n 10 prlimit --as=4294967296 -- \
    "$BIN" --engine generic --ladder-t none --budget 10800 --planted 0 --random 1 "$@" \
    --out "$OUT/$name.jsonl" 2>>"$OUT/$name.log" && st=ok || st="exit $?"
  echo "$(date -u +%FT%TZ) $name $st" >>"$OUT/progress.txt"
}
K=(--kinds kummer --controls tower)
run m4N12-v40 --p 29686813949953 "${K[@]}" --m 4 --t-min 3 --max-t-m4 3
run m2N20-v48 --p proth:48 "${K[@]}" --m 2 --t-min 10 --max-t 10
run m2N20-v56 --p proth:56 "${K[@]}" --m 2 --t-min 10 --max-t 10
run m3N15-v40 --p 29686813949953 "${K[@]}" --m 3 --t-min 5 --max-t-m3 5
run m4N12-P-224 --p P-224 "${K[@]}" --m 4 --t-min 3 --max-t-m4 3 --verify
run m2N20-P-224 --p P-224 "${K[@]}" --m 2 --t-min 10 --max-t 10 --verify
run m2N20-P-192 --p P-192 --kinds dickson-torus --controls tower --m 2 --t-min 10 --max-t 10 --verify
