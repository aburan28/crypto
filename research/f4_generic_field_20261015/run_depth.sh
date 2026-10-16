#!/usr/bin/env bash
# Measured cells (PROTOCOL.md §5): Kummer at v2(p - 1) = 40, 48, 56 (u64
# builder, backend auto) and at NIST P-224 (wide builder, Fp224); the Dickson
# torus at P-192 (FpP192) and P-256 (FpWide<4>). Two random targets per cell,
# plus one planted target at t <= 4 in a separate small run.
# Usage: run_depth.sh BIN OUTDIR
set -euo pipefail
BIN=$1; OUT=$2
mkdir -p "$OUT"
run() { # name args...
  local name=$1; shift
  echo "$(date -u +%FT%TZ) $name start" >>"$OUT/progress.txt"
  RAYON_NUM_THREADS=2 nice -n 10 prlimit --as=4294967296 -- \
    "$BIN" --engine generic --ladder-t none --budget 5400 "$@" \
    --out "$OUT/$name.jsonl" 2>>"$OUT/$name.log" && st=ok || st="exit $?"
  echo "$(date -u +%FT%TZ) $name $st" >>"$OUT/progress.txt"
}
K=(--kinds kummer --controls tower --t-min 1)
# small planted runs (planted_ok and the verifier's planted check)
for p in 29686813949953 proth:48 proth:56; do
  run "planted-${p//:/-}" --p $p "${K[@]}" --m 2,3 --max-t 4 --max-t-m3 3 --planted 1 --random 0
done
run planted-P-224 --p P-224 "${K[@]}" --m 2,3 --max-t 4 --max-t-m3 3 --planted 1 --random 0 --verify
run planted-P-192 --p P-192 --kinds dickson-torus --controls tower --t-min 1 --m 2,3 --max-t 4 --max-t-m3 3 --planted 1 --random 0 --verify
# the measured cells
for p in 29686813949953 proth:48 proth:56; do
  run "kummer-${p//:/-}" --p $p "${K[@]}" --m 2,3,4 --max-t 9 --max-t-m3 4 --max-t-m4 2 --planted 0 --random 2
done
run kummer-P-224 --p P-224 "${K[@]}" --m 2,3,4 --max-t 9 --max-t-m3 4 --max-t-m4 2 --planted 0 --random 2 --verify
run torus-P-192 --p P-192 --kinds dickson-torus --controls tower,null --t-min 1 --m 2,3 --max-t 9 --max-t-m3 4 --planted 0 --random 2 --verify
run torus-P-256 --p P-256 --kinds dickson-torus --controls tower --t-min 1 --m 2,3 --max-t 8 --max-t-m3 3 --planted 0 --random 2 --verify
# optional (PROTOCOL §5.3): m = 2, N = 20 at v2 = 40, one target
run optional-N20-v40 --p 29686813949953 "${K[@]}" --m 2 --t-min 10 --max-t 10 --planted 0 --random 1
