#!/usr/bin/env bash
# Pre-registered speed arms (PROTOCOL.md §4): interleaved, one core, one thread.
# Usage: run_speed.sh BIN OUTDIR REPS
set -euo pipefail
BIN=$1; OUT=$2; REPS=${3:-7}
mkdir -p "$OUT"
common=(--planted 0 --random 1 --targets 0 --ladder-t none --budget 3600)
# cell id | arm | args
cells=(
  "C1|u32|--engine tower --p 2013265921 --kinds dickson --m 2 --controls tower --t-min 8 --max-t 8"
  "C1|fp32|--engine generic --backend fp32 --p 2013265921 --kinds dickson --m 2 --controls tower --t-min 8 --max-t 8"
  "C2|u32|--engine tower --p 2046820351 --kinds dickson-torus --m 2 --controls tower --t-min 8 --max-t 8"
  "C2|fp32|--engine generic --backend fp32 --p 2046820351 --kinds dickson-torus --m 2 --controls tower --t-min 8 --max-t 8"
  "C3|u32|--engine tower --p 3221225473 --kinds kummer --m 2 --controls tower --t-min 8 --max-t 8"
  "C3|fp32|--engine generic --backend fp32 --p 3221225473 --kinds kummer --m 2 --controls tower --t-min 8 --max-t 8"
  "C4|u32|--engine tower --p 3221225473 --kinds kummer --m 3 --controls tower --t-min 4 --max-t-m3 4"
  "C4|fp32|--engine generic --backend fp32 --p 3221225473 --kinds kummer --m 3 --controls tower --t-min 4 --max-t-m3 4"
  "C5|fp64|--engine generic --backend fp64 --p 29686813949953 --kinds kummer --m 2 --controls tower --t-min 8 --max-t 8"
  "C5|fp52|--engine generic --backend fp52 --p 29686813949953 --kinds kummer --m 2 --controls tower --t-min 8 --max-t 8"
  "C6|fp64|--engine generic --backend fp64 --p proth:56 --kinds kummer --m 2 --controls tower --t-min 8 --max-t 8"
  "C7|fp224|--engine generic --backend fp224 --p P-224 --kinds kummer --m 2 --controls tower --t-min 7 --max-t 7"
  "C8|p192|--engine generic --backend p192 --p P-192 --kinds dickson-torus --m 2 --controls tower --t-min 7 --max-t 7"
)
n=${#cells[@]}
for rep in $(seq 1 "$REPS"); do
  # odd repetitions run the list forward, even ones backward, so each
  # cell's arms stay adjacent and alternate which goes first
  for k in $(seq 0 $((n - 1))); do
    if (( rep % 2 == 1 )); then i=$k; else i=$((n - 1 - k)); fi
    IFS='|' read -r cell arm args <<<"${cells[$i]}"
    f="$OUT/$cell.$arm.rep$rep.jsonl"
    # shellcheck disable=SC2086
    RAYON_NUM_THREADS=1 taskset -c 3 nice -n 0 prlimit --as=4294967296 -- \
      "$BIN" $args "${common[@]}" --out "$f" 2>>"$OUT/stderr.log"
    echo "$(date -u +%FT%TZ) rep=$rep cell=$cell arm=$arm load=$(cut -d' ' -f1 /proc/loadavg)" >>"$OUT/progress.txt"
  done
done
