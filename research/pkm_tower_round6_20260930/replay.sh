#!/bin/bash
# The identity check of note §15.2: the engine with the full-rank exit re-runs
# every committed system of rounds 2 and 3 with the flags they were run with,
# and `compare_exit.py` requires every field of every row, and every step of
# every trace, to come out as they wrote it, but for the multiply-adds and the
# residual-row count, which may only fall. A second pass runs a subset with the
# exit off, which must reproduce the rows exactly, multiply-adds included.
# From the repository root:
#
#     cargo build --release --example pkm_tower_pilot
#     research/pkm_tower_round6_20260930/replay.sh
#
# Each run writes replay/<name>.jsonl and .log (exit on) or off/<name>.jsonl
# and .log (exit off), named after the file it repeats, and progress.txt
# records when each run started and ended. A run that ended is not repeated
# when the script is started again.
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:-$here/../../target/release/examples/pkm_tower_pilot}
mkdir -p "$here/replay" "$here/off"
ulimit -v 14000000
P0=786433
P1=2013265921
run() {
  dir=$1
  name=$2
  shift 2
  if grep -q " end $dir/$name exit" "$here/progress.txt" 2>/dev/null; then
    return
  fi
  echo "$(date -u +%FT%TZ) start $dir/$name: $*" >> "$here/progress.txt"
  "$bin" --engine tower "$@" --out "$here/$dir/$name.jsonl" 2> "$here/$dir/$name.log"
  # Read the status before anything else runs.
  rc=$?
  echo "$(date -u +%FT%TZ) end $dir/$name exit $rc" >> "$here/progress.txt"
}
xv() {
  dir=$1
  shift
  run "$dir" XV1-m2-tower-null "$@" --kinds kummer,dickson,isogeny --m 2 --controls tower,null --max-t 9 --planted 2 --random 2 --planted-max-t 7 --budget 900 --ladder-t none
  run "$dir" XV2-m3-tower "$@" --kinds kummer,dickson,isogeny --m 3 --controls tower --max-t-m3 4 --planted 2 --random 2 --planted-max-t 3 --budget 900 --ladder-t none
  run "$dir" XV3-m4-tower-random "$@" --kinds kummer,isogeny --m 4 --controls tower --max-t-m4 3 --planted 0 --random 2 --budget 900 --ladder-t none
  run "$dir" XV3b-m4-tower-planted-stop64 "$@" --kinds kummer,isogeny --m 4 --controls tower --max-t-m4 3 --planted 2 --random 0 --stop-below 64 --budget 900 --ladder-t none
  run "$dir" XV4-ladder "$@" --kinds kummer,isogeny --m none --planted 2 --random 0 --budget 300 --ladder-t 4,6 --ladder-g 1,2,4,8,12
  run "$dir" XV5-kummer-m2-p32 "$@" --p 3221225473 --kinds kummer --m 2 --controls tower --max-t 8 --planted 2 --random 2 --planted-max-t 7 --budget 900 --ladder-t none
}
common=(--budget 7200 --max-nnz 2000000000 --ladder-t none --planted 0 --random 2 --trace)
d1=(--p $P1 --kinds kummer --controls tower --planted 0 --random 1 --ladder-t none --trace --m 4 --t-min 3 --max-t-m4 3 --budget 600)
m4=("${common[@]}" --p $P1 --kinds kummer --m 4 --controls tower --t-min 2 --max-t-m4 3)

# The exit off first: the refactor must change nothing at all.
xv off --no-full-rank-exit
run off D1-kummer-m4-p1-N12-trace --no-full-rank-exit "${d1[@]}"
run off M4-kummer-m4-p1 --no-full-rank-exit "${m4[@]}"

# The exit on (the default): round 2's cross-check, cells, confirmations and
# D1, as round 3 replayed them (its replay_round2.sh), then round 3's cells.
xv replay
run replay K1-kummer-m2-p1 "${common[@]}" --p $P1 --kinds kummer --m 2 --controls tower --t-min 6 --max-t 11
run replay K1n-kummer-m2-p1-null "${common[@]}" --p $P1 --kinds kummer --m 2 --controls null --t-min 6 --max-t 10
run replay K0-kummer-m2-p0-N20 "${common[@]}" --p $P0 --kinds kummer --m 2 --controls tower --t-min 10 --max-t 10 --stop-below 8
run replay I0-isogeny-m2-p0-N20 "${common[@]}" --p $P0 --kinds isogeny --m 2 --controls tower --t-min 10 --max-t 10 --stop-below 8
run replay M4-kummer-m4-p1 "${m4[@]}"
run replay M3-kummer-m3-p1 "${common[@]}" --p $P1 --kinds kummer --m 3 --controls tower --t-min 3 --max-t-m3 5
run replay C-I0-isogeny-m2-p0-N20-t10 "${common[@]}" --p $P0 --kinds isogeny --m 2 --controls tower --t-min 10 --max-t 10 --cap 5 --stop-below 8
run replay C-K0-kummer-m2-p0-N20-t10 "${common[@]}" --p $P0 --kinds kummer --m 2 --controls tower --t-min 10 --max-t 10 --cap 5 --stop-below 8
run replay C-K1-kummer-m2-p1-t10 "${common[@]}" --p $P1 --kinds kummer --m 2 --controls tower --t-min 10 --max-t 10 --cap 5
run replay C-K1-kummer-m2-p1-t11 "${common[@]}" --p $P1 --kinds kummer --m 2 --controls tower --t-min 11 --max-t 11 --cap 5
run replay C-K1n-kummer-m2-p1-null-t10 "${common[@]}" --p $P1 --kinds kummer --m 2 --controls null --t-min 10 --max-t 10 --cap 5
run replay D1-kummer-m4-p1-N12-trace "${d1[@]}"
m4b=(--p $P1 --kinds kummer --m 4 --controls tower --t-min 4 --max-t-m4 4 --planted 0 --random 2 --ladder-t none --budget 14400 --max-nnz 2000000000 --trace)
run replay C-M4b-kummer-m4-p1-N16 "${m4b[@]}" --cap 6
run replay M4b-kummer-m4-p1-N16 "${m4b[@]}"
echo "$(date -u +%FT%TZ) replay done" >> "$here/progress.txt"
