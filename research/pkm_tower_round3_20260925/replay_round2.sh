#!/bin/bash
# The identity check of note §12.2: the engine with the compact basis re-runs
# every system of round 2 and of its cross-check that finished there, with
# the same flags, and `compare_builds.py` requires every deterministic field
# of every row, and every step of every trace, to come out as round 2 wrote
# it. Sizes round 2 did not finish (m = 4 at N = 16, m = 3 at N = 18, the D2
# diagnostic) are left out: they are round 3's cells. From the repository
# root:
#
#     cargo build --release --example pkm_tower_pilot
#     research/pkm_tower_round3_20260925/replay_round2.sh
#
# Each run writes replay/<name>.jsonl and replay/<name>.log, named after the
# round-2 file it repeats, and progress.txt records start and end times.
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:-$here/../../target/release/examples/pkm_tower_pilot}
cd "$here/replay"
ulimit -v 14000000
P0=786433
P1=2013265921
run() {
  name=$1
  shift
  echo "$(date -u +%FT%TZ) start $name: $*" >> ../progress.txt
  "$bin" --engine tower "$@" --out "$name.jsonl" 2> "$name.log"
  rc=$?
  echo "$(date -u +%FT%TZ) end $name exit $rc" >> ../progress.txt
}
# The cross-check cells (round-2 README, "Cross-check").
run XV1-m2-tower-null --kinds kummer,dickson,isogeny --m 2 --controls tower,null --max-t 9 --planted 2 --random 2 --planted-max-t 7 --budget 900 --ladder-t none
run XV2-m3-tower --kinds kummer,dickson,isogeny --m 3 --controls tower --max-t-m3 4 --planted 2 --random 2 --planted-max-t 3 --budget 900 --ladder-t none
run XV3-m4-tower-random --kinds kummer,isogeny --m 4 --controls tower --max-t-m4 3 --planted 0 --random 2 --budget 900 --ladder-t none
run XV3b-m4-tower-planted-stop64 --kinds kummer,isogeny --m 4 --controls tower --max-t-m4 3 --planted 2 --random 0 --stop-below 64 --budget 900 --ladder-t none
run XV4-ladder --kinds kummer,isogeny --m none --planted 2 --random 0 --budget 300 --ladder-t 4,6 --ladder-g 1,2,4,8,12
run XV5-kummer-m2-p32 --p 3221225473 --kinds kummer --m 2 --controls tower --max-t 8 --planted 2 --random 2 --planted-max-t 7 --budget 900 --ladder-t none
# Round 2's cells (its run_round2.sh), to the sizes they finished, with traces.
common=(--budget 7200 --max-nnz 2000000000 --ladder-t none --planted 0 --random 2 --trace)
run K1-kummer-m2-p1 "${common[@]}" --p $P1 --kinds kummer --m 2 --controls tower --t-min 6 --max-t 11
run K1n-kummer-m2-p1-null "${common[@]}" --p $P1 --kinds kummer --m 2 --controls null --t-min 6 --max-t 10
run K0-kummer-m2-p0-N20 "${common[@]}" --p $P0 --kinds kummer --m 2 --controls tower --t-min 10 --max-t 10 --stop-below 8
run I0-isogeny-m2-p0-N20 "${common[@]}" --p $P0 --kinds isogeny --m 2 --controls tower --t-min 10 --max-t 10 --stop-below 8
run M4-kummer-m4-p1 "${common[@]}" --p $P1 --kinds kummer --m 4 --controls tower --t-min 2 --max-t-m4 3
run M3-kummer-m3-p1 "${common[@]}" --p $P1 --kinds kummer --m 3 --controls tower --t-min 3 --max-t-m3 5
# Its confirmations (confirm_round2.sh, confirm_cells.txt).
run C-I0-isogeny-m2-p0-N20-t10 "${common[@]}" --p $P0 --kinds isogeny --m 2 --controls tower --t-min 10 --max-t 10 --cap 5 --stop-below 8
run C-K0-kummer-m2-p0-N20-t10 "${common[@]}" --p $P0 --kinds kummer --m 2 --controls tower --t-min 10 --max-t 10 --cap 5 --stop-below 8
run C-K1-kummer-m2-p1-t10 "${common[@]}" --p $P1 --kinds kummer --m 2 --controls tower --t-min 10 --max-t 10 --cap 5
run C-K1-kummer-m2-p1-t11 "${common[@]}" --p $P1 --kinds kummer --m 2 --controls tower --t-min 11 --max-t 11 --cap 5
run C-K1n-kummer-m2-p1-null-t10 "${common[@]}" --p $P1 --kinds kummer --m 2 --controls null --t-min 10 --max-t 10 --cap 5
# Its D1 diagnostic, whose trace round 2 kept.
run D1-kummer-m4-p1-N12-trace --p $P1 --kinds kummer --controls tower --planted 0 --random 1 --ladder-t none --trace --m 4 --t-min 3 --max-t-m4 3 --budget 600
echo "$(date -u +%FT%TZ) replay done" >> ../progress.txt
