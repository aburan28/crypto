#!/bin/bash
# Round 2 of the PKM tower measurement (note §11.4): the cells in the order
# the note fixes, one system at a time, with the example built from this
# commit. Usage, from the repository root:
#
#     cargo build --release --example pkm_tower_pilot
#     research/pkm_tower_round2_20260925/run_round2.sh
#
# Each cell writes runs/<name>.jsonl (one row per system) and runs/<name>.log
# (the example's stderr), and progress.txt records start and end times.
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:-$here/../../target/release/examples/pkm_tower_pilot}
cd "$here/runs"
ulimit -v 14000000
P0=786433
P1=2013265921
common=(--engine tower --budget 7200 --max-nnz 2000000000 --ladder-t none --planted 0 --random 2)
run() {
  name=$1
  shift
  echo "$(date -u +%FT%TZ) start $name: $*" >> ../progress.txt
  "$bin" "${common[@]}" "$@" --out "$name.jsonl" 2> "$name.log"
  # Read the status before anything else runs: `$(date)` would reset it.
  # (Fixed after round 2 ran; its progress.txt logs every exit as 0.)
  rc=$?
  echo "$(date -u +%FT%TZ) end $name exit $rc" >> ../progress.txt
}
run K1-kummer-m2-p1 --p $P1 --kinds kummer --m 2 --controls tower --t-min 6 --max-t 11
run K1n-kummer-m2-p1-null --p $P1 --kinds kummer --m 2 --controls null --t-min 6 --max-t 10
run K0-kummer-m2-p0-N20 --p $P0 --kinds kummer --m 2 --controls tower --t-min 10 --max-t 10 --stop-below 8
run I0-isogeny-m2-p0-N20 --p $P0 --kinds isogeny --m 2 --controls tower --t-min 10 --max-t 10 --stop-below 8
run M4-kummer-m4-p1 --p $P1 --kinds kummer --m 4 --controls tower --t-min 2 --max-t-m4 4
run M3-kummer-m3-p1 --p $P1 --kinds kummer --m 3 --controls tower --t-min 3 --max-t-m3 6
echo "$(date -u +%FT%TZ) all done" >> ../progress.txt
