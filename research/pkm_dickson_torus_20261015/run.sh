#!/bin/bash
# The runs of PROTOCOL.md §4, in order, one at a time, single-threaded, under
# a 2.5 GB address-space cap. Usage (from anywhere):
#
#     BIN=/path/to/pkm_tower_pilot research/pkm_dickson_torus_20261015/run.sh
#
# Each run writes runs/<name>.jsonl (rows) and runs/<name>.log (stderr);
# progress.txt records start, end and exit status. No run starts after the
# deadline of PROTOCOL.md §4 (150 minutes after R0 started).
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:?set BIN to the pkm_tower_pilot binary}
cd "$here/runs"
ulimit -v 2621440
export RAYON_NUM_THREADS=1
P0=786433
P1=2013265921
P2=3221225473
T24=14680063
T31=2046820351
common=(--engine tower --budget 7200 --max-nnz 2000000000 --ladder-t none --planted 0 --random 2)
start=$(date +%s)
deadline=$((start + 150 * 60))
run() {
  name=$1
  shift
  if [ "$(date +%s)" -gt "$deadline" ]; then
    echo "$(date -u +%FT%TZ) not run (deadline) $name: $*" >> ../progress.txt
    return
  fi
  echo "$(date -u +%FT%TZ) start $name: $*" >> ../progress.txt
  nice -n 10 "$bin" "${common[@]}" "$@" --out "$name.jsonl" 2> "$name.log"
  rc=$?
  echo "$(date -u +%FT%TZ) end $name exit $rc" >> ../progress.txt
}
TOR=(--kinds dickson-torus)
CH=(--m 3,4 --t-min 2 --max-t-m3 4 --max-t-m4 3)
run R0-isogeny-m2-p0-poscontrol --p $P0 --kinds isogeny --m 2 --controls tower --t-min 1 --max-t 5
run R1-torus-m2-T31 --p $T31 "${TOR[@]}" --m 2 --controls tower --t-min 1 --max-t 10
run R2-torus-m2-T24 --p $T24 "${TOR[@]}" --m 2 --controls tower --t-min 1 --max-t 10
run R3-split-dickson-m2-p1 --p $P1 --kinds dickson --m 2 --controls tower --t-min 1 --max-t 10
run R4-torus-m2-T31-null --p $T31 "${TOR[@]}" --m 2 --controls null --t-min 1 --max-t 10
run R5-torus-m2-T24-null --p $T24 "${TOR[@]}" --m 2 --controls null --t-min 1 --max-t 10
run R6-torus-m34-T31 --p $T31 "${TOR[@]}" "${CH[@]}" --controls tower
run R7-torus-m34-T24 --p $T24 "${TOR[@]}" "${CH[@]}" --controls tower
run R8-split-dickson-m34-p1 --p $P1 --kinds dickson "${CH[@]}" --controls tower
run R9-torus-m34-T31-null --p $T31 "${TOR[@]}" "${CH[@]}" --controls null --null-chain
run R10-torus-m34-T24-null --p $T24 "${TOR[@]}" "${CH[@]}" --controls null --null-chain
run R11-kummer-p2-deep --p $P2 --kinds kummer --m 2,3 --controls tower --t-min 1 --max-t 9 --max-t-m3 4
run R11n-kummer-m2-p2-deep-null --p $P2 --kinds kummer --m 2 --controls null --t-min 1 --max-t 9
M15=(--m 3 --t-min 5 --max-t-m3 5 --controls tower)
run R12-torus-m3-N15-T31-t0 --p $T31 "${TOR[@]}" "${M15[@]}" --targets 0
run R13-torus-m3-N15-T24-t0 --p $T24 "${TOR[@]}" "${M15[@]}" --targets 0
run R14-split-dickson-m3-N15-p1-t0 --p $P1 --kinds dickson "${M15[@]}" --targets 0
run R15a-torus-m3-N15-T31-t1 --p $T31 "${TOR[@]}" "${M15[@]}" --targets 1
run R15b-torus-m3-N15-T24-t1 --p $T24 "${TOR[@]}" "${M15[@]}" --targets 1
echo "$(date -u +%FT%TZ) all done" >> ../progress.txt
