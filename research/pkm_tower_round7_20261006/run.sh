#!/bin/bash
# Round 7 of the PKM tower oracle (note section 16): Kummer, m = 3, p1, N = 18,
# on f4_fp_tower, both of round 2's targets. From the repository root:
#
#     cargo build --release --example pkm_tower_pilot
#     research/pkm_tower_round7_20261006/run.sh stage1
#     research/pkm_tower_round7_20261006/run.sh stage2 C      # C by section 16.3's rule
#     research/pkm_tower_round7_20261006/run.sh confirm C D T # target T finished with D
#
# Each run writes runs/<name>.jsonl and runs/<name>.log (the example's stderr:
# the step trace, and a `tower stop in step` line if a stop cut a step off).
# progress.txt records when each run started and ended and its exit status,
# read before anything else runs. A run that ended is not repeated.
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:-$here/../../target/release/examples/pkm_tower_pilot}
mkdir -p "$here/runs"
# The 14 GB address-space cap of rounds 2-6, in KiB.
ulimit -v 14000000
cell=(--engine tower --p 2013265921 --kinds kummer --m 3 --controls tower
      --t-min 6 --max-t-m3 6 --planted 0 --random 2 --ladder-t none --trace)
run() {
  name=$1
  shift
  if grep -q " end $name exit" "$here/progress.txt" 2>/dev/null; then
    return
  fi
  echo "$(date -u +%FT%TZ) start $name: $*" >> "$here/progress.txt"
  "$bin" "$@" --out "$here/runs/$name.jsonl" 2> "$here/runs/$name.log"
  rc=$?
  echo "$(date -u +%FT%TZ) end $name exit $rc" >> "$here/progress.txt"
}
case "${1:-}" in
  stage1)
    # Round 2's flags at t = 6, target 0: the identity check and the sizing.
    run M3b-sizing-t0 "${cell[@]}" --targets 0 --budget 7200 --max-nnz 2000000000
    ;;
  stage2)
    c=$2
    run M3b-t0 "${cell[@]}" --targets 0 --budget 43200 --max-dense "$c"
    run M3b-t1 "${cell[@]}" --targets 1 --budget 43200 --max-dense "$c"
    ;;
  confirm)
    c=$2
    cap=$(($3 - 1))
    t=$4
    run "C-M3b-t$t-cap$cap" "${cell[@]}" --targets "$t" --budget 43200 --max-dense "$c" --cap "$cap"
    ;;
  *)
    echo "usage: $0 stage1 | stage2 C | confirm C D T" >&2
    exit 2
    ;;
esac
