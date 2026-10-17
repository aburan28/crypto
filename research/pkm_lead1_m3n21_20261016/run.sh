#!/bin/bash
# PKM lead 1, m = 3 at N = 21 (PROTOCOL.md): Kummer tower, p1 = 15*2^27 + 1,
# f4_fp_tower, both of the cell's random targets. From the repository root:
#
#     cargo build --release --example pkm_tower_pilot --example pkm_tower_check
#     research/pkm_lead1_m3n21_20261016/run.sh host
#     research/pkm_lead1_m3n21_20261016/run.sh identity
#     research/pkm_lead1_m3n21_20261016/run.sh main 0      # then: main 1
#     research/pkm_lead1_m3n21_20261016/run.sh confirm 0 D # target 0 finished with D
#     research/pkm_lead1_m3n21_20261016/run.sh verify
#
# SHARE=k divides the memory caps by k when k processes run at once on the
# same host (PROTOCOL.md "Memory caps"). Each run writes runs/<name>.jsonl and
# runs/<name>.log (the example's stderr: the step trace). progress.txt records
# when each run started and ended and its exit status; a run that ended is not
# repeated.
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:-$here/../../target/release/examples/pkm_tower_pilot}
chk=${CHK:-$here/../../target/release/examples/pkm_tower_check}
share=${SHARE:-1}
mkdir -p "$here/runs"

# Memory caps from the host, fixed by PROTOCOL.md:
#   A = 0.92 * MemTotal / SHARE           (address-space cap, ulimit -v)
#   C = largest multiple of 1e8 with 4 * 1.45 * C <= A - 4e9   (--max-dense)
mem_kib=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)
a_bytes=$(( mem_kib * 1024 * 92 / 100 / share ))
c_entries=$(( (a_bytes - 4000000000) * 100 / 580 / 100000000 * 100000000 ))
a_kib=$(( a_bytes / 1024 ))

cell=(--engine tower --p 2013265921 --kinds kummer --m 3 --controls tower
      --t-min 7 --max-t-m3 7 --planted 0 --random 2 --ladder-t none --trace)
ident=(--engine tower --p 2013265921 --kinds kummer --m 3 --controls tower
       --t-min 5 --max-t-m3 5 --planted 0 --random 2 --ladder-t none --trace)
budget=345600

run() {
  name=$1
  shift
  if grep -q " end $name exit" "$here/progress.txt" 2>/dev/null; then
    return
  fi
  echo "$(date -u +%FT%TZ) start $name (A=${a_kib}KiB C=$c_entries share=$share): $*" >> "$here/progress.txt"
  ( ulimit -v "$a_kib"; "$bin" "$@" --out "$here/runs/$name.jsonl" 2> "$here/runs/$name.log" )
  # Read the status before anything else runs.
  rc=$?
  echo "$(date -u +%FT%TZ) end $name exit $rc" >> "$here/progress.txt"
}

case "${1:-}" in
  host)
    {
      echo "date: $(date -u +%FT%TZ)"
      echo "commit: $(git -C "$here" rev-parse HEAD)"
      echo "binary sha256: $(sha256sum "$bin" | cut -d' ' -f1)"
      echo "rustc: $(rustc --version)"
      echo "MemTotal KiB: $mem_kib; A KiB: $a_kib; C entries: $c_entries; SHARE: $share"
      echo "RAYON_NUM_THREADS: ${RAYON_NUM_THREADS:-unset (all cores)}"
      lscpu
    } > "$here/runs/host-$(hostname).txt"
    ;;
  identity)
    # One identity run per host: the measuring host repeats it.
    run "identity-N15-t0-$(hostname)" "${ident[@]}" --targets 0 --budget 7200
    ;;
  main)
    t=$2
    run "M3-N21-t$t" "${cell[@]}" --targets "$t" --budget "$budget" --max-dense "$c_entries"
    ;;
  confirm)
    t=$2
    cap=$(($3 - 1))
    run "C-M3-N21-t$t-cap$cap" "${cell[@]}" --targets "$t" --budget "$budget" --max-dense "$c_entries" --cap "$cap"
    ;;
  verify)
    "$chk" verify "$here"/runs/*.jsonl | tee "$here/verify.txt"
    ;;
  *)
    echo "usage: $0 host | identity | main T | confirm T D | verify" >&2
    exit 2
    ;;
esac
