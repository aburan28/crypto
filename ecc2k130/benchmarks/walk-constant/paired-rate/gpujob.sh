#!/bin/bash
# The paired rate WALK-CONSTANT.md section 11 is missing: the sigma kernel
# against the table kernel under the extended cycle rule (rule v2), built from
# one tree, verified, and benched alternating on one card in one session.
# Runs inside nvidia/cuda:13.3.1-devel-ubuntu24.04 with the ecc2k130 tree at
# /work and an output directory at $RESULTS (default /results).  Launchers,
# all running this same script on one RTX PRO 6000:
#   modal run --detach modal_job.py --job benchmarks/walk-constant/paired-rate/gpujob.sh --out DIR
#   python3 runpod_job.py --job benchmarks/walk-constant/paired-rate/gpujob.sh --out DIR
#   python3 aws/bench_job.py --job benchmarks/walk-constant/paired-rate/gpujob.sh --out DIR
# then: python3 benchmarks/walk-constant/paired-rate/summarize.py DIR
#
# The knob set is the one behind section 6's three paired sessions
# (benchmarks/table-walk/gpujob.sh, the audited 256 x 2 preset), so the new
# ratio is comparable with the old kernel's 0.870-0.916 row for row; the two
# geometries are its automatic worker count and the campaign's 385,024.
set -uo pipefail
cd /work || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq make g++ python3 >/dev/null 2>&1
R=${RESULTS:-/results}; mkdir -p "$R"
REPS=${REPS:-6}
{
  nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,power.limit,memory.total --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
  grep -n "cycle-rule=" aws/protocol.py
} | tee "$R/host.txt"

F='ARCH=-gencode arch=compute_120,code=sm_120'
K="BATCH=16 THREADS=256 MINBLOCKS=2 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256 PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1"

fail=0
build() {
  local name=$1; shift
  echo "=== build $name: $*"
  # shellcheck disable=SC2086
  make -B -s ecc2k130 "$F" $K "$@" > "$R/build-$name.log" 2>&1
  grep -A2 "eccPacked131.*walk" "$R/build-$name.log" | grep -E "registers|spill" | head -2 | tee "$R/build-$name.txt"
  if [ ! -x ecc2k130 ]; then echo "BUILD FAILED $name" | tee -a "$R/failures.txt"; tail -20 "$R/build-$name.log"; return 1; fi
  mv ecc2k130 "ecc2k130-$name"
}
build sigma WALK_TABLE=0 || fail=1
build table WALK_TABLE=1 || fail=1

# The device's cycle rule against the host reference's: the unit probe, then
# 300 device reports of each binary re-walked by the host reference.
echo "=== unit test: device table-walk primitives, rule v2 histories"
make test-table-walk-cuda "$F" $K > "$R/test-table-walk-cuda.log" 2>&1 || { echo "test-table-walk-cuda FAILED" | tee -a "$R/failures.txt"; fail=1; }
tail -8 "$R/test-table-walk-cuda.log"
verify() {
  local b=$1 status=0
  [ -x "ecc2k130-$b" ] || return 1
  echo "=== verify $b"
  ./ecc2k130-$b --curve 131 --packed --dp-weight 48 --dp-cap 262144 --steps 96 --launches 6 \
      --verify 300 --run-id 7 > "$R/verify-$b.log" 2>&1 || status=$?
  grep -E "MISMATCH|finished|table walk|resident|OVERFLOW" "$R/verify-$b.log" | tee "$R/verify-$b.txt" || true
  grep -q "(300 verified against the reference, 0 dropped)" "$R/verify-$b.log" || status=1
  [ "$status" = 0 ] || echo "VERIFY FAILED $b" | tee -a "$R/failures.txt"
  return "$status"
}
verify sigma || fail=1
verify table || fail=1

# Rate: alternating sigma and table, REPS rounds per geometry.  The order
# flips every round so neither binary always runs on a warmer card.
sample() {
  local b=$1; shift
  [ -x "ecc2k130-$b" ] || return
  ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 32 --verify 0 "$@" 2>&1 \
      | grep -E "finished|resident" | tr '\n' ' '
  nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader
}
for geometry in auto 385024; do
  args=(); [ "$geometry" = auto ] || args=(--threads "$geometry")
  for rep in $(seq 1 "$REPS"); do
    order="sigma table"; [ $((rep % 2)) = 0 ] && order="table sigma"
    for b in $order; do
      echo "=== bench $geometry $b rep $rep" | tee -a "$R/bench.txt"
      sample "$b" "${args[@]}" | tee -a "$R/bench.txt"
    done
  done
done
echo "=== done"
exit "$fail"
