#!/bin/bash
# GPU job behind ROOFLINE.md: one run settles the two open pipe rates of the
# roofline and times the two integer cuts it pointed at.  Runs inside
# nvidia/cuda:13.3.1-devel-ubuntu24.04 with the ecc2k130 tree at /work and an
# output directory at /results (or $RESULTS), on one RTX PRO 6000:
#   modal run modal_job.py --job benchmarks/roofline/gpujob.sh --out DIR
#   python3 runpod_job.py --job benchmarks/roofline/gpujob.sh --out DIR \
#       --extra benchmarks/roofline --extra roofline.py
#   python3 aws/bench_job.py --job benchmarks/roofline/gpujob.sh --out DIR \
#       --extra benchmarks/roofline --extra roofline.py
# (Modal ships the whole tree; the other two ship Makefile, src, include,
# codegen and generated only, so the probe and roofline.py go in --extra.)
#
#  1. pipes.cu: LOP3/SHF/IMAD/FFMA alone and mixed (does the walk's ALU work
#     share the FMA pipe's issue?) and CLMAD in the walk's patterns (what does
#     the carry-less unit sustain?), with the SASS of every timed loop audited
#     by sass_loops.py so no stream the compiler rewrote can report a rate.
#  2. Six binaries, each the 20 B/s build of ONE-BLOCK-GEOMETRY.md plus what
#     its name says: ref, sqtab (PACKED_SQUARE_TABLE=1), invpoly
#     (PACKED_INV_POLY=1: eight inlined links), invpoly2 (PACKED_INV_POLY=2:
#     one out-of-line copy, 1,664 fewer instructions and 10 fewer registers
#     than invpoly for the same work), both and both2 (the square table with
#     each), and fused2: both2 in the one-pass kernel (TABLE_FUSED=1), an
#     exploratory arm -- the least ALU work of any build (1,218 per update),
#     but the fused kernel alone measured +0.2% where the model said +1.0%.
#     roofline.py prices each one with this container's toolchain, so the
#     prediction and the rate share a compiler.
#  3. Verification: 300 device reports re-walked by the host reference per
#     binary, at one forced common walk count, and the distinguished-point
#     sets must be identical across binaries -- the arms change how the field
#     is computed, not what is computed.
#  4. Rate: REPS rounds, binaries alternating, SM clock and power sampled.
#  5. If Nsight Compute can read the counters here: per-pipe utilisation of the
#     ref and both2 kernels, the hardware's own answer to the roofline's.
set -uo pipefail
cd /work || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq make g++ python3 >/dev/null 2>&1
R=${RESULTS:-/results}; mkdir -p "$R"
REPS=${REPS:-5}
{
  nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
} | tee "$R/host.txt"
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
SM="sm_${CAP}"
fail=0

# 1. Pipe probe.
echo "=== pipes ($SM)"
( cd benchmarks/roofline &&
  nvcc -O3 -std=c++17 -arch="$SM" -Xptxas -v pipes.cu -o /tmp/pipes > "$R/pipes-build.log" 2>&1 &&
  nvcc -O3 -std=c++17 -arch="$SM" -cubin pipes.cu -o /tmp/pipes.cubin >> "$R/pipes-build.log" 2>&1 &&
  nvdisasm -c /tmp/pipes.cubin > "$R/pipes.sass" &&
  python3 sass_loops.py "$R/pipes.sass" > "$R/pipes-sass-loops.jsonl" ) || { echo "PIPES BUILD FAILED" | tee -a "$R/failures.txt"; fail=1; }
[ -x /tmp/pipes ] && { /tmp/pipes 16384 3 2>&1 | tee "$R/pipes.jsonl" || fail=1; }

# 2. Binaries.
build() {
  local name=$1 target=$2; shift 2
  echo "=== build $name: $target $*"
  make -s "$target" "$@" "PRO6000_ARCH=-gencode arch=compute_${CAP},code=${SM}" > "$R/build-$name.log" 2>&1
  grep -A2 "eccPacked131.*walk" "$R/build-$name.log" | grep -E "registers|spill" | head -2 | tee "$R/build-$name.txt"
  if [ ! -x ecc2k130 ]; then echo "BUILD FAILED $name" | tee -a "$R/failures.txt"; tail -20 "$R/build-$name.log"; return 1; fi
  mv ecc2k130 "ecc2k130-$name"
  python3 roofline.py --target "$target" --gpu rtx-pro-6000 --arch "$SM" "${@/#/--make-var=}" --json "$R/roofline-$name.json" \
      > "$R/roofline-$name.txt" 2>&1 || echo "roofline.py failed for $name" | tee -a "$R/failures.txt"
}
VARIANTS="ref sqtab invpoly invpoly2 both both2 fused2"
build ref      gpu-rtx-pro6000-20b || fail=1
build sqtab    gpu-rtx-pro6000-20b PACKED_SQUARE_TABLE=1 || fail=1
build invpoly  gpu-rtx-pro6000-20b PACKED_INV_POLY=1 || fail=1
build invpoly2 gpu-rtx-pro6000-20b PACKED_INV_POLY=2 || fail=1
build both     gpu-rtx-pro6000-20b PACKED_SQUARE_TABLE=1 PACKED_INV_POLY=1 || fail=1
build both2    gpu-rtx-pro6000-20b PACKED_SQUARE_TABLE=1 PACKED_INV_POLY=2 || fail=1
build fused2   gpu-preset "KNOBS=TABLE_FUSED=1 TABLE_PIPE_SELECT=0 PACKED_SQUARE_TABLE=1 PACKED_INV_POLY=2" || fail=1

# 3. Verification at the reference's automatic walk count.
SMS=$(./ecc2k130-ref --curve 131 --packed --bench --steps 1 --launches 1 --verify 0 2>&1 | grep -oE '[0-9]+ SMs' | grep -oE '[0-9]+' | head -1)
T512=$(( ${SMS:-188} * 512 ))
echo "sms: ${SMS:-unknown}, verify threads $T512" | tee -a "$R/host.txt"
verify() {
  local b=$1 status=0
  [ -x "ecc2k130-$b" ] || return 1
  echo "=== verify $b (threads $T512)"
  ./ecc2k130-$b --curve 131 --packed --threads "$T512" --dp-weight 48 --dp-cap 262144 \
      --steps 96 --launches 6 --verify 300 --run-id 7 --dp-file "$R/dp-$b.bin" \
      > "$R/verify-$b.log" 2>&1 || status=$?
  grep -E "MISMATCH|finished|packed chains|resident|registers|persist|OVERFLOW|square table|polynomial inversion|shared bytes" \
      "$R/verify-$b.log" | tee "$R/verify-$b.txt" || true
  return "$status"
}
for b in $VARIANTS; do verify "$b" || fail=1; done
python3 - "$R" $VARIANTS <<'PY' | tee "$R/dp-identity.txt" || fail=1
import hashlib, os, sys
r = sys.argv[1]
ref = None
bad = 0
for name in sys.argv[2:]:
    path = os.path.join(r, "dp-%s.bin" % name)
    if not os.path.exists(path):
        print("%-20s missing" % name)
        bad = 1
        continue
    data = open(path, "rb").read()
    recs = sorted(data[i:i+32] for i in range(0, len(data) - len(data) % 32, 32))
    digest = hashlib.sha256(b"".join(recs)).hexdigest()
    if ref is None: ref = digest
    if digest != ref: bad = 1
    print("%-20s %8d records  sha256 %s  %s" % (name, len(recs), digest[:16], "IDENTICAL to ref" if digest == ref else "DIFFERS from ref"))
sys.exit(bad)
PY

# 4. Rate: alternating binaries, REPS rounds, automatic worker count.
sample() {
  local b=$1
  [ -x "ecc2k130-$b" ] || return
  ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 64 --verify 0 2>&1 \
      | grep -E "finished|resident" | tr '\n' ' '
  nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader
}
for rep in $(seq 1 "$REPS"); do
  for b in $VARIANTS; do
    echo "=== bench $b rep $rep" | tee -a "$R/bench.txt"
    sample "$b" | tee -a "$R/bench.txt"
  done
done

# 5. Hardware counters, where the host lets a container read them.  Failure
# here (ERR_NVGPUCTRPERM on most rented hosts) is recorded, not fatal.
NCU=$(command -v ncu || ls /usr/local/cuda/bin/ncu /opt/nvidia/nsight-compute/*/ncu 2>/dev/null | head -1)
if [ -n "$NCU" ]; then
  # Every per-pipe counter this ncu knows for the part, so no metric name is
  # guessed: an unknown name would fail the whole profile.
  "$NCU" --query-metrics > "$R/ncu-query.txt" 2>&1
  METRICS=$(grep -oE '^(sm__inst_executed_pipe_[a-z0-9_]+|sm__pipe_[a-z0-9_]+_cycles_active|smsp__issue_active)\b' \
      "$R/ncu-query.txt" | sort -u | sed 's/$/.avg.pct_of_peak_sustained_active/' | paste -sd, -)
  for b in ref both2; do
    [ -x "ecc2k130-$b" ] || continue
    echo "=== ncu $b"
    "$NCU" --kernel-name regex:walk --launch-skip 2 --launch-count 1 --print-details all \
        --section SpeedOfLight --section ComputeWorkloadAnalysis --section InstructionStats --section Occupancy \
        ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 4 --verify 0 > "$R/ncu-$b.txt" 2>&1 \
        || echo "ncu $b sections: exit $? (see ncu-$b.txt)" | tee -a "$R/host.txt"
    [ -n "$METRICS" ] && { "$NCU" --kernel-name regex:walk --launch-skip 2 --launch-count 1 --metrics "$METRICS,sm__cycles_elapsed.avg.per_second" \
        ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 4 --verify 0 > "$R/ncu-pipes-$b.txt" 2>&1 \
        || echo "ncu $b pipes: exit $? (see ncu-pipes-$b.txt)" | tee -a "$R/host.txt"; }
  done
else
  echo "ncu: not in this image" | tee -a "$R/host.txt"
fi
echo "=== done"
exit "$fail"
