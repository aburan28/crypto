#!/bin/bash
# Headline protocol (benchmarks/sigma-fused/HEADLINE-PROTOCOL.md) for the
# Frobenius nibble table: control is the confirmed fused preset, candidate adds
# PACKED_SIGMA_TABLE=1 at THREADS=512 MINBLOCKS=1.  Same gates, warmups, A/A
# and A/B panels, 64 launches per sample, same native summarizer.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1
apt-get install -y -qq make g++ >/dev/null 2>&1

R=${RESULTS:-/results}
mkdir -p "$R"
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCH="-gencode arch=compute_${CAP},code=sm_${CAP}"
VERIFY_THREADS=96256
BENCH_THREADS=385024
fail=0

if [ "$GPU_NAME" != "NVIDIA RTX PRO 6000 Blackwell Server Edition" ] || [ "$CAP" != 120 ]; then
  echo "frozen hardware mismatch: $GPU_NAME sm_$CAP" >&2
  exit 2
fi
nvcc --version | grep -q 'release 13\.3, V13\.3\.73' || {
  echo "frozen compiler mismatch" >&2
  nvcc --version >&2
  exit 2
}

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
  echo "verify threads: $VERIFY_THREADS; benchmark threads: $BENCH_THREADS"
} | tee "$R/host.txt"

sha256sum Makefile include/packedkernels.cuh include/packedengine.cuh \
  include/packed131.h include/packedsigma131.h src/main.cu \
  src/testpackedcuda.cu benchmarks/sigma-table/gpujob-headline.sh \
  benchmarks/sigma-fused/corpus_identity.cpp benchmarks/sigma-fused/summarize.cpp \
  > "$R/source-files.sha256"

g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-fused/corpus_identity.cpp -o /tmp/corpus-identity || exit 1
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-fused/summarize.cpp -o /tmp/sigma-fused-summarize || exit 1

COMMON=(
  BATCH=16 SIGMA_FUSED=1
  PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
  PACKED_INLINE_POLY=3 WITNESS=0 WALK_TABLE=0
)

# control: the confirmed fused preset; candidate: the same with the Frobenius
# nibble table, one 512-thread block per SM (SIGMA-TABLE.md).
# HEADLINE_PERSIST=1 measures the recommended configuration instead: the
# candidate adds PACKED_L2_PERSIST=1 and runs at one wave (96,256 workers),
# its best population, against the preset at four waves, its best.
HEADLINE_PERSIST=${HEADLINE_PERSIST:-0}
knobs() {
  case "$1" in
    control)   echo "THREADS=256 MINBLOCKS=2 PACKED_SIGMA_TABLE=0" ;;
    candidate) echo "THREADS=512 MINBLOCKS=1 PACKED_SIGMA_TABLE=1 PACKED_L2_PERSIST=$HEADLINE_PERSIST" ;;
  esac
}
threadsOf() { if [ "$HEADLINE_PERSIST" = 1 ] && [ "$1" = candidate ]; then echo 96256; else echo "$BENCH_THREADS"; fi; }
tableOf() { case "$1" in control) echo 0 ;; candidate) echo 1 ;; esac; }
build() {
  local name=$1
  echo "=== build $name $(knobs "$name")"
  # shellcheck disable=SC2046
  make -s -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" $(knobs "$name") \
    > "$R/build-$name.log" 2>&1 || return 1
  grep -E "Function properties for.*walk|registers|spill|stack frame" \
    "$R/build-$name.log" | tail -12 | tee "$R/build-$name.txt" || true
  [ -x ecc2k130 ] || return 1
  mv ecc2k130 "ecc2k130-$name"
  echo "=== arithmetic $name"
  # shellcheck disable=SC2046
  make -s test-packed-cuda ARCH="$ARCH" "${COMMON[@]}" $(knobs "$name") > "$R/arithmetic-$name.log" 2>&1 || return 1
  grep -E "PASS|mismatch|sigma table" "$R/arithmetic-$name.log" | tee "$R/arithmetic-$name.txt"
  grep -qx "packed arithmetic sigma table: $(tableOf "$name")" "$R/arithmetic-$name.log" || return 1
  echo "=== integration $name"
  python3 codegen/testpackedclient.py "./ecc2k130-$name" > "$R/integration-$name.log" 2>&1 || return 1
  grep -E "PASS|Error|assert" "$R/integration-$name.log" | tee "$R/integration-$name.txt"
}

build control || fail=1
build candidate || fail=1
if [ "$fail" != 0 ]; then
  tail -60 "$R"/build-*.log
  exit "$fail"
fi
sha256sum ecc2k130-control ecc2k130-candidate > "$R/binary-sha256.txt"

verify() {
  local name=$1 status=0
  ./ecc2k130-$name --curve 131 --packed --threads "$VERIFY_THREADS" \
    --dp-weight 48 --dp-cap 262144 --steps 95 --launches 7 \
    --verify 300 --run-id 23 --dp-file "$R/dp-$name.bin" \
    > "$R/verify-$name.log" 2>&1 || status=$?
  grep -E "MISMATCH|OVERFLOW|finished|packed sigma|resident|registers" \
    "$R/verify-$name.log" | tee "$R/verify-$name.txt" || true
  grep -qx "packed sigma fused: 1" "$R/verify-$name.log" || status=1
  grep -q "^packed sigma table: $(tableOf "$name") " "$R/verify-$name.log" || status=1
  grep -qx "packed witness: 0" "$R/verify-$name.log" || status=1
  grep -Eq "^backend cuda-packed131: $VERIFY_THREADS threads x 16 slots x 1 lanes = $((VERIFY_THREADS*16)) walks, dp weight 48, 95 steps per launch$" \
    "$R/verify-$name.log" || status=1
  grep -Eq "\(300 verified against the reference, 0 dropped\)" "$R/verify-$name.log" || status=1
  ! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-$name.log" || status=1
  return "$status"
}

verify control || fail=1
verify candidate || fail=1
/tmp/corpus-identity "$R/dp-control.bin" "$R/dp-candidate.bin" \
  | tee "$R/corpus-identity.txt" || fail=1

if [ "$fail" != 0 ]; then
  echo "correctness preflight failed; timing suppressed" | tee "$R/preflight.txt"
  exit "$fail"
fi
echo "PASS: replay and sorted corpus identity" | tee "$R/preflight.txt"

printf 'phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 pair=$2 order=$3 variant=$4 binary log rc count rate digest state
  binary=$variant
  case "$variant" in
    control|control_a|control_b) binary=control ;;
    candidate) binary=candidate ;;
    *) return 2 ;;
  esac
  log="$R/${phase}-${pair}-${order}-${variant}.log"
  local threads; threads=$(threadsOf "$binary")
  ./ecc2k130-$binary --curve 131 --packed --threads "$threads" \
    --bench --steps 1024 --launches 64 --verify 0 > "$log" 2>&1
  rc=$?
  count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  grep -qx "packed sigma fused: 1" "$log" || rc=1
  grep -q "^packed sigma table: $(tableOf "$binary") " "$log" || rc=1
  grep -qx "packed witness: 0" "$log" || rc=1
  grep -Eq "^backend cuda-packed131: $threads threads x 16 slots x 1 lanes = $((threads*16)) walks, dp weight 0, 1024 steps per launch$" \
    "$log" || rc=1
  if [ "$count" != 1 ] || ! awk -v r="$rate" 'BEGIN{exit !(r+0>0)}'; then rc=1; fi
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu \
    --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$phase" "$pair" "$order" "$variant" "$rate" "$digest" "$state" \
    | tee -a "$R/samples.tsv"
  return "$rc"
}

# Four excluded warmups, two per binary.
sample warmup 0 1 control  || fail=1
sample warmup 0 2 candidate || fail=1
sample warmup 0 3 candidate || fail=1
sample warmup 0 4 control  || fail=1

# Five A/A pairs, alternating order.
for pair in 1 2 3 4 5; do
  if [ $((pair % 2)) = 1 ]; then
    sample aa "$pair" 1 control_a || fail=1
    sample aa "$pair" 2 control_b || fail=1
  else
    sample aa "$pair" 1 control_b || fail=1
    sample aa "$pair" 2 control_a || fail=1
  fi
done

# Five A/B pairs, alternating order.
for pair in 1 2 3 4 5; do
  if [ $((pair % 2)) = 1 ]; then
    sample ab "$pair" 1 control || fail=1
    sample ab "$pair" 2 candidate || fail=1
  else
    sample ab "$pair" 1 candidate || fail=1
    sample ab "$pair" 2 control || fail=1
  fi
done

if [ "$fail" != 0 ]; then
  echo "timing row failed" | tee -a "$R/failures.txt"
  exit "$fail"
fi
/tmp/sigma-fused-summarize "$R/samples.tsv" "$R/preflight.txt" "$R/result.json" \
  | tee "$R/summary.txt" || exit 1
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo "=== done"
