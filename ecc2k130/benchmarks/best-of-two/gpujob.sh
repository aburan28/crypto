#!/bin/bash
# SIGMA_SIGN (BEST-OF-TWO.md): the campaign's sigma preset (aws/build.sh knobs,
# WITNESS=0) against the same build stepping to R + sigma^j(R) or R - sigma^j(R)
# by bit 4 of the normal-basis weight.  Two different walks, so no corpus
# identity: each binary is checked against its own reference instead.
#
#   1. correctness: packed arithmetic, client integration (replay, restart,
#      resume), 300 device reports re-walked by the reference, and 16 x 4
#      planted GF(2^41) logs on the host build of each walk;
#   2. timing: four warmups, five A/A pairs, five A/B pairs, 64 launches of
#      1024 steps at the campaign's 385,024 workers, native summarizer.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1
apt-get install -y -qq make g++ python3 >/dev/null 2>&1

R=${RESULTS:-/results}
mkdir -p "$R"
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCH="-gencode arch=compute_${CAP},code=sm_${CAP}"
VERIFY_THREADS=96256
BENCH_THREADS=385024
fail=0

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
  echo "host cpus: $(nproc)"
  echo "verify threads: $VERIFY_THREADS; benchmark threads: $BENCH_THREADS"
} | tee "$R/host.txt"

sha256sum Makefile include/packedkernels.cuh include/packedengine.cuh include/ref.h \
  include/solver.h include/walk.h src/main.cu benchmarks/best-of-two/gpujob.sh \
  benchmarks/sigma-fused/summarize.cpp > "$R/source-files.sha256"
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-fused/summarize.cpp -o /tmp/sigma-fused-summarize || exit 1

# The campaign knob set (aws/build.sh), CLMAD on for sm_120.
COMMON=(
  BATCH=16 THREADS=256 MINBLOCKS=2 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1
  PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
  WALK_TABLE=0 TABLE_BRANCHES=8 WITNESS=0
)
# CANDIDATE_SIGN: 1 = bit 4 of the weight, 2 = parity of popcount(x AND sigma^j x)
CANDIDATE_SIGN=${CANDIDATE_SIGN:-2}
signOf() { case "$1" in control) echo 0 ;; candidate) echo "$CANDIDATE_SIGN" ;; esac; }

build() {
  local name=$1 sign
  sign=$(signOf "$name")
  echo "=== build $name SIGMA_SIGN=$sign"
  make -s -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" SIGMA_SIGN="$sign" \
    > "$R/build-$name.log" 2>&1 || return 1
  grep -E "Function properties for.*walk|registers|spill|stack frame" \
    "$R/build-$name.log" | tail -12 | tee "$R/build-$name.txt" || true
  [ -x ecc2k130 ] || return 1
  mv ecc2k130 "ecc2k130-$name"
  make -s -B ecc2k130-cpu "${COMMON[@]}" SIGMA_SIGN="$sign" > "$R/build-cpu-$name.log" 2>&1 || return 1
  mv ecc2k130-cpu "ecc2k130-cpu-$name"
  echo "=== arithmetic $name"
  make -s test-packed-cuda ARCH="$ARCH" "${COMMON[@]}" SIGMA_SIGN="$sign" \
    > "$R/arithmetic-$name.log" 2>&1 || return 1
  grep -E "PASS|mismatch" "$R/arithmetic-$name.log" | tee "$R/arithmetic-$name.txt"
  echo "=== integration $name"
  python3 codegen/testpackedclient.py "./ecc2k130-$name" > "$R/integration-$name.log" 2>&1 || return 1
  grep -E "PASS|Error|assert" "$R/integration-$name.log" | tee "$R/integration-$name.txt"
  echo "=== host self-test $name (GF(2^23))"
  "./ecc2k130-cpu-$name" --test --curve 23 > "$R/selftest-$name.log" 2>&1 || return 1
  tail -n 2 "$R/selftest-$name.log" | tee "$R/selftest-$name.txt"
}

build control || fail=1
build candidate || fail=1
if [ "$fail" != 0 ]; then
  tail -60 "$R"/build-*.log
  exit "$fail"
fi
sha256sum ecc2k130-control ecc2k130-candidate ecc2k130-cpu-control ecc2k130-cpu-candidate \
  > "$R/binary-sha256.txt"

verify() {
  local name=$1 status=0
  ./ecc2k130-$name --curve 131 --packed --threads "$VERIFY_THREADS" \
    --dp-weight 48 --dp-cap 262144 --steps 95 --launches 7 \
    --verify 300 --run-id 23 --dp-file "$R/dp-$name.bin" \
    > "$R/verify-$name.log" 2>&1 || status=$?
  grep -E "MISMATCH|OVERFLOW|finished|packed sigma|resident|registers" \
    "$R/verify-$name.log" | tee "$R/verify-$name.txt" || true
  grep -qx "packed sigma sign: $(signOf "$name")" "$R/verify-$name.log" || status=1
  grep -qx "packed witness: 0" "$R/verify-$name.log" || status=1
  grep -Eq "\(300 verified against the reference, 0 dropped\)" "$R/verify-$name.log" || status=1
  ! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-$name.log" || status=1
  return "$status"
}
verify control || fail=1
verify candidate || fail=1

# Planted GF(2^41) logs on the host build: 16 instances x 4 run ids per walk.
printf 'variant\tinstance\trunId\tsolved\titerations\twalks\n' > "$R/planted.tsv"
planted() {
  local name=$1 inst rid log ok it walks status=0
  for inst in $(seq 0 15); do
    for rid in 1 2 3 4; do
      log="$R/planted-$name-$inst-$rid.log"
      timeout 900 "./ecc2k130-cpu-$name" --curve 41 --instance "$inst" --ref \
        --threads "$(nproc)" --dp-weight 8 --run-id "$rid" > "$log" 2>&1
      ok=0
      grep -q "matches the planted discrete log: yes" "$log" && ok=1
      it=$(sed -nE 's/.*solved after ([0-9]+) iterations of ([0-9]+) walks.*/\1/p' "$log")
      walks=$(sed -nE 's/.*solved after ([0-9]+) iterations of ([0-9]+) walks.*/\2/p' "$log")
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$inst" "$rid" "$ok" "${it:-}" "${walks:-}" >> "$R/planted.tsv"
      [ "$ok" = 1 ] || status=1
    done
  done
  echo "planted $name: $(awk -v n="$name" '$1==n && $4==1' "$R/planted.tsv" | wc -l)/64 recovered" \
    | tee -a "$R/planted.txt"
  return "$status"
}
planted control || fail=1
planted candidate || fail=1

if [ "$fail" != 0 ]; then
  echo "correctness preflight failed; timing suppressed" | tee "$R/preflight.txt"
  exit "$fail"
fi
echo "PASS: arithmetic, integration, host self-test, 300/300 replay and 64/64 planted logs per walk" \
  | tee "$R/preflight.txt"

printf 'phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 pair=$2 order=$3 variant=$4 binary log rc count rate digest state
  case "$variant" in
    control|control_a|control_b) binary=control ;;
    candidate) binary=candidate ;;
    *) return 2 ;;
  esac
  log="$R/${phase}-${pair}-${order}-${variant}.log"
  ./ecc2k130-$binary --curve 131 --packed --threads "$BENCH_THREADS" \
    --bench --steps 1024 --launches 64 --verify 0 > "$log" 2>&1
  rc=$?
  count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  grep -qx "packed sigma sign: $(signOf "$binary")" "$log" || rc=1
  grep -qx "packed witness: 0" "$log" || rc=1
  grep -Eq "^backend cuda-packed131: $BENCH_THREADS threads x 16 slots x 1 lanes = $((BENCH_THREADS*16)) walks, dp weight 0, 1024 steps per launch$" \
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

sample warmup 0 1 control  || fail=1
sample warmup 0 2 candidate || fail=1
sample warmup 0 3 candidate || fail=1
sample warmup 0 4 control  || fail=1
for pair in 1 2 3 4 5; do
  if [ $((pair % 2)) = 1 ]; then
    sample aa "$pair" 1 control_a || fail=1
    sample aa "$pair" 2 control_b || fail=1
  else
    sample aa "$pair" 1 control_b || fail=1
    sample aa "$pair" 2 control_a || fail=1
  fi
done
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
sha256sum "$R"/* > "$R/artifact-files.sha256" 2>/dev/null
echo "=== done"
