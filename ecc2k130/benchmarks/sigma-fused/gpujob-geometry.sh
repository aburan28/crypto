#!/bin/bash
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1
apt-get install -y -qq make g++ >/dev/null 2>&1
R=${RESULTS:-/results}; mkdir -p "$R"
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
if [ "$GPU_NAME" != "NVIDIA RTX PRO 6000 Blackwell Server Edition" ] || [ "$CAP" != 120 ]; then
  echo "frozen hardware mismatch: $GPU_NAME sm_$CAP" >&2; exit 2
fi
nvcc --version | grep -q 'release 13\.3, V13\.3\.73' || { nvcc --version >&2; exit 2; }
ARCH="-gencode arch=compute_120,code=sm_120"
fail=0
ARMS="b16t256 b32t256 b32t512 b64t256 b64t512"

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
} | tee "$R/host.txt"
sha256sum Makefile include/packedkernels.cuh include/packedengine.cuh include/packed131.h \
  benchmarks/sigma-fused/GEOMETRY-PROTOCOL.md benchmarks/sigma-fused/gpujob-geometry.sh \
  benchmarks/sigma-fused/corpus_identity.cpp benchmarks/sigma-fused/geometry_summarize.cpp \
  > "$R/source-files.sha256"
g++ -O2 -std=c++17 -Wall -Wextra -Werror benchmarks/sigma-fused/corpus_identity.cpp -o /tmp/corpus-identity || exit 1
g++ -O2 -std=c++17 -Wall -Wextra -Werror benchmarks/sigma-fused/geometry_summarize.cpp -o /tmp/geometry-summarize || exit 1

COMMON=(
  PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3
  PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1
  PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
  PACKED_INLINE_POLY=3 SIGMA_FUSED=1 WITNESS=0 WALK_TABLE=0
)
build() {
  local arm=$1 batch=$2 block=$3 minblocks=$4
  make -s -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" BATCH="$batch" THREADS="$block" MINBLOCKS="$minblocks" \
    > "$R/build-$arm.log" 2>&1 || return 1
  [ -x ecc2k130 ] || return 1
  mv ecc2k130 "ecc2k130-$arm"
}
build b16t256 16 256 2 || fail=1
build b32t256 32 256 2 || fail=1
build b32t512 32 512 1 || fail=1
build b64t256 64 256 2 || fail=1
build b64t512 64 512 1 || fail=1
if [ "$fail" != 0 ]; then tail -40 "$R"/build-*.log; exit "$fail"; fi
sha256sum ecc2k130-* > "$R/binary-sha256.txt"

verify() {
  local arm=$1 batch=$2 block=$3 blocks=$4 threads=$5 status=0
  ./ecc2k130-$arm --curve 131 --packed --threads "$threads" --dp-weight 48 --dp-cap 262144 \
    --steps 95 --launches 7 --verify 300 --run-id 31 --dp-file "$R/dp-$arm.bin" \
    > "$R/verify-$arm.log" 2>&1 || status=$?
  grep -qx "packed sigma fused: 1" "$R/verify-$arm.log" || status=1
  grep -qx "packed witness: 0" "$R/verify-$arm.log" || status=1
  grep -qx "packed launch bounds: $block threads, $blocks min blocks" "$R/verify-$arm.log" || status=1
  grep -Eq "^backend cuda-packed131: $threads threads x $batch slots x 1 lanes = 1540096 walks, dp weight 48, 95 steps per launch$" "$R/verify-$arm.log" || status=1
  grep -Eq "\(300 verified against the reference, 0 dropped\)" "$R/verify-$arm.log" || status=1
  ! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-$arm.log" || status=1
  return "$status"
}
verify b16t256 16 256 2 96256 || fail=1
verify b32t256 32 256 2 48128 || fail=1
verify b32t512 32 512 1 48128 || fail=1
verify b64t256 64 256 2 24064 || fail=1
verify b64t512 64 512 1 24064 || fail=1
/tmp/corpus-identity --canonical-out "$R/corpus-sorted.bin" \
  "$R/dp-b16t256.bin" "$R/dp-b32t256.bin" "$R/dp-b32t512.bin" \
  "$R/dp-b64t256.bin" "$R/dp-b64t512.bin" | tee "$R/corpus-identity.txt" || fail=1
if [ "$fail" != 0 ]; then echo "correctness preflight failed" | tee "$R/preflight.txt"; exit "$fail"; fi
sha256sum "$R/corpus-sorted.bin" > "$R/corpus-sorted.sha256"
wc -c "$R/corpus-sorted.bin" > "$R/corpus-sorted.bytes"
rm -f "$R"/dp-*.bin
echo "PASS: replay and sorted v1 corpus identity" | tee "$R/preflight.txt"

printf 'phase\tround\torder\tvariant\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 round=$2 order=$3 arm=$4 batch=$5 block=$6 blocks=$7 threads=$8 launches=$9 log rc count rate digest state
  log="$R/${phase}-${round}-${order}-${arm}.log"
  ./ecc2k130-$arm --curve 131 --packed --threads "$threads" --bench --steps 1024 --launches "$launches" --verify 0 > "$log" 2>&1
  rc=$?; count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  grep -qx "packed sigma fused: 1" "$log" || rc=1
  grep -qx "packed witness: 0" "$log" || rc=1
  grep -qx "packed launch bounds: $block threads, $blocks min blocks" "$log" || rc=1
  grep -Eq "^backend cuda-packed131: $threads threads x $batch slots x 1 lanes = 6160384 walks, dp weight 0, 1024 steps per launch$" "$log" || rc=1
  grep -Eq "\(0 verified against the reference, 0 dropped\)" "$log" || rc=1
  ! grep -Eq "MISMATCH|OVERFLOW" "$log" || rc=1
  if [ "$count" != 1 ] || ! awk -v r="$rate" 'BEGIN{exit !(r+0>0)}'; then rc=1; fi
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$phase" "$round" "$order" "$arm" "$rate" "$digest" "$state" | tee -a "$R/samples.tsv"
  return "$rc"
}
sample warmup 0 1 b16t256 16 256 2 385024 16 || fail=1
sample warmup 0 2 b32t256 32 256 2 192512 16 || fail=1
sample warmup 0 3 b32t512 32 512 1 192512 16 || fail=1
sample warmup 0 4 b64t256 64 256 2 96256 16 || fail=1
sample warmup 0 5 b64t512 64 512 1 96256 16 || fail=1
for round in 1 3; do
  sample screen "$round" 1 b16t256 16 256 2 385024 32 || fail=1
  sample screen "$round" 2 b32t256 32 256 2 192512 32 || fail=1
  sample screen "$round" 3 b32t512 32 512 1 192512 32 || fail=1
  sample screen "$round" 4 b64t256 64 256 2 96256 32 || fail=1
  sample screen "$round" 5 b64t512 64 512 1 96256 32 || fail=1
done
sample screen 2 1 b64t512 64 512 1 96256 32 || fail=1
sample screen 2 2 b64t256 64 256 2 96256 32 || fail=1
sample screen 2 3 b32t512 32 512 1 192512 32 || fail=1
sample screen 2 4 b32t256 32 256 2 192512 32 || fail=1
sample screen 2 5 b16t256 16 256 2 385024 32 || fail=1
if [ "$fail" != 0 ]; then exit "$fail"; fi
/tmp/geometry-summarize "$R/samples.tsv" "$R/preflight.txt" "$R/result.json" | tee "$R/summary.txt" || exit 1
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo "=== done"
