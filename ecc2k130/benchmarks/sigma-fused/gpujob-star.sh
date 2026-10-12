#!/bin/bash
# Frozen native-only B16 fused-sigma one-knob star.  This script intentionally
# contains no Python dependency.  STAR-PROTOCOL.md defines the admission gate.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1
apt-get install -y -qq make g++ >/dev/null 2>&1

R=${RESULTS:-/results}
mkdir -p "$R"
fail=0

GPU_ROWS=$(nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader)
if [ "$(printf '%s\n' "$GPU_ROWS" | wc -l | tr -d ' ')" != 1 ]; then
  echo "frozen hardware requires exactly one visible GPU" >&2
  exit 2
fi
GPU_NAME=$(printf '%s\n' "$GPU_ROWS" | cut -d, -f1)
CAP=$(printf '%s\n' "$GPU_ROWS" | cut -d, -f2 | tr -d ' .')
if [ "$GPU_NAME" != "NVIDIA RTX PRO 6000 Blackwell Server Edition" ] || [ "$CAP" != 120 ]; then
  echo "frozen hardware mismatch: $GPU_NAME sm_$CAP" >&2
  exit 2
fi
if ! nvcc --version | grep -q 'release 13\.3, V13\.3\.73'; then
  echo "frozen compiler mismatch" >&2
  nvcc --version >&2
  exit 2
fi
SOURCE=${SOURCE_REV:-}
if ! printf '%s\n' "$SOURCE" | grep -Eq '^[0-9a-f]{40}$'; then
  echo "SOURCE_REV must be one exact clean 40-hex commit" >&2
  exit 2
fi
if [ -f SOURCE_REV ] && [ "$(cat SOURCE_REV)" != "$SOURCE" ]; then
  echo "SOURCE_REV file/environment mismatch" >&2
  exit 2
fi
if [ -d .git ] || [ -f .git ]; then
  GIT_HEAD=$(git rev-parse HEAD 2>/dev/null || true)
  [ "$GIT_HEAD" = "$SOURCE" ] || { echo "Git HEAD/SOURCE_REV mismatch" >&2; exit 2; }
  git diff --quiet -- . || { echo "tracked source is dirty" >&2; exit 2; }
  git diff --cached --quiet -- . || { echo "index is dirty" >&2; exit 2; }
fi

ARCH="-gencode arch=compute_120,code=sm_120"
VERIFY_THREADS=96256
BENCH_THREADS=385024
UPDATES=201863462912

g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-fused/star_log_check.cpp -o /tmp/star-log-check || exit 1
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-fused/star_summarize.cpp -o /tmp/star-summarize || exit 1
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-fused/corpus_identity.cpp -o /tmp/corpus-identity || exit 1
ARMS=()
while IFS=$'\t' read -r arm _; do ARMS+=("$arm"); done < <(/tmp/star-log-check list)
if [ "${#ARMS[@]}" != 11 ] || [ "${ARMS[0]}" != baseline ]; then
  echo "native arm contract did not return the frozen eleven-arm star" >&2
  exit 2
fi
CANDIDATES=("${ARMS[@]:1}")

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap \
    --format=csv,noheader
  nvcc --version | tail -2
  echo "source: $SOURCE"
  echo "geometry: B16/T256/minBlocks2; verify threads $VERIFY_THREADS; benchmark threads $BENCH_THREADS"
  echo "work: $UPDATES complete scalar updates per timing row"
  echo "protocol: benchmarks/sigma-fused/STAR-PROTOCOL.md"
} | tee "$R/host.txt"

{
  find Makefile include generated src -type f -print
  printf '%s\n' \
    benchmarks/sigma-fused/STAR-PROTOCOL.md \
    benchmarks/sigma-fused/gpujob-star.sh \
    benchmarks/sigma-fused/star_arms.hpp \
    benchmarks/sigma-fused/star_log_check.cpp \
    benchmarks/sigma-fused/star_summarize.cpp \
    benchmarks/sigma-fused/test_star_native.sh \
    benchmarks/sigma-fused/corpus_identity.cpp
} | LC_ALL=C sort -u | xargs sha256sum > "$R/source-files.sha256" || exit 1

COMMON=(
  BATCH=16 THREADS=256 MINBLOCKS=2
  PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_INLINE_POLY=3 PACKED_CLMAD=1
  PACKED_CLMAD_SQUARE=0 PACKED_KARAT3=0 PACKED_WEIGHTED_PREFIX=2
  PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 PACKED_TOP_CLMAD=0
  PACKED_STATE_TILE=256 PACKED_ADD_COMBINE=0 PACKED_TOP_HOIST=0
  PACKED_ONB_INV=0 PACKED_ALU_SQR=0 PACKED_SQUARE_TABLE=0
  PACKED_CHAINS=1 PACKED_SLOT_PREFETCH=0 PACKED_SLOT_PIPELINE=0
  SIGMA_FUSED=1 WITNESS=0 WALK_TABLE=0 PHASE_PROFILE=0
)

printf 'arm\tdeclaredOneKnob\tfullVariedFlags\n' > "$R/arm-flags.tsv"
build() {
  local arm=$1 flagLine declared
  local -a varied
  flagLine=$(/tmp/star-log-check flags "$arm") || return 1
  read -r -a varied <<< "$flagLine"
  declared=$(/tmp/star-log-check list | awk -F '\t' -v arm="$arm" '$1 == arm {print $2}')
  printf '%s\t%s\t%s\n' "$arm" "$declared" "$flagLine" | tee -a "$R/arm-flags.tsv"
  {
    printf 'make -s -B ecc2k130 ARCH=%q' "$ARCH"
    printf ' %q' "${COMMON[@]}" "${varied[@]}"
    printf '\n'
  } > "$R/build-command-$arm.txt"
  echo "=== build $arm: ${declared:-baseline}"
  make -s -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" "${varied[@]}" \
    > "$R/build-$arm.log" 2>&1 || return 1
  [ -x ecc2k130 ] || return 1
  grep -E "Function properties for.*walk|Used [0-9]+ registers|spill|stack frame" \
    "$R/build-$arm.log" | tail -20 > "$R/build-$arm-resources.txt" || true
  mv ecc2k130 "ecc2k130-$arm"
}

for arm in "${ARMS[@]}"; do build "$arm" || fail=1; done
if [ "$fail" != 0 ]; then
  echo "one or more builds failed; correctness and timing suppressed" | tee "$R/failures.txt"
  tail -80 "$R"/build-*.log
  exit "$fail"
fi
: > "$R/binary-sha256.txt"
for arm in "${ARMS[@]}"; do sha256sum "ecc2k130-$arm" >> "$R/binary-sha256.txt"; done

printf 'arm\tregisters\tlocalBytesPerThread\tstaticSharedBytesPerBlock\tl2WindowBytes\tl2FieldBytes\tl2CapBytes\n' \
  > "$R/resources.tsv"
verify() {
  local arm=$1 log="$R/verify-$1.log" status=0
  ./"ecc2k130-$arm" --curve 131 --packed --threads "$VERIFY_THREADS" \
    --dp-weight 48 --dp-cap 262144 --steps 95 --launches 7 \
    --verify 300 --run-id 37 --dp-file "$R/dp-$arm.bin" \
    > "$log" 2>&1 || status=$?
  if [ "$status" = 0 ]; then
    /tmp/star-log-check verify "$arm" "$log" >> "$R/resources.tsv" || status=1
  fi
  grep -E "MISMATCH|OVERFLOW|finished|packed kernel:|packed launch bounds:|packed L2 persist window:" \
    "$log" > "$R/verify-$arm-summary.txt" || true
  return "$status"
}
for arm in "${ARMS[@]}"; do verify "$arm" || fail=1; done
if [ "$fail" = 0 ]; then
  corpusArgs=()
  for arm in "${ARMS[@]}"; do corpusArgs+=("$R/dp-$arm.bin"); done
  /tmp/corpus-identity --canonical-out "$R/corpus-sorted.bin" "${corpusArgs[@]}" \
    | tee "$R/corpus-identity.txt" || fail=1
fi
if [ "$fail" != 0 ]; then
  echo "correctness preflight failed; timing suppressed" | tee "$R/preflight.txt"
  exit "$fail"
fi
sha256sum "$R/corpus-sorted.bin" > "$R/corpus-sorted.sha256"
wc -c "$R/corpus-sorted.bin" > "$R/corpus-sorted.bytes"
rm -f "$R"/dp-*.bin
echo "PASS: 300/300 replay and sorted v1 corpus identity across 11 arms" | tee "$R/preflight.txt"

printf 'phase\tcomparison\tpair\torder\tvariant\tbinary\trateMps\tupdates\tlogSha256\tgpuState\n' \
  > "$R/samples.tsv"
sample() {
  local phase=$1 comparison=$2 pair=$3 order=$4 variant=$5 arm=$6
  local log="$R/${phase}-${comparison}-${pair}-${order}-${variant}-${arm}.log"
  local rc=0 rate digest state
  ./"ecc2k130-$arm" --curve 131 --packed --threads "$BENCH_THREADS" \
    --bench --steps 1024 --launches 32 --verify 0 > "$log" 2>&1 || rc=$?
  if [ "$rc" = 0 ]; then
    rate=$(/tmp/star-log-check bench "$arm" "$log") || rc=1
  fi
  if [ "$rc" != 0 ]; then
    echo "invalid timing row: $phase $comparison $pair $variant $arm" | tee -a "$R/failures.txt"
    return 1
  fi
  digest=$(sha256sum "$log" | cut -d' ' -f1)
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu \
    --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$phase" "$comparison" "$pair" "$order" "$variant" "$arm" "$rate" \
    "$UPDATES" "$digest" "$state" | tee -a "$R/samples.tsv"
}

# Every excluded warmup is the same work as a retained row.
for arm in "${ARMS[@]}"; do
  sample warmup "$arm" 0 1 warmup "$arm" || fail=1
done

# Five matched A/A pairs establish the session's baseline noise envelope.
for pair in 1 2 3 4 5; do
  if [ $((pair % 2)) = 1 ]; then
    sample aa baseline-aa "$pair" 1 a baseline || fail=1
    sample aa baseline-aa "$pair" 2 b baseline || fail=1
  else
    sample aa baseline-aa "$pair" 1 b baseline || fail=1
    sample aa baseline-aa "$pair" 2 a baseline || fail=1
  fi
done

# Three paired screens per candidate.  Candidate traversal is forward,
# reverse, then half-rotated; A/B position alternates by arm and round.
for round in 1 2 3; do
  orderArms=()
  if [ "$round" = 1 ]; then
    orderArms=("${CANDIDATES[@]}")
  elif [ "$round" = 2 ]; then
    for ((i=${#CANDIDATES[@]}-1; i>=0; --i)); do orderArms+=("${CANDIDATES[i]}"); done
  else
    orderArms=("${CANDIDATES[@]:5}" "${CANDIDATES[@]:0:5}")
  fi
  for arm in "${orderArms[@]}"; do
    ordinal=0
    for ((i=0; i<${#CANDIDATES[@]}; ++i)); do
      [ "${CANDIDATES[i]}" = "$arm" ] && ordinal=$i
    done
    if [ $(((ordinal + round) % 2)) = 0 ]; then
      sample screen "$arm" "$round" 1 baseline baseline || fail=1
      sample screen "$arm" "$round" 2 candidate "$arm" || fail=1
    else
      sample screen "$arm" "$round" 1 candidate "$arm" || fail=1
      sample screen "$arm" "$round" 2 baseline baseline || fail=1
    fi
  done
done

if [ "$fail" != 0 ]; then
  echo "timing panel invalid; summarization suppressed" | tee -a "$R/failures.txt"
  exit "$fail"
fi
/tmp/star-summarize "$R/samples.tsv" "$R/preflight.txt" "$R/result.json" \
  | tee "$R/summary.txt" || exit 1
find "$R" -maxdepth 1 -type f ! -name artifact-files.sha256 -print0 \
  | LC_ALL=C sort -z | xargs -0 sha256sum > "$R/artifact-files.sha256"
echo "=== done"
