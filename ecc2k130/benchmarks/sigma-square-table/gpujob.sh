#!/usr/bin/env bash
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export LC_ALL=C CUDA_DISABLE_PTX_JIT=1 CUDA_FORCE_PTX_JIT=0 CUDA_FORCE_JIT=0
R=${RESULTS:-/results}
mkdir -p "$R"

ARCH='-gencode arch=compute_120,code=sm_120'
VERIFY_THREADS=96256
BENCH_THREADS=385024
PARTIAL_THREADS=513
EXPECTED_UPDATES=201863462912
RUN_ID=43
fail=0

if [[ ! ${SOURCE_REV:-} =~ ^[0-9a-f]{40}$ ]]; then
  echo "SOURCE_REV must be one clean 40-hex commit" | tee "$R/failures.txt"
  exit 2
fi
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
if [ "$GPU_NAME" != 'NVIDIA RTX PRO 6000 Blackwell Server Edition' ] || [ "$CAP" != 120 ] ||
   ! nvcc --version | grep -q 'release 13\.3, V13\.3\.73'; then
  echo "frozen hardware/compiler mismatch: $GPU_NAME sm_$CAP" | tee "$R/failures.txt"
  nvcc --version | tee -a "$R/failures.txt"
  exit 2
fi
{
  nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version
  uname -a
  echo "source=$SOURCE_REV"
  echo "geometry=B16/T256/minBlocks2 WITNESS0; verifyThreads=$VERIFY_THREADS; benchThreads=$BENCH_THREADS; partialThreads=$PARTIAL_THREADS"
  echo "updatesPerTimingRow=$EXPECTED_UPDATES"
} > "$R/host.txt" 2>&1

{
  find Makefile include generated src -type f -print0
  printf '%s\0' \
    benchmarks/sigma-square-table/PROTOCOL.md \
    benchmarks/sigma-square-table/gpujob.sh \
    benchmarks/sigma-square-table/log_check.cpp \
    benchmarks/sigma-square-table/summarize.cpp \
    benchmarks/sigma-square-table/compile_audit.cpp \
    benchmarks/sigma-square-table/replay.cpp \
    benchmarks/sigma-square-table/audit.cpp \
    benchmarks/sigma-square-table/run-native.sh \
    benchmarks/sigma-fused/corpus_identity.cpp
} | sort -z | xargs -0 sha256sum > "$R/source-files.sha256" || exit 1
sha256sum -c "$R/source-files.sha256" > "$R/source-files.check" 2>&1 || exit 1

g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-square-table/log_check.cpp -o /tmp/sigma-square-log-check || exit 1
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-square-table/summarize.cpp -o /tmp/sigma-square-summarize || exit 1
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-square-table/compile_audit.cpp -o /tmp/sigma-square-compile-audit || exit 1
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-fused/corpus_identity.cpp -o /tmp/corpus-identity || exit 1
/tmp/sigma-square-log-check --self-test || exit 1
/tmp/sigma-square-summarize --self-test || exit 1
make -s test-sigma-square-table-native CXX=g++ > "$R/native-field-replay.log" 2>&1 || exit 1

COMMON=(
  BATCH=16 THREADS=256 MINBLOCKS=2 WITNESS=0 STREAM_KARAT=0 SMEM_SPILL=0 GLOBAL_CG=0
  PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_INLINE_POLY=3 PACKED_CLMAD=1
  PACKED_CLMAD_SQUARE=0 PACKED_KARAT3=0 PACKED_WEIGHTED_PREFIX=2
  PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 PACKED_TOP_CLMAD=0
  PACKED_STATE_TILE=256 PACKED_ADD_COMBINE=0 PACKED_PAIR_ILP=0
  UNROLL_SLOTS=1 PACKED_SLOT_PREFETCH=0 PACKED_SLOT_PIPELINE=0
  SIGMA_FUSED=1 SIGMA_FUSED_LATE_Y=0 PACKED_L2_PERSIST=0 PACKED_TOP_HOIST=0
  PACKED_ONB_INV=0 PACKED_FROM_REDUCED=0 PACKED_CLMUL_FLAT=0
  PACKED_PAIR_CLMUL=0 PACKED_ALU_SQR=0 PACKED_ALU_SQUARE=0
  PACKED_SQUARE_TABLE=0 PACKED_INV_POLY=0 PACKED_CHAINS=1
  WALK_TABLE=0 PHASE_PROFILE=0 TABLE_GLOBAL_HINTS=0
)

build_arm() {
  local name=$1 table=$2
  echo "=== build $name SIGMA_SQUARE_TABLE=$table" | tee -a "$R/progress.txt"
  {
    printf 'make -s -B ecc2k130 ARCH=%q' "$ARCH"
    printf ' %q' "${COMMON[@]}" "SIGMA_SQUARE_TABLE=$table"
    printf '\n'
  } > "$R/build-command-$name.txt"
  make -s -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" SIGMA_SQUARE_TABLE="$table" \
    > "$R/build-$name.log" 2>&1 || return 1
  [ -x ecc2k130 ] || return 1
  mv ecc2k130 "client-$name"
  cuobjdump --dump-resource-usage "client-$name" > "$R/resource-$name.txt" || return 1
}
build_arm control 0 || fail=1
build_arm candidate 1 || fail=1
if [ "$fail" != 0 ]; then
  tail -80 "$R"/build-*.log | tee -a "$R/failures.txt"
  exit 1
fi
sha256sum client-control client-candidate > "$R/binary-sha256.txt"

# Bind the running binaries to the same exact static SASS/resource facts as the
# no-GPU admission artifact, without retaining two 65 MiB disassemblies.
mkdir -p /tmp/sigma-runtime-compile
cp client-control /tmp/sigma-runtime-compile/control
cp client-candidate /tmp/sigma-runtime-compile/candidate
cp "$R/build-control.log" /tmp/sigma-runtime-compile/control-build.log
cp "$R/build-candidate.log" /tmp/sigma-runtime-compile/candidate-build.log
cp "$R/resource-control.txt" /tmp/sigma-runtime-compile/control-resources.txt
cp "$R/resource-candidate.txt" /tmp/sigma-runtime-compile/candidate-resources.txt
nvcc --version > /tmp/sigma-runtime-compile/nvcc-version.txt
cp "$R/build-command-control.txt" /tmp/sigma-runtime-compile/common-command.txt
cuobjdump --dump-sass client-control > /tmp/sigma-runtime-compile/control.sass
cuobjdump --dump-sass client-candidate > /tmp/sigma-runtime-compile/candidate.sass
/tmp/sigma-square-compile-audit /tmp/sigma-runtime-compile "$SOURCE_REV" \
  > "$R/runtime-compile-audit.json" || exit 1
sha256sum /tmp/sigma-runtime-compile/control.sass /tmp/sigma-runtime-compile/candidate.sass \
  > "$R/runtime-sass.sha256"

for arm in control candidate; do
  table=0; [ "$arm" = candidate ] && table=1
  for target in test-packed-cuda test-packed-storage-cuda test-shared-sigma-cuda; do
    make -s -B "$target" ARCH="$ARCH" "${COMMON[@]}" SIGMA_SQUARE_TABLE="$table" \
      > "$R/gate-$arm-$target.log" 2>&1 || fail=1
  done
done
if [ "$fail" != 0 ]; then
  echo "device arithmetic/storage/shared-sigma gate failed" | tee -a "$R/failures.txt"
  exit 1
fi

check_log() {
  local arm=$1 mode=$2 log=$3 threads=$4 dp=$5 steps=$6 out=$7
  /tmp/sigma-square-log-check "$arm" "$mode" "$log" "$threads" "$dp" "$steps" "$out"
}

checkpoint_run() {
  local arm=$1 launches=$2 checkpoint=$3 label=$4 rc=0 log
  log="$R/checkpoint-$label.log"
  "./client-$arm" --curve 131 --packed --threads "$PARTIAL_THREADS" --steps 16 \
    --launches "$launches" --verify 0 --run-id "$RUN_ID" --bench \
    --checkpoint "$checkpoint" --checkpoint-every 3600 > "$log" 2>&1 || rc=$?
  [ "$rc" = 0 ] || return 1
  check_log "$arm" checkpoint "$log" "$PARTIAL_THREADS" 0 16 "$R/checkpoint-$label.check"
}
checkpoint_run control 2 "$R/control-whole.ck" control-whole || fail=1
checkpoint_run candidate 2 "$R/candidate-whole.ck" candidate-whole || fail=1
checkpoint_run control 1 "$R/control-prefix.ck" control-prefix || fail=1
cp "$R/control-prefix.ck" "$R/control-to-candidate.ck"
checkpoint_run candidate 1 "$R/control-to-candidate.ck" control-to-candidate || fail=1
checkpoint_run candidate 1 "$R/candidate-prefix.ck" candidate-prefix || fail=1
cp "$R/candidate-prefix.ck" "$R/candidate-to-control.ck"
checkpoint_run control 1 "$R/candidate-to-control.ck" candidate-to-control || fail=1
for checkpoint in candidate-whole.ck control-to-candidate.ck candidate-to-control.ck; do
  cmp "$R/control-whole.ck" "$R/$checkpoint" || fail=1
done
sha256sum "$R"/*.ck > "$R/checkpoint-sha256.txt"

verify_arm() {
  local arm=$1 rc=0 log="$R/verify-$1.log"
  "./client-$arm" --curve 131 --packed --threads "$VERIFY_THREADS" \
    --dp-weight 48 --dp-cap 262144 --steps 95 --launches 7 --verify 300 \
    --run-id "$RUN_ID" --dp-file "$R/dp-$arm.bin" > "$log" 2>&1 || rc=$?
  [ "$rc" = 0 ] || return 1
  check_log "$arm" verify "$log" "$VERIFY_THREADS" 48 95 "$R/verify-$arm.check"
}
verify_arm control || fail=1
verify_arm candidate || fail=1
/tmp/corpus-identity --canonical-out "$R/corpus-sorted.bin" \
  "$R/dp-control.bin" "$R/dp-candidate.bin" > "$R/corpus-identity.txt" || fail=1

if [ "$fail" != 0 ]; then
  echo "correctness, occupancy, partial-block, checkpoint or corpus gate failed; timing suppressed" \
    | tee "$R/preflight.txt"
  exit 1
fi
echo "PASS: all correctness and identity gates" | tee "$R/preflight.txt"

printf 'phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 pair=$2 order=$3 variant=$4 arm rateFile log digest state rc=0
  arm=$variant
  case "$variant" in
    control|control_a|control_b) arm=control ;;
    candidate) arm=candidate ;;
    *) return 2 ;;
  esac
  log="$R/${phase}-${pair}-${order}-${variant}.log"
  rateFile="$R/${phase}-${pair}-${order}-${variant}.rate"
  "./client-$arm" --curve 131 --packed --threads "$BENCH_THREADS" --bench \
    --steps 1024 --launches 32 --verify 0 > "$log" 2>&1 || rc=$?
  [ "$rc" = 0 ] || return 1
  check_log "$arm" timing "$log" "$BENCH_THREADS" 0 1024 "$rateFile" || return 1
  rate=$(cat "$rateFile")
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu \
    --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$phase" "$pair" "$order" "$variant" "$rate" "$digest" "$state" \
    | tee -a "$R/samples.tsv"
}

sample warmup 0 1 control || fail=1
sample warmup 0 2 candidate || fail=1
for pair in 1 2 3 4 5; do
  if [ $((pair % 2)) = 1 ]; then
    sample ab "$pair" 1 control || fail=1
    sample ab "$pair" 2 candidate || fail=1
    sample aa "$pair" 1 control_a || fail=1
    sample aa "$pair" 2 control_b || fail=1
  else
    sample aa "$pair" 1 control_b || fail=1
    sample aa "$pair" 2 control_a || fail=1
    sample ab "$pair" 1 candidate || fail=1
    sample ab "$pair" 2 control || fail=1
  fi
done
if [ "$fail" != 0 ]; then
  echo "timing row failed" | tee -a "$R/failures.txt"
  exit 1
fi
/tmp/sigma-square-summarize "$R/samples.tsv" "$R/preflight.txt" "$R/result.json" \
  | tee "$R/summary.txt" || exit 1
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo "=== done"
