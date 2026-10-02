#!/usr/bin/env bash
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export LC_ALL=C
export CUDA_DISABLE_PTX_JIT=1 CUDA_FORCE_PTX_JIT=0 CUDA_FORCE_JIT=0
R=${RESULTS:-/results}
mkdir -p "$R"

ARCH='-gencode arch=compute_120,code=sm_120'
VERIFY_THREADS=96256
BENCH_THREADS=385024
EXPECTED_UPDATES=403726925824
RUN_ID=31
TABLE_SHA=0d385683efbc15afd4f31618ed569ff2b076722881a51091ddf51adca827e047
fail=0

COMMON=(
  BATCH=16 THREADS=512 MINBLOCKS=1 WITNESS=0
  STREAM_KARAT=0 SMEM_SPILL=0 GLOBAL_CG=0
  PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_INLINE_POLY=3 PACKED_CLMAD=1
  PACKED_CLMAD_SQUARE=0 PACKED_KARAT3=0 PACKED_WEIGHTED_PREFIX=2
  PACKED_COMPACT_STATE=1 PACKED_TOP_CLMAD=0 PACKED_STATE_TILE=256
  PACKED_ADD_COMBINE=0 PACKED_PAIR_ILP=0 UNROLL_SLOTS=1
  PACKED_SLOT_PREFETCH=0 PACKED_SLOT_PIPELINE=0 SIGMA_FUSED=1
  PACKED_L2_PERSIST=0 PACKED_TOP_HOIST=0 PACKED_ONB_INV=0
  PACKED_FROM_REDUCED=0 PACKED_CLMUL_FLAT=0 PACKED_PAIR_CLMUL=0
  PACKED_ALU_SQR=0 WALK_TABLE=0 TABLE_BRANCHES=8 TABLE_PIVOT_BYTES=0
  TABLE_PHASE_POPC=0 TABLE_GLOBAL=0 TABLE_ADDEND_GLOBAL=0
  TABLE_TAG_DENOM=0 TABLE_FUSED=0 TABLE_FUSED_PIPE=0 PHASE_PROFILE=0
  TABLE_PIPE_SELECT=0 TABLE_SPLIT_FORWARD=0 TABLE_BATCH_HINTS=0
  CYCLE_FAST2=0 CYCLE_PROFILE=0 PACKED_CHAIN_FIRST=0 PACKED_CHAINS=1
  PACKED_ALU_SQUARE=0 PACKED_SQUARE_TABLE=0 PACKED_INV_POLY=0
)

GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
if [ "$GPU_NAME" != 'NVIDIA RTX PRO 6000 Blackwell Server Edition' ] || [ "$CAP" != 120 ] ||
   ! nvcc --version | grep -q 'release 13\.3, V13\.3\.73'; then
  echo "frozen hardware/compiler mismatch: $GPU_NAME sm_$CAP" | tee "$R/failures.txt"
  nvcc --version >> "$R/failures.txt" 2>&1
  exit 2
fi

{
  nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version
  uname -a
  echo "source=${SOURCE_REV:-unknown}"
  echo "geometry=B16/T512/min1 workers=$BENCH_THREADS updates=$EXPECTED_UPDATES"
  printf 'commonFlags='; printf '%s ' "${COMMON[@]}"; printf '\n'
} > "$R/host.txt" 2>&1

# Rebuild the exact map from repository arithmetic before compiling CUDA.
mkdir -p /tmp/direct-sigma
g++ -O3 -std=c++17 -Wall -Wextra -Werror -DECC_PACKED_PERM_SIGMA=1 \
  benchmarks/direct-sigma-map/synthesize.cpp -o /tmp/direct-sigma/synthesize || exit 1
/tmp/direct-sigma/synthesize /tmp/direct-sigma/table3.h /tmp/direct-sigma/half5.h \
  /tmp/direct-sigma/gpu.cuh > "$R/native-synthesis.txt" || exit 1
if ! cmp -s /tmp/direct-sigma/gpu.cuh include/packeddirectsigma131.cuh ||
   [ "$(sha256sum include/packeddirectsigma131.cuh | awk '{print $1}')" != "$TABLE_SHA" ]; then
  echo 'generated direct-sigma table/header identity failed' | tee -a "$R/failures.txt"
  exit 1
fi

g++ -O2 -std=c++17 -Wall -Wextra -Werror benchmarks/sigma-fused/corpus_identity.cpp \
  -o /tmp/corpus-identity || exit 1
g++ -O2 -std=c++17 -Wall -Wextra -Werror benchmarks/direct-sigma-map/summarize_gpu.cpp \
  -o /tmp/direct-sigma-summarize || exit 1

{
  find Makefile include src generated -type f -print0
  printf '%s\0' benchmarks/direct-sigma-map/PROTOCOL.md \
    benchmarks/direct-sigma-map/GPU-PROTOCOL.md benchmarks/direct-sigma-map/gpujob.sh \
    benchmarks/direct-sigma-map/synthesize.cpp benchmarks/direct-sigma-map/summarize_gpu.cpp \
    benchmarks/sigma-fused/corpus_identity.cpp
} | sort -z | xargs -0 sha256sum > "$R/source-files.sha256"

build_arm() {
  local name=$1 direct=$2 shared=$3
  echo "=== build $name direct=$direct shared=$shared" | tee -a "$R/progress.txt"
  make -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" DIRECT_SIGMA="$direct" \
    PACKED_SHARED_SIGMA="$shared" > "$R/build-$name.log" 2>&1 || return 1
  mv ecc2k130 "client-$name"
  if grep -Eq '[1-9][0-9]* bytes spill (stores|loads)' "$R/build-$name.log"; then
    echo "nonzero ptxas spill for $name" | tee -a "$R/failures.txt"; return 1
  fi
  make -B test-packed-cuda ARCH="$ARCH" "${COMMON[@]}" DIRECT_SIGMA="$direct" \
    PACKED_SHARED_SIGMA="$shared" > "$R/gate-$name-arithmetic.log" 2>&1 || return 1
  make -B test-packed-storage-cuda ARCH="$ARCH" "${COMMON[@]}" DIRECT_SIGMA="$direct" \
    PACKED_SHARED_SIGMA="$shared" > "$R/gate-$name-storage.log" 2>&1 || return 1
  if [ "$shared" = 1 ]; then
    make -B test-shared-sigma-cuda ARCH="$ARCH" "${COMMON[@]}" DIRECT_SIGMA=0 \
      PACKED_SHARED_SIGMA=1 > "$R/gate-$name-shared-sigma.log" 2>&1 || return 1
  fi
  if [ "$direct" = 1 ]; then
    make -B test-direct-sigma-cuda ARCH="$ARCH" "${COMMON[@]}" DIRECT_SIGMA=1 \
      PACKED_SHARED_SIGMA=0 > "$R/gate-$name-direct-sigma.log" 2>&1 || return 1
    grep -qx 'PASS: 5144 GPU direct-sigma cases (1048 basis and 4096 dense)' \
      "$R/gate-$name-direct-sigma.log" || return 1
  fi
}

build_arm control 0 1 || fail=1
build_arm candidate 1 0 || fail=1
if [ "$fail" != 0 ] || [ ! -x client-control ] || [ ! -x client-candidate ]; then
  tail -80 "$R"/build-*.log 2>/dev/null || true
  exit 1
fi
sha256sum client-control client-candidate > "$R/binary-sha256.txt"

# Bind each executable to its single intended knob pair through build commands.
grep -q -- '-DECC_DIRECT_SIGMA=0' "$R/build-control.log" || fail=1
grep -q -- '-DECC_PACKED_SHARED_SIGMA=1' "$R/build-control.log" || fail=1
grep -q -- '-DECC_DIRECT_SIGMA=1' "$R/build-candidate.log" || fail=1
grep -q -- '-DECC_PACKED_SHARED_SIGMA=0' "$R/build-candidate.log" || fail=1

checkpoint_run() {
  local binary=$1 launches=$2 checkpoint=$3 label=$4 rc=0
  "./client-$binary" --curve 131 --packed --threads 512 --steps 16 --launches "$launches" \
    --verify 0 --run-id "$RUN_ID" --bench --checkpoint "$checkpoint" \
    --checkpoint-every 3600 > "$R/checkpoint-$label.log" 2>&1 || rc=$?
  [ "$rc" = 0 ] && ! grep -Eq 'MISMATCH|OVERFLOW|stopping:|collision found|solved' \
    "$R/checkpoint-$label.log"
}
checkpoint_run control 2 "$R/control-whole.ck" control-whole || fail=1
checkpoint_run candidate 2 "$R/candidate-whole.ck" candidate-whole || fail=1
checkpoint_run control 1 "$R/control-prefix.ck" control-prefix || fail=1
cp "$R/control-prefix.ck" "$R/control-to-candidate.ck"
checkpoint_run candidate 1 "$R/control-to-candidate.ck" control-to-candidate || fail=1
checkpoint_run candidate 1 "$R/candidate-prefix.ck" candidate-prefix || fail=1
cp "$R/candidate-prefix.ck" "$R/candidate-to-control.ck"
checkpoint_run control 1 "$R/candidate-to-control.ck" candidate-to-control || fail=1
if ! cmp -s "$R/control-whole.ck" "$R/candidate-whole.ck" ||
   ! cmp -s "$R/control-whole.ck" "$R/control-to-candidate.ck" ||
   ! cmp -s "$R/control-whole.ck" "$R/candidate-to-control.ck"; then
  echo 'cross-binary checkpoint identity failed' | tee -a "$R/failures.txt"; fail=1
fi
sha256sum "$R"/*.ck > "$R/checkpoint-sha256.txt"

verify_arm() {
  local name=$1 direct=$2 shared=$3 log="$R/verify-$1.log" rc=0
  "./client-$name" --curve 131 --packed --threads "$VERIFY_THREADS" \
    --dp-weight 48 --dp-cap 262144 --steps 95 --launches 7 --verify 300 \
    --run-id "$RUN_ID" --dp-file "$R/dp-$name.bin" > "$log" 2>&1 || rc=$?
  [ "$rc" = 0 ] || return 1
  grep -qx "packed direct sigma: $direct" "$log" || return 1
  grep -qx "packed shared sigma: $shared" "$log" || return 1
  grep -qx 'packed sigma fused: 1' "$log" || return 1
  grep -Eq '^device: NVIDIA RTX PRO 6000 Blackwell Server Edition, 188 SMs, 1 block\(s\) of 512 packed threads resident per SM$' "$log" || return 1
  grep -Eq '^packed kernel: [0-9]+ registers/thread, 0 local bytes/thread, (0|1792) shared bytes/block,' "$log" || return 1
  if [ "$direct" = 1 ]; then
    grep -qx 'packed direct sigma shared bytes: 56320' "$log" || return 1
    grep -qx 'packed table walk: 0 (0 branches, 56320 shared bytes)' "$log" || return 1
  else
    grep -qx 'packed table walk: 0 (0 branches, 0 shared bytes)' "$log" || return 1
  fi
  grep -Eq "^backend cuda-packed131: $VERIFY_THREADS threads x 16 slots x 1 lanes = 1540096 walks, dp weight 48, 95 steps per launch$" "$log" || return 1
  grep -Eq '\(300 verified against the reference, 0 dropped\)$' "$log" || return 1
  ! grep -Eq 'MISMATCH|OVERFLOW|stopping:|collision found|solved' "$log"
}
verify_arm control 0 1 || fail=1
verify_arm candidate 1 0 || fail=1
/tmp/corpus-identity "$R/dp-control.bin" "$R/dp-candidate.bin" \
  | tee "$R/corpus-identity.txt" || fail=1

if [ "$fail" != 0 ]; then
  echo 'correctness/resource preflight failed; timing suppressed' | tee "$R/preflight.txt"
  exit 1
fi
echo 'PASS: native/GPU map, replay, checkpoints, resources and sorted corpus identity' \
  | tee "$R/preflight.txt"

printf 'phase\tpair\torder\tvariant\trateMps\titerations\tdropped\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 pair=$2 order=$3 variant=$4 binary=$variant direct shared
  local log="$R/${phase}-${pair}-${order}-${variant}.log" rc=0 count rate iterations dropped digest state
  case "$variant" in
    control|control_a|control_b) binary=control; direct=0; shared=1 ;;
    candidate) direct=1; shared=0 ;;
    *) return 2 ;;
  esac
  "./client-$binary" --curve 131 --packed --threads "$BENCH_THREADS" --bench \
    --steps 1024 --launches 64 --verify 0 --run-id "$RUN_ID" > "$log" 2>&1 || rc=$?
  count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  iterations=$(awk '/M it\/s/ && /iterations/ {for(i=1;i<=NF;i++) if($i=="iterations") v=$(i-1)} END{print v+0}' "$log")
  dropped=$(sed -nE 's/.*\(([0-9]+) verified against the reference, ([0-9]+) dropped\).*/\2/p' "$log" | tail -1)
  grep -qx "packed direct sigma: $direct" "$log" || rc=1
  grep -qx "packed shared sigma: $shared" "$log" || rc=1
  grep -qx 'packed sigma fused: 1' "$log" || rc=1
  grep -Eq '^device: NVIDIA RTX PRO 6000 Blackwell Server Edition, 188 SMs, 1 block\(s\) of 512 packed threads resident per SM$' "$log" || rc=1
  grep -Eq "^backend cuda-packed131: $BENCH_THREADS threads x 16 slots x 1 lanes = 6160384 walks, dp weight 0, 1024 steps per launch$" "$log" || rc=1
  ! grep -Eq 'MISMATCH|OVERFLOW|stopping:|collision found|solved' "$log" || rc=1
  if [ "$count" != 1 ] || [ "$iterations" != "$EXPECTED_UPDATES" ] ||
     [ "${dropped:-1}" != 0 ] || ! awk -v r="${rate:-0}" 'BEGIN{exit !(r+0>0)}'; then rc=1; fi
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$phase" "$pair" "$order" \
    "$variant" "${rate:-0}" "$iterations" "${dropped:-1}" "$digest" "$state" | tee -a "$R/samples.tsv"
  return "$rc"
}

sample warmup 0 1 control || fail=1
sample warmup 0 2 candidate || fail=1
sample warmup 0 3 candidate || fail=1
sample warmup 0 4 control || fail=1
for pair in 1 2 3; do
  if [ $((pair % 2)) = 1 ]; then
    sample aa "$pair" 1 control_a || fail=1; sample aa "$pair" 2 control_b || fail=1
  else
    sample aa "$pair" 1 control_b || fail=1; sample aa "$pair" 2 control_a || fail=1
  fi
done
for pair in 1 2 3 4 5; do
  if [ $((pair % 2)) = 1 ]; then
    sample ab "$pair" 1 control || fail=1; sample ab "$pair" 2 candidate || fail=1
  else
    sample ab "$pair" 1 candidate || fail=1; sample ab "$pair" 2 control || fail=1
  fi
done
if [ "$fail" != 0 ]; then echo 'timing row failed' | tee -a "$R/failures.txt"; exit 1; fi
/tmp/direct-sigma-summarize "$R/samples.tsv" "$R/preflight.txt" "$R/result.json" \
  | tee "$R/summary.txt" || exit 1
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo '=== done'
