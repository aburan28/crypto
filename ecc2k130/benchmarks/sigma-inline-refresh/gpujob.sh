#!/usr/bin/env bash
# Frozen current-source sigma PACKED_INLINE_POLY=0/3 RTX PRO 6000 comparison.
set -uo pipefail

cd "${WORK:-/work}" || exit 1
R=${RESULTS:-/results}
mkdir -p "$R"
export CUDA_DISABLE_PTX_JIT=1
export CUDA_FORCE_PTX_JIT=0
export CUDA_FORCE_JIT=0
export LC_ALL=C

if ! command -v jq >/dev/null 2>&1; then
  export DEBIAN_FRONTEND=noninteractive
  if ! apt-get update -qq > "$R/apt-jq.log" 2>&1 ||
     ! apt-get install -y -qq jq >> "$R/apt-jq.log" 2>&1; then
    echo 'could not install jq for the native/shell result auditor' | tee -a "$R/failures.txt"
    exit 1
  fi
fi

ARCH='-gencode arch=compute_120,code=sm_120'
WORKERS=385024
BATCH=16
STEPS=1024
LAUNCHES=32
RUN_ID=7
EXPECTED_UPDATES=201863462912
VERIFY_WEIGHT=42
VERIFY_STEPS=24
VERIFY_LAUNCHES=6
fail=0

# Spell every walk-affecting or resource-affecting Make variable out.  This does
# not inherit the table-walk gpu-preset and the only per-arm change appended by
# build_arm is PACKED_INLINE_POLY.
COMMON=(
  BATCH=16 THREADS=256 MINBLOCKS=2 WITNESS=0
  STREAM_KARAT=0 SMEM_SPILL=0 GLOBAL_CG=0
  PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_CLMAD_SQUARE=0
  PACKED_KARAT3=0 PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1
  PACKED_SHARED_SIGMA=1 PACKED_TOP_CLMAD=0 PACKED_STATE_TILE=256
  PACKED_ADD_COMBINE=0 PACKED_PAIR_ILP=0 UNROLL_SLOTS=1
  PACKED_SLOT_PREFETCH=0 PACKED_SLOT_PIPELINE=0 PACKED_L2_PERSIST=0
  PACKED_TOP_HOIST=0 PACKED_ONB_INV=0 PACKED_FROM_REDUCED=0
  PACKED_CLMUL_FLAT=0 PACKED_PAIR_CLMUL=0 PACKED_ALU_SQR=0
  WALK_TABLE=0 TABLE_BRANCHES=8 TABLE_PIVOT_BYTES=0 TABLE_PHASE_POPC=0
  TABLE_GLOBAL=0 TABLE_ADDEND_GLOBAL=0 TABLE_TAG_DENOM=0
  TABLE_FUSED=0 TABLE_FUSED_PIPE=0 PHASE_PROFILE=0 TABLE_PIPE_SELECT=0
  TABLE_SPLIT_FORWARD=0 TABLE_BATCH_HINTS=0 CYCLE_FAST2=0 CYCLE_PROFILE=0
  PACKED_CHAIN_FIRST=0 PACKED_CHAINS=1 PACKED_ALU_SQUARE=0
  PACKED_SQUARE_TABLE=0 PACKED_INV_POLY=0
)

{
  echo "source=${SOURCE_REV:-unknown}"
  echo "image=nvidia/cuda:13.3.1-devel-ubuntu24.04"
  echo "arch=sm_120"
  echo "workers=$WORKERS"
  echo "batch=$BATCH"
  echo "blockThreads=256"
  echo "minBlocks=2"
  echo "steps=$STEPS"
  echo "launches=$LAUNCHES"
  echo "updatesPerRankedRow=$EXPECTED_UPDATES"
  echo "runId=$RUN_ID"
  echo "verificationDpWeight=$VERIFY_WEIGHT"
  echo "verificationSteps=$VERIFY_STEPS"
  echo "verificationLaunches=$VERIFY_LAUNCHES"
  printf 'commonFlags='
  printf '%s ' "${COMMON[@]}"
  printf '\n'
} > "$R/config.txt"

{
  nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.max.sm,clocks.max.memory,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version
  uname -a
} > "$R/host-before.txt" 2>&1 || fail=1

if [ "$(wc -l < "$R/host-before.txt")" -lt 3 ] ||
   ! grep -q '^NVIDIA RTX PRO 6000 Blackwell Server Edition,' "$R/host-before.txt" ||
   ! grep -q 'release 13.3, V13.3.73' "$R/host-before.txt"; then
  echo 'GPU or compiler identity gate failed' | tee -a "$R/failures.txt"
  fail=1
fi

# Compilation source plus the exact job/comparator/auditor that interpret it.
{
  find Makefile include src generated -type f -print0
  printf '%s\0' benchmarks/dp_identity.py \
    benchmarks/sigma-inline-refresh/README.md \
    benchmarks/sigma-inline-refresh/gpujob.sh \
    benchmarks/sigma-inline-refresh/audit.sh
} | sort -z | xargs -0 sha256sum > "$R/source-files.sha256"
sha256sum "$R/source-files.sha256" > "$R/source-manifest.sha256"

build_arm() {
  local name=$1 mode=$2 target binary
  echo "=== BUILD $name inline=$mode" | tee -a "$R/progress.txt"
  if ! make -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" PACKED_INLINE_POLY="$mode" \
      > "$R/build-$name.log" 2>&1; then
    echo "client build failed: $name" | tee -a "$R/failures.txt"
    return 1
  fi
  mv ecc2k130 "client-$name"

  for target in test-packed-cuda test-packed-storage-cuda test-shared-sigma-cuda; do
    case "$target" in
      test-packed-cuda) binary=build/test-packed-cuda ;;
      test-packed-storage-cuda) binary=build/test-packed-storage-cuda ;;
      test-shared-sigma-cuda) binary=build/test-shared-sigma-cuda ;;
    esac
    if ! make -B "$target" ARCH="$ARCH" "${COMMON[@]}" PACKED_INLINE_POLY="$mode" \
        > "$R/build-$name-$target.log" 2>&1; then
      echo "$target build/initial gate failed: $name" | tee -a "$R/failures.txt"
      return 1
    fi
    cp "$binary" "${target#test-}-$name"
  done
}

build_arm control 0 || fail=1
build_arm inline3 3 || fail=1

if [ ! -x client-control ] || [ ! -x client-inline3 ]; then
  echo 'one or more client binaries missing; later gates and timing suppressed' | tee -a "$R/failures.txt"
  exit 1
fi

{
  sha256sum client-control client-inline3
  sha256sum packed-cuda-control packed-cuda-inline3
  sha256sum packed-storage-cuda-control packed-storage-cuda-inline3
  sha256sum shared-sigma-cuda-control shared-sigma-cuda-inline3
} > "$R/binary-sha256.txt"

extract_code() {
  local name=$1 dir="$R/code-$1" cubin
  mkdir -p "$dir"
  cuobjdump --dump-resource-usage "client-$name" > "$R/resources-$name.txt" 2>&1 || return 1
  (cd "$dir" && cuobjdump --extract-elf all "/work/client-$name" > extract.log 2>&1) || return 1
  cubin=$(find "$dir" -maxdepth 1 -type f -name '*.cubin' -print)
  if [ "$(printf '%s\n' "$cubin" | sed '/^$/d' | wc -l)" -ne 1 ]; then
    echo "expected one sm_120 cubin for $name" | tee -a "$R/failures.txt"
    return 1
  fi
  nvdisasm -c "$cubin" > "$R/sass-$name.txt" 2> "$R/sass-$name.stderr" || return 1
  sha256sum "$cubin" "$R/sass-$name.txt" "$R/resources-$name.txt" >> "$R/code-sha256.txt"
  gzip -9 "$R/sass-$name.txt"
  sha256sum "$R/sass-$name.txt.gz" >> "$R/code-sha256.txt"
}

extract_code control || fail=1
extract_code inline3 || fail=1

run_gate() {
  local name=$1 kind=$2 binary log="$R/gate-$1-$2.log" rc=0
  case "$kind" in
    arithmetic) binary="packed-cuda-$name" ;;
    storage) binary="packed-storage-cuda-$name" ;;
    shared-sigma) binary="shared-sigma-cuda-$name" ;;
  esac
  timeout 600 "./$binary" > "$log" 2>&1 || rc=$?
  if [ "$rc" -ne 0 ]; then
    echo "$kind gate failed for $name with $rc" | tee -a "$R/failures.txt"
    return 1
  fi
  case "$kind" in
    arithmetic)
      grep -qx 'packed arithmetic native carryless multiply: 1' "$log" &&
      grep -qx 'packed arithmetic weighted prefix: 2' "$log" &&
      grep -qx 'PASS: 3120 GPU Frobenius vectors, every field basis vector for all selected powers plus dense cases' "$log" &&
      grep -qx 'PASS: 2526 GPU polynomial reductions against long division, including ignored upper-word bits and canonical outputs' "$log" &&
      grep -qx 'PASS: 18194 GPU polynomial products, including all 17161 basis pairs' "$log" &&
      grep -qx 'PASS: 18194 GPU paired polynomial products against independent multiplication' "$log" &&
      grep -qx 'PASS: 1157 GPU polynomial squares against independent multiplication and long division' "$log" &&
      grep -qx 'PASS: 6240 GPU paired Frobenius vectors, both inputs against independent routing' "$log"
      ;;
    storage)
      grep -qx 'packed storage compact state: 1' "$log" &&
      grep -qx 'packed storage batch: 16' "$log" &&
      grep -qx 'PASS: 128 GPU storage cases, 297344 records, independent physical images and logical reads with canaries' "$log"
      ;;
    shared-sigma)
      grep -qx 'packed shared sigma probe: 1' "$log" &&
      grep -qx 'PASS: 21 GPU sigma scenarios, 21036 input pairs, global and selected helpers against independent routing' "$log" &&
      grep -qx 'PASS: 114 complete block mask snapshots, 51072 words, output guards and inactive blocks' "$log"
      ;;
  esac || {
    echo "$kind coverage marker failed for $name" | tee -a "$R/failures.txt"
    return 1
  }
}

for name in control inline3; do
  run_gate "$name" arithmetic || fail=1
  run_gate "$name" storage || fail=1
  run_gate "$name" shared-sigma || fail=1
done

checkpoint_run() {
  local name=$1 launches=$2 checkpoint=$3 label=$4 rc=0
  "./client-$name" --curve 131 --packed --threads 8 --steps 16 --launches "$launches" \
    --verify 0 --run-id "$RUN_ID" --bench --checkpoint "$checkpoint" \
    --checkpoint-every 3600 > "$R/checkpoint-$label.log" 2>&1 || rc=$?
  if [ "$rc" -ne 0 ] || grep -Eq 'MISMATCH|OVERFLOW|stopping:|collision found|solved' "$R/checkpoint-$label.log"; then
    return 1
  fi
}

# Byte-identical uninterrupted and split/resumed checkpoints in both directions.
checkpoint_run control 2 "$R/control-whole.ck" control-whole || fail=1
checkpoint_run inline3 2 "$R/inline3-whole.ck" inline3-whole || fail=1
checkpoint_run control 1 "$R/control-prefix.ck" control-prefix || fail=1
cp "$R/control-prefix.ck" "$R/control-to-inline3.ck"
checkpoint_run inline3 1 "$R/control-to-inline3.ck" control-to-inline3 || fail=1
checkpoint_run inline3 1 "$R/inline3-prefix.ck" inline3-prefix || fail=1
cp "$R/inline3-prefix.ck" "$R/inline3-to-control.ck"
checkpoint_run control 1 "$R/inline3-to-control.ck" inline3-to-control || fail=1
if ! cmp -s "$R/control-whole.ck" "$R/inline3-whole.ck" ||
   ! cmp -s "$R/control-whole.ck" "$R/control-to-inline3.ck" ||
   ! cmp -s "$R/control-whole.ck" "$R/inline3-to-control.ck"; then
  echo 'cross-binary checkpoint identity failed' | tee -a "$R/failures.txt"
  fail=1
fi
sha256sum "$R"/*.ck > "$R/checkpoint-sha256.txt"

verify_arm() {
  local name=$1 mode=$2 log="$R/verify-$1.log" rc=0
  "./client-$name" --curve 131 --packed --threads "$WORKERS" --steps "$VERIFY_STEPS" \
    --launches "$VERIFY_LAUNCHES" --dp-weight "$VERIFY_WEIGHT" --dp-cap 2000000 \
    --verify 300 --run-id "$RUN_ID" --dp-file "$R/dp-$name.bin" > "$log" 2>&1 || rc=$?
  if [ "$rc" -ne 0 ] ||
     ! grep -qx "packed inline polynomial: $mode" "$log" ||
     ! grep -qx 'packed shared sigma: 1' "$log" ||
     ! grep -qx 'packed table walk: 0 (0 branches, 0 shared bytes)' "$log" ||
     ! grep -Eq "^backend cuda-packed131: $WORKERS threads x 16 slots x 1 lanes = 6160384 walks, dp weight $VERIFY_WEIGHT, $VERIFY_STEPS steps per launch$" "$log" ||
     ! grep -Eq '\(300 verified against the reference, 0 dropped\)$' "$log" ||
     grep -Eq 'MISMATCH|OVERFLOW|stopping:|collision found|verified \[k\]P|solved' "$log" ||
     [ ! -s "$R/dp-$name.bin" ]; then
    echo "reference replay/corpus gate failed for $name" | tee -a "$R/failures.txt"
    return 1
  fi
}

verify_arm control 0 || fail=1
verify_arm inline3 3 || fail=1
python3 benchmarks/dp_identity.py "$R" control inline3 > "$R/dp-identity.txt" 2>&1 || fail=1
if ! grep -Eq '^inline3 +[1-9][0-9]* records .* IDENTICAL to control$' "$R/dp-identity.txt"; then
  echo 'format-aware sorted DP identity failed or corpus empty' | tee -a "$R/failures.txt"
  fail=1
fi

if [ "$fail" -ne 0 ]; then
  echo 'pre-timing gate failed; warmups and ranked timing suppressed' | tee -a "$R/failures.txt"
  bash benchmarks/sigma-inline-refresh/audit.sh "$R" --allow-incomplete || true
  exit 1
fi

printf 'phase\tpair\torder\tvariant\trateMps\titerations\tdropped\tlogSha256\tgpuState\n' > "$R/samples.tsv"

check_timing_markers() {
  local log=$1 mode=$2
  grep -qx "packed inline polynomial: $mode" "$log" &&
  grep -qx 'packed denominator cache: 1' "$log" &&
  grep -qx 'packed multiply by value: 1' "$log" &&
  grep -qx 'packed Frobenius network: 3' "$log" &&
  grep -qx 'packed polynomial chain: 1' "$log" &&
  grep -qx 'packed polynomial state: 1' "$log" &&
  grep -qx 'packed unrolled inversion: 1' "$log" &&
  grep -qx 'packed paired products: 1' "$log" &&
  grep -qx 'packed direct reduction: 1' "$log" &&
  grep -qx 'packed generated product: 1' "$log" &&
  grep -qx 'packed native carryless multiply: 1' "$log" &&
  grep -qx 'packed compact state: 1' "$log" &&
  grep -qx 'packed shared sigma: 1' "$log" &&
  grep -qx 'packed state tile: 256' "$log" &&
  grep -qx 'packed table walk: 0 (0 branches, 0 shared bytes)' "$log" &&
  grep -Eq "^backend cuda-packed131: $WORKERS threads x 16 slots x 1 lanes = 6160384 walks, dp weight 0, $STEPS steps per launch$" "$log" &&
  ! grep -Eq 'MISMATCH|OVERFLOW|stopping:|collision found|solved' "$log"
}

sample() {
  local phase=$1 pair=$2 order=$3 name=$4 mode=$5
  local log="$R/${phase}-${pair}-${order}-${name}.log" rc=0 count rate iterations dropped digest state
  "./client-$name" --curve 131 --packed --threads "$WORKERS" --steps "$STEPS" \
    --launches "$LAUNCHES" --verify 0 --run-id "$RUN_ID" --bench > "$log" 2>&1 || rc=$?
  count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  iterations=$(awk '/M it\/s/ && /iterations/ {
    for (i=1;i<=NF;i++) if ($i=="iterations") v=$(i-1)
  } END {print v+0}' "$log")
  dropped=$(sed -nE 's/.*\(([0-9]+) verified against the reference, ([0-9]+) dropped\).*/\2/p' "$log" | tail -1)
  check_timing_markers "$log" "$mode" || rc=1
  if [ "$count" -ne 1 ] || [ "$iterations" -ne "$EXPECTED_UPDATES" ] ||
     [ "${dropped:-1}" -ne 0 ] || ! awk -v r="${rate:-0}" 'BEGIN { exit !(r+0 > 0) }'; then
    rc=1
  fi
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu \
    --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$phase" "$pair" "$order" "$name" "${rate:-0}" "$iterations" \
    "${dropped:-1}" "$digest" "$state" | tee -a "$R/samples.tsv"
  return "$rc"
}

# Exactly two full-work warmups, excluded by the auditor.
sample warmup 0 1 control 0 || fail=1
sample warmup 0 2 inline3 3 || fail=1

# Exactly five equal-work pairs with alternating order AB, BA, AB, BA, AB.
for pair in 1 2 3 4 5; do
  if [ $((pair % 2)) -eq 1 ]; then
    sample ranked "$pair" 1 control 0 || fail=1
    sample ranked "$pair" 2 inline3 3 || fail=1
  else
    sample ranked "$pair" 1 inline3 3 || fail=1
    sample ranked "$pair" 2 control 0 || fail=1
  fi
done

nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.current.sm,clocks.current.memory,power.limit,temperature.gpu --format=csv,noheader \
  > "$R/host-after.txt" 2>&1 || fail=1

bash benchmarks/sigma-inline-refresh/audit.sh "$R" || fail=1
exit "$fail"
