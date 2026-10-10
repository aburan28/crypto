#!/usr/bin/env bash
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export LC_ALL=C CUDA_DISABLE_PTX_JIT=1 CUDA_FORCE_PTX_JIT=0 CUDA_FORCE_JIT=0
R=${RESULTS:-/results}
mkdir -p "$R"

ARCH='-gencode arch=compute_120,code=sm_120'
VERIFY_THREADS=96256
BENCH_THREADS=385024
EXPECTED_UPDATES=403726925824
RUN_ID=33
fail=0

COMMON=(
  BATCH=16 WITNESS=0 STREAM_KARAT=0 SMEM_SPILL=0 GLOBAL_CG=0
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
  PACKED_PAIR_CLMUL=0 PACKED_ALU_SQR=0 WALK_TABLE=0 TABLE_BRANCHES=8
  TABLE_PIVOT_BYTES=0 TABLE_PHASE_POPC=0 TABLE_GLOBAL=0
  TABLE_ADDEND_GLOBAL=0 TABLE_TAG_DENOM=0 TABLE_FUSED=0
  TABLE_FUSED_PIPE=0 PHASE_PROFILE=0 TABLE_PIPE_SELECT=0
  TABLE_SPLIT_FORWARD=0 TABLE_BATCH_HINTS=0 CYCLE_FAST2=0
  CYCLE_PROFILE=0 PACKED_CHAIN_FIRST=0 PACKED_CHAINS=1
  PACKED_ALU_SQUARE=0 PACKED_SQUARE_TABLE=0 PACKED_INV_POLY=0
)

GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
if [ "$GPU_NAME" != 'NVIDIA RTX PRO 6000 Blackwell Server Edition' ] || [ "$CAP" != 120 ] ||
   ! nvcc --version | grep -q 'release 13\.3, V13\.3\.73'; then
  echo "frozen hardware/compiler mismatch: $GPU_NAME sm_$CAP" | tee "$R/failures.txt"
  exit 2
fi
{
  nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version
  uname -a
  echo "source=${SOURCE_REV:-unknown}"
  echo "geometry=t256/min2 versus t512/min1; batch=16 workers=$BENCH_THREADS updates=$EXPECTED_UPDATES"
  printf 'commonFlags='; printf '%s ' "${COMMON[@]}"; printf '\n'
} > "$R/host.txt" 2>&1

g++ -O2 -std=c++17 -Wall -Wextra -Werror benchmarks/sigma-fused/corpus_identity.cpp \
  -o /tmp/corpus-identity || exit 1
g++ -O2 -std=c++17 -Wall -Wextra -Werror benchmarks/sigma-fused/b16-geometry/summarize.cpp \
  -o /tmp/geometry-summarize || exit 1
{
  find Makefile include src generated -type f -print0
  printf '%s\0' benchmarks/sigma-fused/b16-geometry/PROTOCOL.md \
    benchmarks/sigma-fused/b16-geometry/gpujob.sh \
    benchmarks/sigma-fused/b16-geometry/summarize.cpp \
    benchmarks/sigma-fused/b16-geometry/audit.cpp benchmarks/sigma-fused/corpus_identity.cpp
} | sort -z | xargs -0 sha256sum > "$R/source-files.sha256"
sha256sum -c "$R/source-files.sha256" > "$R/source-files.check" 2>&1 || exit 1

build_arm() {
  local name=$1 threads=$2 minblocks=$3
  echo "=== build $name T$threads min$minblocks" | tee -a "$R/progress.txt"
  make -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" THREADS="$threads" MINBLOCKS="$minblocks" \
    > "$R/build-$name.log" 2>&1 || return 1
  mv ecc2k130 "client-$name"
  grep -q -- "-DECC_THREADS=$threads" "$R/build-$name.log" || return 1
  grep -q -- "-DECC_MINBLOCKS=$minblocks" "$R/build-$name.log" || return 1
  grep -q -- '-DECC_PACKED_INLINE_POLY=3' "$R/build-$name.log" || return 1
  grep -q -- '-DECC_PACKED_SHARED_SIGMA=1' "$R/build-$name.log" || return 1
  grep -q -- '-DECC_SIGMA_FUSED=1' "$R/build-$name.log" || return 1
  grep -q -- '-DECC_SIGMA_FUSED_LATE_Y=0' "$R/build-$name.log" || return 1
  grep -q -- '-DECC_WITNESS=0' "$R/build-$name.log" || return 1
  if ! awk '
    /Function properties for _ZN12eccPacked1314walkE10WalkParamsIjEPj/ {
      getline
      if ($0 == "    0 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads") ok++
    }
    END { exit ok == 1 ? 0 : 1 }
  ' "$R/build-$name.log"; then return 1; fi
  for target in test-packed-cuda test-packed-storage-cuda test-shared-sigma-cuda; do
    make -B "$target" ARCH="$ARCH" "${COMMON[@]}" THREADS="$threads" MINBLOCKS="$minblocks" \
      > "$R/gate-$name-$target.log" 2>&1 || return 1
  done
}
build_arm t256 256 2 || fail=1
build_arm t512 512 1 || fail=1
if [ "$fail" != 0 ] || [ ! -x client-t256 ] || [ ! -x client-t512 ]; then exit 1; fi
sha256sum client-t256 client-t512 > "$R/binary-sha256.txt"
sha256sum -c "$R/binary-sha256.txt" > "$R/binary-sha256.check" 2>&1 || exit 1

checkpoint_run() {
  local binary=$1 launches=$2 checkpoint=$3 label=$4 rc=0
  "./client-$binary" --curve 131 --packed --threads 512 --steps 16 --launches "$launches" \
    --verify 0 --run-id "$RUN_ID" --bench --checkpoint "$checkpoint" \
    --checkpoint-every 3600 > "$R/checkpoint-$label.log" 2>&1 || rc=$?
  [ "$rc" = 0 ] && ! grep -Eq 'MISMATCH|OVERFLOW|stopping:|collision found|solved' \
    "$R/checkpoint-$label.log"
}
checkpoint_run t256 2 "$R/t256-whole.ck" t256-whole || fail=1
checkpoint_run t512 2 "$R/t512-whole.ck" t512-whole || fail=1
checkpoint_run t256 1 "$R/t256-prefix.ck" t256-prefix || fail=1
cp "$R/t256-prefix.ck" "$R/t256-to-t512.ck"
checkpoint_run t512 1 "$R/t256-to-t512.ck" t256-to-t512 || fail=1
checkpoint_run t512 1 "$R/t512-prefix.ck" t512-prefix || fail=1
cp "$R/t512-prefix.ck" "$R/t512-to-t256.ck"
checkpoint_run t256 1 "$R/t512-to-t256.ck" t512-to-t256 || fail=1
if ! cmp -s "$R/t256-whole.ck" "$R/t512-whole.ck" ||
   ! cmp -s "$R/t256-whole.ck" "$R/t256-to-t512.ck" ||
   ! cmp -s "$R/t256-whole.ck" "$R/t512-to-t256.ck"; then fail=1; fi
sha256sum "$R"/*.ck > "$R/checkpoint-sha256.txt"

verify_arm() {
  local name=$1 expectedBlocks=$2 expectedThreads=$3 log="$R/verify-$1.log" rc=0
  "./client-$name" --curve 131 --packed --dp-weight 48 --dp-cap 262144 \
    --steps 95 --launches 7 --verify 300 --run-id "$RUN_ID" \
    --dp-file "$R/dp-$name.bin" > "$log" 2>&1 || rc=$?
  [ "$rc" = 0 ] || return 1
  grep -qx 'packed sigma fused: 1' "$log" || return 1
  grep -qx 'packed sigma fused late y: 0' "$log" || return 1
  grep -qx 'packed shared sigma: 1' "$log" || return 1
  grep -qx "packed launch bounds: $expectedThreads threads, $expectedBlocks min blocks" "$log" || return 1
  grep -Eq "^device: NVIDIA RTX PRO 6000 Blackwell Server Edition, 188 SMs, $expectedBlocks block\(s\) of $expectedThreads packed threads resident per SM$" "$log" || return 1
  grep -Eq '^packed kernel: [0-9]+ registers/thread, 0 local bytes/thread, 1792 shared bytes/block,' "$log" || return 1
  grep -Eq "^backend cuda-packed131: $VERIFY_THREADS threads x 16 slots x 1 lanes = 1540096 walks, dp weight 48, 95 steps per launch$" "$log" || return 1
  grep -Eq '\(300 verified against the reference, 0 dropped\)$' "$log" || return 1
  ! grep -Eq 'MISMATCH|OVERFLOW|stopping:|collision found|solved' "$log"
}
verify_arm t256 2 256 || fail=1
verify_arm t512 1 512 || fail=1
/tmp/corpus-identity --canonical-out "$R/corpus-sorted.bin" "$R/dp-t256.bin" "$R/dp-t512.bin" \
  | tee "$R/corpus-identity.txt" || fail=1
if [ "$fail" != 0 ]; then echo 'preflight failed; timing suppressed' | tee "$R/preflight.txt"; exit 1; fi
sha256sum "$R/corpus-sorted.bin" > "$R/corpus-sorted.sha256"
wc -c "$R/corpus-sorted.bin" > "$R/corpus-sorted.bytes"
rm -f "$R/dp-t256.bin" "$R/dp-t512.bin"
echo 'PASS: arithmetic, storage, shared sigma, resources, occupancy, replay, checkpoints and corpus identity' \
  | tee "$R/preflight.txt"

printf 'phase\tpair\torder\tvariant\trateMps\titerations\tdropped\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 pair=$2 order=$3 variant=$4 binary rc=0
  local log="$R/${phase}-${pair}-${order}-${variant}.log" count rate iterations dropped digest state
  case "$variant" in
    control|control_a|control_b) binary=t256 ;;
    candidate) binary=t512 ;;
    *) return 2 ;;
  esac
  "./client-$binary" --curve 131 --packed --threads "$BENCH_THREADS" --bench \
    --steps 1024 --launches 64 --verify 0 --run-id "$RUN_ID" > "$log" 2>&1 || rc=$?
  count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  iterations=$(awk '/M it\/s/ && /iterations/ {for(i=1;i<=NF;i++) if($i=="iterations") v=$(i-1)} END{print v+0}' "$log")
  dropped=$(sed -nE 's/.*\(([0-9]+) verified against the reference, ([0-9]+) dropped\).*/\2/p' "$log" | tail -1)
  grep -qx 'packed sigma fused: 1' "$log" || rc=1
  grep -qx 'packed sigma fused late y: 0' "$log" || rc=1
  grep -qx 'packed shared sigma: 1' "$log" || rc=1
  case "$binary" in
    t256) grep -qx 'packed launch bounds: 256 threads, 2 min blocks' "$log" || rc=1 ;;
    t512) grep -qx 'packed launch bounds: 512 threads, 1 min blocks' "$log" || rc=1 ;;
  esac
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
for pair in 1 2 3 4 5; do
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
awk -F '\t' -v root="$R" 'NR > 1 {
  print $8 "  " root "/" $1 "-" $2 "-" $3 "-" $4 ".log"
}' "$R/samples.tsv" > "$R/sample-files.sha256"
sha256sum -c "$R/sample-files.sha256" > "$R/sample-files.check" 2>&1 || exit 1
/tmp/geometry-summarize "$R/samples.tsv" "$R/preflight.txt" "$R/result.json" \
  | tee "$R/summary.txt" || exit 1
g++ -O2 -std=c++17 -Wall -Wextra -Werror benchmarks/sigma-fused/b16-geometry/audit.cpp \
  -o /tmp/geometry-audit || exit 1
/tmp/geometry-audit "$R" > "$R/native-audit.json" || exit 1
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo '=== done'
