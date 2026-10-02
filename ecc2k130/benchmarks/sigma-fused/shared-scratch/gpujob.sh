#!/bin/bash
# Frozen four-arm fused-sigma shared-scratch panel. Benchmark only.
set -euo pipefail
cd "${WORK:-/work}"
export LC_ALL=C CUDA_DISABLE_PTX_JIT=1 CUDA_FORCE_PTX_JIT=0 CUDA_FORCE_JIT=0
R=${RESULTS:-/results}
mkdir -p "$R"

mapfile -t GPU_NAMES < <(nvidia-smi --query-gpu=name --format=csv,noheader)
[ "${#GPU_NAMES[@]}" = 1 ]
GPU_NAME=${GPU_NAMES[0]}
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | tr -d '. ')
[ "$GPU_NAME" = 'NVIDIA RTX PRO 6000 Blackwell Server Edition' ] && [ "$CAP" = 120 ]
nvcc --version | grep -q 'release 13\.3, V13\.3\.73'
[[ ${SOURCE_REV:-} =~ ^[0-9a-f]{40}$ ]]
if [ -f SOURCE_REV ]; then [ "$(cat SOURCE_REV)" = "$SOURCE_REV" ]; fi
if [ -e .git ]; then
    [ "$(git rev-parse HEAD)" = "$SOURCE_REV" ]
    [ -z "$(git status --porcelain -- .)" ]
fi

ARCH='-gencode arch=compute_120,code=sm_120'
VERIFY_THREADS=96256
BENCH_THREADS=385024
UPDATES=201863462912
RUN_ID=41
ARMS='control cache2 cache3 cache4'
COMMON=(
    BATCH=16 THREADS=256 MINBLOCKS=2 WITNESS=0 WALK_TABLE=0
    STREAM_KARAT=0 SMEM_SPILL=0 GLOBAL_CG=0
    PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
    PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
    PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
    PACKED_GENERATED_PRODUCT=1 PACKED_INLINE_POLY=3 PACKED_CLMAD=1
    PACKED_CLMAD_SQUARE=0 PACKED_KARAT3=0 PACKED_WEIGHTED_PREFIX=2
    PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 PACKED_TOP_CLMAD=0
    PACKED_STATE_TILE=256 PACKED_ADD_COMBINE=0 PACKED_PAIR_ILP=0
    UNROLL_SLOTS=1 PACKED_SLOT_PREFETCH=0 PACKED_SLOT_PIPELINE=0
    SIGMA_FUSED=1 SIGMA_FUSED_LATE_Y=0 PACKED_L2_PERSIST=0
    PACKED_TOP_HOIST=0 PACKED_ONB_INV=0 PACKED_FROM_REDUCED=0
    PACKED_CLMUL_FLAT=0 PACKED_PAIR_CLMUL=0 PACKED_ALU_SQR=0
    PACKED_ALU_SQUARE=0 PACKED_SQUARE_TABLE=0 PACKED_INV_POLY=0
    TABLE_PHASE_POPC=0 TABLE_GLOBAL=0 TABLE_ADDEND_GLOBAL=0
    TABLE_TAG_DENOM=0 TABLE_FUSED=0 TABLE_FUSED_PIPE=0
    PHASE_PROFILE=0 TABLE_PIPE_SELECT=0 TABLE_SPLIT_FORWARD=0
    TABLE_BATCH_HINTS=0 TABLE_BLOCK_HINTS=0 TABLE_GLOBAL_HINTS=0
    CYCLE_FAST2=0 CYCLE_PROFILE=0 PACKED_CHAIN_FIRST=0 PACKED_CHAINS=1
)

slots() {
    case "$1" in
        control) echo 0 ;;
        cache2) echo 2 ;;
        cache3) echo 3 ;;
        cache4) echo 4 ;;
        *) return 1 ;;
    esac
}
shared_bytes() {
    echo $((1792 + 10240 * $(slots "$1")))
}

{
    nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
    nvcc --version
    echo "source: $SOURCE_REV"
    echo "geometry: B16/T256/minBlocks2; workers $BENCH_THREADS; updates $UPDATES"
    printf 'common flags:'; printf ' %s' "${COMMON[@]}"; printf '\n'
    for arm in $ARMS; do echo "$arm: SIGMA_FUSED_SHARED_SLOTS=$(slots "$arm")"; done
} | tee "$R/host.txt"

g++ -O2 -std=c++17 -Wall -Wextra -Werror \
    benchmarks/sigma-fused/shared-scratch/model.cpp -o /tmp/scratch-model
g++ -O2 -std=c++17 -Wall -Wextra -Werror -Wno-unknown-pragmas \
    benchmarks/sigma-fused/shared-scratch/test_native.cpp -o /tmp/scratch-native
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
    benchmarks/sigma-fused/shared-scratch/summarize.cpp -o /tmp/scratch-summarize
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
    benchmarks/sigma-fused/shared-scratch/log_check.cpp -o /tmp/scratch-log-check
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
    benchmarks/sigma-fused/corpus_identity.cpp -o /tmp/corpus-identity
{
    /tmp/scratch-model --self-test
    /tmp/scratch-native
    /tmp/scratch-summarize --self-test
    /tmp/scratch-log-check --self-test
} > "$R/native-controls.log" 2>&1

{
    find Makefile include generated src -type f -print
    printf '%s\n' benchmarks/sigma-fused/corpus_identity.cpp \
        benchmarks/sigma-fused/shared-scratch/PROPOSAL.md \
        benchmarks/sigma-fused/shared-scratch/PROTOCOL.md \
        benchmarks/sigma-fused/shared-scratch/IMPLEMENTATION.md \
        benchmarks/sigma-fused/shared-scratch/model.cpp \
        benchmarks/sigma-fused/shared-scratch/model.json \
        benchmarks/sigma-fused/shared-scratch/test_native.cpp \
        benchmarks/sigma-fused/shared-scratch/compile_audit.cpp \
        benchmarks/sigma-fused/shared-scratch/compilecheck.sh \
        benchmarks/sigma-fused/shared-scratch/gpujob.sh \
        benchmarks/sigma-fused/shared-scratch/summarize.cpp \
        benchmarks/sigma-fused/shared-scratch/log_check.cpp
} | LC_ALL=C sort -u | xargs sha256sum > "$R/source-files.sha256"

build() {
    local arm=$1 cached
    cached=$(slots "$arm")
    make -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" \
        SIGMA_FUSED_SHARED_SLOTS="$cached" > "$R/build-$arm.log" 2>&1
    grep -q -- "-DECC_SIGMA_FUSED_SHARED_SLOTS=$cached" "$R/build-$arm.log"
    awk '
      /Function properties for _ZN12eccPacked1314walkE10WalkParamsIjEPj/ {
        getline; if ($0 == "    0 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads") ok++
      }
      END { exit ok == 1 ? 0 : 1 }
    ' "$R/build-$arm.log"
    mv ecc2k130 "ecc2k130-$arm"
    cp "ecc2k130-$arm" "$R/"
}
for arm in $ARMS; do build "$arm"; done
sha256sum "$R"/ecc2k130-* > "$R/binary-sha256.txt"

# Existing arithmetic/layout controls use the control build. Scratch-specific
# device controls compile and run each candidate against the production walk.
for target in test-packed-cuda test-packed-storage-cuda test-shared-sigma-cuda; do
    make -B "$target" ARCH="$ARCH" "${COMMON[@]}" SIGMA_FUSED_SHARED_SLOTS=0 \
        > "$R/gate-$target.log" 2>&1
done
for cached in 2 3 4; do
    make -B "build/test-sigma-fused-shared-scratch-cuda-$cached" \
        ARCH="$ARCH" "${COMMON[@]}" > "$R/device-build-$cached.log" 2>&1
    "./build/test-sigma-fused-shared-scratch-cuda-$cached" \
        > "$R/device-$cached.log" 2>&1
    expected=$((1792 + 10240 * cached))
    grep -Eq "^shared scratch resources: slots $cached, registers ([1-9][0-9]?|1[01][0-9]|12[0-8]), local 0, static shared $expected, reserved shared [0-9]+, shared/SM [0-9]+, active blocks/SM 2$" \
        "$R/device-$cached.log"
    grep -Eq '^PASS: 8 device scenarios, 73856 compact/shared field records, [0-9]+ tagged denominators, 21 blocks; canaries and partial blocks intact$' \
        "$R/device-$cached.log"
done

markers() {
    local log=$1 arm=$2 workers=$3 dp=$4 steps=$5
    local cached expectedShared
    cached=$(slots "$arm")
    expectedShared=$(shared_bytes "$arm")
    grep -qx 'packed sigma fused: 1' "$log" || return 1
    grep -qx 'packed sigma fused late y: 0' "$log" || return 1
    grep -qx "packed sigma fused shared slots: $cached" "$log" || return 1
    grep -qx 'packed shared sigma: 1' "$log" || return 1
    grep -qx 'packed compact state: 1' "$log" || return 1
    grep -qx 'packed launch bounds: 256 threads, 2 min blocks' "$log" || return 1
    grep -Eq '^device: NVIDIA RTX PRO 6000 Blackwell Server Edition, 188 SMs, 2 block\(s\) of 256 packed threads resident per SM$' "$log" || return 1
    grep -Eq "^packed kernel: ([1-9][0-9]?|1[01][0-9]|12[0-8]) registers/thread, 0 local bytes/thread, $expectedShared shared bytes/block, single-product multiplier$" "$log" || return 1
    grep -qx "backend cuda-packed131: $workers threads x 16 slots x 1 lanes = $((workers * 16)) walks, dp weight $dp, $steps steps per launch" "$log" || return 1
    if grep -Eq 'MISMATCH|OVERFLOW|unusable|collision found|solved|stopping:' "$log"; then
        return 1
    fi
    if [ -f "$R/resource-$arm.txt" ]; then
        diff "$R/resource-$arm.txt" <(grep '^packed kernel:' "$log") || return 1
    fi
}

verify() {
    local arm=$1 log="$R/verify-$1.log"
    "./ecc2k130-$arm" --curve 131 --packed --threads "$VERIFY_THREADS" \
        --dp-weight 48 --dp-cap 262144 --steps 95 --launches 7 \
        --verify 300 --run-id "$RUN_ID" --dp-file "$R/dp-$arm.bin" \
        > "$log" 2>&1
    markers "$log" "$arm" "$VERIFY_THREADS" 48 95
    /tmp/scratch-log-check "$log" $((VERIFY_THREADS * 16 * 95 * 7)) 300 -1
}
for arm in $ARMS; do verify "$arm"; done
for arm in $ARMS; do grep '^packed kernel:' "$R/verify-$arm.log" > "$R/resource-$arm.txt"; done
/tmp/corpus-identity --canonical-out "$R/corpus-sorted.bin" \
    "$R/dp-control.bin" "$R/dp-cache2.bin" "$R/dp-cache3.bin" "$R/dp-cache4.bin" \
    > "$R/corpus-identity.txt"

# One launch at each frozen boundary length. Every arm must leave exactly the
# same checkpoint; pchain and denominators are intentionally absent from it.
for boundary in 1 2 3 7 16 95; do
    expected=$((513 * 16 * boundary))
    for arm in $ARMS; do
        log="$R/boundary-$boundary-$arm.log"
        "./ecc2k130-$arm" --curve 131 --packed --threads 513 --bench \
            --steps "$boundary" --launches 1 --verify 0 --run-id "$RUN_ID" \
            --checkpoint "$R/boundary-$boundary-$arm.ckpt" --checkpoint-every 3600 \
            > "$log" 2>&1
        markers "$log" "$arm" 513 0 "$boundary"
        /tmp/scratch-log-check "$log" "$expected" 0 -1
    done
    for arm in cache2 cache3 cache4; do
        cmp "$R/boundary-$boundary-control.ckpt" "$R/boundary-$boundary-$arm.ckpt"
    done
done

# Four-launch prefixes and every three-launch continuation at a partial block.
for arm in $ARMS; do
    "./ecc2k130-$arm" --curve 131 --packed --threads 513 --dp-weight 48 \
        --dp-cap 262144 --steps 95 --launches 4 --verify 300 --run-id "$RUN_ID" \
        --checkpoint "$R/prefix-$arm.ckpt" --dp-file "$R/prefix-$arm.bin" \
        > "$R/prefix-$arm.log" 2>&1
    markers "$R/prefix-$arm.log" "$arm" 513 48 95
    /tmp/scratch-log-check "$R/prefix-$arm.log" $((513 * 16 * 95 * 4)) 300 -1
done
for arm in cache2 cache3 cache4; do
    cmp "$R/prefix-control.ckpt" "$R/prefix-$arm.ckpt"
done
/tmp/corpus-identity "$R/prefix-control.bin" "$R/prefix-cache2.bin" \
    "$R/prefix-cache3.bin" "$R/prefix-cache4.bin" > "$R/prefix-identity.txt"
for prefix in $ARMS; do
    for arm in $ARMS; do
        label="$prefix-to-$arm"
        cp "$R/prefix-$prefix.ckpt" "$R/$label.ckpt"
        "./ecc2k130-$arm" --curve 131 --packed --threads 513 --dp-weight 48 \
            --dp-cap 262144 --steps 95 --launches 3 --verify 0 --run-id "$RUN_ID" \
            --checkpoint "$R/$label.ckpt" --dp-file "$R/$label.bin" \
            > "$R/$label.log" 2>&1
        markers "$R/$label.log" "$arm" 513 48 95
        /tmp/scratch-log-check "$R/$label.log" $((513 * 16 * 95 * 3)) 0 380
    done
done
for prefix in $ARMS; do
    for arm in $ARMS; do
        label="$prefix-to-$arm"
        [ "$label" = control-to-control ] && continue
        cmp "$R/control-to-control.ckpt" "$R/$label.ckpt"
    done
done
continuations=()
for prefix in $ARMS; do for arm in $ARMS; do continuations+=("$R/$prefix-to-$arm.bin"); done; done
/tmp/corpus-identity "${continuations[@]}" > "$R/continuation-identity.txt"

echo 'PASS native model/helper, exact sm120 resources, three device scratch controls, arithmetic/storage/shared-sigma gates, 300/300 replay, six boundary lengths, partial cross-arm checkpoints and sorted corpus identity' \
    > "$R/preflight.txt"

printf 'phase\tpair\torder\tvariant\trateMps\titerations\tdropped\tlogSha256\tgpuState\n' \
    > "$R/samples.tsv"
sample() {
    local phase=$1 pair=$2 order=$3 variant=$4 arm
    local log rc=0 count rate iterations dropped digest state
    arm=$variant
    case "$variant" in control_a|control_b) arm=control ;; esac
    log="$R/$phase-$pair-$order-$variant.log"
    "./ecc2k130-$arm" --curve 131 --packed --threads "$BENCH_THREADS" --bench \
        --steps 1024 --launches 32 --verify 0 > "$log" 2>&1 || rc=$?
    markers "$log" "$arm" "$BENCH_THREADS" 0 1024 || rc=1
    count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
    rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
    iterations=$(sed -nE 's/.*[[:space:]]([0-9]+) iterations.*/\1/p' "$log" | tail -1)
    dropped=$(sed -nE 's/.*\(([0-9]+) verified against the reference, ([0-9]+) dropped\).*/\2/p' "$log" | tail -1)
    [ "$count" = 1 ] && [ "$iterations" = "$UPDATES" ] && [ "${dropped:-1}" = 0 ] || rc=1
    grep -Eq '\(0 verified against the reference, 0 dropped\)' "$log" || rc=1
    digest=$(sha256sum "$log" | awk '{print $1}')
    state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu \
        --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$phase" "$pair" "$order" "$variant" "$rate" "$iterations" \
        "${dropped:-1}" "$digest" "$state" | tee -a "$R/samples.tsv"
    return "$rc"
}
sample warmup 0 1 control
sample warmup 0 2 cache2
sample warmup 0 3 cache3
sample warmup 0 4 cache4
for pair in 1 2 3 4 5; do
    if [ $((pair % 2)) = 1 ]; then a=control_a; b=control_b; else a=control_b; b=control_a; fi
    sample aa "$pair" 1 "$a"; sample aa "$pair" 2 "$b"
done
sample screen 1 1 control; sample screen 1 2 cache2; sample screen 1 3 cache3; sample screen 1 4 cache4
sample screen 2 1 cache4; sample screen 2 2 cache3; sample screen 2 3 cache2; sample screen 2 4 control
sample screen 3 1 cache2; sample screen 3 2 control; sample screen 3 3 cache4; sample screen 3 4 cache3
sample screen 4 1 cache3; sample screen 4 2 cache4; sample screen 4 3 control; sample screen 4 4 cache2
sample screen 5 1 control; sample screen 5 2 cache2; sample screen 5 3 cache3; sample screen 5 4 cache4

awk -F '\t' -v root="$R" 'NR>1 {print $8 "  " root "/" $1 "-" $2 "-" $3 "-" $4 ".log"}' \
    "$R/samples.tsv" > "$R/sample-files.sha256"
sha256sum -c "$R/sample-files.sha256" > "$R/sample-files.check" 2>&1
preflight_sha=$(sha256sum "$R/preflight.txt" | awk '{print $1}')
/tmp/scratch-summarize "$R" "$preflight_sha" "$R/result.json"

for file in "$R"/*; do
    case "$file" in "$R/artifact-files.sha256"|"$R/job.log"|"$R/exit-code") continue ;; esac
    [ ! -f "$file" ] || sha256sum "$file"
done > /tmp/shared-scratch-artifact-files.sha256
mv /tmp/shared-scratch-artifact-files.sha256 "$R/artifact-files.sha256"
echo '=== done'
