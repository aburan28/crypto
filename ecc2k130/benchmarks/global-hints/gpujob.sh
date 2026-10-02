#!/bin/bash
# Frozen full-cost v3 control versus GPU-wide queue; no search/solver job.
set -euo pipefail
cd "${WORK:-/work}"
R=${RESULTS:-/results}
mkdir -p "$R"
mapfile -t GPU_NAMES < <(nvidia-smi --query-gpu=name --format=csv,noheader)
[ "${#GPU_NAMES[@]}" = 1 ]
GPU_NAME=${GPU_NAMES[0]}
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | tr -d '. ')
[ "$GPU_NAME" = 'NVIDIA RTX PRO 6000 Blackwell Server Edition' ] && [ "$CAP" = 120 ]
[[ ${SOURCE_REV:-} =~ ^[0-9a-f]{40}$ ]]
if [ -f SOURCE_REV ]; then [ "$(cat SOURCE_REV)" = "$SOURCE_REV" ]; fi
if [ -e .git ]; then
    [ "$(git rev-parse HEAD)" = "$SOURCE_REV" ]
    [ -z "$(git status --porcelain -- .)" ]
fi
nvcc --version | grep -q 'release 13\.3, V13\.3\.73'
ARCH='PRO6000_ARCH=-gencode arch=compute_120,code=sm_120'
COMMON='BATCH=16 THREADS=512 MINBLOCKS=1 WITNESS=0 GLOBAL_CG=0 TABLE_FUSED=0 TABLE_GLOBAL=0 TABLE_ADDEND_GLOBAL=0 TABLE_PHASE_POPC=0 CYCLE_PROFILE=0 PHASE_PROFILE=0 PACKED_CHAINS=1 PACKED_TOP_HOIST=0 PACKED_ONB_INV=0 PACKED_SLOT_PREFETCH=0 PACKED_SLOT_PIPELINE=0 PACKED_CLMAD_SQUARE=0 PACKED_KARAT3=0 PACKED_TOP_CLMAD=0 PACKED_ADD_COMBINE=0 PACKED_ALU_SQR=0 SIGMA_FUSED=0 SIGMA_FUSED_LATE_Y=0'
VERIFY_THREADS=96256
BENCH_THREADS=385024
{
    nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
    nvcc --version
    echo "source: $SOURCE_REV"
    echo "control: $COMMON TABLE_BLOCK_HINTS=1 TABLE_GLOBAL_HINTS=0"
    echo "candidate: $COMMON TABLE_BLOCK_HINTS=0 TABLE_GLOBAL_HINTS=1"
} | tee "$R/host.txt"
sha256sum Makefile include/packed*.h include/packed*.cuh include/tablewalk.h include/cycleanchor_body.h \
    include/ref.h include/tablev3replay.h src/main.cu src/tablev3replay.cpp \
    benchmarks/global-hints/* benchmarks/block-both2-confirm5/corpus_identity.cpp > "$R/source-files.sha256"
g++ -O2 -std=c++17 -Wall -Wextra -Werror benchmarks/global-hints/test_queue.cpp -o /tmp/test-queue
g++ -O2 -std=c++17 -Wall -Wextra -Werror benchmarks/global-hints/summarize.cpp -o /tmp/global-summarize
g++ -O2 -std=c++17 -Wall -Wextra -Werror benchmarks/block-both2-confirm5/corpus_identity.cpp -o /tmp/corpus-identity
{
    /tmp/test-queue
    /tmp/global-summarize --self-test
    /tmp/corpus-identity --self-test
    make CXX=g++ test-packed test-table-v3-replay
} > "$R/native-controls.log" 2>&1
build() {
    local name=$1 extra=$2
    make -s gpu-preset "$ARCH" KNOBS="$COMMON $extra" > "$R/build-$name.log" 2>&1
    mv ecc2k130 "ecc2k130-$name"
    cp "ecc2k130-$name" "$R/"
}
build control 'TABLE_BLOCK_HINTS=1 TABLE_GLOBAL_HINTS=0'
build candidate 'TABLE_BLOCK_HINTS=0 TABLE_GLOBAL_HINTS=1'
make -s test-global-hints-preset "$ARCH" KNOBS="$COMMON TABLE_BLOCK_HINTS=0 TABLE_GLOBAL_HINTS=1" \
    > "$R/device-controls.log" 2>&1
cp build/test-global-hints-cuda /tmp/test-queue /tmp/global-summarize /tmp/corpus-identity build/table-v3-replay "$R/"
sha256sum "$R"/ecc2k130-* "$R"/test-global-hints-cuda "$R"/test-queue "$R"/global-summarize \
    "$R"/corpus-identity "$R"/table-v3-replay > "$R/binary-sha256.txt"
markers() {
    local log=$1 name=$2 workers=$3 global=0 block=1 dp=${4:-48} steps=${5:-95}
    local blob window
    if [ "$name" = candidate ]; then global=1; block=0; fi
    grep -qx "packed table GPU-wide hints: $global" "$log" || return 1
    grep -qx "packed table block hints: $block, queue 512" "$log" || return 1
    grep -qx 'packed table split forward: 1' "$log" || return 1
    grep -qx 'packed table batch hints: 1' "$log" || return 1
    grep -qx 'packed cycle fast2: 1' "$log" || return 1
    grep -qx 'packed witness: 0' "$log" || return 1
    grep -qx 'packed square table: 1' "$log" || return 1
    grep -qx 'packed polynomial inversion: 2' "$log" || return 1
    grep -qx 'packed table walk: 1 (8 branches, 57052 shared bytes)' "$log" || return 1
    grep -qx 'packed launch bounds: 512 threads, 1 min blocks' "$log" || return 1
    grep -qx "backend cuda-packed131: $workers threads x 16 slots x 1 lanes = $((workers*16)) walks, dp weight $dp, $steps steps per launch" "$log" || return 1
    grep -qx 'packed hint scheduling: table fused 0, pipe select 1, chain first 1, inline polynomial 3, phase profile 0, cycle profile 0' "$log" || return 1
    local feature
    for feature in 'denominator cache: 1' 'multiply by value: 1' 'Frobenius network: 3' \
        'polynomial chain: 1' 'polynomial state: 1' 'unrolled inversion: 1' 'paired products: 1' \
        'pair ilp: 1' 'pair clmul: 0' 'clmul flat: 0' 'from reduced: 1' 'slot unroll: 1' 'chains: 1' \
        'L2 persist: 1' 'direct reduction: 1' 'generated product: 1' 'native carryless multiply: 1' \
        'weighted prefix: 2' 'compact state: 1' 'shared sigma: 1' 'state tile: 256' 'alu square: 1' \
        'table global: 0' 'table addend global: 0' 'top hoist: 0' 'onb inv: 0' 'slot prefetch: 0' \
        'slot pipeline: 0' 'native carryless square: 0' 'three-limb Karatsuba: 0' 'top clmad: 0' \
        'add combine: 0' 'alu onb square: 0' 'table phase popc: 0' 'sigma fused: 0' \
        'sigma fused late y: 0' 'profile ranges: 0'; do
        grep -qx "packed $feature" "$log" || return 1
    done
    grep -qx 'packed table pivot bytes: 1, table shared bytes 57052' "$log" || return 1
    blob=$((((workers+255)/256)*16*4352*3))
    window=$blob; if [ "$window" -gt 83886080 ]; then window=83886080; fi
    grep -qx "packed L2 persist window: $window of $blob field bytes, cap 83886080" "$log" || return 1
    grep -Eq '^packed kernel: [0-9]+ registers/thread, [0-9]+ local bytes/thread, [0-9]+ shared bytes/block, single-product multiplier$' "$log" || return 1
    if [ "$global" = 1 ]; then
        grep -qx "packed GPU-wide hint queue: $((workers*16)) entries, 188 resolver blocks of 128 threads" "$log" || return 1
        grep -Eq '^packed hint select kernel: [0-9]+ registers/thread,' "$log" || return 1
        grep -Eq '^packed hint resolve kernel: [0-9]+ registers/thread,' "$log" || return 1
    fi
    if [ -f "$R/resource-$name.txt" ]; then
        diff "$R/resource-$name.txt" <(grep -E '^packed (kernel:|hint (select|resolve) kernel:)' "$log") || return 1
    fi
    ! grep -Eq 'MISMATCH|OVERFLOW|unusable|collision:' "$log" || return 1
}
identity() {
    local a=$1 b=$2 label=$3
    /tmp/corpus-identity "$a" "$b" /tmp/sorted-a /tmp/sorted-b > "$R/identity-$label.txt"
    [ -s /tmp/sorted-a ] && [ -s /tmp/sorted-b ]
    cmp /tmp/sorted-a /tmp/sorted-b
    sha256sum /tmp/sorted-a /tmp/sorted-b >> "$R/identity-$label.txt"
}
verify() {
    local name=$1 workers=$2 label=$3 log="$R/verify-$3-$1.log"
    "./ecc2k130-$name" --curve 131 --packed --threads "$workers" --dp-weight 48 --dp-cap 262144 \
        --steps 95 --launches 7 --verify 300 --run-id 7 --dp-file "$R/dp-$label-$name.bin" > "$log" 2>&1
    markers "$log" "$name" "$workers"
    grep -Eq '\(300 verified against the reference, 0 dropped\)' "$log" || return 1
}
verify control "$VERIFY_THREADS" full
verify candidate "$VERIFY_THREADS" full
for name in control candidate; do
    grep -E '^packed (kernel:|hint (select|resolve) kernel:)' "$R/verify-full-$name.log" > "$R/resource-$name.txt"
done
identity "$R/dp-full-control.bin" "$R/dp-full-candidate.bin" full
build/table-v3-replay --run-id 7 --dp-weight 48 --samples 300 --min-nonzero 299 \
    "$R/dp-full-control.bin" "$R/dp-full-candidate.bin" > "$R/spread-replay.json"
for workers in 511 513; do
    verify control "$workers" "partial-$workers"
    verify candidate "$workers" "partial-$workers"
    identity "$R/dp-partial-$workers-control.bin" "$R/dp-partial-$workers-candidate.bin" "partial-$workers"
done
# Prefix from either arm, continuation by either arm: identical persistent
# state after odd launch boundaries and host reseeding, across both directions.
for name in control candidate; do
    "./ecc2k130-$name" --curve 131 --packed --threads 513 --dp-weight 48 --dp-cap 262144 \
        --steps 95 --launches 4 --verify 300 --run-id 7 --checkpoint "$R/prefix-$name.ckpt" \
        --dp-file "$R/prefix-$name.bin" > "$R/prefix-$name.log" 2>&1
    markers "$R/prefix-$name.log" "$name" 513
done
cmp "$R/prefix-control.ckpt" "$R/prefix-candidate.ckpt"
identity "$R/prefix-control.bin" "$R/prefix-candidate.bin" prefix
for prefix in control candidate; do
    for name in control candidate; do
        label="$prefix-to-$name"
        cp "$R/prefix-$prefix.ckpt" "$R/$label.ckpt"
        "./ecc2k130-$name" --curve 131 --packed --threads 513 --dp-weight 48 --dp-cap 262144 \
            --steps 95 --launches 3 --verify 0 --run-id 7 --checkpoint "$R/$label.ckpt" \
            --dp-file "$R/$label.bin" > "$R/$label.log" 2>&1
        markers "$R/$label.log" "$name" 513
        grep -q ' at iteration 380$' "$R/$label.log"
    done
done
for label in control-to-candidate candidate-to-control candidate-to-candidate; do
    cmp "$R/control-to-control.ckpt" "$R/$label.ckpt"
    identity "$R/control-to-control.bin" "$R/$label.bin" "$label"
done
echo 'PASS native/device queue controls, odd300 replay, spread299 nonzero, full/partial corpus identity and bidirectional checkpoints' > "$R/preflight.txt"
printf 'phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
    local phase=$1 pair=$2 order=$3 variant=$4 name=control
    local log rc=0 rate count final digest state
    [ "$variant" != candidate ] || name=candidate
    log="$R/$phase-$pair-$order-$variant.log"
    "./ecc2k130-$name" --curve 131 --packed --threads "$BENCH_THREADS" --bench \
        --steps 1024 --launches 32 --verify 0 > "$log" 2>&1 || rc=$?
    markers "$log" "$name" "$BENCH_THREADS" 0 1024 || rc=1
    count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
    rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
    final=$(sed -nE 's/.*[[:space:]]([0-9]+) iterations.*/\1/p' "$log" | tail -1)
    [ "$count" = 1 ] && [ "$final" = 201863462912 ] || rc=1
    grep -Eq '\(0 verified against the reference, 0 dropped\)' "$log" || rc=1 || return 1
    digest=$(sha256sum "$log" | awk '{print $1}')
    state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$phase" "$pair" "$order" "$variant" "$rate" "$digest" "$state" | tee -a "$R/samples.tsv"
    return "$rc"
}
sample warmup 0 1 control
sample warmup 0 2 candidate
for phase in aa ab; do
    for pair in 1 2 3 4 5; do
        if [ "$phase" = aa ]; then a=aa1; b=aa2; else a=control; b=candidate; fi
        if [ $((pair%2)) = 0 ]; then t=$a; a=$b; b=$t; fi
        sample "$phase" "$pair" 1 "$a"
        sample "$phase" "$pair" 2 "$b"
    done
done
preflight_sha=$(sha256sum "$R/preflight.txt" | awk '{print $1}')
/tmp/global-summarize "$R/samples.tsv" "$R/preflight.txt" "$preflight_sha" "$R/result.json"
# job.log is still being written by Modal; exit-code is added after return.
# The immutable raw archive and an independent post-run manifest bind those.
for file in "$R"/*; do
    case "$file" in "$R/artifact-files.sha256"|"$R/job.log"|"$R/exit-code") continue ;; esac
    [ ! -f "$file" ] || sha256sum "$file"
done > /tmp/global-artifact-files.sha256
mv /tmp/global-artifact-files.sha256 "$R/artifact-files.sha256"
echo '=== done'
