#!/usr/bin/env bash
# The evidence runs declared in PROTOCOL.md: A/A then interleaved A/B on the
# six frozen §23 parameter files, both key paths, each process under
# isolated_bench on one reserved CPU.  Thin orchestration only: every figure
# comes from the `ic price` reports this writes, and nothing is overwritten.
#
#   BASE=<ic built at 7424539b> CAND=<ic built at a98ba95a> \
#   BENCH=<isolated_bench> [MAX_OTHER_CPU=0.3] bash run.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
runs="$here/runs"
sizes="k0n41 k1n47 k0n53 k0n57 k1n59 k0n61"
: "${BASE:?}" "${CAND:?}" "${BENCH:?}"

price() { # phase size path arm round
    local phase=$1 size=$2 path=$3 arm=$4 round=$5 bin simd
    bin=$BASE; [ "$arm" = candidate ] && bin=$CAND
    simd=1; [ "$path" = portable ] && simd=0
    local out="$runs/$phase/$size-$path-$arm-r$round.price.json"
    [ -e "$out" ] && return 0
    mkdir -p "$runs/$phase"
    "$BENCH" run --wait --cpus 3 --max-other-cpu "${MAX_OTHER_CPU:-0.1}" --label "$phase/$size-$path-$arm-r$round" \
        --out "$runs/$phase/isolation.jsonl" -- \
        env KIC_SCAN_SIMD=$simd RAYON_NUM_THREADS=1 "$bin" price \
        --params "$root/research/ic_single_target_20260930/runs/$size/T01.params.json" \
        --repeats 1 --repeats-fast 1 --json --out "$out" > /dev/null
}

for round in 1 2 3; do
    for size in $sizes; do
        for path in portable avx512; do
            price aa "$size" "$path" baseline "$round"
            price aa "$size" "$path" baseline "$((round + 3))"
        done
    done
done
for round in 1 2 3 4 5; do
    for size in $sizes; do
        for path in portable avx512; do
            price ab "$size" "$path" baseline "$round"
            price ab "$size" "$path" candidate "$round"
        done
    done
done
echo "runs complete"
