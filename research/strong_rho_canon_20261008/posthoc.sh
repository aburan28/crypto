#!/usr/bin/env bash
# The two M3 follow-ups that ran after the declared runs, exactly as they ran.
# Not declared in PROTOCOL.md; README.md reports them apart.
#
#   callgrind  One Callgrind run per arm, portable n61 (Valgrind hides
#              AVX-512), one price round.
#   screen     Five interleaved rounds of three arms on portable n61:
#              baseline, item1 (23ed8d3a alone, which runs the same
#              index-calculus source as the baseline) and candidate.
#
#   BIN=<dir holding ic-{baseline,item1,candidate}> BENCH=<isolated_bench> \
#     CG=<dir for the raw Callgrind files> bash posthoc.sh callgrind|screen
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
runs="$here/runs/m3"
slug=icv1-f2m61-t158598901-ab42b6c5
params="$root/research/ic_single_target_20260930/runs/k0n61/T01.params.json"
: "${BIN:?}"

case ${1:?callgrind or screen} in
callgrind)
    : "${CG:?}"
    d="$runs/callgrind"; mkdir -p "$d"
    for arm in baseline candidate; do
        [ -e "$d/$slug-portable-$arm.price.json" ] && continue
        KIC_SCAN_SIMD=0 RAYON_NUM_THREADS=1 valgrind --tool=callgrind \
            --callgrind-out-file="$CG/m3-$arm.callgrind" "$BIN/ic-$arm" price \
            --params "$params" --repeats 1 --repeats-fast 1 --price-rounds 1 --json \
            --out "$d/$slug-portable-$arm.price.json" > /dev/null 2> "$d/$slug-portable-$arm.valgrind.log"
        callgrind_annotate --threshold=99 "$CG/m3-$arm.callgrind" \
            | sed "s#$CG/##g; s#$BIN/##g; s#$root/##g" | head -45 > "$d/$slug-portable-$arm.annotate.txt"
    done ;;
screen)
    : "${BENCH:?}"
    d="$runs/posthoc"; mkdir -p "$d"
    for round in 1 2 3 4 5; do
        for arm in baseline item1 candidate; do
            out="$d/$slug-portable-$arm-r$round.price.json"
            [ -e "$out" ] && continue
            "$BENCH" run --wait --cpus 3 --max-other-cpu "${MAX_OTHER_CPU:-0.3}" \
                --label "m3/posthoc/$slug-portable-$arm-r$round" --out "$d/isolation.jsonl" -- \
                env KIC_SCAN_SIMD=0 RAYON_NUM_THREADS=1 "$BIN/ic-$arm" price --params "$params" \
                --repeats 1 --repeats-fast 1 --json --out "$out" > /dev/null
        done
    done ;;
*) echo "usage: posthoc.sh callgrind|screen" >&2; exit 2 ;;
esac
