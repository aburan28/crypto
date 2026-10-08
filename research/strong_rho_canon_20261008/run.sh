#!/usr/bin/env bash
# The evidence runs PROTOCOL.md declares: M1 (m = 83 reference rate),
# M2 (narrow strong rho fixtures) and M3 (index-calculus pipeline), each an
# A/A pass then interleaved A/B rounds, every process under isolated_bench
# on one reserved CPU.  Thin orchestration: every figure is in the files
# this writes, and nothing is overwritten.
#
#   BIN=<dir holding {strong_rho_step_rate,koblitz_rho_fixture,ic}-{baseline,candidate}>
#   BENCH=<isolated_bench> [MAX_OTHER_CPU=0.3] bash run.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
runs="$here/runs"
: "${BIN:?}" "${BENCH:?}"
# The six §23 curves: slug = n:a for the fixture, and the §23 run directory.
narrow="icv1-f2m41-tm2308219-7f48b14a=41:0:k0n41 icv1-f2m47-t22705043-f4e44623=47:1:k1n47
        icv1-f2m53-tm56619371-dac20a85=53:0:k0n53 icv1-f2m57-tm747311035-c1f545af=57:0:k0n57
        icv1-f2m59-tm943548413-98844ecc=59:1:k1n59 icv1-f2m61-t158598901-ab42b6c5=61:0:k0n61"
wide=icv1-f2m83-tm6151469093347-debefd74

bench() { # dir label out-file command...
    local dir=$1 label=$2 out=$3; shift 3
    [ -e "$out" ] && return 0
    mkdir -p "$dir"
    "$BENCH" run --wait --cpus 3 --max-other-cpu "${MAX_OTHER_CPU:-0.1}" --label "$label" \
        --out "$dir/isolation.jsonl" -- env RAYON_NUM_THREADS=1 "$@" > "$out.tmp"
    mv "$out.tmp" "$out"
}

m1() { # phase arm round
    local d="$runs/m1/$1"
    bench "$d" "m1/$1/$wide-$2-r$3" "$d/$wide-$2-r$3.json" \
        "$BIN/strong_rho_step_rate-$2" 262144
}

m2() { # phase slug=n:a:dir arm round
    local slug=${2%%=*} spec=${2#*=} d="$runs/m2/$1"
    local n=${spec%%:*} rest=${spec#*:}
    local a=${rest%%:*}
    bench "$d" "m2/$1/$slug-$3-r$4" "$d/$slug-$3-r$4.jsonl" \
        "$BIN/koblitz_rho_fixture-$3" "$n" "$a" signed_frobenius 4 strong
}

m3() { # phase slug=n:a:dir path arm round
    local slug=${2%%=*} spec=${2#*=} d="$runs/m3/$1" simd=1
    local dir=${spec##*:}
    [ "$3" = portable ] && simd=0
    local out="$d/$slug-$3-$4-r$5.price.json"
    [ -e "$out" ] && return 0
    mkdir -p "$d"
    "$BENCH" run --wait --cpus 3 --max-other-cpu "${MAX_OTHER_CPU:-0.1}" \
        --label "m3/$1/$slug-$3-$4-r$5" --out "$d/isolation.jsonl" -- \
        env KIC_SCAN_SIMD=$simd RAYON_NUM_THREADS=1 "$BIN/ic-$4" price \
        --params "$root/research/ic_single_target_20260930/runs/$dir/T01.params.json" \
        --repeats 1 --repeats-fast 1 --json --out "$out" > /dev/null
}

for round in 1 2 3; do
    m1 aa baseline "$round"; m1 aa baseline "$((round + 3))"
    for s in $narrow; do m2 aa "$s" baseline "$round"; m2 aa "$s" baseline "$((round + 3))"; done
    for s in $narrow; do
        case $s in *k0n53|*k0n61) ;; *) continue ;; esac
        for path in portable avx512; do
            m3 aa "$s" "$path" baseline "$round"; m3 aa "$s" "$path" baseline "$((round + 3))"
        done
    done
done
for round in 1 2 3 4 5; do
    m1 ab baseline "$round"; m1 ab candidate "$round"
    for s in $narrow; do m2 ab "$s" baseline "$round"; m2 ab "$s" candidate "$round"; done
    for s in $narrow; do
        case $s in *k0n53|*k0n61) ;; *) continue ;; esac
        for path in portable avx512; do
            m3 ab "$s" "$path" baseline "$round"; m3 ab "$s" "$path" candidate "$round"
        done
    done
done
echo "runs complete"
