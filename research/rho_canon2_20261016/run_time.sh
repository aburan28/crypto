#!/usr/bin/env bash
# Native time, as PROTOCOL.md declares: koblitz_rho_fixture n 0
# signed_frobenius 4 strong for n = 53 and 61, an A/A pass (the baseline six
# times, in 3 rounds of two) and then 5 interleaved A/B rounds, every process
# under isolated_bench on one reserved CPU with RAYON_NUM_THREADS=1.  Nothing
# is overwritten.
#
#   BIN=<dir with koblitz_rho_fixture-{baseline,candidate}> BENCH=<isolated_bench> \
#     [MAX_OTHER_CPU=0.3] bash run_time.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
runs="$here/runs/time"
: "${BIN:?}" "${BENCH:?}"
curves="icv1-f2m53-tm56619371-dac20a85=53 icv1-f2m61-t158598901-ab42b6c5=61"

one() { # phase slug=n arm round
    local slug=${2%%=*} n=${2#*=} d="$runs/$1"
    local out="$d/$slug-$3-r$4.jsonl"
    [ -e "$out" ] && return 0
    mkdir -p "$d"
    "$BENCH" run --wait --cpus 3 --max-other-cpu "${MAX_OTHER_CPU:-0.1}" --label "time/$1/$slug-$3-r$4" \
        --out "$d/isolation.jsonl" -- env RAYON_NUM_THREADS=1 \
        "$BIN/koblitz_rho_fixture-$3" "$n" 0 signed_frobenius 4 strong > "$out.tmp"
    mv "$out.tmp" "$out"
}

for round in 1 2 3; do
    for c in $curves; do one aa "$c" baseline "$round"; one aa "$c" baseline "$((round + 3))"; done
done
for round in 1 2 3 4 5; do
    for c in $curves; do one ab "$c" baseline "$round"; one ab "$c" candidate "$round"; done
done
echo "runs complete"
