#!/bin/sh
# Thin native-binary orchestration for the preregistered one-target panel.
# Run from the repository root after building the three release examples.
set -eu

note=research/notes/ecc2k130/n37_native_one_target_20261003
points=research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b02.points.jsonl
first=${1:-0}
last=${2:-31}
rep_first=${3:-0}
rep_last=${4:-4}
ic=target/release/examples/n37_native_m6_one_target
rho=target/release/examples/koblitz_rho_batch_ks_strong_online
replay=target/release/examples/n37_native_m6_one_target_replay
mkdir -p "$note/raw"

run_ic() {
    prefix=$1
    index=$2
    if [ -f "$prefix.ic.json" ]; then
        return
    fi
    if [ -e "$prefix.ic.stdout" ] || [ -e "$prefix.ic.stderr" ]; then
        echo "partial IC artifact at $prefix; preserve and inspect it" >&2
        exit 1
    fi
    /usr/bin/time -p env RAYON_NUM_THREADS=1 "$ic" "$index" "$prefix.ic.json" \
        >"$prefix.ic.stdout" 2>"$prefix.ic.stderr"
}

run_rho() {
    prefix=$1
    seed=$2
    public_point=$3
    if [ -f "$prefix.rho.jsonl" ]; then
        return
    fi
    if [ -e "$prefix.rho.stderr" ]; then
        echo "partial rho artifact at $prefix; preserve and inspect it" >&2
        exit 1
    fi
    /usr/bin/time -p env RAYON_NUM_THREADS=1 KIC_RHO_RUNG=3 \
        KIC_RHO_LANES=32 KIC_RHO_DP_BITS=8 \
        KIC_RHO_BATCH_CORPUS=compact-disjoint-cold-v2-n37-L1024-b02-20261001 \
        KIC_RHO_TARGET_POINT="$public_point" \
        "$rho" 37 0 signed_frobenius 1 "$seed" \
        >"$prefix.rho.jsonl" 2>"$prefix.rho.stderr"
}

i=$first
while [ "$i" -le "$last" ]; do
    line=$((i + 1))
    public_point=$(sed -n "${line}p" "$points" | jq -r 'join(",")')
    if [ -z "$public_point" ]; then
        echo "missing public point $i" >&2
        exit 1
    fi
    j=$rep_first
    while [ "$j" -le "$rep_last" ]; do
        name=$(printf 'q%02d_r%d' "$i" "$j")
        prefix="$note/raw/$name"
        seed=$((2026100110102 + 1000 * i + j))
        if [ $(((i + j) % 2)) -eq 0 ]; then
            run_ic "$prefix" "$i"
            run_rho "$prefix" "$seed" "$public_point"
        else
            run_rho "$prefix" "$seed" "$public_point"
            run_ic "$prefix" "$i"
        fi
        if [ ! -f "$prefix.replay.json" ]; then
            "$replay" "$prefix.ic.json" "$prefix.rho.jsonl" \
                "$prefix.replay.json" >"$prefix.replay.stdout" 2>"$prefix.replay.stderr"
        fi
        echo "$name complete"
        j=$((j + 1))
    done
    i=$((i + 1))
done
