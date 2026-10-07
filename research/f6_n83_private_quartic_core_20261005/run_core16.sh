#!/bin/sh
# Amendment 1: frozen k=16 planted and ordinary core reduction.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_private_quartic_core_20261005"
binary="$root/target/bench_bins/f6_n83_private_quartic_core16_probe"
export RAYON_NUM_THREADS=1
printf 'mode\tk\texit\n' > "$out/core16_status.tsv"
run_one() {
    mode=$1
    label=$2
    if [ "$mode" = planted ]; then
        if gtimeout -k 5s 300s "$binary" planted 16 > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
            code=0
        else
            code=$?
        fi
    else
        if gtimeout -k 5s 300s "$binary" ordinary 0 16 core16 > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
            code=0
        else
            code=$?
        fi
    fi
    printf '%s\t16\t%s\n' "$mode" "$code" >> "$out/core16_status.tsv"
    printf '%s k=16 exit=%s\n' "$mode" "$code"
    return "$code"
}
run_one planted core16_planted || exit 1
run_one ordinary core16_ordinary
