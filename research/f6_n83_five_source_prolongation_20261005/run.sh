#!/bin/sh
# Frozen k sequence. Stop after a column/RSS cap, timeout, or nonzero exit.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_five_source_prolongation_20261005"
binary="$root/target/bench_bins/f6_n83_five_source_prolongation_probe"
export RAYON_NUM_THREADS=1
printf 'mode\toffset\tk\texit\n' > "$out/status.tsv"
run_one() {
    mode=$1
    offset=$2
    k=$3
    label=$4
    if [ "$mode" = planted ]; then
        if gtimeout -k 5s 120s "$binary" planted "$k" > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
            code=0
        else
            code=$?
        fi
    else
        if gtimeout -k 5s 120s "$binary" ordinary "$offset" "$k" > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
            code=0
        else
            code=$?
        fi
    fi
    printf '%s\t%s\t%s\t%s\n' "$mode" "$offset" "$k" "$code" >> "$out/status.tsv"
    printf '%s %s k=%s exit=%s\n' "$mode" "$offset" "$k" "$code"
    return "$code"
}
for k in 1 4 8 16; do
    run_one planted none "$k" "planted_$k" || break
    run_one ordinary 0 "$k" "ordinary_0_$k" || break
    status=$(jq -r 'select(.phase == "root_reduction") | .status' "$out/ordinary_0_$k.jsonl")
    if [ "$status" != reduced ]; then
        printf 'stopping at k=%s status=%s\n' "$k" "$status"
        break
    fi
done
