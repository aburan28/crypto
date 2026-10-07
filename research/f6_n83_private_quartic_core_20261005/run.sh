#!/bin/sh
# Frozen k=16,90 panel. Stop on a failed or incomplete certificate.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_private_quartic_core_20261005"
binary="$root/target/bench_bins/f6_n83_private_quartic_core_probe"
export RAYON_NUM_THREADS=1
printf 'mode\tk\texit\n' > "$out/status.tsv"
run_one() {
    mode=$1
    k=$2
    label=$3
    if [ "$mode" = planted ]; then
        if gtimeout -k 5s 300s "$binary" planted "$k" > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
            code=0
        else
            code=$?
        fi
    else
        if gtimeout -k 5s 300s "$binary" ordinary 0 "$k" > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
            code=0
        else
            code=$?
        fi
    fi
    printf '%s\t%s\t%s\n' "$mode" "$k" "$code" >> "$out/status.tsv"
    printf '%s k=%s exit=%s\n' "$mode" "$k" "$code"
    return "$code"
}
for k in 16 90; do
    run_one planted "$k" "planted_$k" || break
    run_one ordinary "$k" "ordinary_0_$k" || break
    certificate_status=$(jq -r 'select(.phase == "private_certificate") | .status' "$out/ordinary_0_$k.jsonl")
    if [ "$certificate_status" != complete ]; then
        printf 'stopping at k=%s certificate=%s\n' "$k" "$certificate_status"
        break
    fi
done
