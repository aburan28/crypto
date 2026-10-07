#!/bin/sh
# One bounded native process per planted control and ordinary torsion offset.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_fivesum_d18_20261005"
binary=${F6_FIVESUM_PROBE:-"$root/target/bench_bins/f6_n83_fivesum_d18_probe"}
printf 'label\texit\n' > "$out/status.tsv"
run_one() {
    label=$1
    shift
    if gtimeout -k 5s 120s "$binary" "$@" > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
        code=0
    else
        code=$?
    fi
    printf '%s\t%s\n' "$label" "$code" >> "$out/status.tsv"
    printf '%s exit=%s\n' "$label" "$code"
}
run_one planted planted
run_one ordinary_0 ordinary 0
run_one ordinary_1 ordinary 1
run_one ordinary_2 ordinary 2
run_one ordinary_3 ordinary 3
