#!/bin/sh
# Same n83 public T001 and exact full-base query, baseline/final/final/baseline.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_compact_xmap_20261005"
bins="$root/target/bench_bins"
printf 'arm\texit\n' > "$out/direct_status.tsv"
run_one() {
    arm=$1
    binary=$2
    label=$3
    if gtimeout -k 5s 120s "$binary" > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
        code=0
    else
        code=$?
    fi
    printf '%s\t%s\n' "$arm" "$code" >> "$out/direct_status.tsv"
    printf '%s exit=%s\n' "$arm" "$code"
}
run_one baseline "$bins/f6_pmull_full_candidate" direct_1_baseline
run_one final "$bins/f6_xmap_full_candidate" direct_2_final
run_one final "$bins/f6_xmap_full_candidate" direct_3_final
run_one baseline "$bins/f6_pmull_full_candidate" direct_4_baseline
