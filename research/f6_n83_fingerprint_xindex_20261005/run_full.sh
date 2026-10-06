#!/bin/sh
# Only run after the preregistered dimension-8/10 screen passes.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_fingerprint_xindex_20261005"
bins="$root/target/bench_bins"
export RAYON_NUM_THREADS=1
printf 'order\tarm\texit\n' > "$out/full_status.tsv"
run_one() {
    order=$1
    arm=$2
    binary=$3
    label=$4
    if gtimeout -k 5s 120s "$binary" > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
        code=0
    else
        code=$?
    fi
    printf '%s\t%s\t%s\n' "$order" "$arm" "$code" >> "$out/full_status.tsv"
    printf '%s %s exit=%s\n' "$order" "$arm" "$code"
}
run_one 1 baseline "$bins/f6_xmap_full_candidate" full_1_baseline
run_one 2 candidate "$bins/f6_fingerprint_full_candidate" full_2_candidate
run_one 3 candidate "$bins/f6_fingerprint_full_candidate" full_3_candidate
run_one 4 baseline "$bins/f6_xmap_full_candidate" full_4_baseline
run_one 5 baseline "$bins/f6_xmap_full_candidate" full_5_baseline
run_one 6 candidate "$bins/f6_fingerprint_full_candidate" full_6_candidate
run_one 7 candidate "$bins/f6_fingerprint_full_candidate" full_7_candidate
run_one 8 baseline "$bins/f6_xmap_full_candidate" full_8_baseline
run_one 9 baseline "$bins/f6_xmap_full_candidate" full_9_baseline
run_one 10 candidate "$bins/f6_fingerprint_full_candidate" full_10_candidate
