#!/bin/sh
# Frozen binaries. Never compile or modify a binary between paired arms.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_flat_xindex_20261005"
bins="$root/target/bench_bins"
printf 'order\tarm\tprobe\texit\n' > "$out/status.tsv"
run_one() {
    order=$1
    arm=$2
    probe=$3
    binary=$4
    label=$5
    if gtimeout -k 5s 120s "$binary" > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
        code=0
    else
        code=$?
    fi
    printf '%s\t%s\t%s\t%s\n' "$order" "$arm" "$probe" "$code" >> "$out/status.tsv"
    printf '%s %s %s exit=%s\n' "$order" "$arm" "$probe" "$code"
}
export RAYON_NUM_THREADS=1
run_one 1 baseline small "$bins/f6_xmap_index_candidate" small_1_baseline
run_one 2 candidate small "$bins/f6_flat_index_candidate" small_2_candidate
run_one 3 candidate small "$bins/f6_flat_index_candidate" small_3_candidate
run_one 4 baseline small "$bins/f6_xmap_index_candidate" small_4_baseline
run_one 1 baseline full "$bins/f6_xmap_full_candidate" full_1_baseline
run_one 2 candidate full "$bins/f6_flat_full_candidate" full_2_candidate
run_one 3 candidate full "$bins/f6_flat_full_candidate" full_3_candidate
run_one 4 baseline full "$bins/f6_xmap_full_candidate" full_4_baseline
run_one 5 baseline full "$bins/f6_xmap_full_candidate" full_5_baseline
run_one 6 candidate full "$bins/f6_flat_full_candidate" full_6_candidate
run_one 7 candidate full "$bins/f6_flat_full_candidate" full_7_candidate
run_one 8 baseline full "$bins/f6_xmap_full_candidate" full_8_baseline
run_one 9 baseline full "$bins/f6_xmap_full_candidate" full_9_baseline
run_one 10 candidate full "$bins/f6_flat_full_candidate" full_10_candidate
