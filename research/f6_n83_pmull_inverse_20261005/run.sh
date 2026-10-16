#!/bin/sh
# Frozen native Rust binaries, no builds between full-base arms.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_pmull_inverse_20261005"
bins="$root/target/bench_bins"
printf 'arm\tprobe\texit\n' > "$out/status.tsv"
run_one() {
    arm=$1
    probe=$2
    binary=$3
    label=$4
    if "$binary" > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
        code=0
    else
        code=$?
    fi
    printf '%s\t%s\t%s\n' "$arm" "$probe" "$code" >> "$out/status.tsv"
    printf '%s %s exit=%s\n' "$arm" "$probe" "$code"
}
run_one baseline small "$bins/f6_index_index_candidate" baseline_small
run_one candidate small "$bins/f6_inv_index_candidate" candidate_small
run_one baseline full "$bins/f6_index_full_candidate" full_1_baseline
run_one candidate full "$bins/f6_inv_full_candidate" full_2_candidate
run_one candidate full "$bins/f6_inv_full_candidate" full_3_candidate
run_one baseline full "$bins/f6_index_full_candidate" full_4_baseline
