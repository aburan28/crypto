#!/bin/sh
# Run only the frozen native Rust probe binaries; do not rebuild between arms.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_batch_size_20261004"
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
run_one baseline small "$bins/f6_wide_n83_index_probe" baseline_small
run_one candidate small "$bins/f6_wide_n83_index_probe_candidate" candidate_small
run_one baseline full "$bins/f6_wide_n83_full_index_probe" full_1_baseline
run_one candidate full "$bins/f6_wide_n83_full_index_probe_candidate" full_2_candidate
run_one candidate full "$bins/f6_wide_n83_full_index_probe_candidate" full_3_candidate
run_one baseline full "$bins/f6_wide_n83_full_index_probe" full_4_baseline
