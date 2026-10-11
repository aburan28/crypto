#!/bin/sh
# Frozen order and binaries for PAIRED_REPLAY_PROTOCOL.md.
set -u

out_dir=/Volumes/SSD990/crypto/worktrees/f6-n83-reduction-20261004/research/f6_n83_reduction_20261004
baseline=/private/tmp/f6_n83_reduction_baseline_full
candidate=/private/tmp/f6_n83_reduction_candidate_full

: > "$out_dir/replay_status.tsv"
run_one() {
    number=$1
    arm=$2
    binary=$3
    "$binary" > "$out_dir/replay_${number}_${arm}.jsonl" 2> "$out_dir/replay_${number}_${arm}.stderr"
    status=$?
    printf '%s\t%s\t%s\n' "$number" "$arm" "$status" >> "$out_dir/replay_status.tsv"
}

run_one 1 baseline "$baseline"
run_one 2 candidate "$candidate"
run_one 3 candidate "$candidate"
run_one 4 baseline "$baseline"
run_one 5 baseline "$baseline"
run_one 6 candidate "$candidate"
