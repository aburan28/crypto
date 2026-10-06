#!/bin/sh
# Frozen small probes: baseline, candidate, candidate, baseline.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_compact_xmap_20261005"
bins="$root/target/bench_bins"
printf 'arm\texit\n' > "$out/small_replay_status.tsv"
run_one() {
    arm=$1
    binary=$2
    label=$3
    if gtimeout -k 5s 120s "$binary" > "$out/$label.jsonl" 2> "$out/$label.stderr.txt"; then
        code=0
    else
        code=$?
    fi
    printf '%s\t%s\n' "$arm" "$code" >> "$out/small_replay_status.tsv"
    printf '%s exit=%s\n' "$arm" "$code"
}
run_one baseline "$bins/f6_compact_index_candidate" small_replay_1_baseline
run_one candidate "$bins/f6_xmap_index_candidate" small_replay_2_candidate
run_one candidate "$bins/f6_xmap_index_candidate" small_replay_3_candidate
run_one baseline "$bins/f6_compact_index_candidate" small_replay_4_baseline
