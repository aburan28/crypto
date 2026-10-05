#!/bin/sh
# Run from repository root against the two frozen planted-witness binaries.
set -u
out=research/f6_n83_valid_coordinates_20261005
baseline=target/bench_bins/f6_wide_witness_lazy
candidate=target/bench_bins/f6_wide_witness_indexed
: > "$out/witness_status.tsv"
failed=0
for pair in 1 2 3 4 5; do
    if [ $((pair % 2)) -eq 1 ]; then
        order='baseline candidate'
    else
        order='candidate baseline'
    fi
    for arm in $order; do
        if [ "$arm" = baseline ]; then
            bin=$baseline
        else
            bin=$candidate
        fi
        env RAYON_NUM_THREADS=1 gtimeout -k 5s 120s "$bin" \
            > "$out/witness_${pair}_${arm}.jsonl" \
            2> "$out/witness_${pair}_${arm}.stderr.txt"
        code=$?
        printf '%s\t%s\t%s\n' "$pair" "$arm" "$code" >> "$out/witness_status.tsv"
        if [ "$code" -ne 0 ]; then
            failed=1
        fi
    done
done
exit "$failed"
