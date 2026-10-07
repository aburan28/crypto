#!/bin/sh
# Deterministic interleaving of the already-published baseline and selection.
set -eu

if [ "$#" -ne 2 ]; then
    echo "usage: $0 <n:41|53> <release-example-binary-directory>" >&2
    exit 2
fi
n=$1
bin_dir=$2
study=experiments/koblitz-base-size-cold-panel-20261004
selection=$study/PILOT_ANALYSIS.json
baseline=$(jq -er --argjson n "$n" '.cells[] | select(.n==$n) | .baseline_K' "$selection")
selected=$(jq -er --argjson n "$n" '.cells[] | select(.n==$n) | .selected_smaller_K' "$selection")
pair=$study/run_holdout_pair.sh

for repeat in 1 2 3 4 5 6; do
    if [ $((repeat % 2)) -eq 1 ]; then
        "$pair" "$n" "$baseline" "$repeat" "$bin_dir"
        "$pair" "$n" "$selected" "$repeat" "$bin_dir"
    else
        "$pair" "$n" "$selected" "$repeat" "$bin_dir"
        "$pair" "$n" "$baseline" "$repeat" "$bin_dir"
    fi
done
