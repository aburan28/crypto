#!/bin/sh
# End-to-end CLI regression using one fixed public point and two base sizes.
set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 <release-example-binary-directory>" >&2
    exit 2
fi
bin_dir=$1
study=experiments/koblitz-failed-target-exit-20261004
target=$study/input/n41_public_q.jsonl
temporary=$(mktemp -d)
trap 'rm -rf "$temporary"' EXIT HUP INT TERM
timeout_bin=${TIMEOUT_BIN:-timeout}

run_case() {
    columns=$1
    output=$temporary/k$columns
    mkdir -p "$output"
    if KIC_DUMP_BASE="$output/base.jsonl" KIC_DUMP_RANK="$output/rank.jsonl" \
        "$timeout_bin" 30 "$bin_dir/koblitz_orbit_dlp_fast_online" \
        "construct:41:0:$columns" "$target" 410041 "$output/target.jsonl" \
        > "$output/summary.jsonl" 2> "$output/stderr"; then
        status=0
    else
        status=$?
    fi
    [ -s "$output/base.jsonl" ]
    [ -s "$output/rank.jsonl" ]
    [ -s "$output/summary.jsonl" ]
    [ -s "$output/target.jsonl" ]
    if [ "$columns" -eq 20 ]; then
        [ "$status" -eq 1 ]
        jq -e '.rank==20 and .targets==1 and .targets_solved==0 and .targets_failed==1' "$output/summary.jsonl" > /dev/null
        jq -e '.exit_code==1 and .group_verified==null and .recovered_scalar==null' "$output/target.jsonl" > /dev/null
    else
        [ "$status" -eq 0 ]
        jq -e '.rank==85 and .targets==1 and .targets_solved==1 and .targets_failed==0' "$output/summary.jsonl" > /dev/null
        jq -e '.exit_code==0 and .group_verified==true and (.recovered_scalar|type=="number")' "$output/target.jsonl" > /dev/null
    fi
}

run_case 20
run_case 85
echo "PASS: failed target exits 1 after preserving evidence; solved target exits 0"
