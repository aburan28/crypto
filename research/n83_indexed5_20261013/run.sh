#!/bin/sh
set -eu

study_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$study_dir/../.." && pwd)
producer=${1:?pass compiled release producer binary}
phase=${2:?control, comparison or primary}
case "$phase" in
  control) arm=k512_indexed5; columns=512; limit=900 ;;
  comparison) arm=k1798_indexed5; columns=1798; limit=1800 ;;
  primary) arm=k2500_indexed5; columns=2500; limit=2400 ;;
  *) echo "invalid phase" >&2; exit 2 ;;
esac

expected_source=42f21942e5e85874dfe9017f720697ab07f4cb0546a5c7d27e448a5a61fc448b
expected_binary=7a0946cd6bfb419a4d3fa60c6d8dfd4a72197021f6c0d60e8339390d5f54136d
expected_input=21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8
expected_wrapper=1b1b175e9e06e1b27943a213b419b751ebaed109a553961529a644aea4907ec0
wrapper=$repo_dir/research/n83_load75_20261013/run_with_usage.py
source_hash=$(shasum -a 256 "$repo_dir/examples/koblitz_orbit_dlp_fast.rs" | awk '{print $1}')
binary_hash=$(shasum -a 256 "$producer" | awk '{print $1}')
input_hash=$(shasum -a 256 "$study_dir/primary_target.jsonl" | awk '{print $1}')
wrapper_hash=$(shasum -a 256 "$wrapper" | awk '{print $1}')
[ "$source_hash" = "$expected_source" ] || { echo "source hash mismatch" >&2; exit 2; }
[ "$binary_hash" = "$expected_binary" ] || { echo "binary hash mismatch" >&2; exit 2; }
[ "$input_hash" = "$expected_input" ] || { echo "target hash mismatch" >&2; exit 2; }
[ "$wrapper_hash" = "$expected_wrapper" ] || { echo "wrapper hash mismatch" >&2; exit 2; }

run_dir=$study_dir/runs/$phase/$arm
[ ! -e "$run_dir" ] || { echo "run directory exists: $run_dir" >&2; exit 2; }
mkdir -p "$run_dir"
printf 'arm\t%s\ncolumns\t%s\nsource_sha256\t%s\nbinary_sha256\t%s\ninput_sha256\t%s\nwrapper_sha256\t%s\n' \
  "$arm" "$columns" "$source_hash" "$binary_hash" "$input_hash" "$wrapper_hash" > "$run_dir/launch.tsv"
date -u +%Y-%m-%dT%H:%M:%SZ > "$run_dir/started_utc.txt"
set +e
python3 "$wrapper" --receipt "$run_dir/usage.json" --limit-seconds "$limit" -- \
  env -i PATH=/usr/bin:/bin:/opt/homebrew/bin KIC_PACKED_ROOT_TABLE=1 \
    KIC_PACKED_ROOT_TABLE_DENSE=0 KIC_PACKED_ROOT_TABLE_LOAD75=0 \
    KIC_STATE28=1 KIC_INDEXED_ROOT5=1 KIC_ROOT_BLOOM=0 KIC_SPARSE_CANONICAL=0 \
    KIC_RANK_THREADS=12 KIC_TARGET_THREADS=1 KIC_DUMP_BASE="$run_dir/base_dump.jsonl" \
    "$producer" "construct:83:1:$columns" "$study_dir/primary_target.jsonl" 29 \
    "$run_dir/targets_result.jsonl" > "$run_dir/stdout.jsonl" 2> "$run_dir/stderr.log"
rc_n83=$?
set -e
printf '%s\n' "$rc_n83" > "$run_dir/exit_status.txt"
date -u +%Y-%m-%dT%H:%M:%SZ > "$run_dir/finished_utc.txt"
echo "$arm status=$rc_n83"
