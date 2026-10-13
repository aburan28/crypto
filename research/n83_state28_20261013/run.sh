#!/bin/sh
set -eu

study_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$study_dir/../.." && pwd)
producer=${1:?pass compiled release producer binary}
phase=${2:?control or primary}
case "$phase" in control|primary) ;; *) echo "invalid phase" >&2; exit 2 ;; esac

expected_source=5f8e8ec971e159322f54639c9d8415024e0197360e36b3272e6f9e3ad4628596
expected_binary=13c223fdff4b5341b71277bed68f2bcfaf5d0f21013220bb9c75866f1fff2e15
expected_input=21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8
expected_wrapper=1b1b175e9e06e1b27943a213b419b751ebaed109a553961529a644aea4907ec0
source_hash=$(shasum -a 256 "$repo_dir/examples/koblitz_orbit_dlp_fast.rs" | awk '{print $1}')
binary_hash=$(shasum -a 256 "$producer" | awk '{print $1}')
input_hash=$(shasum -a 256 "$study_dir/primary_target.jsonl" | awk '{print $1}')
wrapper=$repo_dir/research/n83_load75_20261013/run_with_usage.py
wrapper_hash=$(shasum -a 256 "$wrapper" | awk '{print $1}')
[ "$source_hash" = "$expected_source" ] || { echo "source hash mismatch" >&2; exit 2; }
[ "$binary_hash" = "$expected_binary" ] || { echo "binary hash mismatch" >&2; exit 2; }
[ "$input_hash" = "$expected_input" ] || { echo "target hash mismatch" >&2; exit 2; }
[ "$wrapper_hash" = "$expected_wrapper" ] || { echo "wrapper hash mismatch" >&2; exit 2; }

if [ "$phase" = control ]; then
  arms='k512_wide:512:0:1:0 k512_packed:512:1:1:0'
  limit=300
else
  arms='k2202_packed75:2202:1:0:1 k1798_packed_dense:1798:1:1:0'
  limit=1800
fi

for entry in $arms; do
  arm=${entry%%:*}; rest=${entry#*:}
  columns=${rest%%:*}; rest=${rest#*:}
  packed_state=${rest%%:*}; rest=${rest#*:}
  dense=${rest%%:*}; load75=${rest#*:}
  run_dir=$study_dir/runs/$phase/$arm
  [ ! -e "$run_dir" ] || { echo "run directory exists: $run_dir" >&2; exit 2; }
  mkdir -p "$run_dir"
  printf 'arm\t%s\ncolumns\t%s\nstate28\t%s\ndense\t%s\nload75\t%s\nsource_sha256\t%s\nbinary_sha256\t%s\ninput_sha256\t%s\n' \
    "$arm" "$columns" "$packed_state" "$dense" "$load75" "$source_hash" "$binary_hash" "$input_hash" > "$run_dir/launch.tsv"
  date -u +%Y-%m-%dT%H:%M:%SZ > "$run_dir/started_utc.txt"
  set +e
  python3 "$wrapper" --receipt "$run_dir/usage.json" --limit-seconds "$limit" -- \
    env -i PATH=/usr/bin:/bin:/opt/homebrew/bin KIC_PACKED_ROOT_TABLE=1 \
      KIC_PACKED_ROOT_TABLE_DENSE="$dense" KIC_PACKED_ROOT_TABLE_LOAD75="$load75" \
      KIC_STATE28="$packed_state" KIC_ROOT_BLOOM=0 KIC_SPARSE_CANONICAL=0 \
      KIC_RANK_THREADS=12 KIC_TARGET_THREADS=1 KIC_DUMP_BASE="$run_dir/base_dump.jsonl" \
      "$producer" "construct:83:1:$columns" "$study_dir/primary_target.jsonl" 29 \
      "$run_dir/targets_result.jsonl" > "$run_dir/stdout.jsonl" 2> "$run_dir/stderr.log"
  rc_n83=$?
  set -e
  printf '%s\n' "$rc_n83" > "$run_dir/exit_status.txt"
  date -u +%Y-%m-%dT%H:%M:%SZ > "$run_dir/finished_utc.txt"
  echo "$arm status=$rc_n83"
done
