#!/bin/sh
set -eu

study_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$study_dir/../.." && pwd)
producer=${1:?pass compiled release producer binary}
expected_source=d4bf0f34d9d57eb26ed795a610be4baa469b3835455d69ff902ef46a20036fdd
expected_binary=fc4b62da8e5dc7980c3e2a27e774166d594b5661641833252e8ec8ae6894fa40
expected_wrapper=1b1b175e9e06e1b27943a213b419b751ebaed109a553961529a644aea4907ec0
expected_points=e9ed7ddb0aeafdbfc7354d67279a8712572eea9d2aee3896f99633bdc2e937cc
source_hash=$(shasum -a 256 "$repo_dir/examples/koblitz_orbit_dlp_fast.rs" | awk '{print $1}')
binary_hash=$(shasum -a 256 "$producer" | awk '{print $1}')
wrapper_hash=$(shasum -a 256 "$study_dir/run_with_usage.py" | awk '{print $1}')
points_hash=$(shasum -a 256 "$study_dir/targets.jsonl" | awk '{print $1}')
[ "$source_hash" = "$expected_source" ] || { echo "source hash mismatch" >&2; exit 2; }
[ "$binary_hash" = "$expected_binary" ] || { echo "binary hash mismatch" >&2; exit 2; }
[ "$wrapper_hash" = "$expected_wrapper" ] || { echo "resource wrapper hash mismatch" >&2; exit 2; }
[ "$points_hash" = "$expected_points" ] || { echo "points hash mismatch" >&2; exit 2; }

for entry in k2202_load75:2202:0:1 k2000_load75:2000:0:1 k1798_dense:1798:1:0; do
  arm=${entry%%:*}
  rest=${entry#*:}
  columns=${rest%%:*}
  rest=${rest#*:}
  dense=${rest%%:*}
  load75=${rest#*:}
  run_dir=$study_dir/runs/panel/$arm
  [ ! -e "$run_dir" ] || { echo "run directory exists: $run_dir" >&2; exit 2; }
  mkdir -p "$run_dir"
  printf 'phase\tpanel\narm\t%s\ncolumns\t%s\ndense\t%s\nload75\t%s\nsource_sha256\t%s\nbinary_sha256\t%s\nwrapper_sha256\t%s\npoints_sha256\t%s\n' \
    "$arm" "$columns" "$dense" "$load75" "$source_hash" "$binary_hash" "$wrapper_hash" "$points_hash" > "$run_dir/launch.tsv"
  date -u +%Y-%m-%dT%H:%M:%SZ > "$run_dir/started_utc.txt"
  set +e
  python3 "$study_dir/run_with_usage.py" --receipt "$run_dir/usage.json" --limit-seconds 900 -- \
    env KIC_PACKED_ROOT_TABLE=1 KIC_PACKED_ROOT_TABLE_DENSE="$dense" KIC_PACKED_ROOT_TABLE_LOAD75="$load75" \
      KIC_RANK_THREADS=12 KIC_TARGET_THREADS=1 KIC_DUMP_BASE="$run_dir/base_dump.jsonl" \
      "$producer" "construct:83:1:$columns" "$study_dir/targets.jsonl" 29 \
      "$run_dir/targets_result.jsonl" > "$run_dir/stdout.jsonl" 2> "$run_dir/stderr.log"
  rc_n83=$?
  set -e
  printf '%s\n' "$rc_n83" > "$run_dir/exit_status.txt"
  date -u +%Y-%m-%dT%H:%M:%SZ > "$run_dir/finished_utc.txt"
  echo "$arm status=$rc_n83"
done
