#!/bin/sh
set -eu

study_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$study_dir/../.." && pwd)
producer=${1:?pass compiled release producer binary}
phase=${2:?primary or panel}
case "$phase" in primary|panel) ;; *) echo "invalid phase" >&2; exit 2;; esac

expected_source=d4bf0f34d9d57eb26ed795a610be4baa469b3835455d69ff902ef46a20036fdd
expected_binary=fc4b62da8e5dc7980c3e2a27e774166d594b5661641833252e8ec8ae6894fa40
expected_points=e9ed7ddb0aeafdbfc7354d67279a8712572eea9d2aee3896f99633bdc2e937cc
source_hash=$(shasum -a 256 "$repo_dir/examples/koblitz_orbit_dlp_fast.rs" | awk '{print $1}')
binary_hash=$(shasum -a 256 "$producer" | awk '{print $1}')
points_hash=$(shasum -a 256 "$study_dir/targets.jsonl" | awk '{print $1}')
[ "$source_hash" = "$expected_source" ] || { echo "source hash mismatch" >&2; exit 2; }
[ "$binary_hash" = "$expected_binary" ] || { echo "binary hash mismatch" >&2; exit 2; }
[ "$points_hash" = "$expected_points" ] || { echo "points hash mismatch" >&2; exit 2; }

if [ "$phase" = primary ]; then
  input=$study_dir/primary_target.jsonl
  [ "$(shasum -a 256 "$input" | awk '{print $1}')" = "$(head -1 "$study_dir/targets.jsonl" | shasum -a 256 | awk '{print $1}')" ] || exit 2
  arms='k1798_dense:1798:1:0 k2000_load75:2000:0:1 k2202_load75:2202:0:1'
else
  input=$study_dir/targets.jsonl
  arms='k2202_load75:2202:0:1 k2000_load75:2000:0:1 k1798_dense:1798:1:0'
fi

for entry in $arms; do
  arm=${entry%%:*}
  rest=${entry#*:}
  columns=${rest%%:*}
  rest=${rest#*:}
  dense=${rest%%:*}
  load75=${rest#*:}
  run_dir=$study_dir/runs/$phase/$arm
  [ ! -e "$run_dir" ] || { echo "run directory exists: $run_dir" >&2; exit 2; }
  mkdir -p "$run_dir"
  printf 'phase\t%s\narm\t%s\ncolumns\t%s\ndense\t%s\nload75\t%s\nsource_sha256\t%s\nbinary_sha256\t%s\npoints_sha256\t%s\n' \
    "$phase" "$arm" "$columns" "$dense" "$load75" "$source_hash" "$binary_hash" "$points_hash" > "$run_dir/launch.tsv"
  date -u +%Y-%m-%dT%H:%M:%SZ > "$run_dir/started_utc.txt"
  set +e
  if [ "$phase" = panel ]; then
    /usr/bin/time -l /opt/homebrew/bin/timeout 900 \
      env KIC_PACKED_ROOT_TABLE=1 KIC_PACKED_ROOT_TABLE_DENSE="$dense" KIC_PACKED_ROOT_TABLE_LOAD75="$load75" \
      KIC_RANK_THREADS=12 KIC_TARGET_THREADS=1 KIC_DUMP_BASE="$run_dir/base_dump.jsonl" \
      "$producer" "construct:83:1:$columns" "$input" 29 "$run_dir/targets_result.jsonl" \
      > "$run_dir/stdout.jsonl" 2> "$run_dir/stderr.log"
  else
    /usr/bin/time -l /opt/homebrew/bin/timeout 900 \
      env KIC_PACKED_ROOT_TABLE=1 KIC_PACKED_ROOT_TABLE_DENSE="$dense" KIC_PACKED_ROOT_TABLE_LOAD75="$load75" \
      KIC_RANK_THREADS=12 KIC_TARGET_THREADS=1 \
      "$producer" "construct:83:1:$columns" "$input" 29 "$run_dir/targets_result.jsonl" \
      > "$run_dir/stdout.jsonl" 2> "$run_dir/stderr.log"
  fi
  status=$?
  set -e
  printf '%s\n' "$status" > "$run_dir/exit_status.txt"
  date -u +%Y-%m-%dT%H:%M:%SZ > "$run_dir/finished_utc.txt"
  echo "$arm status=$status"
done
