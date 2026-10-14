#!/bin/sh
set -eu

study_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
producer=${1:?pass the compiled release producer binary}
expected=caa723d4e2894b8573c04f6ec2ec474ce4fb41eb8a2c452b793b3d2d286aa351
actual=$(shasum -a 256 "$producer" | awk '{print $1}')
[ "$actual" = "$expected" ] || { echo "producer binary hash mismatch" >&2; exit 2; }
target_hash=$(shasum -a 256 "$study_dir/primary_target.jsonl" | awk '{print $1}')
[ "$target_hash" = "$(head -1 "$study_dir/targets.jsonl" | shasum -a 256 | awk '{print $1}')" ] || {
  echo "primary target does not match frozen panel" >&2; exit 2;
}

run_arm() {
  arm=$1
  columns=$2
  dense=$3
  dir=$study_dir/runs/primary/$arm
  if [ -e "$dir" ]; then
    echo "run directory already exists: $dir" >&2
    exit 2
  fi
  mkdir -p "$dir"
  printf 'arm\t%s\ncolumns\t%s\ndense\t%s\nproducer_sha256\t%s\nprimary_target_sha256\t%s\n' \
    "$arm" "$columns" "$dense" "$actual" "$target_hash" > "$dir/launch.tsv"
  date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/started_utc.txt"
  set +e
  /usr/bin/time -l /opt/homebrew/bin/timeout 900 \
    env KIC_PACKED_ROOT_TABLE=1 KIC_PACKED_ROOT_TABLE_DENSE="$dense" \
    KIC_RANK_THREADS=12 KIC_TARGET_THREADS=1 \
    "$producer" "construct:83:1:$columns" "$study_dir/primary_target.jsonl" 29 \
    "$dir/targets_result.jsonl" > "$dir/stdout.jsonl" 2> "$dir/stderr.log"
  status=$?
  set -e
  printf '%s\n' "$status" > "$dir/exit_status.txt"
  date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/finished_utc.txt"
  echo "$arm status=$status"
}

run_arm k1271_regular 1271 0
run_arm k1271_dense 1271 1
run_arm k1500_dense 1500 1
run_arm k1500_regular 1500 0
