#!/bin/sh
set -eu

study_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
producer=${1:?pass the compiled release producer binary}
class=${2:?pass primary or panel}
case "$class" in
  primary) input=$study_dir/primary_target.jsonl; expected_input=21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8 ;;
  panel) input=$study_dir/targets.jsonl; expected_input=e9ed7ddb0aeafdbfc7354d67279a8712572eea9d2aee3896f99633bdc2e937cc ;;
  *) echo "class must be primary or panel" >&2; exit 2 ;;
esac
expected=caa723d4e2894b8573c04f6ec2ec474ce4fb41eb8a2c452b793b3d2d286aa351
actual=$(shasum -a 256 "$producer" | awk '{print $1}')
[ "$actual" = "$expected" ] || { echo "producer binary hash mismatch" >&2; exit 2; }
dir=$study_dir/runs/$class/k1798_dense
if [ -e "$dir" ]; then
  echo "run directory already exists: $dir" >&2
  exit 2
fi
mkdir -p "$dir"
input_hash=$(shasum -a 256 "$input" | awk '{print $1}')
[ "$input_hash" = "$expected_input" ] || { echo "input hash mismatch" >&2; exit 2; }
printf 'arm\tk1798_dense\ncolumns\t1798\ndense\t1\nproducer_sha256\t%s\ninput_sha256\t%s\n' \
  "$actual" "$input_hash" > "$dir/launch.tsv"
date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/started_utc.txt"
set +e
if [ "$class" = panel ]; then
  /usr/bin/time -l /opt/homebrew/bin/timeout 900 \
    env KIC_PACKED_ROOT_TABLE=1 KIC_PACKED_ROOT_TABLE_DENSE=1 \
    KIC_RANK_THREADS=12 KIC_TARGET_THREADS=1 KIC_DUMP_BASE="$dir/base_dump.jsonl" \
    "$producer" construct:83:1:1798 "$input" 29 \
    "$dir/targets_result.jsonl" > "$dir/stdout.jsonl" 2> "$dir/stderr.log"
else
  /usr/bin/time -l /opt/homebrew/bin/timeout 900 \
    env KIC_PACKED_ROOT_TABLE=1 KIC_PACKED_ROOT_TABLE_DENSE=1 \
    KIC_RANK_THREADS=12 KIC_TARGET_THREADS=1 \
    "$producer" construct:83:1:1798 "$input" 29 \
    "$dir/targets_result.jsonl" > "$dir/stdout.jsonl" 2> "$dir/stderr.log"
fi
status=$?
set -e
printf '%s\n' "$status" > "$dir/exit_status.txt"
date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/finished_utc.txt"
echo "k1798_dense $class status=$status"
