#!/bin/sh
set -eu

study_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
worker=${1:?pass the compiled release rho worker}
expected=90557e0c14f9fb652bb4a90c8c0a5719eeeefae26bf0f236f34a18ffa42f516c
actual=$(shasum -a 256 "$worker" | awk '{print $1}')
[ "$actual" = "$expected" ] || { echo "rho binary hash mismatch" >&2; exit 2; }
config_hash=$(shasum -a 256 "$study_dir/rho_config.json" | awk '{print $1}')
[ "$config_hash" = cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b ] || {
  echo "rho config hash mismatch" >&2; exit 2;
}
dir=$study_dir/runs/rho_primary
[ ! -e "$dir" ] || { echo "rho directory already exists" >&2; exit 2; }
mkdir -p "$dir"
printf 'worker_sha256\t%s\nconfig_sha256\t%s\nseed\t2026101301\n' \
  "$actual" "$config_hash" > "$dir/launch.tsv"
date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/started_utc.txt"
set +e
/usr/bin/time -l /opt/homebrew/bin/timeout 60 \
  "$worker" "$study_dir/primary_target.jsonl" 2026101301 \
  "$study_dir/rho_config.json" "$dir/receipt.json" \
  > "$dir/stdout.log" 2> "$dir/stderr.log"
status=$?
set -e
printf '%s\n' "$status" > "$dir/exit_status.txt"
date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/finished_utc.txt"
echo "rho primary status=$status"
