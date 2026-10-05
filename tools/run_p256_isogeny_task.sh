#!/bin/sh
# Build and execute one content-addressed P-256 walk task locally.
# This script never contacts TaskQ, Cairn, AWS, S3, or any other network.
set -eu

usage() {
    cat <<'EOF'
usage: tools/run_p256_isogeny_task.sh [OUT_DIR] [STEPS] [SOURCE_COMMIT]

Defaults:
  OUT_DIR       ./p256-isogeny-task-<UTC timestamp>
  STEPS         4096
  SOURCE_COMMIT current Git HEAD

The output directory must not already exist. The script builds the release
binary, creates a task manifest, runs it, independently replays the result,
and derives an offline coordination-only Cairn receipt. It does not submit or
publish anything.
EOF
}

case "${1:-}" in
    -h|--help)
        usage
        exit 0
        ;;
esac

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
out_dir=${1:-"$PWD/p256-isogeny-task-$(date -u +%Y%m%dT%H%M%SZ)"}
steps=${2:-4096}
source_commit=${3:-$(git -C "$repo_root" rev-parse HEAD)}

if [ -e "$out_dir" ]; then
    printf 'p256-isogeny-task: refusing existing output path %s\n' "$out_dir" >&2
    exit 1
fi

case "$source_commit" in
    *[!0-9a-f]*)
        printf 'p256-isogeny-task: SOURCE_COMMIT must be exactly 40 lowercase hex characters\n' >&2
        exit 1
        ;;
esac
if [ "${#source_commit}" -ne 40 ]; then
    printf 'p256-isogeny-task: SOURCE_COMMIT must be exactly 40 lowercase hex characters\n' >&2
    exit 1
fi

mkdir -p "$out_dir"
run_id=local-$(date -u +%Y%m%dT%H%M%SZ)
task_path=$out_dir/input-task.json
receipt_path=$out_dir/cairn-receipt.json

cargo build --manifest-path "$repo_root/Cargo.toml" --locked --release --bin p256_isogeny_task
binary=$repo_root/target/release/p256_isogeny_task

"$binary" manifest \
    --run-id "$run_id" \
    --task-id prefix-0 \
    --source-commit "$source_commit" \
    --steps "$steps" \
    --output "$task_path"

"$binary" run --task "$task_path" --out "$out_dir/result"
"$binary" verify --task "$task_path" --dir "$out_dir/result"
"$binary" cairn-receipt \
    --task "$task_path" \
    --dir "$out_dir/result" \
    --output "$receipt_path"

printf 'p256-isogeny-task: verified result in %s\n' "$out_dir"
printf 'p256-isogeny-task: no queue or network submission was made\n'
