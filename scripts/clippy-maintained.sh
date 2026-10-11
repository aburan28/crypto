#!/usr/bin/env bash
# Strict lint for maintained code. Frozen measurement source remains byte-bound
# and is compiled by rust-all-targets.yml, without current-style lint rewrites.
set -euo pipefail

repository_root=$(cd -- "$(dirname -- "$0")/.." && pwd)
cd -- "$repository_root"
metadata_file=$(mktemp "${TMPDIR:-/tmp}/clippy-metadata.XXXXXX")
target_file=$(mktemp "${TMPDIR:-/tmp}/clippy-targets.XXXXXX")
trap 'rm -f "$metadata_file" "$target_file"' EXIT

cargo metadata --no-deps --format-version 1 > "$metadata_file"
jq -r --arg root "$repository_root" '
  .packages[] | select(.name == "crypto") | .targets[]
  | select(.kind | index("example"))
  | select(.src_path | startswith($root + "/examples/"))
  | .name
' "$metadata_file" > "$target_file"

selected_examples=()
while IFS= read -r target; do
    selected_examples+=(--example "$target")
done < "$target_file"
if [ ${#selected_examples[@]} -eq 0 ]; then
    echo "no maintained examples found in Cargo metadata" >&2
    exit 1
fi
if [ "${1:-}" = "--print-targets" ]; then
    cat "$target_file"
    exit 0
fi
cargo clippy --keep-going --lib --bins --tests --benches "${selected_examples[@]}" -- -D warnings
