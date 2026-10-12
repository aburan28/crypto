#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo 'usage: run.sh discovery|holdout OUTPUT_DIR [DISCOVERY_BUNDLE]' >&2
  exit 2
fi
phase=$1
output_dir=$2
discovery_bundle=${3:-}
if [[ "$phase" != discovery && "$phase" != holdout ]]; then
  echo 'phase must be discovery or holdout' >&2
  exit 2
fi
if [[ "$phase" == holdout && -z "$discovery_bundle" ]]; then
  echo 'holdout requires a qualified discovery bundle' >&2
  exit 2
fi
if [[ -e "$output_dir" ]]; then
  echo 'refuse to overwrite an attempt directory' >&2
  exit 2
fi

study_dir=$(cd "$(dirname "$0")" && pwd)
repo_root=$(git -C "$study_dir" rev-parse --show-toplevel)
mkdir -p "$output_dir"
output_dir=$(cd "$output_dir" && pwd)
for name in Cargo.toml Cargo.lock PROTOCOL.md README.md protocol.json worker.rs run.sh; do
  cp "$study_dir/$name" "$output_dir/$name"
done
git -C "$repo_root" rev-parse HEAD > "$output_dir/git-head.txt"
rustc --version --verbose > "$output_dir/rustc.txt"
uname -a > "$output_dir/host.txt"

python3 "$repo_root/tools/isolated_bench.py" busy -- \
  cargo test --offline --locked --release --manifest-path "$study_dir/Cargo.toml" -- --test-threads=1 \
  > "$output_dir/test.stdout" 2> "$output_dir/test.stderr"
python3 "$repo_root/tools/isolated_bench.py" busy -- \
  cargo build --offline --locked --release --manifest-path "$study_dir/Cargo.toml" \
  > "$output_dir/build.stdout" 2> "$output_dir/build.stderr"
cp "$study_dir/target/release/boolean-pivot-screen" "$output_dir/boolean-pivot-screen"

if [[ "$phase" == discovery ]]; then
  python3 "$repo_root/tools/isolated_bench.py" busy -- \
    "$output_dir/boolean-pivot-screen" --run discovery "$output_dir/result.json" \
    > "$output_dir/worker.stdout" 2> "$output_dir/worker.stderr"
else
  python3 "$repo_root/tools/isolated_bench.py" busy -- \
    "$output_dir/boolean-pivot-screen" --run holdout "$output_dir/result.json" "$discovery_bundle" \
    > "$output_dir/worker.stdout" 2> "$output_dir/worker.stderr"
fi
"$output_dir/boolean-pivot-screen" --seal "$output_dir"
"$output_dir/boolean-pivot-screen" --verify-bundle "$output_dir"
