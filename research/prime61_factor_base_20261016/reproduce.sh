#!/bin/sh
# Thin build/run orchestration; all research work runs in the native driver.
set -eu
root=$(git rev-parse --show-toplevel)
cd "$root"
out=${1:?usage: reproduce.sh NEW_OUTPUT_DIRECTORY}
mkdir "$out"
cargo build --release --no-default-features --example prime61_factor_base
timeout 300 target/release/examples/prime61_factor_base 20 100000 1 16 > "$out/control.jsonl" 2> "$out/control.stderr"
timeout 300 target/release/examples/prime61_factor_base 61 100000 12 3 > "$out/prime61.jsonl" 2> "$out/prime61.stderr"
timeout 300 target/release/examples/prime61_factor_base 61 100000 12 3 > "$out/replay-prime61.jsonl" 2> "$out/replay-prime61.stderr"
cmp "$out/prime61.jsonl" "$out/replay-prime61.jsonl"
cmp research/prime61_factor_base_20261016/results/control.jsonl "$out/control.jsonl"
cmp research/prime61_factor_base_20261016/results/prime61.jsonl "$out/prime61.jsonl"
sha256sum "$out"/*.jsonl > "$out/records.sha256"
