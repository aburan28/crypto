#!/bin/bash
# Native research code; shell only controls fresh paths, logs and immutable chunks.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
out=${1:?usage: reproduce.sh NEW_OUTPUT_DIRECTORY}
mkdir "$out"
spec=research/prime61_disk_walk_20261017/spec.json
worker=target/release/examples/isogeny_disk_walk
cargo build --release --no-default-features --example isogeny_disk_walk --example prime61_isogeny_experiment
cargo test --release --no-default-features --example isogeny_disk_walk --example prime61_isogeny_experiment
"$worker" init "$spec" "$out/walk.db" > "$out/init.json"
for expanded in 1 8 32; do
    timeout 300 "$worker" advance "$spec" "$out/walk.db" "$expanded" > "$out/expand-$expanded.jsonl" 2> "$out/expand-$expanded.stderr"
    if [ "$expanded" = 1 ]; then measure=16; else measure=64; fi
    timeout 300 "$worker" measure "$spec" "$out/walk.db" "$measure" > "$out/measure-$expanded.json" 2> "$out/measure-$expanded.stderr"
    timeout 300 "$worker" verify "$spec" "$out/walk.db" 16 > "$out/verify-$expanded.jsonl" 2> "$out/verify-$expanded.stderr"
    "$worker" snapshot "$spec" "$out/walk.db" "$out/checkpoint-$expanded.db" > "$out/snapshot-$expanded.json"
    "$worker" dump "$spec" "$out/checkpoint-$expanded.db" > "$out/checkpoint-$expanded.jsonl"
    gzip -n -c "$out/checkpoint-$expanded.db" > "$out/checkpoint-$expanded.db.gz"
    sha256sum "$out/checkpoint-$expanded.db" "$out/checkpoint-$expanded.db.gz" "$out/checkpoint-$expanded.jsonl" > "$out/checkpoint-$expanded.sha256"
done
"$worker" status "$spec" "$out/walk.db" > "$out/final-status.json"
# Upload immutable checkpoint-*.db.gz plus spec and hash manifest through the
# connected Drive upload_file action. Read back metadata and restore/download
# the bytes before considering a remote checkpoint verified. Drive credentials
# are intentionally not embedded here. Close the local writer before archiving.
