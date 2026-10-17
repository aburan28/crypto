#!/bin/bash
# Thin orchestration only; all arithmetic, collection and summaries are native.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
out=${1:?usage: reproduce.sh NEW_OUTPUT_DIRECTORY}
mkdir "$out"
cargo build --release --no-default-features --bin isogeny_walk --example prime61_isogeny_experiment
cargo test --release --no-default-features --example prime61_isogeny_experiment
curve=(--curve custom --name prime61-isogeny-source --p 2305843009213693951 --a 2305843009213693948 --b 111 --order 2305843010818082053 --gx 1 --gy 1796364306280342617)
timeout 300 target/release/isogeny_walk walk "${curve[@]}" --primes 3,7,11,13,19,23,53,61 --max-curves 64 --max-depth 1 --audit-points 2 --threads 1 --no-class-audits --out "$out/walk" > "$out/walk.stdout" 2> "$out/walk.stderr"
target/release/isogeny_walk verify "${curve[@]}" --dir "$out/walk" > "$out/verify.json" 2> "$out/verify.stderr"
timeout 300 target/release/examples/prime61_isogeny_experiment "$out/walk/isogeny_routes.json" 100000 12 > "$out/panel.jsonl" 2> "$out/panel.stderr"
timeout 300 target/release/examples/prime61_isogeny_experiment "$out/walk/isogeny_routes.json" 100000 12 > "$out/replay.jsonl" 2> "$out/replay.stderr"
cmp "$out/panel.jsonl" "$out/replay.jsonl"
cmp research/prime61_isogeny_20261016/results/factor-base-l61-identified.jsonl "$out/panel.jsonl"
json_rlibs=(target/release/deps/libserde_json-*.rlib)
crypto_rlibs=(target/release/deps/libcrypto_lib-*.rlib)
rustc --edition=2021 research/prime61_isogeny_20261016/summary.rs -L dependency=target/release/deps --extern serde_json="${json_rlibs[0]}" -o "$out/summary"
"$out/summary" "$out/panel.jsonl" > "$out/summary.json"
cmp research/prime61_isogeny_20261016/results/summary.json "$out/summary.json"
rustc --edition=2021 research/prime61_isogeny_20261016/resources.rs -L dependency=target/release/deps --extern crypto_lib="${crypto_rlibs[0]}" -o "$out/resources"
"$out/resources" 68162026 > "$out/memory.json"
sha256sum "$out"/*.jsonl "$out/walk/isogeny_routes.json" > "$out/records.sha256"
# Separate structural audit, when native PARI/GP is installed:
# timeout 300 gp -fq research/prime61_isogeny_20261016/class_number.gp
# timeout 300 gp -fq research/prime61_isogeny_20261016/class_order_certificate.gp
# The first command previously timed out during bnfcertify; retain its exit.
