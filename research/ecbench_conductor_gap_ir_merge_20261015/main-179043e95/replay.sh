#!/usr/bin/env bash
set -euo pipefail
cd /workspace/crypto
out=/workspace/pr1850-followup/latest-main-20261015
export CARGO_HOME=/workspace/toolchains/cargo
export RUSTUP_HOME=/workspace/toolchains/rustup
export CARGO_TARGET_DIR=/workspace/toolchains/volcano-ir/test-target
export CARGO_PROFILE_RELEASE_CODEGEN_UNITS=256
cp -a /workspace/pr1850-followup/replay "$out/replay"
replay_study="$out/replay/research/ecbench_conductor_gap_ir_20261014"
python3 tools/isolated_bench.py busy -- /workspace/toolchains/cargo/bin/cargo run --locked --release --example conductor_gap_ir -j3 -- analyze "$replay_study" "$out/replay/research/ecbench_conductor_gap_20261012"
cmp research/ecbench_conductor_gap_ir_20261014/analysis.json "$replay_study/analysis.json"
cmp research/ecbench_conductor_gap_ir_20261014/analysis.md "$replay_study/analysis.md"
sha256sum "$replay_study/analysis.json" "$replay_study/analysis.md" > "$out/native-replay.sha256"
