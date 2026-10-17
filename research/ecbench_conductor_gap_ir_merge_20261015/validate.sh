#!/usr/bin/env bash
# Validation only. No measurement or frozen input is changed.
set -u
cd /workspace/crypto || exit 2
out=/workspace/pr1850-followup
export CARGO_HOME=/workspace/toolchains/cargo
export RUSTUP_HOME=/workspace/toolchains/rustup
export CARGO_TARGET_DIR=/workspace/toolchains/volcano-ir/test-target
export CARGO_PROFILE_RELEASE_CODEGEN_UNITS=256
cargo=/workspace/toolchains/cargo/bin/cargo
date -u +%FT%TZ > "$out/validation.started"
git rev-parse HEAD MERGE_HEAD > "$out/validation.revisions"
"$cargo" --version > "$out/toolchain.txt"
: > "$out/validation.status"
run_check() {
  local label="$1"
  shift
  printf '%q ' "$@" > "$out/$label.command"
  printf '\n' >> "$out/$label.command"
  "$@" > "$out/$label.log" 2>&1
  local result=$?
  printf '%s=%s\n' "$label" "$result" >> "$out/validation.status"
}
run_check library-tests python3 tools/isolated_bench.py busy -- "$cargo" test --locked --release --lib -j3 -- --test-threads=1
run_check driver-tests python3 tools/isolated_bench.py busy -- "$cargo" test --locked --release --example conductor_gap_ir -j3
run_check cargo-check python3 tools/isolated_bench.py busy -- "$cargo" check --locked --release --lib -j3
date -u +%FT%TZ > "$out/validation.finished"
