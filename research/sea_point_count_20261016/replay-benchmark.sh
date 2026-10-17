#!/usr/bin/env bash
# Thin build/run orchestration; measurement and analysis execute in Rust.
set -euo pipefail
export RUSTUP_TOOLCHAIN=${RUSTUP_TOOLCHAIN:-1.99.1}
study=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$study/../.." && pwd)
baseline=${BASELINE_WORKTREE:-/tmp/sea-original-replay}
changed=${CANDIDATE_WORKTREE:-/tmp/sea-candidate-replay}
out=${BENCH_OUTPUT:-$study/evidence/replay-$(date -u +%Y%m%dT%H%M%SZ)}
ref=cf2745214d83e7bbf737f07d758e486d6af8fa7f
candidate=44f8d17f1c5b8988a57a16e9dddc8e5ce099d325
native=research/sea_point_count_20261016/native
prepare() {
  local tree=$1 revision=$2
  if [ ! -e "$tree" ]; then
    git -C "$repo" worktree add --detach --no-checkout "$tree" "$revision"
    git -C "$tree" sparse-checkout set isogeny_algos
    git -C "$tree" checkout "$revision"
  fi
  test "$(git -C "$tree" rev-parse HEAD)" = "$revision"
  git -C "$tree" diff --exit-code HEAD -- isogeny_algos
  mkdir -p "$tree/$native/src/bin"
  cp "$study/native/Cargo.toml" "$study/native/Cargo.lock" "$tree/$native/"
  cp "$study/native/src/bin/measure.rs" "$tree/$native/src/bin/"
  flock /tmp/crypto-bench.lock cargo build --locked --release --manifest-path "$tree/$native/Cargo.toml" --bin measure
}
# Both timed arms stay frozen even if current CM generation is corrected later.
prepare "$baseline" "$ref"
prepare "$changed" "$candidate"
mkdir -p "$out/phi-cache"
flock /tmp/crypto-bench.lock cargo build --locked --release --manifest-path "$study/native/Cargo.toml" --bin reserve --bin analyze
flock /tmp/crypto-bench.lock env SEA_PHI_CACHE_DIR="$out/phi-cache" \
  "$baseline/$native/target/release/measure" sea192_cached > "$out/cache-preparation.json"
sha256sum "$changed/$native/src/bin/measure.rs" "$baseline/$native/src/bin/measure.rs" \
  "$changed/$native/target/release/measure" "$baseline/$native/target/release/measure" \
  "$out/phi-cache/"* > "$out/hashes.txt"
"$study/native/target/release/reserve" "$baseline/$native/target/release/measure" \
  "$changed/$native/target/release/measure" "$out" | tee "$out/suite.txt"
"$study/native/target/release/analyze" "$out" "$repo/docs/curves/registry.json"
