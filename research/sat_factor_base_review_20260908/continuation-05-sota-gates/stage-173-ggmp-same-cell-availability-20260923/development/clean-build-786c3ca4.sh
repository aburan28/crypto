#!/bin/zsh
set -eu

repo=/Volumes/SSD990/crypto-f4-single-target-stage171
root=/Volumes/SSD990/koblitz-ggmp-build24-786c3ca4
commit=786c3ca41d5d4643143c671595e4c1bfede4561b

mkdir -p "$root/source" "$root/build" "$root/bin" "$root/source/.cargo"
/usr/bin/git -C "$repo" archive --format=tar --output="$root/source.tar" "$commit" \
  Cargo.toml \
  research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock \
  src \
  examples/koblitz_public_factor_base_discovery.rs \
  docs/ic/calibration.json
/usr/bin/tar -xf "$root/source.tar" -C "$root/source"
/bin/cp \
  "$root/source/research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock" \
  "$root/source/Cargo.lock"
/bin/cp \
  /Volumes/SSD990/koblitz-native-f4-build20-080d3c2b/source/.cargo/config.toml \
  "$root/source/.cargo/config.toml"
cd "$root/source"
CARGO_HOME=/Volumes/SSD990/koblitz-native-f4-build20-080d3c2b/cargo-home \
CARGO_INCREMENTAL=0 \
RAYON_NUM_THREADS=1 \
/opt/homebrew/bin/cargo build --release --locked --offline --jobs 1 \
  --target-dir "$root/build" \
  --example koblitz_public_factor_base_discovery
/bin/cp \
  "$root/build/release/examples/koblitz_public_factor_base_discovery" \
  "$root/bin/koblitz_public_factor_base_discovery"
