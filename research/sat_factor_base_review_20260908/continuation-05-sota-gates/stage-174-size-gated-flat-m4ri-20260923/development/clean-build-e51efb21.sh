#!/bin/zsh
set -eu

repo=/Volumes/SSD990/crypto-kic-stage174-native-f4
root=/Volumes/SSD990/koblitz-native-f4-build23-e51efb21
commit=e51efb219edd4511a08d545de21184928e5e771c

mkdir -p "$root/source" "$root/build" "$root/bin" "$root/tmp"
/usr/bin/git -C "$repo" archive --format=tar --output="$root/source.tar" "$commit" \
  Cargo.toml \
  research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock \
  src \
  examples/koblitz_pdp_export.rs \
  examples/koblitz_pdp_backend.rs \
  docs/ic/calibration.json
/usr/bin/tar -xf "$root/source.tar" -C "$root/source"
/bin/cp \
  "$root/source/research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock" \
  "$root/source/Cargo.lock"
mkdir -p "$root/source/.cargo"
/bin/cp \
  /Volumes/SSD990/koblitz-native-f4-build20-080d3c2b/source/.cargo/config.toml \
  "$root/source/.cargo/config.toml"
cd "$root/source"
TMPDIR="$root/tmp" \
CARGO_HOME=/Volumes/SSD990/koblitz-native-f4-build20-080d3c2b/cargo-home \
CARGO_INCREMENTAL=0 \
RAYON_NUM_THREADS=1 \
/opt/homebrew/bin/cargo build --release --locked --offline --jobs 1 \
  --target-dir "$root/build" \
  --example koblitz_pdp_export \
  --example koblitz_pdp_backend
/bin/cp "$root/build/release/examples/koblitz_pdp_export" "$root/bin/koblitz_pdp_export"
/bin/cp "$root/build/release/examples/koblitz_pdp_backend" "$root/bin/koblitz_pdp_backend"
