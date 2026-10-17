#!/usr/bin/env bash
# Reproduce from repository root. Paths below are relative to that root.
set -euo pipefail
report=research/cm_volcano_crt_20261016
cargo test --release --manifest-path isogeny_algos/Cargo.toml --lib -- --nocapture > "$report/records/library-tests.log" 2>&1
cargo test --release --manifest-path isogeny_algos/Cargo.toml --test api --test sutherland > "$report/records/integration-tests.log" 2>&1
cargo fmt --manifest-path isogeny_algos/Cargo.toml --check > "$report/records/format.log" 2>&1
cargo clippy --release --manifest-path isogeny_algos/Cargo.toml --all-targets -- -D warnings > "$report/records/clippy-all-targets.log" 2>&1
cargo build --release --manifest-path isogeny_algos/Cargo.toml --lib
mapfile -t rlibs < <(find isogeny_algos/target/release/deps -maxdepth 1 -name 'libisogeny_algos-*.rlib' -print)
# A clean target directory has one rlib. Choose the newest when replaying with other profiles.
rlib=$(ls -t "${rlibs[@]}" | head -1)
rustc --edition 2021 -O "$report/validate.rs" --extern "isogeny_algos=$rlib" -L dependency=isogeny_algos/target/release/deps -o /tmp/cm-modpoly-validate
/tmp/cm-modpoly-validate > "$report/records/acceptance.jsonl" 2> "$report/records/acceptance.log"
# On this host pkg-config is absent; these are the installed OpenSSL locations.
OPENSSL_LIB_DIR=/usr/lib/x86_64-linux-gnu OPENSSL_INCLUDE_DIR=/usr/include CARGO_BUILD_JOBS=1 cargo test --release --lib > "$report/records/root-library-tests.log" 2>&1
# Preserve log content; normalize terminal blank lines for git diff --check.
perl -0pi -e 's/\n+\z/\n/' "$report"/records/*.log
git diff --check
