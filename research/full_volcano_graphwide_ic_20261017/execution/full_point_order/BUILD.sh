#!/usr/bin/env bash
# Native standalone diagnostic. Explicit matching dependencies, fresh output.
set -euo pipefail
rustc_bin=${1:?rustc}
deps=${2:?dependency directory}
crypto_rlib=${3:?matching crypto rlib}
serde_rlib=${4:?matching serde rlib}
sqlite_rlib=${5:?matching SQLite rlib}
build_out=${6:?fresh output directory}
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
mkdir -- "$build_out"
common=(--edition=2021 -D warnings -C opt-level=3 -C strip=symbols
  --extern "crypto_lib=$crypto_rlib" --extern "serde_json=$serde_rlib"
  --extern "rusqlite=$sqlite_rlib" -L "dependency=$deps")
"$rustc_bin" --version > "$build_out/toolchain.txt"
sha256sum "$crypto_rlib" "$serde_rlib" "$sqlite_rlib" > "$build_out/linked-rlibs.sha256"
sha256sum "$script_dir/diagnostic.rs" "$script_dir/PROTOCOL.md" > "$build_out/source.sha256"
"$rustc_bin" "${common[@]}" --test "$script_dir/diagnostic.rs" -o "$build_out/full-point-tests"
"$build_out/full-point-tests" --test-threads=1 > "$build_out/tests.log" 2>&1
"$rustc_bin" "${common[@]}" "$script_dir/diagnostic.rs" -o "$build_out/full-point"
"$rustc_bin" "${common[@]}" "$script_dir/summarize.rs" -o "$build_out/full-point-summary"
sha256sum "$build_out/full-point" > "$build_out/binary.sha256"
