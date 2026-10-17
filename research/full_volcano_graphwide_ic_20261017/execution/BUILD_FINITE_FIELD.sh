#!/usr/bin/env bash
# Match all rlibs to the supplied dependency directory; do not guess suffixes.
# Compile/test only. Actual measured runs must wait for other heavy jobs to end.
set -euo pipefail
rustc_bin=${1:?rustc executable}
deps=${2:?dependency directory}
crypto_rlib=${3:?matching crypto_lib rlib}
serde_rlib=${4:?matching serde_json rlib}
sqlite_rlib=${5:?matching rusqlite rlib}
out_dir=${6:?fresh output directory}
task_tmp=${7:?linker temporary directory}
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
mkdir -- "$out_dir"
mkdir -p -- "$task_tmp"
common=(--edition=2021 -C opt-level=3 -C strip=symbols
  --extern "crypto_lib=$crypto_rlib" --extern "serde_json=$serde_rlib"
  --extern "rusqlite=$sqlite_rlib" -L "dependency=$deps")
"$rustc_bin" --version > "$out_dir/toolchain.txt"
sha256sum -- "$crypto_rlib" "$serde_rlib" "$sqlite_rlib" > "$out_dir/linked-rlibs.sha256"
TMPDIR="$task_tmp" "$rustc_bin" "${common[@]}" --test \
  "$script_dir/finite_field_k1.rs" -o "$out_dir/finite-field-k1-tests"
"$out_dir/finite-field-k1-tests" --test-threads=1 > "$out_dir/tests.log" 2>&1
TMPDIR="$task_tmp" "$rustc_bin" "${common[@]}" \
  "$script_dir/finite_field_k1.rs" -o "$out_dir/finite-field-k1"
sha256sum -- "$out_dir/finite-field-k1" > "$out_dir/binary.sha256"
