#!/usr/bin/env bash
# Exact local build selectors, retained for native preparation provenance.
set -euo pipefail
repo=${1:?repository directory}
target_dir=${2:?cargo debug dependency directory}
rustc_bin=${3:?rustc executable}
out_bin=${4:?output executable}
task_tmp=${5:?linker temporary directory}
mkdir -p -- "$task_tmp"
# Build the root library first; these suffixes belong to the recorded build.
TMPDIR="$task_tmp" "$rustc_bin" --edition=2021 -C strip=symbols \
  "$repo/research/full_volcano_graphwide_ic_20261017/inventory.rs" \
  --extern "crypto_lib=$target_dir/libcrypto_lib-038cc0b52d457b73.rlib" \
  --extern "serde_json=$target_dir/libserde_json-103d5c7fd4432996.rlib" \
  --extern "rusqlite=$target_dir/librusqlite-36fcccce347a3989.rlib" \
  -L "dependency=$target_dir" -o "$out_bin"
