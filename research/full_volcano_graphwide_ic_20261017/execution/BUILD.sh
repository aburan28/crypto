#!/usr/bin/env bash
# Build the repository with cargo first. Supply the matching .rlibs from that
# build explicitly; never guess between incompatible cached dependencies.
set -euo pipefail
rustc_bin=${1:?rustc executable}
deps=${2:?dependency directory}
crypto_rlib=${3:?matching crypto_lib rlib}
serde_rlib=${4:?matching serde_json rlib}
sqlite_rlib=${5:?matching rusqlite rlib}
out_dir=${6:?output directory}
task_tmp=${7:?linker temporary directory}
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
mkdir -p -- "$out_dir" "$task_tmp"
for target in driver summarize checkpoint match_profiles; do
  TMPDIR="$task_tmp" "$rustc_bin" --edition=2021 -C strip=symbols \
    "$script_dir/$target.rs" \
    --extern "crypto_lib=$crypto_rlib" --extern "serde_json=$serde_rlib" \
    --extern "rusqlite=$sqlite_rlib" -L "dependency=$deps" \
    -o "$out_dir/volcano-ic-$target"
done
