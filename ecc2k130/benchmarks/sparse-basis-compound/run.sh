#!/bin/sh
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$here/../.." && pwd)
build="$root/build/sparse-basis-compound"
cxx=${CXX:-clang++}

mkdir -p "$build/tmp"
export TMPDIR="$build/tmp"

"$cxx" -O3 -std=c++17 -Wall -Wextra -Werror -Wno-unknown-pragmas \
  -DECC_PACKED_PERM_SIGMA=1 -DECC_PACKED_FROM_REDUCED=1 \
  -DECC_PACKED_DIRECT_REDUCE=1 -DECC_PACKED_INLINE_POLY=3 \
  "$here/audit.cpp" -o "$build/audit"

"$build/audit" > "$here/result.json"

shasum -a 256 "$here/PROTOCOL.md" "$here/audit.cpp" "$here/run.sh" \
  "$here/result.json" > "$here/files.sha256"

cat "$here/result.json"
echo "PASS: native sparse-basis compound static screen"
