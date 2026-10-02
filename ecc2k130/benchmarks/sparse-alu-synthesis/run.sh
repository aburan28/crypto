#!/bin/sh
set -eu

here=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH='' cd -- "$here/../.." && pwd)
build="$root/build/sparse-alu-synthesis"
cxx=${CXX:-clang++}

mkdir -p "$build/tmp"
export TMPDIR="$build/tmp"

"$cxx" -O3 -std=c++17 -Wall -Wextra -Werror -Wno-unknown-pragmas \
  -DECC_PACKED_PERM_SIGMA=1 -DECC_PACKED_FROM_REDUCED=1 \
  -DECC_PACKED_DIRECT_REDUCE=1 -DECC_PACKED_INLINE_POLY=3 \
  "$here/check.cpp" -o "$build/check"

"$build/check" > "$here/result.json"

(cd "$root" && shasum -a 256 \
  benchmarks/sparse-alu-synthesis/PROTOCOL.md \
  benchmarks/sparse-alu-synthesis/RESULTS.md \
  benchmarks/sparse-alu-synthesis/check.cpp \
  benchmarks/sparse-alu-synthesis/run.sh \
  benchmarks/sparse-alu-synthesis/result.json) > "$here/MANIFEST.sha256"

cat "$here/result.json"
echo "PASS: native sparse ALU synthesis screen"
