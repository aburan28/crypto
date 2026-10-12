#!/bin/sh
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$here/../.." && pwd)
build="$root/build/direct-sigma-map"
cxx=${CXX:-clang++}
mkdir -p "$build/tmp"
export TMPDIR="$build/tmp"

common="-O3 -std=c++17 -DECC_PACKED_PERM_SIGMA=1"

"$cxx" $common -Wall -Wextra -Werror "$here/synthesize.cpp" -o "$build/synthesize"
(cd "$root" && "$build/synthesize" \
  benchmarks/direct-sigma-map/direct_sigma_table3.generated.h \
  benchmarks/direct-sigma-map/direct_sigma_half5.generated.h) | tee "$build/synthesis.txt"

"$cxx" $common -Wall -Wextra -Werror "$here/verify_generated.cpp" -o "$build/verify-generated"
"$build/verify-generated" | tee "$build/verification.txt"

"$cxx" -O2 -std=c++17 -Wall -Wextra -Werror "$here/count_assembly.cpp" \
  -o "$build/count-assembly"

for reduced in 0 1; do
  "$cxx" $common -S -mllvm -inline-threshold=1000000 \
    -DECC_PACKED_FROM_REDUCED=$reduced "$here/static_compile.cpp" \
    -o "$build/static-inline-reduced-$reduced.s"
  "$build/count-assembly" "$build/static-inline-reduced-$reduced.s" \
    > "$build/counts-inline-reduced-$reduced.txt"
  cat "$build/counts-inline-reduced-$reduced.txt"
  "$cxx" $common -c -mllvm -inline-threshold=1000000 \
    -DECC_PACKED_FROM_REDUCED=$reduced "$here/static_compile.cpp" \
    -o "$build/static-inline-reduced-$reduced.o"
  nm -nm "$build/static-inline-reduced-$reduced.o" \
    > "$build/symbols-inline-reduced-$reduced.txt"
  size -m "$build/static-inline-reduced-$reduced.o" \
    > "$build/size-inline-reduced-$reduced.txt"
done

"$cxx" $common -S -mllvm -inline-threshold=1000000 "$here/static_compile5.cpp" \
  -o "$build/static-half5.s"
"$build/count-assembly" "$build/static-half5.s" > "$build/counts-half5.txt"
cat "$build/counts-half5.txt"
"$cxx" $common -c -mllvm -inline-threshold=1000000 "$here/static_compile5.cpp" \
  -o "$build/static-half5.o"
nm -nm "$build/static-half5.o" > "$build/symbols-half5.txt"
size -m "$build/static-half5.o" > "$build/size-half5.txt"

shasum -a 256 "$here"/*.cpp "$here"/*.generated.h "$build"/*.txt \
  > "$build/files.sha256"
echo "PASS: native direct-sigma static screen"
