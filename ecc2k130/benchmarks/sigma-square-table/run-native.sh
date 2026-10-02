#!/bin/sh
set -eu

here=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH='' cd -- "$here/../.." && pwd)
build="$here/.build"
cxx=${CXX:-clang++}

mkdir -p "$build"
flags='-O3 -std=c++17 -Wall -Wextra -Werror -Wno-unknown-pragmas -DECC_PACKED_DIRECT_REDUCE=1 -DECC_PACKED_GENERATED_PRODUCT=1 -DECC_PACKED_ALU_SQUARE=0'

# shellcheck disable=SC2086
"$cxx" $flags -DECC_SIGMA_SQUARE_TABLE=0 "$here/replay.cpp" -o "$build/replay-control"
# shellcheck disable=SC2086
"$cxx" $flags -DECC_SIGMA_SQUARE_TABLE=1 "$here/replay.cpp" -o "$build/replay-candidate"
"$build/replay-control" "$build/control.bin" > "$build/control.json"
"$build/replay-candidate" "$build/candidate.bin" > "$build/candidate.json"
cmp "$build/control.bin" "$build/candidate.bin"

guard_flags='-DECC_PACKED_POLY_STATE=1 -DECC_PACKED_CACHE_DENOM=1 -DECC_PACKED_POLY_CHAIN=1 -DECC_PACKED_PAIR_PRODUCTS=1 -DECC_PACKED_WEIGHTED_PREFIX=2 -DECC_PACKED_PERM_SIGMA=3 -DECC_PACKED_CHAINS=1 -DECC_TABLE_FUSED=0 -DECC_TABLE_TAG_DENOM=0 -DECC_PACKED_SLOT_PIPELINE=0 -DECC_PACKED_SLOT_PREFETCH=0 -DECC_PHASE_PROFILE=0 -DECC_PACKED_CHAIN_FIRST=0'
# shellcheck disable=SC2086
"$cxx" -E -x c++ $guard_flags -DECC_SIGMA_FUSED=1 -DECC_WALK_TABLE=0 \
  -DECC_SIGMA_SQUARE_TABLE=1 -DECC_PACKED_SQUARE_TABLE=0 -DECC_PACKED_ALU_SQUARE=0 \
  "$here/guard-probe.cu" > "$build/guard-valid.i"
# shellcheck disable=SC2086
if "$cxx" -E -x c++ $guard_flags -DECC_SIGMA_FUSED=0 -DECC_WALK_TABLE=0 \
  -DECC_SIGMA_SQUARE_TABLE=1 -DECC_PACKED_SQUARE_TABLE=0 -DECC_PACKED_ALU_SQUARE=0 \
  "$here/guard-probe.cu" > "$build/guard-no-fused.i" 2> "$build/guard-no-fused.err"; then
  echo "missing SIGMA_FUSED guard" >&2; exit 1
fi
grep -q 'ECC_SIGMA_SQUARE_TABLE requires the polynomial-state sigma-fused walk' \
  "$build/guard-no-fused.err"
# shellcheck disable=SC2086
if "$cxx" -E -x c++ $guard_flags -DECC_SIGMA_FUSED=1 -DECC_WALK_TABLE=0 \
  -DECC_SIGMA_SQUARE_TABLE=1 -DECC_PACKED_SQUARE_TABLE=0 -DECC_PACKED_ALU_SQUARE=1 \
  "$here/guard-probe.cu" > "$build/guard-alu.i" 2> "$build/guard-alu.err"; then
  echo "missing ALU-square combination guard" >&2; exit 1
fi
grep -q 'keep PACKED_SQUARE_TABLE and PACKED_ALU_SQUARE off' "$build/guard-alu.err"
# shellcheck disable=SC2086
if "$cxx" -E -x c++ $guard_flags -DECC_SIGMA_FUSED=1 -DECC_WALK_TABLE=0 \
  -DECC_SIGMA_SQUARE_TABLE=1 -DECC_PACKED_SQUARE_TABLE=1 -DECC_PACKED_ALU_SQUARE=0 \
  "$here/guard-probe.cu" > "$build/guard-table.i" 2> "$build/guard-table.err"; then
  echo "missing table-walk square combination guard" >&2; exit 1
fi
grep -q 'keep PACKED_SQUARE_TABLE and PACKED_ALU_SQUARE off' "$build/guard-table.err"

"$cxx" -O3 -std=c++17 -Wall -Wextra -Werror "$here/log_check.cpp" \
  -o "$build/log-check"
"$cxx" -O3 -std=c++17 -Wall -Wextra -Werror "$here/summarize.cpp" \
  -o "$build/summarize"
"$cxx" -O3 -std=c++17 -Wall -Wextra -Werror "$here/postrun_audit.cpp" \
  -o "$build/postrun-audit"
"$build/log-check" --self-test > "$build/log-check-self-test.txt"
"$build/summarize" --self-test > "$build/summarize-self-test.txt"
"$build/postrun-audit" --self-test > "$build/postrun-audit-self-test.txt"

"$cxx" -O3 -std=c++17 -Wall -Wextra -Werror -Wno-unknown-pragmas \
  -DECC_PACKED_DIRECT_REDUCE=1 -DECC_PACKED_GENERATED_PRODUCT=1 \
  "$here/audit.cpp" -o "$build/audit"
"$build/audit" "$build/control.bin" "$build/candidate.bin" "$root" > "$here/result.json"

(cd "$here" && shasum -a 256 .build/control.bin .build/candidate.bin) > "$here/streams.sha256"
(cd "$root" && shasum -a 256 \
  Makefile \
  include/kernel.h \
  include/packed131.h \
  include/packedkernels.cuh \
  include/packedengine.cuh \
  benchmarks/sigma-square-table/PROTOCOL.md \
  benchmarks/sigma-square-table/RESULTS.md \
  benchmarks/sigma-square-table/replay.cpp \
  benchmarks/sigma-square-table/audit.cpp \
  benchmarks/sigma-square-table/compile_audit.cpp \
  benchmarks/sigma-square-table/log_check.cpp \
  benchmarks/sigma-square-table/summarize.cpp \
  benchmarks/sigma-square-table/postrun_audit.cpp \
  benchmarks/sigma-square-table/gpujob.sh \
  benchmarks/sigma-square-table/guard-probe.cu \
  benchmarks/sigma-square-table/run-native.sh \
  benchmarks/sigma-square-table/attempt-1-compile-failure.json \
  benchmarks/sigma-square-table/attempt-2-native-portability-failure.json \
  benchmarks/sigma-square-table/compile-artifact.json \
  benchmarks/sigma-square-table/compile-artifact-core.json \
  benchmarks/sigma-square-table/compile-files.sha256 \
  benchmarks/sigma-square-table/compile-files-core.sha256 \
  benchmarks/sigma-square-table/compile-result.json \
  benchmarks/sigma-square-table/compile-result-core.json \
  benchmarks/sigma-square-table/independent-compile-review.json \
  benchmarks/sigma-square-table/result.json \
  benchmarks/sigma-square-table/streams.sha256) > "$here/MANIFEST.sha256"

cat "$here/result.json"
echo "PASS: sigma square table native replay and independent audit"
