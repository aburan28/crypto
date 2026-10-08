#!/bin/bash
# Compile-only CUDA 13.3/sm_120 resource producer. No GPU or driver required.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
OUT=${1:?usage: compilecheck.sh OUTPUT_DIR}
[ ! -e "$OUT" ] || { echo "output already exists: $OUT" >&2; exit 2; }
mkdir -p "$OUT"
cd "$ROOT"

nvcc --version | tee "$OUT/nvcc-version.txt"
grep -q 'release 13\.3, V13\.3\.73' "$OUT/nvcc-version.txt"
SOURCE_REV=${SOURCE_REV:-}
if git -C "$ROOT" rev-parse HEAD >/dev/null 2>&1; then
    GIT_REV=$(git -C "$ROOT" rev-parse HEAD)
    [ -z "$SOURCE_REV" ] || [ "$SOURCE_REV" = "$GIT_REV" ]
    SOURCE_REV=$GIT_REV
    [ -z "$(git -C "$ROOT" status --porcelain -- .)" ]
fi
[[ $SOURCE_REV =~ ^[0-9a-f]{40}$ ]]
printf '%s\n' "$SOURCE_REV" > "$OUT/source-rev.txt"

for spec in control:0 cache2:2 cache3:3 cache4:4; do
    name=${spec%%:*}
    slots=${spec##*:}
    make -B gpu-rtx-pro6000-sigma-fused \
        PRO6000_ARCH='-gencode arch=compute_120,code=sm_120' \
        SIGMA_FUSED_SHARED_SLOTS="$slots" > "$OUT/build-$name.log" 2>&1
    mv ecc2k130 "$OUT/ecc2k130-$name"
    cuobjdump --dump-resource-usage "$OUT/ecc2k130-$name" \
        > "$OUT/resources-$name.txt" 2>&1
done

make -B compile-sigma-fused-shared-scratch-cuda \
    PRO6000_ARCH='-gencode arch=compute_120,code=sm_120' \
    > "$OUT/device-build.log" 2>&1
cp build/test-sigma-fused-shared-scratch-cuda-2 \
   build/test-sigma-fused-shared-scratch-cuda-0 \
   build/test-sigma-fused-shared-scratch-cuda-3 \
   build/test-sigma-fused-shared-scratch-cuda-4 "$OUT/"

g++ -O2 -std=c++17 -Wall -Wextra -Werror \
    benchmarks/sigma-fused/shared-scratch/compile_audit.cpp \
    -o "$OUT/compile-audit"
"$OUT/compile-audit" --self-test > "$OUT/compile-audit-self-test.txt"
"$OUT/compile-audit" "$OUT" "$OUT/result.json" | tee "$OUT/compile-audit.txt"

sha256sum Makefile include/packed131.h include/packedcompactstate.cuh \
    include/packedengine.cuh include/packedkernels.cuh include/packedsigmascratch.h \
    src/main.cu src/testsigmasharedscratchcuda.cu \
    benchmarks/sigma-fused/shared-scratch/PROPOSAL.md \
    benchmarks/sigma-fused/shared-scratch/PROTOCOL.md \
    benchmarks/sigma-fused/shared-scratch/IMPLEMENTATION.md \
    benchmarks/sigma-fused/shared-scratch/ATTEMPTS.md \
    benchmarks/sigma-fused/shared-scratch/model.cpp \
    benchmarks/sigma-fused/shared-scratch/model.json \
    benchmarks/sigma-fused/shared-scratch/test_native.cpp \
    benchmarks/sigma-fused/shared-scratch/compile_audit.cpp \
    benchmarks/sigma-fused/shared-scratch/source_audit.cpp \
    benchmarks/sigma-fused/shared-scratch/summarize.cpp \
    benchmarks/sigma-fused/shared-scratch/log_check.cpp \
    benchmarks/sigma-fused/shared-scratch/gpujob.sh \
    benchmarks/sigma-fused/shared-scratch/compilecheck.sh \
    > "$OUT/source-files.sha256"
sha256sum "$OUT"/ecc2k130-* "$OUT"/test-sigma-fused-shared-scratch-cuda-* \
    "$OUT/compile-audit" > "$OUT/binary-files.sha256"
echo '=== compile-only gate complete'
