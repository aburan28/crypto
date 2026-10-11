#!/bin/bash
# gpu/ecc2k on one GPU: native clmad product (default device build) against
# the software product (-DF2M_CLMAD=0). Both self-test against the host
# implementation; then field and rho benchmarks alternate REPS times.
set -euo pipefail
cd "${WORK:-/work}/benchmarks/ecc2k-clmad/tree"
R=${RESULTS:-/results}; mkdir -p "$R"
ARCH=${ARCH:-sm_120}; REPS=${REPS:-3}
{ nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,power.limit,compute_cap --format=csv,noheader
  nvcc --version | tail -1; echo "source: ${SOURCE_REV:-unknown}"; echo "arch $ARCH reps $REPS"; } | tee "$R/host.txt"
sha256sum f2m.cuh koblitz.cuh rho2k.cuh kernels2k.cuh bench2k.cu > "$R/source-files.sha256"
python3 ecref2k.py params --curve ecc2k95 > curve_ecc2k95.h
python3 ecref2k.py vectors --curve ecc2k95 > vec_ecc2k95.h
nvcc -O3 -std=c++17 -arch=$ARCH -lineinfo -DGPU_ECC2K_CURVE_HEADER='"curve_ecc2k95.h"' -Xptxas -v -o bench_native bench2k.cu 2> "$R/ptxas-native.txt"
nvcc -O3 -std=c++17 -arch=$ARCH -lineinfo -DF2M_CLMAD=0 -DGPU_ECC2K_CURVE_HEADER='"curve_ecc2k95.h"' -Xptxas -v -o bench_software bench2k.cu 2> "$R/ptxas-software.txt"
cuobjdump -sass bench_native | grep -c CLMAD > "$R/clmad-count-native.txt" || true
cuobjdump -sass bench_software | grep -c CLMAD > "$R/clmad-count-software.txt" || true
for v in native software; do echo "=== selftest $v"; ./bench_$v selftest 2>&1 | tee "$R/selftest-$v.txt"; done
for i in $(seq 1 "$REPS"); do for v in software native; do
  echo "=== field $v rep $i"; ./bench_$v field 2>&1 | tee "$R/field-$v-$i.txt"
  echo "=== rho $v rep $i"; ./bench_$v rho --iters 2048 2>&1 | tee "$R/rho-$v-$i.txt"
done; done
echo JOB_OK
