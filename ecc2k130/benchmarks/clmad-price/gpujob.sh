#!/bin/bash
# Build and run the carry-less instruction-rate probe on whatever GPU the
# container has (TWO-CHAINS.md section 6: does the CLMAD rate track the FP64
# pipe across SKUs?).  Same container contract as benchmarks/two-chains/gpujob.sh:
#   modal run modal_job.py --job benchmarks/clmad-price/gpujob.sh --out DIR --gpu B200
set -uo pipefail
cd /work || exit 1
R=${RESULTS:-/results}; mkdir -p "$R"
{
  nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,power.limit,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-unknown}"
} | tee "$R/host.txt"
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
echo "=== build probe for sm_$CAP"
nvcc -O3 -std=c++17 -arch=sm_$CAP -Xptxas -v benchmarks/clmad-price/probe.cu -o probe > "$R/build.log" 2>&1 || { tail -20 "$R/build.log"; exit 1; }
grep -E "registers|spill" "$R/build.log" | head -4
for pass in 1 2 3; do
  echo "=== probe pass $pass"
  ./probe 20000 3 2>&1 | tee -a "$R/probe.txt"
  nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader | tee -a "$R/probe.txt"
done
echo "=== done"
