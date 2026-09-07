#!/bin/bash
# ptx_stats_kangaroo.sh -- static analysis of the kangaroo kernels, no GPU.
#
# Reports ptxas register / stack / occupancy figures for Hopper and both
# Blackwell targets.  These are the numbers quoted in README.md.
#
# Shares the throwaway CUDA_PATH that ../ecc/ptx_stats.sh builds:
#   ../ecc/ptx_stats.sh --setup
#   ./ptx_stats_kangaroo.sh
set -u
CUDA_FAKE="${CUDA_FAKE:-/tmp/cudafake}"
CLANG="${CLANG:-clang++}"
PTXAS="${PTXAS:-$CUDA_FAKE/bin/ptxas}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${WORK:-$(mktemp -d)}"

[ -x "$PTXAS" ] || { echo "no ptxas at $PTXAS -- run ../ecc/ptx_stats.sh --setup"; exit 1; }
cd "$HERE"
make -C ../ecc -s curve_secp256k1.h >/dev/null 2>&1 || true

CFLAGS="--cuda-device-only --cuda-path=$CUDA_FAKE --cuda-gpu-arch=sm_90
        -Wno-unknown-cuda-version -O3 -std=c++17 -I$HERE -I$HERE/../ecc"

cat > "$WORK/inst.cu" <<'EOF'
#include "kernels_kangaroo.cuh"
template __global__ void k_kang_walk<8>(kg_ctx, uint32_t);
template __global__ void k_kang_walk_lowmem<8>(kg_ctx, uint32_t);
EOF

echo "=== ptxas: registers, local stack, occupancy (block = 128) ==="
for mb in 0 4; do
  $CLANG $CFLAGS -DKG_MIN_BLOCKS=$mb -S -o "$WORK/i_$mb.ptx" "$WORK/inst.cu" 2>/dev/null
  for arch in sm_90 sm_100 sm_120; do
    sed -e 's/^\.version .*/.version 8.7/' -e "s/^\.target sm_90.*/.target $arch/" \
        "$WORK/i_$mb.ptx" > "$WORK/i_${mb}_$arch.ptx"
    "$PTXAS" -arch=$arch -O3 -v "$WORK/i_${mb}_$arch.ptx" -o /dev/null 2>&1 |
      awk -v a="$arch" -v mb="$mb" '
        /Compiling entry/ { split($0, f, "\x27"); name = f[2] }
        /Used [0-9]+ registers/ {
          for (i = 1; i <= NF; i++) if ($i == "registers,") regs = $(i-1)
          stack = 0
          for (i = 1; i <= NF; i++) if ($i == "cumulative") stack = $(i-2)
          thr = int(65536 / regs); res = int(thr / 128) * 128
          printf "  %-7s min_blocks=%s %-24s %3d regs %5d B stack %5d thr/SM\n",
                 a, mb, substr(name, 1, 24), regs, stack, res
        }'
  done
done
echo
echo "(threads/SM assumes 64K registers/SM; consumer parts cap SM occupancy at"
echo " 1536 threads rather than 2048 -- ./bench prints the real limit.)"
