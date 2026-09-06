#!/bin/bash
# ptx_stats.sh -- static analysis of the kernels without a GPU.
#
# Counts PTX instructions per field multiply and reports the ptxas register /
# local-stack / occupancy figures for Hopper and both Blackwell targets.
# These are the numbers quoted in OPTIMIZATION_BLACKWELL.md.
#
# Needs a CUDA-capable clang (18+) and a ptxas new enough to know sm_100.
# Neither has to come from a full CUDA install; the pip wheels are enough:
#
#   pip install nvidia-cuda-nvcc-cu12 nvidia-cuda-runtime-cu12 \
#               nvidia-cuda-cccl-cu12 nvidia-curand-cu12
#   ./ptx_stats.sh --setup      # assemble a CUDA_PATH from those wheels
#   ./ptx_stats.sh
#
# clang 18 emits PTX for sm_90 at the newest; since the kernels use no
# architecture-specific instructions, retargeting the .target directive and
# letting ptxas compile for sm_100 / sm_120 gives valid Blackwell figures.
set -u

CUDA_FAKE="${CUDA_FAKE:-/tmp/cudafake}"
CLANG="${CLANG:-clang++}"
PTXAS="${PTXAS:-$CUDA_FAKE/bin/ptxas}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${WORK:-$(mktemp -d)}"

setup() {
    local nv
    nv=$(python3 -c 'import site;print(site.getsitepackages()[0])')/nvidia
    [ -d "$nv/cuda_nvcc" ] || { echo "install the nvidia-cuda-*-cu12 wheels first"; exit 1; }
    rm -rf "$CUDA_FAKE"
    mkdir -p "$CUDA_FAKE/bin" "$CUDA_FAKE/include/crt" "$CUDA_FAKE/nvvm"
    cp -r "$nv/cuda_nvcc/nvvm/libdevice" "$CUDA_FAKE/nvvm/"
    cp "$nv/cuda_nvcc/bin/ptxas" "$CUDA_FAKE/bin/"
    cp -r "$nv/cuda_runtime/include/"* "$CUDA_FAKE/include/"
    cp -r "$nv/cuda_nvcc/include/crt/"* "$CUDA_FAKE/include/crt/" 2>/dev/null
    cp -r "$nv/cuda_cccl/include/"* "$CUDA_FAKE/include/" 2>/dev/null
    cp -r "$nv/curand/include/"* "$CUDA_FAKE/include/" 2>/dev/null
    echo "CUDA Version 12.3.0" > "$CUDA_FAKE/version.txt"
    echo "assembled $CUDA_FAKE"
}

[ "${1:-}" = "--setup" ] && { setup; exit 0; }
[ -x "$PTXAS" ] || { echo "no ptxas at $PTXAS -- run '$0 --setup'"; exit 1; }

cd "$HERE"
make -s curve_secp256k1.h >/dev/null 2>&1 || true

CFLAGS="--cuda-device-only --cuda-path=$CUDA_FAKE --cuda-gpu-arch=sm_90
        -Wno-unknown-cuda-version -O3 -std=c++17 -I$HERE
        -DGPU_ECC_CURVE_HEADER=\"curve_secp256k1.h\""

# ---------------------------------------------------------------- #
# 1. instructions per field multiply, both reductions, both codegen paths
# ---------------------------------------------------------------- #
cat > "$WORK/count.cu" <<'EOF'
#include "fp256.cuh"
template<int N> __global__ void kmul(fp256 *o, const fp256 *i) {
    fp256 a = i[0], b = i[1];
#pragma unroll
    for (int k = 0; k < N; k++) a = Fp::mul(a, b);
    o[threadIdx.x] = a;
}
template __global__ void kmul<1>(fp256*, const fp256*);
template __global__ void kmul<11>(fp256*, const fp256*);
EOF

echo "=== PTX instructions per field multiply (sm_90, clang) ==="
for fast in 1 0; do
  for ptx in 0 1; do
    $CLANG $CFLAGS -DFP_FAST=$fast -DFP_PTX=$ptx -S \
        -o "$WORK/c_${fast}_${ptx}.ptx" "$WORK/count.cu" 2>/dev/null || continue
    python3 - "$WORK/c_${fast}_${ptx}.ptx" "$fast" "$ptx" <<'PY'
import re, sys
path, fast, ptx = sys.argv[1], sys.argv[2], sys.argv[3]
txt = open(path).read()
sizes = []
for m in re.finditer(r'\.visible \.entry (\S+?)\(', txt):
    i = txt.index('{', m.end()); j = i; d = 0
    while True:
        if txt[j] == '{': d += 1
        elif txt[j] == '}':
            d -= 1
            if d == 0: break
        j += 1
    ops = [l.strip() for l in txt[i:j].split('\n')
           if l.strip() and not l.strip().startswith(('//', '.', '$', '{', '}'))]
    sizes.append(len(ops))
sizes.sort()
red = "secp256k1 special" if fast == '1' else "Montgomery CIOS "
gen = "inline PTX" if ptx == '1' else "portable  "
print(f"  {red}  {gen}  {(sizes[1]-sizes[0])/10:6.1f}")
PY
  done
done

# ---------------------------------------------------------------- #
# 2. registers / stack / occupancy per architecture
# ---------------------------------------------------------------- #
cat > "$WORK/inst.cu" <<'EOF'
#include "kernels.cuh"
template __global__ void k_rho_walk<8>(rho_ctx, uint32_t);
template __global__ void k_rho_walk_lowmem<8>(rho_ctx, uint32_t);
EOF

echo
echo "=== ptxas: registers, local stack, occupancy (block = 128) ==="
for mb in 0 3 4; do
  $CLANG $CFLAGS -DRHO_MIN_BLOCKS=$mb -S -o "$WORK/i_$mb.ptx" "$WORK/inst.cu" 2>/dev/null
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
          thr = int(65536 / regs); blocks = int(thr / 128); res = blocks * 128
          printf "  %-8s min_blocks=%s %-24s %3d regs %5d B stack  %4d thr/SM (%4.1f%%)\n",
                 a, mb, substr(name, 1, 24), regs, stack, res, res * 100 / 2048
        }'
  done
done
echo
echo "(occupancy assumes 2048 threads/SM and 64K registers/SM)"
