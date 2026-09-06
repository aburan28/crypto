#!/bin/bash
# ptx_stats2k.sh -- static analysis of the ECC2K kernels without a GPU.
#
# Reports instructions per binary-field operation and the ptxas register /
# stack / occupancy figures for Hopper and both Blackwell targets.  These
# are the numbers quoted in README.md.
#
# Shares the throwaway CUDA_PATH that gpu/ecc/ptx_stats.sh builds:
#
#   ../ecc/ptx_stats.sh --setup
#   ./ptx_stats2k.sh
set -u

CUDA_FAKE="${CUDA_FAKE:-/tmp/cudafake}"
CLANG="${CLANG:-clang++}"
PTXAS="${PTXAS:-$CUDA_FAKE/bin/ptxas}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${WORK:-$(mktemp -d)}"

[ -x "$PTXAS" ] || { echo "no ptxas at $PTXAS -- run ../ecc/ptx_stats.sh --setup"; exit 1; }
cd "$HERE"
make -s curve_ecc2k95.h >/dev/null 2>&1 || true

CFLAGS="--cuda-device-only --cuda-path=$CUDA_FAKE --cuda-gpu-arch=sm_90
        -Wno-unknown-cuda-version -O3 -std=c++17 -I$HERE
        -DGPU_ECC2K_CURVE_HEADER=\"curve_ecc2k95.h\""

cat > "$WORK/count.cu" <<'EOF'
#include "f2m.cuh"
template<int N> __global__ void kmul(f2e *o, const f2e *i) {
    f2e a = i[0], b = i[1];
#pragma unroll
    for (int k = 0; k < N; k++) a = F2::mul(a, b);
    o[threadIdx.x] = a;
}
template<int N> __global__ void ksqr(f2e *o, const f2e *i) {
    f2e a = i[0];
#pragma unroll
    for (int k = 0; k < N; k++) a = F2::sqr(a);
    o[threadIdx.x] = a;
}
template __global__ void kmul<1>(f2e*, const f2e*);
template __global__ void kmul<11>(f2e*, const f2e*);
template __global__ void ksqr<1>(f2e*, const f2e*);
template __global__ void ksqr<11>(f2e*, const f2e*);
EOF

echo "=== PTX instructions per field operation (sm_90, clang) ==="
$CLANG $CFLAGS -S -o "$WORK/c.ptx" "$WORK/count.cu" 2>/dev/null
python3 - "$WORK/c.ptx" <<'PY'
import re, sys
txt = open(sys.argv[1]).read()
funcs = {}
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
    funcs[m.group(1)] = ops
def pick(tag):
    ks = sorted([n for n in funcs if tag in n], key=lambda n: len(funcs[n]))
    return (len(funcs[ks[1]]) - len(funcs[ks[0]])) / 10.0, funcs[ks[1]]
mu, ops = pick('kmul')
sq, _ = pick('ksqr')
print("  multiply  %6.1f" % mu)
print("  squaring  %6.1f   (%.2f of a multiply)" % (sq, sq / mu))
w = sum(1 for o in ops if o.startswith('mul.wide'))
print("  %d widening multiplies per field multiply (9 carry-less 32x32 products)"
      % (w / 11))
PY

cat > "$WORK/inst.cu" <<'EOF'
#include "kernels2k.cuh"
template __global__ void k2k_rho_walk<8>(rho2k_ctx, uint32_t);
template __global__ void k2k_rho_walk_lowmem<8>(rho2k_ctx, uint32_t);
EOF

echo
echo "=== ptxas: registers, local stack, occupancy (block = 128) ==="
for mb in 0 4; do
  $CLANG $CFLAGS -DR2K_MIN_BLOCKS=$mb -S -o "$WORK/i_$mb.ptx" "$WORK/inst.cu" 2>/dev/null
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
          printf "  %-7s min_blocks=%s %-26s %3d regs %5d B stack %5d thr/SM\n",
                 a, mb, substr(name, 1, 26), regs, stack, res
        }'
  done
done
echo
echo "(threads/SM assumes 64K registers/SM; consumer parts cap SM occupancy at"
echo " 1536 threads rather than 2048 -- ./bench2k prints the real limit.)"
