// gpu/fes/fes.metal -- Metal compute shader for the Gray-code FES search, the
// Apple-silicon counterpart of fes.cu.  Same algorithm as fes.cuh (see there
// for the derivation); MSL has different builtins so the kernel is
// self-contained.  Compile-checked on macOS with
//   xcrun -sdk macosx metal -c fes.metal -o fes.air
// One thread per sub-cube (fixed top n-low_bits variables).
#include <metal_stdlib>
using namespace metal;

// Coefficient of x_i x_j from the lower-triangular flat array.
static inline ulong fes_quad(device const ulong *quad_tri, int i, int j) {
  int a = i > j ? i : j;
  int b = i > j ? j : i;
  return quad_tri[a * (a + 1) / 2 + b];
}

static inline ulong fes_eval(ulong cst, device const ulong *lin,
                             device const ulong *quad_tri, int n, ulong x) {
  ulong v = cst;
  for (int i = 0; i < n; ++i) {
    if ((x >> i) & 1ul) {
      v ^= lin[i];
      for (int j = 0; j <= i; ++j) {
        if ((x >> j) & 1ul) v ^= quad_tri[i * (i + 1) / 2 + j];
      }
    }
  }
  return v;
}

kernel void fes_kernel(device const ulong *lin [[buffer(0)]],
                       device const ulong *quad_tri [[buffer(1)]],
                       constant ulong &cst [[buffer(2)]],
                       constant int &n [[buffer(3)]],
                       constant int &low_bits [[buffer(4)]],
                       device ulong *out [[buffer(5)]],
                       device atomic_uint *count [[buffer(6)]],
                       constant int &max_out [[buffer(7)]],
                       uint gid [[thread_position_in_grid]]) {
  const ulong prefix = (ulong)gid;
  const ulong nprefix = 1ul << (uint)(n - low_bits);
  if (prefix >= nprefix) return;

  const ulong base = prefix << (uint)low_bits;
  thread ulong df[40];
  for (int i = 0; i < n; ++i) {
    ulong d = lin[i] ^ fes_quad(quad_tri, i, i);
    for (int j = 0; j < n; ++j) {
      if (j != i && ((base >> j) & 1ul)) d ^= fes_quad(quad_tri, i, j);
    }
    df[i] = d;
  }
  ulong val = fes_eval(cst, lin, quad_tri, n, base);
  if (val == 0ul) {
    uint idx = atomic_fetch_add_explicit(count, 1u, memory_order_relaxed);
    if ((int)idx < max_out) out[idx] = base;
  }
  const ulong steps = (low_bits >= 63) ? 0ul : (1ul << (uint)low_bits);
  for (ulong s = 1ul; s < steps; ++s) {
    int i1 = ctz(s);
    val ^= df[i1];
    for (int j = 0; j < n; ++j) {
      if (j != i1) df[j] ^= fes_quad(quad_tri, i1, j);
    }
    if (val == 0ul) {
      ulong x = base | (s ^ (s >> 1));
      uint idx = atomic_fetch_add_explicit(count, 1u, memory_order_relaxed);
      if ((int)idx < max_out) out[idx] = x;
    }
  }
}
