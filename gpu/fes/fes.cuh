// gpu/fes/fes.cuh -- Gray-code fast exhaustive search (FES) for quadratic
// Boolean systems, written to compile BOTH as host C++ (g++/clang++) and as
// CUDA (nvcc), so the algorithm is verified on a GPU-less machine and only the
// launch is device-specific.  This is the parallelizable core of a binary
// index-calculus decomposition: the 2^n cube of assignments splits into
// independent sub-cubes, one per thread.
//
// System (packed, matching fref.py):
//   value(x) = CONST ^ (+)_{i: x_i=1} LIN[i] ^ (+)_{i>=j: x_i=x_j=1} QUAD(i,j)
// where each coefficient is a u64 whose bit e is that term's coefficient in
// equation e, and QUAD is stored as a lower-triangular flat array
//   QUAD_TRI[i*(i+1)/2 + j] = coefficient of x_i x_j   (i >= j).
// A point x is a solution iff value(x) == 0.
//
// The walk maintains, at the current point, the first derivative df[i] of the
// whole system along each variable i (df[i] = f(x ^ e_i) ^ f(x)).  Flipping
// variable i costs one XOR into the value plus n XORs to update the other
// derivatives by the (constant) second derivatives QUAD(i,j).  That is O(n)
// per point; the classic libfes O(1)-amortized refinement layers a second
// Gray code over the df updates and is a drop-in optimization of this loop.
#pragma once

#include <stdint.h>

#if defined(__CUDACC__)
#define FES_HD __host__ __device__
#else
#define FES_HD static inline
#endif

// Lower-triangular index for the unordered pair {i, j}.
FES_HD int fes_tri_index(int i, int j) {
  int a = i > j ? i : j;
  int b = i > j ? j : i;
  return a * (a + 1) / 2 + b;
}

FES_HD uint64_t fes_quad(const uint64_t *quad_tri, int i, int j) {
  return quad_tri[fes_tri_index(i, j)];
}

// Direct evaluation, used to verify the walk and to initialize a sub-cube.
FES_HD uint64_t fes_eval(uint64_t cst, const uint64_t *lin,
                         const uint64_t *quad_tri, int n, uint64_t x) {
  uint64_t v = cst;
  for (int i = 0; i < n; ++i) {
    if ((x >> i) & 1ULL) {
      v ^= lin[i];
      for (int j = 0; j <= i; ++j) {
        if ((x >> j) & 1ULL) {
          v ^= quad_tri[i * (i + 1) / 2 + j];
        }
      }
    }
  }
  return v;
}

// Reflected binary Gray code of s.
FES_HD uint64_t fes_gray(uint64_t s) { return s ^ (s >> 1); }

// Trailing-zero count for a nonzero 64-bit word.
FES_HD int fes_ctz(uint64_t s) {
#if defined(__CUDACC__)
  return __ffsll((unsigned long long)s) - 1;
#elif defined(__GNUC__)
  return __builtin_ctzll((unsigned long long)s);
#else
  int c = 0;
  while (((s >> c) & 1ULL) == 0ULL) ++c;
  return c;
#endif
}

// Enumerate one sub-cube: the top (n - low_bits) variables are fixed to
// `prefix` (bit p of prefix sets variable low_bits + p), and the low
// `low_bits` variables run over all 2^low_bits assignments in Gray-code order.
// Solutions x (full n-bit assignments) are written to `out` up to `max_out`;
// `*found` is the true count even if it exceeds `max_out`.
FES_HD void fes_search_subcube(uint64_t cst, const uint64_t *lin,
                               const uint64_t *quad_tri, int n, int low_bits,
                               uint64_t prefix, uint64_t *out, int *found,
                               int max_out, uint64_t *df_scratch) {
  const uint64_t base = prefix << low_bits;
  // Initialize df[] and the value at the sub-cube origin (low bits = 0).
  uint64_t *df = df_scratch; // length n, caller-provided (no dynamic alloc)
  for (int i = 0; i < n; ++i) {
    uint64_t d = lin[i] ^ fes_quad(quad_tri, i, i);
    for (int j = 0; j < n; ++j) {
      if (j != i && ((base >> j) & 1ULL)) {
        d ^= fes_quad(quad_tri, i, j);
      }
    }
    df[i] = d;
  }
  uint64_t val = fes_eval(cst, lin, quad_tri, n, base);
  if (val == 0) {
    if (*found < max_out) out[*found] = base;
    ++*found;
  }
  const uint64_t steps = (low_bits >= 63) ? 0 : (1ULL << low_bits);
  for (uint64_t s = 1; s < steps; ++s) {
    int i1 = fes_ctz(s); // flips a low variable, 0 <= i1 < low_bits
    val ^= df[i1];
    // x_{i1} toggled: refresh the other derivatives by the second derivative.
    for (int j = 0; j < n; ++j) {
      if (j != i1) df[j] ^= fes_quad(quad_tri, i1, j);
    }
    if (val == 0) {
      uint64_t x = base | fes_gray(s);
      if (*found < max_out) out[*found] = x;
      ++*found;
    }
  }
}
