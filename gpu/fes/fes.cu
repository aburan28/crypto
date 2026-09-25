// gpu/fes/fes.cu -- CUDA launcher for the Gray-code FES search.  One thread per
// sub-cube (a fixed assignment of the top n-low_bits variables); each thread
// runs the same fes_search_subcube the host test verifies.  Build with
// `make cuda` (nvcc); this file is the device-specific launch only, the
// algorithm lives in fes.cuh and is checked by `make test` without a GPU.
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#include "fes.cuh"
#include "fes_io.hpp"
#include "system.h"

// Upper bound on variables, so per-thread scratch is a fixed-size local array.
#ifndef FES_MAX_VARS
#define FES_MAX_VARS 40
#endif
// Per-sub-cube solution cap before spilling to the global counter only.
#ifndef FES_LOCAL_CAP
#define FES_LOCAL_CAP 64
#endif

__global__ void fes_kernel(uint64_t cst, const uint64_t *lin,
                           const uint64_t *quad_tri, int n, int low_bits,
                           uint64_t *gout, unsigned int *gcount, int max_out) {
  const uint64_t prefix = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
  const uint64_t nprefix = 1ULL << (n - low_bits);
  if (prefix >= nprefix) return;

  uint64_t local[FES_LOCAL_CAP];
  uint64_t df[FES_MAX_VARS];
  int found = 0;
  fes_search_subcube(cst, lin, quad_tri, n, low_bits, prefix, local, &found,
                     FES_LOCAL_CAP, df);
  for (int i = 0; i < found && i < FES_LOCAL_CAP; ++i) {
    unsigned int idx = atomicAdd(gcount, 1u);
    if ((int)idx < max_out) gout[idx] = local[i];
  }
}

// Launch the kernel over the whole cube of a system and copy back the
// solutions. Returns the device count (>= sols.size() if it exceeded max_out),
// or a negative CUDA error code.
static long long run_on_device(const FesSystem &sys, int low_bits,
                               std::vector<uint64_t> &sols) {
  const int quad_len = sys.n * (sys.n + 1) / 2;
  const int max_out = 1 << 20;
  uint64_t *d_lin = nullptr, *d_quad = nullptr, *d_out = nullptr;
  unsigned int *d_count = nullptr;
#define TRY(call)                                                              \
  do {                                                                        \
    if ((call) != cudaSuccess) {                                             \
      fprintf(stderr, "CUDA error at %s:%d\n", __FILE__, __LINE__);          \
      return -1;                                                             \
    }                                                                        \
  } while (0)
  TRY(cudaMalloc(&d_lin, (sys.n ? sys.n : 1) * sizeof(uint64_t)));
  TRY(cudaMalloc(&d_quad, (quad_len ? quad_len : 1) * sizeof(uint64_t)));
  TRY(cudaMalloc(&d_out, max_out * sizeof(uint64_t)));
  TRY(cudaMalloc(&d_count, sizeof(unsigned int)));
  if (sys.n)
    TRY(cudaMemcpy(d_lin, sys.lin.data(), sys.n * sizeof(uint64_t), cudaMemcpyHostToDevice));
  if (quad_len)
    TRY(cudaMemcpy(d_quad, sys.quad_tri.data(), quad_len * sizeof(uint64_t), cudaMemcpyHostToDevice));
  TRY(cudaMemset(d_count, 0, sizeof(unsigned int)));

  const uint64_t nprefix = 1ULL << (sys.n - low_bits);
  const int threads = 256;
  const int blocks = (int)((nprefix + threads - 1) / threads);
  fes_kernel<<<blocks, threads>>>(sys.cst, d_lin, d_quad, sys.n, low_bits, d_out,
                                  d_count, max_out);
  TRY(cudaGetLastError());
  TRY(cudaDeviceSynchronize());
  unsigned int count = 0;
  TRY(cudaMemcpy(&count, d_count, sizeof(unsigned int), cudaMemcpyDeviceToHost));
  unsigned int listed = count < (unsigned)max_out ? count : (unsigned)max_out;
  sols.resize(listed);
  if (listed)
    TRY(cudaMemcpy(sols.data(), d_out, listed * sizeof(uint64_t), cudaMemcpyDeviceToHost));
  cudaFree(d_lin);
  cudaFree(d_quad);
  cudaFree(d_out);
  cudaFree(d_count);
#undef TRY
  return (long long)count;
}

int main(int argc, char **argv) {
  // Worker mode: `fes_cuda --in <system-file>` reads a system in the shared
  // contract, runs the kernel, and prints SOLUTIONS — this is what the Rust
  // binary drives. With no --in, it runs the compiled-in self-test.
  const char *in_path = nullptr;
  int low_bits = -1;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--in") && i + 1 < argc) in_path = argv[++i];
    else if (!strcmp(argv[i], "--low-bits") && i + 1 < argc) low_bits = atoi(argv[++i]);
  }

  if (in_path) {
    FesSystem sys;
    if (!fes_read_system(in_path, sys)) {
      fprintf(stderr, "fes_cuda: could not parse %s\n", in_path);
      return 2;
    }
    int lb = (low_bits >= 1 && low_bits <= sys.n) ? low_bits : sys.n;
    std::vector<uint64_t> sols;
    long long count = run_on_device(sys, lb, sols);
    if (count < 0) return 2;
    fes_write_solutions(sols);
    return 0;
  }

  // Self-test on the compiled-in system.h.
  FesSystem sys;
  sys.n = FES_N;
  sys.m = FES_M;
  sys.cst = FES_CONST;
  sys.lin.assign(FES_LIN, FES_LIN + FES_N);
  sys.quad_tri.assign(FES_QUAD_TRI, FES_QUAD_TRI + FES_N * (FES_N + 1) / 2);
  int lb = FES_N > 10 ? FES_N - 10 : FES_N;
  std::vector<uint64_t> sols;
  long long count = run_on_device(sys, lb, sols);
  if (count < 0) return 2;
  printf("solutions=%lld low_bits=%d\n", count, lb);
  int ok = ((int)count == FES_NUM_SOLUTIONS);
  printf(ok ? "MATCH: %d expected\n" : "MISMATCH: expected %d\n", FES_NUM_SOLUTIONS);
  return ok ? 0 : 1;
}
