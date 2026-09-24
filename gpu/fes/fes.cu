// gpu/fes/fes.cu -- CUDA launcher for the Gray-code FES search.  One thread per
// sub-cube (a fixed assignment of the top n-low_bits variables); each thread
// runs the same fes_search_subcube the host test verifies.  Build with
// `make cuda` (nvcc); this file is the device-specific launch only, the
// algorithm lives in fes.cuh and is checked by `make test` without a GPU.
#include <cstdint>
#include <cstdio>
#include <vector>

#include "fes.cuh"
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

#define CUDA_OK(call)                                                           \
  do {                                                                          \
    cudaError_t e = (call);                                                     \
    if (e != cudaSuccess) {                                                     \
      fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e),        \
              __FILE__, __LINE__);                                              \
      return 2;                                                                 \
    }                                                                           \
  } while (0)

int main(int argc, char **argv) {
  int low_bits = FES_N > 10 ? FES_N - 10 : FES_N; // ~1024 threads by default
  if (argc > 1) low_bits = atoi(argv[1]);
  if (low_bits < 1 || low_bits > FES_N) low_bits = FES_N;

  const int quad_len = FES_N * (FES_N + 1) / 2;
  const int max_out = 1 << 16;

  uint64_t *d_lin = nullptr, *d_quad = nullptr, *d_out = nullptr;
  unsigned int *d_count = nullptr;
  CUDA_OK(cudaMalloc(&d_lin, FES_N * sizeof(uint64_t)));
  CUDA_OK(cudaMalloc(&d_quad, quad_len * sizeof(uint64_t)));
  CUDA_OK(cudaMalloc(&d_out, max_out * sizeof(uint64_t)));
  CUDA_OK(cudaMalloc(&d_count, sizeof(unsigned int)));
  CUDA_OK(cudaMemcpy(d_lin, FES_LIN, FES_N * sizeof(uint64_t), cudaMemcpyHostToDevice));
  CUDA_OK(cudaMemcpy(d_quad, FES_QUAD_TRI, quad_len * sizeof(uint64_t), cudaMemcpyHostToDevice));
  CUDA_OK(cudaMemset(d_count, 0, sizeof(unsigned int)));

  const uint64_t nprefix = 1ULL << (FES_N - low_bits);
  const int threads = 256;
  const int blocks = (int)((nprefix + threads - 1) / threads);
  fes_kernel<<<blocks, threads>>>(FES_CONST, d_lin, d_quad, FES_N, low_bits,
                                  d_out, d_count, max_out);
  CUDA_OK(cudaGetLastError());
  CUDA_OK(cudaDeviceSynchronize());

  unsigned int count = 0;
  CUDA_OK(cudaMemcpy(&count, d_count, sizeof(unsigned int), cudaMemcpyDeviceToHost));
  std::vector<uint64_t> out(count < (unsigned)max_out ? count : max_out);
  if (!out.empty())
    CUDA_OK(cudaMemcpy(out.data(), d_out, out.size() * sizeof(uint64_t), cudaMemcpyDeviceToHost));

  printf("threads=%llu low_bits=%d solutions=%u\n", (unsigned long long)nprefix,
         low_bits, count);
  // Cross-check against the generator's known set.
  int ok = ((int)count == FES_NUM_SOLUTIONS);
  printf(ok ? "MATCH: %d expected\n" : "MISMATCH: expected %d\n", FES_NUM_SOLUTIONS);

  cudaFree(d_lin);
  cudaFree(d_quad);
  cudaFree(d_out);
  cudaFree(d_count);
  return ok ? 0 : 1;
}
