#pragma once
#include "macaulay.cuh"
#ifndef MAC_THREADS
#define MAC_THREADS 256
#endif

/* One block per matrix. */
__global__ void rref_batch_kernel(uint32_t* a, int rows, int cols, int* piv,
                                  int* rank) {
    extern __shared__ uint32_t smem[];
    const int inst = blockIdx.x;
    uint32_t* mine = a + (size_t)inst * rows * cols;
    int* mypiv = piv + (size_t)inst * rows;
    int r = mac_rref_block<MAC_THREADS>(mine, rows, cols, mypiv, smem);
    if (threadIdx.x == 0) rank[inst] = r;
}

/* Montgomery conversion, flat over the whole batch. */
__global__ void to_mont_kernel(uint32_t* a, size_t n) {
    for (size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x; i < n;
         i += (size_t)gridDim.x * blockDim.x) {
        a[i] = fp_to_mont(a[i]);
    }
}
__global__ void from_mont_kernel(uint32_t* a, size_t n) {
    for (size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x; i < n;
         i += (size_t)gridDim.x * blockDim.x) {
        a[i] = fp_from_mont(a[i]);
    }
}

