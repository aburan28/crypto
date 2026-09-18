#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../include/curveparams.h"
#include "../include/packed131.h"

#if !ECC_WALK_HALVING
#error "build with -DECC_WALK_HALVING=1"
#endif

using eccPacked131::P131;
using R = Ref<CfgF131>;

static void checked(cudaError_t e) {
    if (e != cudaSuccess) {
        std::fprintf(stderr, "CUDA point-halving test: %s\n", cudaGetErrorString(e));
        std::exit(1);
    }
}

static P131 pack(const R::Elem &a) {
    return P131{{uint32_t(a.v[0]), uint32_t(a.v[0] >> 32),
                 uint32_t(a.v[1]), uint32_t(a.v[1] >> 32), uint32_t(a.v[2])}};
}

__global__ void halveProbe(const P131 *x, const P131 *y, P131 *hx, P131 *hy, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) eccPacked131::pointHalf131(x[i], y[i], hx + i, hy + i);
}

int main() {
    const R::Point p = R::make(R::fromLimbs(eccF131::PX), R::fromLimbs(eccF131::PY));
    std::vector<P131> x, y, wantX, wantY;
    for (int i = 1; i <= 512; ++i) {
        const R::Point q = R::scalarMul(p, u192_from(i));
        const R::Point h = R::half(q);
        if (!R::eq(R::dbl(h), q) || R::trace(h.x)) {
            std::fprintf(stderr, "host point-halving oracle failed at %d\n", i);
            return 1;
        }
        x.push_back(pack(q.x)); y.push_back(pack(q.y));
        wantX.push_back(pack(h.x)); wantY.push_back(pack(h.y));
    }
    P131 *dx, *dy, *dhx, *dhy;
    const size_t bytes = x.size() * sizeof(P131);
    checked(cudaMalloc(&dx, bytes)); checked(cudaMalloc(&dy, bytes));
    checked(cudaMalloc(&dhx, bytes)); checked(cudaMalloc(&dhy, bytes));
    checked(cudaMemcpy(dx, x.data(), bytes, cudaMemcpyHostToDevice));
    checked(cudaMemcpy(dy, y.data(), bytes, cudaMemcpyHostToDevice));
    halveProbe<<<2, 256>>>(dx, dy, dhx, dhy, int(x.size()));
    checked(cudaGetLastError());
    std::vector<P131> gotX(x.size()), gotY(x.size());
    checked(cudaMemcpy(gotX.data(), dhx, bytes, cudaMemcpyDeviceToHost));
    checked(cudaMemcpy(gotY.data(), dhy, bytes, cudaMemcpyDeviceToHost));
    int bad = 0;
    for (size_t i = 0; i < x.size(); ++i)
        for (int w = 0; w < 5; ++w)
            bad += gotX[i].v[w] != wantX[i].v[w] || gotY[i].v[w] != wantY[i].v[w];
    cudaFree(dx); cudaFree(dy); cudaFree(dhx); cudaFree(dhy);
    std::printf("point-halving CUDA probe: %zu subgroup points, mismatched words %d\n",
                x.size(), bad);
    std::printf("%s\n", bad ? "FAIL" : "PASS");
    return bad ? 1 : 0;
}
