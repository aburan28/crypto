#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../include/curveparams.h"
#include "../include/packed131.h"

#if !ECC_WALK_HALVING
#error "build with -DECC_WALK_HALVING=1"
#endif
#ifndef ECC_HALVING_POLY_STATE
#define ECC_HALVING_POLY_STATE 0
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
static bool same(P131 a, P131 b) {
    for (int i = 0; i < 5; ++i) if (a.v[i] != b.v[i]) return false;
    return true;
}

__global__ void halveProbe(const P131 *x, const P131 *y, P131 *hx, P131 *hy, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
#if ECC_HALVING_POLY_STATE
        eccPacked131::pointHalfPolynomial131(x[i], y[i], hx + i, hy + i);
#else
        eccPacked131::pointHalf131(x[i], y[i], hx + i, hy + i);
#endif
    }
}
__global__ __launch_bounds__(ECC_THREADS, ECC_MINBLOCKS)
void halveBench(P131 *x, P131 *y, int n, int steps) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    P131 px = x[i], py = y[i];
    for (int step = 0; step < steps; ++step) {
#if ECC_HALVING_POLY_STATE
        eccPacked131::pointHalfPolynomial131(px, py, &px, &py);
#else
        eccPacked131::pointHalf131(px, py, &px, &py);
#endif
    }
    x[i] = px; y[i] = py;
}

int main() {
    const R::Point p = R::make(R::fromLimbs(eccF131::PX), R::fromLimbs(eccF131::PY));
    const U192 ell = u192_from_dec(eccF131::ELL_DEC);
    const U192 inv2 = mod_inv(u192_from(2), ell);
    std::vector<P131> x, y, wantX, wantY;
    int badRoot = 0, badSqrt = 0, badHost = 0;
    int correctTrace[2] = {0, 0}, otherTrace[2] = {0, 0};
    for (int i = 1; i <= 512; ++i) {
        const R::Point q = R::scalarMul(p, u192_from(i));
        const R::Point h = R::half(q);
        const R::Point expectedHalf = R::scalarMul(q, inv2);
        badHost += !R::eq(h, expectedHalf);
        R::Elem lam = R::halfTrace(q.x);
        R::Elem cx = R::sigma(R::add(R::add(q.y, q.x), R::mul(lam, q.x)), 130);
        R::Point c0 = R::make(cx, R::mul(cx, R::add(lam, cx)));
        lam = R::add(lam, R::one());
        cx = R::add(cx, R::sigma(q.x, 130));
        R::Point c1 = R::make(cx, R::mul(cx, R::add(lam, cx)));
        const R::Point other = R::eq(c0, expectedHalf) ? c1 : c0;
        correctTrace[R::trace(R::mul(expectedHalf.y, R::inv(expectedHalf.x)))]++;
        otherTrace[R::trace(R::mul(other.y, R::inv(other.x)))]++;
        if (!R::eq(R::dbl(h), q) || R::trace(h.x)) {
            std::fprintf(stderr, "host point-halving oracle failed at %d\n", i);
            return 1;
        }
        const P131 qx = pack(q.x);
        const P131 root = eccPacked131::halvingQuadraticRoot131(qx);
        const P131 rootCheck = eccPacked131::add131(eccPacked131::sqr131(root), root);
        const P131 sqroot = eccPacked131::halvingSqrt131(qx);
        badRoot += !same(rootCheck, qx);
        badSqrt += !same(eccPacked131::sqr131(sqroot), qx);
        const P131 qy = pack(q.y), hx = pack(h.x), hy = pack(h.y);
#if ECC_HALVING_POLY_STATE
        x.push_back(eccPacked131::toPolynomial131(qx));
        y.push_back(eccPacked131::toPolynomial131(qy));
        wantX.push_back(eccPacked131::toPolynomial131(hx));
        wantY.push_back(eccPacked131::toPolynomial131(hy));
#else
        x.push_back(qx); y.push_back(qy);
        wantX.push_back(hx); wantY.push_back(hy);
#endif
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
    std::printf("point-halving CUDA probe: %zu subgroup points, root failures %d, sqrt failures %d, host subgroup failures %d, mismatched words %d\n",
                x.size(), badRoot, badSqrt, badHost, bad);
    std::printf("  trace(y/x): subgroup [%d,%d], other [%d,%d]\n",
                correctTrace[0], correctTrace[1], otherTrace[0], otherTrace[1]);
    std::printf("  polynomial state: %d\n", ECC_HALVING_POLY_STATE);
    std::printf("%s\n", bad ? "FAIL" : "PASS");
    if (badRoot || badSqrt || badHost || bad) return 1;

    const int benchN = 188 * 512, steps = 4096;
    std::vector<P131> bx(benchN), by(benchN);
    for (int i = 0; i < benchN; ++i) {
        bx[i] = x[i % x.size()]; by[i] = y[i % y.size()];
    }
    P131 *bdx, *bdy;
    checked(cudaMalloc(&bdx, size_t(benchN) * sizeof(P131)));
    checked(cudaMalloc(&bdy, size_t(benchN) * sizeof(P131)));
    checked(cudaMemcpy(bdx, bx.data(), size_t(benchN) * sizeof(P131), cudaMemcpyHostToDevice));
    checked(cudaMemcpy(bdy, by.data(), size_t(benchN) * sizeof(P131), cudaMemcpyHostToDevice));
    halveBench<<<(benchN + ECC_THREADS - 1) / ECC_THREADS, ECC_THREADS>>>(
        bdx, bdy, benchN, 4);
    checked(cudaDeviceSynchronize());
    cudaEvent_t begin, end; checked(cudaEventCreate(&begin)); checked(cudaEventCreate(&end));
    for (int rep = 0; rep < 3; ++rep) {
        checked(cudaEventRecord(begin));
        halveBench<<<(benchN + ECC_THREADS - 1) / ECC_THREADS, ECC_THREADS>>>(
            bdx, bdy, benchN, steps);
        checked(cudaEventRecord(end)); checked(cudaEventSynchronize(end));
        float ms = 0; checked(cudaEventElapsedTime(&ms, begin, end));
        std::printf("  halving raw %.3f B/s (%g ms)\n",
                    double(benchN) * steps / (double(ms) * 1e6), ms);
    }
    cudaFree(bdx); cudaFree(bdy); cudaEventDestroy(begin); cudaEventDestroy(end);
    return 0;
}
