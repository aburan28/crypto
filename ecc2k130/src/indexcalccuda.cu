// GPU pair-enumeration oracle for three-summand index calculus on ECC2K-130.
//
// The algorithm is the null object of RESEARCH_ECC2K130_RR_SOLVER_PANEL.md §7:
// for a factor base F, scan unordered pairs, add, and test whether the remainder
// against a target has Hamming weight within the base.  It cannot beat the
// product-law floor m·2^131; this binary measures that oracle on one G7e
// RTX PRO 6000 in the same packed type-II ONB arithmetic the rho client uses.
//
//   make indexcalc-cuda
//   ./build/indexcalc-cuda --self-test
//   ./build/indexcalc-cuda --weight 2 --planted 32 --search-generator --cpu --table --json out.json

#include <cuda_runtime.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "../include/curveparams.h"
#include "../include/packed131.h"

using eccPacked131::P131;
using R = Ref<CfgF131>;

static void die(const char *msg) {
    std::fprintf(stderr, "indexcalc-cuda: %s\n", msg);
    std::exit(1);
}
static void checked(cudaError_t e) {
    if (e != cudaSuccess) {
        std::fprintf(stderr, "indexcalc-cuda: %s\n", cudaGetErrorString(e));
        std::exit(1);
    }
}

ECC_HD bool isZero131(P131 a) {
    return !(a.v[0] | a.v[1] | a.v[2] | a.v[3] | a.v[4]);
}
ECC_HD int hw131(P131 a) {
#ifdef __CUDA_ARCH__
    return __popc(a.v[0]) + __popc(a.v[1]) + __popc(a.v[2]) + __popc(a.v[3]) + __popc(a.v[4]);
#else
    return __builtin_popcount(a.v[0]) + __builtin_popcount(a.v[1]) + __builtin_popcount(a.v[2])
         + __builtin_popcount(a.v[3]) + __builtin_popcount(a.v[4]);
#endif
}
ECC_HD bool same131(P131 a, P131 b) {
    return !((a.v[0] ^ b.v[0]) | (a.v[1] ^ b.v[1]) | (a.v[2] ^ b.v[2])
             | (a.v[3] ^ b.v[3]) | (a.v[4] ^ b.v[4]));
}
ECC_HD bool affineAdd(P131 x1, P131 y1, P131 x2, P131 y2, P131 *x3, P131 *y3) {
    P131 d = eccPacked131::add131(x1, x2);
    if (isZero131(d)) return false;
    P131 lam = eccPacked131::mul131(eccPacked131::add131(y1, y2), eccPacked131::inv131(d));
    P131 x = eccPacked131::add131(eccPacked131::add131(eccPacked131::sqr131(lam), lam), d);
    P131 y = eccPacked131::add131(eccPacked131::add131(eccPacked131::mul131(lam, eccPacked131::add131(x1, x)), x), y1);
    *x3 = x;
    *y3 = y;
    return true;
}

struct Hit {
    int i, j, k;
    P131 xt, yt;
};

__global__ __launch_bounds__(128, 1)
void searchRows(const P131 *X, const P131 *Y, int B, P131 rx, P131 ry, int maxW,
                Hit *hits, int maxHits, unsigned long long *nHits,
                unsigned long long *nPairs) {
    const int i = (int)blockIdx.x;
    if (i >= B) return;
    P131 xi = X[i], yi = Y[i];
    unsigned long long local = 0;
    for (int j = i + 1 + (int)threadIdx.x; j < B; j += (int)blockDim.x) {
        local++;
        P131 sx, sy;
        if (!affineAdd(xi, yi, X[j], Y[j], &sx, &sy)) continue;
        P131 tx, ty;
        P131 nsy = eccPacked131::add131(sx, sy);
        if (!affineAdd(rx, ry, sx, nsy, &tx, &ty)) continue;
        if (hw131(tx) > maxW) continue;
        if (same131(tx, xi) || same131(tx, X[j])) continue;
        unsigned long long slot = atomicAdd(nHits, 1ULL);
        if ((int)slot < maxHits) {
            hits[slot].i = i;
            hits[slot].j = j;
            hits[slot].k = -1;
            hits[slot].xt = tx;
            hits[slot].yt = ty;
        }
    }
    atomicAdd(nPairs, local);
}

__global__ __launch_bounds__(128, 1)
void addBench(P131 *x, P131 *y, int n, int steps) {
    const int i = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    if (i >= n - 1) return;
    P131 ax = x[i], ay = y[i], bx = x[i + 1], by = y[i + 1];
    for (int s = 0; s < steps; s++) {
        P131 nx, ny;
        if (!affineAdd(ax, ay, bx, by, &nx, &ny)) break;
        ax = nx;
        ay = ny;
    }
    x[i] = ax;
    y[i] = ay;
}

struct PairSum {
    P131 x;
    int i, j;
};

ECC_HD int cmp131(P131 a, P131 b) {
    for (int w = 4; w >= 0; --w) {
        if (a.v[w] < b.v[w]) return -1;
        if (a.v[w] > b.v[w]) return 1;
    }
    return 0;
}

__global__ __launch_bounds__(128, 1)
void fillPairs(const P131 *X, const P131 *Y, int B, PairSum *tab,
               unsigned long long *nStored) {
    const int i = (int)blockIdx.x;
    if (i >= B) return;
    P131 xi = X[i], yi = Y[i];
    for (int j = i + 1 + (int)threadIdx.x; j < B; j += (int)blockDim.x) {
        P131 sx, sy;
        if (!affineAdd(xi, yi, X[j], Y[j], &sx, &sy)) continue;
        unsigned long long slot = atomicAdd(nStored, 1ULL);
        tab[slot].x = sx;
        tab[slot].i = i;
        tab[slot].j = j;
    }
}

ECC_HD int lowerBoundX(const PairSum *tab, int n, P131 x) {
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (cmp131(tab[mid].x, x) < 0) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

__global__ __launch_bounds__(128, 1)
void probeTable(const P131 *X, const P131 *Y, int B, const PairSum *tab, int nTab,
                P131 rx, P131 ry, Hit *hits, int maxHits, unsigned long long *nHits,
                unsigned long long *nAdds, int steps) {
    const int nWork = B * steps;
    unsigned long long localAdds = 0;
    for (int idx = (int)(blockIdx.x * blockDim.x + threadIdx.x);
         idx < nWork;
         idx += (int)(blockDim.x * gridDim.x)) {
        const int k = idx % B;
        const int s = idx / B;
        P131 tx, ty;
        localAdds++;
        P131 npy = eccPacked131::add131(X[k], Y[k]);
        if (!affineAdd(rx, ry, X[k], npy, &tx, &ty)) continue;
        int p = lowerBoundX(tab, nTab, tx);
        while (p < nTab && same131(tab[p].x, tx)) {
            int i = tab[p].i, j = tab[p].j;
            // Same-x as i or j is P or -P (char-2); those triples cancel.
            if (i != k && j != k && !same131(X[k], X[i]) && !same131(X[k], X[j])) {
                P131 sx, sy;
                if (affineAdd(X[i], Y[i], X[j], Y[j], &sx, &sy)
                    && same131(sx, tx) && same131(sy, ty)) {
                    if (s == 0) {
                        unsigned long long slot = atomicAdd(nHits, 1ULL);
                        if ((int)slot < maxHits) {
                            hits[slot].i = i;
                            hits[slot].j = j;
                            hits[slot].k = k;
                            hits[slot].xt = tx;
                            hits[slot].yt = ty;
                        }
                    }
                    break;
                }
            }
            p++;
        }
    }
    atomicAdd(nAdds, localAdds);
}

static P131 pack(const R::Elem &a) {
    return P131{{uint32_t(a.v[0]), uint32_t(a.v[0] >> 32),
                 uint32_t(a.v[1]), uint32_t(a.v[1] >> 32), uint32_t(a.v[2])}};
}
static R::Elem unpack(P131 a) {
    unsigned long long v[3] = {
        (unsigned long long)a.v[0] | ((unsigned long long)a.v[1] << 32),
        (unsigned long long)a.v[2] | ((unsigned long long)a.v[3] << 32),
        a.v[4]
    };
    return R::fromLimbs(v);
}
static R::Point unpackPt(P131 x, P131 y) {
    return R::make(unpack(x), unpack(y));
}

static bool pointFromX(const R::Elem &x, R::Point *out) {
    if (R::isZero(x)) return false;
    R::Elem xx = R::sqr(x);
    if (R::isZero(xx)) return false;
    R::Elem c = R::add(x, R::inv(xx));
    if (R::trace(c)) return false;
    R::Point p = R::make(x, R::mul(x, R::halfTrace(c)));
    if (!R::onCurve(p)) return false;
    *out = p;
    return true;
}

static void recBits(int start, int left, R::Elem cur, std::vector<R::Elem> *out) {
    if (left == 0) {
        out->push_back(cur);
        return;
    }
    for (int i = start; i <= R::M - left; ++i) {
        R::Elem nxt = cur;
        R::setBit(nxt, i);
        recBits(i + 1, left - 1, nxt, out);
    }
}

struct Base {
    std::vector<P131> x, y;
    std::vector<R::Point> pts;
    int weight = 0;
    int abscissae = 0;
};

static Base buildBase(int weight, bool oddOrder) {
    // oddOrder keeps the necessary trace-zero abscissae (even ONB weight).
    // It is not a [r]P certificate: that scalar mul is out of scope for the
    // oracle-throughput measurement and is what made a weight-2 build stall.
    Base b;
    b.weight = weight;
    std::vector<R::Elem> xs;
    for (int w = 1; w <= weight; ++w) recBits(0, w, R::zero(), &xs);
    for (const R::Elem &x : xs) {
        if (oddOrder && (R::weight(x) & 1)) continue;
        R::Point p;
        if (!pointFromX(x, &p)) continue;
        b.abscissae++;
        b.pts.push_back(p);
        b.x.push_back(pack(p.x));
        b.y.push_back(pack(p.y));
        R::Point n = R::neg(p);
        b.pts.push_back(n);
        b.x.push_back(pack(n.x));
        b.y.push_back(pack(n.y));
    }
    return b;
}

struct SearchResult {
    unsigned long long hits = 0;
    unsigned long long pairs = 0;
    double seconds = 0;
    std::vector<Hit> stored;
    bool verified = true;
};

static SearchResult runSearch(const Base &b, P131 rx, P131 ry, int maxHits, bool onGpu) {
    SearchResult out;
    const int B = (int)b.x.size();
    if (B < 3) die("factor base too small");
    if (onGpu) {
        P131 *dx, *dy;
        Hit *dh;
        unsigned long long *dHits, *dPairs;
        checked(cudaMalloc(&dx, b.x.size() * sizeof(P131)));
        checked(cudaMalloc(&dy, b.y.size() * sizeof(P131)));
        checked(cudaMalloc(&dh, maxHits * sizeof(Hit)));
        checked(cudaMalloc(&dHits, sizeof(unsigned long long)));
        checked(cudaMalloc(&dPairs, sizeof(unsigned long long)));
        checked(cudaMemcpy(dx, b.x.data(), b.x.size() * sizeof(P131), cudaMemcpyHostToDevice));
        checked(cudaMemcpy(dy, b.y.data(), b.y.size() * sizeof(P131), cudaMemcpyHostToDevice));
        checked(cudaMemset(dHits, 0, sizeof(unsigned long long)));
        checked(cudaMemset(dPairs, 0, sizeof(unsigned long long)));
        cudaEvent_t start, stop;
        checked(cudaEventCreate(&start));
        checked(cudaEventCreate(&stop));
        checked(cudaDeviceSynchronize());
        checked(cudaEventRecord(start));
        searchRows<<<B, 128>>>(dx, dy, B, rx, ry, b.weight, dh, maxHits, dHits, dPairs);
        checked(cudaEventRecord(stop));
        checked(cudaEventSynchronize(stop));
        float ms = 0;
        checked(cudaEventElapsedTime(&ms, start, stop));
        out.seconds = ms / 1000.0;
        cudaEventDestroy(start);
        cudaEventDestroy(stop);
        checked(cudaMemcpy(&out.hits, dHits, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
        checked(cudaMemcpy(&out.pairs, dPairs, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
        int keep = (int)std::min(out.hits, (unsigned long long)maxHits);
        out.stored.resize(keep);
        if (keep) checked(cudaMemcpy(out.stored.data(), dh, keep * sizeof(Hit), cudaMemcpyDeviceToHost));
        cudaFree(dx); cudaFree(dy); cudaFree(dh); cudaFree(dHits); cudaFree(dPairs);
    } else {
        auto t0 = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 4)
        for (int i = 0; i < B; ++i) {
            unsigned long long localPairs = 0;
            for (int j = i + 1; j < B; ++j) {
                localPairs++;
                P131 sx, sy;
                if (!affineAdd(b.x[i], b.y[i], b.x[j], b.y[j], &sx, &sy)) continue;
                P131 tx, ty;
                if (!affineAdd(rx, ry, sx, eccPacked131::add131(sx, sy), &tx, &ty)) continue;
                if (hw131(tx) > b.weight) continue;
                if (same131(tx, b.x[i]) || same131(tx, b.x[j])) continue;
                Hit h;
                h.i = i; h.j = j; h.k = -1; h.xt = tx; h.yt = ty;
#pragma omp critical
                {
                    out.hits++;
                    if ((int)out.stored.size() < maxHits) out.stored.push_back(h);
                }
            }
#pragma omp atomic
            out.pairs += localPairs;
        }
        out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    R::Point target = unpackPt(rx, ry);
    for (const Hit &h : out.stored) {
        if (h.i < 0 || h.j <= h.i || h.j >= B) { out.verified = false; continue; }
        R::Point s = R::addPt(b.pts[h.i], b.pts[h.j]);
        R::Point t = unpackPt(h.xt, h.yt);
        if (h.k < 0) {
            R::Point sum = R::addPt(s, t);
            if (!R::eq(sum, target) || !R::onCurve(t)) out.verified = false;
        } else {
            if (h.k >= B || h.k == h.i || h.k == h.j) { out.verified = false; continue; }
            if (!R::eq(s, t) || !R::onCurve(t)) out.verified = false;
            if (!R::eq(R::addPt(s, b.pts[h.k]), target)) out.verified = false;
        }
    }
    return out;
}

static bool lessPair(const PairSum &a, const PairSum &b) {
    int c = cmp131(a.x, b.x);
    if (c) return c < 0;
    if (a.i != b.i) return a.i < b.i;
    return a.j < b.j;
}

struct PairTable {
    PairSum *dx = 0;
    P131 *dX = 0, *dY = 0;
    int B = 0;
    int nStored = 0;
    size_t bytes = 0;
    double buildSeconds = 0;
    double sortSeconds = 0;
};

static void freeTable(PairTable *t) {
    if (t->dx) cudaFree(t->dx);
    if (t->dX) cudaFree(t->dX);
    if (t->dY) cudaFree(t->dY);
    t->dx = 0; t->dX = 0; t->dY = 0;
}

static PairTable buildTable(const Base &b) {
    PairTable t;
    t.B = (int)b.x.size();
    unsigned long long cap = (unsigned long long)t.B * (unsigned long long)(t.B - 1) / 2;
    t.bytes = (size_t)cap * sizeof(PairSum);
    checked(cudaMalloc(&t.dX, b.x.size() * sizeof(P131)));
    checked(cudaMalloc(&t.dY, b.y.size() * sizeof(P131)));
    checked(cudaMalloc(&t.dx, t.bytes));
    checked(cudaMemcpy(t.dX, b.x.data(), b.x.size() * sizeof(P131), cudaMemcpyHostToDevice));
    checked(cudaMemcpy(t.dY, b.y.data(), b.y.size() * sizeof(P131), cudaMemcpyHostToDevice));
    unsigned long long *dCount;
    checked(cudaMalloc(&dCount, sizeof(unsigned long long)));
    checked(cudaMemset(dCount, 0, sizeof(unsigned long long)));
    cudaEvent_t start, stop;
    checked(cudaEventCreate(&start));
    checked(cudaEventCreate(&stop));
    checked(cudaDeviceSynchronize());
    checked(cudaEventRecord(start));
    fillPairs<<<t.B, 128>>>(t.dX, t.dY, t.B, t.dx, dCount);
    checked(cudaEventRecord(stop));
    checked(cudaEventSynchronize(stop));
    float ms = 0;
    checked(cudaEventElapsedTime(&ms, start, stop));
    t.buildSeconds = ms / 1000.0;
    unsigned long long stored = 0;
    checked(cudaMemcpy(&stored, dCount, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
    cudaFree(dCount);
    if (stored == 0 || stored > cap) die("pair table fill produced an empty or overrun buffer");
    t.nStored = (int)stored;
    std::vector<PairSum> host(stored);
    checked(cudaMemcpy(host.data(), t.dx, stored * sizeof(PairSum), cudaMemcpyDeviceToHost));
    auto t0 = std::chrono::steady_clock::now();
    std::sort(host.begin(), host.end(), lessPair);
    t.sortSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    checked(cudaMemcpy(t.dx, host.data(), stored * sizeof(PairSum), cudaMemcpyHostToDevice));
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    return t;
}

static SearchResult runTableProbe(const Base &b, const PairTable &tab, P131 rx, P131 ry,
                                  int maxHits, int steps) {
    SearchResult out;
    Hit *dh;
    unsigned long long *dHits, *dAdds;
    checked(cudaMalloc(&dh, maxHits * sizeof(Hit)));
    checked(cudaMalloc(&dHits, sizeof(unsigned long long)));
    checked(cudaMalloc(&dAdds, sizeof(unsigned long long)));
    checked(cudaMemset(dHits, 0, sizeof(unsigned long long)));
    checked(cudaMemset(dAdds, 0, sizeof(unsigned long long)));
    int threads = 128;
    int nWork = tab.B * steps;
    int blocks = (nWork + threads - 1) / threads;
    if (blocks < 1) blocks = 1;
    cudaEvent_t start, stop;
    checked(cudaEventCreate(&start));
    checked(cudaEventCreate(&stop));
    checked(cudaDeviceSynchronize());
    checked(cudaEventRecord(start));
    probeTable<<<blocks, threads>>>(tab.dX, tab.dY, tab.B, tab.dx, tab.nStored,
                                    rx, ry, dh, maxHits, dHits, dAdds, steps);
    checked(cudaEventRecord(stop));
    checked(cudaEventSynchronize(stop));
    float ms = 0;
    checked(cudaEventElapsedTime(&ms, start, stop));
    out.seconds = ms / 1000.0;
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    checked(cudaMemcpy(&out.hits, dHits, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
    checked(cudaMemcpy(&out.pairs, dAdds, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
    int keep = (int)std::min(out.hits, (unsigned long long)maxHits);
    out.stored.resize(keep);
    if (keep) checked(cudaMemcpy(out.stored.data(), dh, keep * sizeof(Hit), cudaMemcpyDeviceToHost));
    cudaFree(dh); cudaFree(dHits); cudaFree(dAdds);
    R::Point target = unpackPt(rx, ry);
    const int B = tab.B;
    for (const Hit &h : out.stored) {
        if (h.i < 0 || h.j <= h.i || h.j >= B || h.k < 0 || h.k >= B) {
            out.verified = false;
            continue;
        }
        if (h.k == h.i || h.k == h.j) { out.verified = false; continue; }
        R::Point s = R::addPt(b.pts[h.i], b.pts[h.j]);
        R::Point t = unpackPt(h.xt, h.yt);
        if (!R::eq(s, t) || !R::onCurve(t)) out.verified = false;
        if (!R::eq(R::addPt(s, b.pts[h.k]), target)) out.verified = false;
    }
    return out;
}

static int selfTest() {
    int bad = 0;
    const R::Point g = R::make(R::fromLimbs(eccF131::PX), R::fromLimbs(eccF131::PY));
    for (int i = 1; i <= 64; ++i) {
        R::Point a = R::scalarMul(g, u192_from(i * 3 + 1));
        R::Point b = R::scalarMul(g, u192_from(i * 5 + 2));
        R::Point want = R::addPt(a, b);
        P131 x3, y3;
        if (!affineAdd(pack(a.x), pack(a.y), pack(b.x), pack(b.y), &x3, &y3)) {
            bad++;
            continue;
        }
        R::Point got = unpackPt(x3, y3);
        bad += !R::eq(got, want);
    }
    if (bad) {
        std::fprintf(stderr, "affine-add mismatches against Ref: %d\n", bad);
        return 1;
    }
    std::fprintf(stderr, "self-test: packed add vs Ref, then one planted triple\n");
    Base b = buildBase(2, true);
    if (b.pts.size() < 8) {
        std::fprintf(stderr, "weight-2 even-trace base too small: %zu\n", b.pts.size());
        return 1;
    }
    R::Point r = R::addPt(R::addPt(b.pts[0], b.pts[2]), b.pts[4]);
    SearchResult gpu = runSearch(b, pack(r.x), pack(r.y), 64, true);
    if (!gpu.verified || gpu.hits == 0) {
        std::fprintf(stderr, "planted GPU search failed: hits=%llu verified=%d\n",
                     gpu.hits, gpu.verified);
        return 1;
    }
    std::printf("self-test ok  base=%zu abscissae=%d planted_hits=%llu pairs=%llu\n",
                b.pts.size(), b.abscissae, gpu.hits, gpu.pairs);
    PairTable tab = buildTable(b);
    SearchResult tableHit = runTableProbe(b, tab, pack(r.x), pack(r.y), 64, 1);
    if (!tableHit.verified || tableHit.hits == 0) {
        std::fprintf(stderr, "planted GPU pair-table probe failed: hits=%llu verified=%d\n",
                     tableHit.hits, tableHit.verified);
        freeTable(&tab);
        return 1;
    }
    std::printf("self-test table ok  stored=%d planted_hits=%llu adds=%llu\n",
                tab.nStored, tableHit.hits, tableHit.pairs);
    freeTable(&tab);
    return 0;
}

struct TableMeas {
    const PairTable *tab = 0;
    int plantedTried = 0;
    int plantedFound = 0;
    bool plantedOk = true;
    const SearchResult *probe = 0;
    const SearchResult *generator = 0;
};

static void writeJson(FILE *f, const char *gpuName, int sm, const Base &b,
                      int plantedTried, int plantedFound, bool plantedOk,
                      const SearchResult *naturalGpu, const SearchResult *naturalCpu,
                      double benchAdds, double benchSeconds, int benchSteps,
                      const TableMeas *table) {
    std::fprintf(f, "{\n");
    std::fprintf(f, "  \"schema\": \"ecc2k130_g7e_index_calculus/v1\",\n");
    std::fprintf(f, "  \"curve\": \"K_0 / F_2^131\",\n");
    std::fprintf(f, "  \"oracle_class\": \"enumerate_pairs\",\n");
    std::fprintf(f, "  \"class\": \"engineering\",\n");
    std::fprintf(f, "  \"gpu_name\": \"%s\",\n", gpuName);
    std::fprintf(f, "  \"sm\": %d,\n", sm);
    std::fprintf(f, "  \"weight\": %d,\n", b.weight);
    std::fprintf(f, "  \"factor_base_points\": %zu,\n", b.pts.size());
    std::fprintf(f, "  \"factor_base_abscissae\": %d,\n", b.abscissae);
    std::fprintf(f, "  \"products_per_affine_add\": 10,\n");
    std::fprintf(f, "  \"products_per_pair\": 20,\n");
    std::fprintf(f, "  \"planted_tried\": %d,\n", plantedTried);
    std::fprintf(f, "  \"planted_found\": %d,\n", plantedFound);
    std::fprintf(f, "  \"planted_verified\": %s,\n", plantedOk ? "true" : "false");
    if (naturalGpu) {
        std::fprintf(f, "  \"generator_gpu_hits\": %llu,\n", naturalGpu->hits);
        std::fprintf(f, "  \"generator_gpu_pairs\": %llu,\n", naturalGpu->pairs);
        std::fprintf(f, "  \"generator_gpu_seconds\": %.9f,\n", naturalGpu->seconds);
        std::fprintf(f, "  \"generator_gpu_verified\": %s,\n", naturalGpu->verified ? "true" : "false");
        double rate = naturalGpu->seconds > 0 ? naturalGpu->pairs / naturalGpu->seconds : 0;
        std::fprintf(f, "  \"generator_gpu_pairs_per_second\": %.6f,\n", rate);
    }
    if (naturalCpu) {
        std::fprintf(f, "  \"generator_cpu_hits\": %llu,\n", naturalCpu->hits);
        std::fprintf(f, "  \"generator_cpu_pairs\": %llu,\n", naturalCpu->pairs);
        std::fprintf(f, "  \"generator_cpu_seconds\": %.9f,\n", naturalCpu->seconds);
        std::fprintf(f, "  \"generator_cpu_verified\": %s,\n", naturalCpu->verified ? "true" : "false");
        double rate = naturalCpu->seconds > 0 ? naturalCpu->pairs / naturalCpu->seconds : 0;
        std::fprintf(f, "  \"generator_cpu_pairs_per_second\": %.6f,\n", rate);
    }
    std::fprintf(f, "  \"bench_affine_adds\": %.0f,\n", benchAdds);
    std::fprintf(f, "  \"bench_seconds\": %.9f,\n", benchSeconds);
    std::fprintf(f, "  \"bench_steps\": %d,\n", benchSteps);
    double benchRate = benchSeconds > 0 ? benchAdds / benchSeconds : 0;
    if (table && table->tab) {
        std::fprintf(f, "  \"bench_affine_adds_per_second\": %.6f,\n", benchRate);
        const PairTable *tab = table->tab;
        std::fprintf(f, "  \"pair_table\": {\n");
        std::fprintf(f, "    \"class\": \"engineering\",\n");
        std::fprintf(f, "    \"pairs_stored\": %d,\n", tab->nStored);
        std::fprintf(f, "    \"bytes\": %zu,\n", tab->bytes);
        std::fprintf(f, "    \"build_seconds\": %.9f,\n", tab->buildSeconds);
        std::fprintf(f, "    \"sort_seconds\": %.9f,\n", tab->sortSeconds);
        std::fprintf(f, "    \"planted_tried\": %d,\n", table->plantedTried);
        std::fprintf(f, "    \"planted_found\": %d,\n", table->plantedFound);
        std::fprintf(f, "    \"planted_verified\": %s,\n", table->plantedOk ? "true" : "false");
        if (table->probe) {
            double pr = table->probe->seconds > 0
                ? table->probe->pairs / table->probe->seconds : 0;
            std::fprintf(f, "    \"probe_adds\": %llu,\n", table->probe->pairs);
            std::fprintf(f, "    \"probe_seconds\": %.9f,\n", table->probe->seconds);
            std::fprintf(f, "    \"probe_adds_per_second\": %.6f,\n", pr);
        }
        if (table->generator) {
            std::fprintf(f, "    \"generator_hits\": %llu,\n", table->generator->hits);
            std::fprintf(f, "    \"generator_adds\": %llu,\n", table->generator->pairs);
            std::fprintf(f, "    \"generator_seconds\": %.9f,\n", table->generator->seconds);
            std::fprintf(f, "    \"generator_verified\": %s\n",
                         table->generator->verified ? "true" : "false");
        } else {
            std::fprintf(f, "    \"generator_hits\": null\n");
        }
        std::fprintf(f, "  }\n");
    } else {
        std::fprintf(f, "  \"bench_affine_adds_per_second\": %.6f\n", benchRate);
    }
    std::fprintf(f, "}\n");
}

static int argi(int argc, char **argv, const char *name, int fallback) {
    for (int i = 1; i < argc - 1; ++i)
        if (!std::strcmp(argv[i], name)) return std::atoi(argv[i + 1]);
    return fallback;
}
static const char *args(int argc, char **argv, const char *name) {
    for (int i = 1; i < argc - 1; ++i)
        if (!std::strcmp(argv[i], name)) return argv[i + 1];
    return 0;
}
static bool argf(int argc, char **argv, const char *name) {
    for (int i = 1; i < argc; ++i)
        if (!std::strcmp(argv[i], name)) return true;
    return false;
}

int main(int argc, char **argv) {
    if (argf(argc, argv, "--self-test")) return selfTest();
    const int weight = argi(argc, argv, "--weight", 2);
    const int plantedN = argi(argc, argv, "--planted", 8);
    const int benchPoints = argi(argc, argv, "--bench-points", 2048);
    const int benchSteps = argi(argc, argv, "--bench-steps", 32);
    const bool wantCpu = argf(argc, argv, "--cpu");
    const bool wantGen = argf(argc, argv, "--search-generator");
    const bool wantTable = argf(argc, argv, "--table");
    const int tableProbeSteps = argi(argc, argv, "--table-probe-steps", 4096);
    const char *jsonPath = args(argc, argv, "--json");
    if (weight < 1 || weight > 3) die("--weight must be 1, 2 or 3");

    cudaDeviceProp prop;
    checked(cudaGetDeviceProperties(&prop, 0));
    std::printf("gpu %s sm_%d %.0f MiB\n", prop.name, prop.major * 10 + prop.minor,
                prop.totalGlobalMem / (1024.0 * 1024.0));

    std::fprintf(stderr, "building weight-%d ONB Hamming base\n", weight);
    Base b = buildBase(weight, true);
    std::printf("factor base weight<=%d even-trace: %zu points, %d abscissae\n",
                weight, b.pts.size(), b.abscissae);

    int plantedFound = 0;
    bool plantedOk = true;
    const int stride = std::max(2, (int)b.pts.size() / (plantedN + 3));
    for (int t = 0; t < plantedN; ++t) {
        int i = (t * stride) % (int)b.pts.size();
        int j = (i + 2) % (int)b.pts.size();
        int k = (i + 4) % (int)b.pts.size();
        if (i == j || j == k || i == k) { plantedOk = false; continue; }
        R::Point r = R::addPt(R::addPt(b.pts[i], b.pts[j]), b.pts[k]);
        SearchResult got = runSearch(b, pack(r.x), pack(r.y), 32, true);
        bool ok = got.verified && got.hits > 0;
        plantedFound += ok;
        plantedOk = plantedOk && ok;
        std::printf("planted %d indices %d+%d+%d hits=%llu verified=%d\n",
                    t, i, j, k, got.hits, got.verified);
    }

    PairTable tab{};
    TableMeas tableMeas{};
    SearchResult tableProbe{}, tableGen{};
    int tablePlantedFound = 0;
    bool tablePlantedOk = true;
    if (wantTable) {
        std::fprintf(stderr, "building pair-sum table (%d points)\n", (int)b.pts.size());
        tab = buildTable(b);
        std::printf("pair table stored=%d bytes=%zu build=%.6fs sort=%.6fs\n",
                    tab.nStored, tab.bytes, tab.buildSeconds, tab.sortSeconds);
        const int strideT = std::max(2, (int)b.pts.size() / (plantedN + 3));
        for (int t = 0; t < plantedN; ++t) {
            int i = (t * strideT) % (int)b.pts.size();
            int j = (i + 2) % (int)b.pts.size();
            int k = (i + 4) % (int)b.pts.size();
            if (i == j || j == k || i == k) { tablePlantedOk = false; continue; }
            R::Point r = R::addPt(R::addPt(b.pts[i], b.pts[j]), b.pts[k]);
            SearchResult got = runTableProbe(b, tab, pack(r.x), pack(r.y), 32, 1);
            bool ok = got.verified && got.hits > 0;
            tablePlantedFound += ok;
            tablePlantedOk = tablePlantedOk && ok;
            std::printf("table planted %d indices %d+%d+%d hits=%llu verified=%d\n",
                        t, i, j, k, got.hits, got.verified);
        }
        R::Point g = R::make(R::fromLimbs(eccF131::PX), R::fromLimbs(eccF131::PY));
        tableProbe = runTableProbe(b, tab, pack(g.x), pack(g.y), 8, tableProbeSteps);
        std::printf("table probe adds=%llu seconds=%.6f rate=%.3e add/s\n",
                    tableProbe.pairs, tableProbe.seconds,
                    tableProbe.seconds > 0 ? tableProbe.pairs / tableProbe.seconds : 0);
        tableGen = runTableProbe(b, tab, pack(g.x), pack(g.y), 256, 1);
        std::printf("table generator hits=%llu adds=%llu seconds=%.6f verified=%d\n",
                    tableGen.hits, tableGen.pairs, tableGen.seconds, tableGen.verified);
        tableMeas.tab = &tab;
        tableMeas.plantedTried = plantedN;
        tableMeas.plantedFound = tablePlantedFound;
        tableMeas.plantedOk = tablePlantedOk;
        tableMeas.probe = &tableProbe;
        tableMeas.generator = &tableGen;
    }

    SearchResult genGpu{}, genCpu{};
    SearchResult *genGpuPtr = 0, *genCpuPtr = 0;
    if (wantGen) {
        R::Point g = R::make(R::fromLimbs(eccF131::PX), R::fromLimbs(eccF131::PY));
        genGpu = runSearch(b, pack(g.x), pack(g.y), 256, true);
        genGpuPtr = &genGpu;
        std::printf("generator GPU hits=%llu pairs=%llu seconds=%.6f verified=%d\n",
                    genGpu.hits, genGpu.pairs, genGpu.seconds, genGpu.verified);
        if (wantCpu) {
            genCpu = runSearch(b, pack(g.x), pack(g.y), 256, false);
            genCpuPtr = &genCpu;
            std::printf("generator CPU hits=%llu pairs=%llu seconds=%.6f verified=%d\n",
                        genCpu.hits, genCpu.pairs, genCpu.seconds, genCpu.verified);
        }
    }

    double benchAdds = 0, benchSeconds = 0;
    int nBench = std::min(benchPoints, (int)b.pts.size());
    if (nBench > 128 && benchSteps > 0) {
        P131 *dx, *dy;
        checked(cudaMalloc(&dx, nBench * sizeof(P131)));
        checked(cudaMalloc(&dy, nBench * sizeof(P131)));
        checked(cudaMemcpy(dx, b.x.data(), nBench * sizeof(P131), cudaMemcpyHostToDevice));
        checked(cudaMemcpy(dy, b.y.data(), nBench * sizeof(P131), cudaMemcpyHostToDevice));
        int threads = 128;
        int blocks = (nBench + threads - 1) / threads;
        addBench<<<blocks, threads>>>(dx, dy, nBench, 1);
        checked(cudaDeviceSynchronize());
        cudaEvent_t start, stop;
        checked(cudaEventCreate(&start));
        checked(cudaEventCreate(&stop));
        checked(cudaEventRecord(start));
        addBench<<<blocks, threads>>>(dx, dy, nBench, benchSteps);
        checked(cudaEventRecord(stop));
        checked(cudaEventSynchronize(stop));
        float ms = 0;
        checked(cudaEventElapsedTime(&ms, start, stop));
        benchSeconds = ms / 1000.0;
        cudaEventDestroy(start);
        cudaEventDestroy(stop);
        benchAdds = (double)(nBench - 1) * (double)benchSteps;
        std::printf("bench adds=%.0f seconds=%.6f rate=%.3e add/s\n",
                    benchAdds, benchSeconds, benchSeconds > 0 ? benchAdds / benchSeconds : 0);
        cudaFree(dx); cudaFree(dy);
    }

    FILE *jsonOut = stdout;
    if (jsonPath && std::strcmp(jsonPath, "-")) {
        jsonOut = std::fopen(jsonPath, "w");
        if (!jsonOut) {
            std::perror(jsonPath);
            jsonOut = stdout;
        }
    }
    if (jsonPath) {
        writeJson(jsonOut, prop.name, prop.major * 10 + prop.minor, b,
                  plantedN, plantedFound, plantedOk, genGpuPtr, genCpuPtr,
                  benchAdds, benchSeconds, benchSteps,
                  wantTable ? &tableMeas : 0);
        if (jsonOut != stdout) {
            std::fclose(jsonOut);
            std::printf("wrote %s\n", jsonPath);
        }
    }
    if (wantTable) freeTable(&tab);
    return plantedOk && plantedFound == plantedN
        && (!wantTable || (tablePlantedOk && tablePlantedFound == plantedN)) ? 0 : 1;
}
