// One scalar walk per slot, retaining the existing iteration and DP record.
#pragma once
#include "kernel.h"
#include "packed131.h"

namespace eccPacked131 {
#ifndef ECC_PACKED_CACHE_DENOM
#define ECC_PACKED_CACHE_DENOM 0
#endif
#if ECC_PACKED_CACHE_DENOM != 0 && ECC_PACKED_CACHE_DENOM != 1
#error "ECC_PACKED_CACHE_DENOM must be 0 or 1"
#endif
#ifndef ECC_PACKED_POLY_CHAIN
#define ECC_PACKED_POLY_CHAIN 0
#endif
#ifndef ECC_PACKED_PAIR_PRODUCTS
#define ECC_PACKED_PAIR_PRODUCTS 0
#endif
#if ECC_PACKED_PAIR_PRODUCTS && !ECC_PACKED_POLY_CHAIN
#error "ECC_PACKED_PAIR_PRODUCTS requires polynomial chains"
#endif
#if ECC_PACKED_POLY_CHAIN != 0 && ECC_PACKED_POLY_CHAIN != 1
#error "ECC_PACKED_POLY_CHAIN must be 0 or 1"
#endif
#if ECC_PACKED_POLY_CHAIN && !ECC_PACKED_CACHE_DENOM
#error "ECC_PACKED_POLY_CHAIN requires the denominator cache"
#endif
#ifndef ECC_PACKED_POLY_STATE
#define ECC_PACKED_POLY_STATE 0
#endif
#if ECC_PACKED_POLY_STATE != 0 && ECC_PACKED_POLY_STATE != 1
#error "ECC_PACKED_POLY_STATE must be 0 or 1"
#endif
#if ECC_PACKED_POLY_STATE && (!ECC_PACKED_CACHE_DENOM || !ECC_PACKED_POLY_CHAIN)
#error "ECC_PACKED_POLY_STATE requires the denominator cache and polynomial chains"
#endif
static __constant__ P131 orbitX[128], orbitY[128], targetX, targetY;

__device__ __forceinline__ P131 load(const unsigned *p, int slot, int tid, int threads) {
    P131 a;
#pragma unroll
    for (int i = 0; i < 5; ++i) a.v[i] = p[(size_t(slot) * 5 + i) * threads + tid];
    return a;
}
__device__ __forceinline__ void store(unsigned *p, int slot, int tid, int threads, P131 a) {
#pragma unroll
    for (int i = 0; i < 5; ++i) p[(size_t(slot) * 5 + i) * threads + tid] = a.v[i];
}
__device__ __forceinline__ void toLimbs(P131 a, unsigned long long *out) {
    out[0] = a.v[0] | (static_cast<unsigned long long>(a.v[1]) << 32);
    out[1] = a.v[2] | (static_cast<unsigned long long>(a.v[3]) << 32);
    out[2] = a.v[4];
}
__device__ __forceinline__ int weight(P131 a) {
    return __popc(a.v[0]) + __popc(a.v[1]) + __popc(a.v[2]) + __popc(a.v[3]) + __popc(a.v[4]);
}

static __global__ void ECC_BOUNDS init(WalkParams<unsigned> p, bool reseed) {
    const size_t id = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (id >= size_t(p.threads) * ECC_BATCH || (reseed && !p.dead[id])) return;
    const unsigned long long seed = reseed ? p.seed[id] + 1 : eccSeedFor(p.runId, id);
    const unsigned long long c0 = eccPrf(seed, 0), c1 = eccPrf(seed, 1);
    P131 x = targetX, y = targetY;
#pragma unroll 1
    for (int i = 0; i < 128; ++i) {
        if (((i < 64 ? c0 >> (i & 63) : c1 >> (i & 63)) & 1) == 0) continue;
        P131 d = add131(x, orbitX[i]), e = add131(y, orbitY[i]);
        P131 lambda = mul131(e, inv131(d));
        P131 nx = add131(add131(sqr131(lambda), lambda), d);
        y = add131(add131(mul131(lambda, add131(x, nx)), nx), y);
        x = nx;
    }
    const int slot = int(id / p.threads), tid = int(id % p.threads);
#if ECC_PACKED_POLY_STATE
    // Seed construction uses the established normal-basis point arithmetic.
    // Only the persistent coordinate representation changes.
    x = toPolynomial131(x);
    y = toPolynomial131(y);
#endif
    store(p.x, slot, tid, p.threads, x);
    store(p.y, slot, tid, p.threads, y);
    p.seed[id] = seed;
    p.startIter[id] = p.iterBase;
    p.dead[id] = 0;
}

static __global__ void ECC_BOUNDS walk(WalkParams<unsigned> p, unsigned *denominators) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= p.threads) return;
#if ECC_PACKED_POLY_CHAIN && !ECC_PACKED_POLY_STATE
    unsigned *polyDenominators=denominators+size_t(p.threads)*ECC_BATCH*5;
#endif
#if !ECC_PACKED_CACHE_DENOM
    unsigned char js[ECC_BATCH];
#endif
    P131 prod, inv;
#pragma unroll 1
    for (int step = 0; step < p.steps; ++step) {
        const unsigned long long now = p.iterBase + step;
        const bool guard = p.maxIters && now % ECC_GUARD_PERIOD == 0;
#pragma unroll 1
        for (int slot = 0; slot < ECC_BATCH; ++slot) {
            P131 x = load(p.x, slot, tid, p.threads);
#if ECC_PACKED_POLY_STATE
            x = fromPolynomial131(x);
#endif
            const int hw = weight(x);
            const size_t id = size_t(slot) * p.threads + tid;
            if (!p.dead[id]) {
                if (hw <= p.dpWeight) {
                    const unsigned dest = atomicAdd(p.dpCount, 1u);
                    if (dest < p.dpCap) {
                        DpRecord rec;
                        rec.seed = p.seed[id];
                        rec.iters = now - p.startIter[id];
                        toLimbs(x, rec.x);
#if ECC_PACKED_POLY_STATE
                        toLimbs(fromPolynomial131(load(p.y, slot, tid, p.threads)), rec.y);
#else
                        toLimbs(load(p.y, slot, tid, p.threads), rec.y);
#endif
                        p.dp[dest] = rec;
                    }
                    p.dead[id] = 1;
                } else if (guard && now - p.startIter[id] >= p.maxIters) {
                    // Overdue walks need a restart, not a false DP report.
                    p.dead[id] = 1;
                    atomicAdd(p.dpCount + 1, 1u);
                }
            }
            const int j = 3 + ((hw >> 1) & 7);
#if !ECC_PACKED_CACHE_DENOM
            js[slot] = j;
#endif
            P131 d = add131(x, sigma131(x, j));
#if ECC_PACKED_POLY_CHAIN
            P131 dp = toPolynomial131(d);
            prod = slot == 0 ? dp : mulPolynomial131(prod, dp);
#if !ECC_PACKED_POLY_STATE
            store(polyDenominators, slot, tid, p.threads, dp);
#endif
#else
            prod = slot == 0 ? d : mul131(prod, d);
#endif
            if (slot + 1 < ECC_BATCH) store(p.pchain, slot, tid, p.threads, prod);
#if ECC_PACKED_CACHE_DENOM
            // The high 29 bits are unused by the field. Keep the jump index
            // beside the denominator, eliminating the separate local array.
#if ECC_PACKED_POLY_STATE
            dp.v[4] |= unsigned(j - 3) << 3;
            store(denominators, slot, tid, p.threads, dp);
#else
            d.v[4] |= unsigned(j - 3) << 3;
            store(denominators, slot, tid, p.threads, d);
#endif
#endif
        }
#if ECC_PACKED_POLY_CHAIN
        inv = toPolynomial131(inv131(fromPolynomial131(prod)));
#else
        inv = inv131(prod);
#endif
#pragma unroll 1
        for (int slot = ECC_BATCH - 1; slot >= 0; --slot) {
            P131 x = load(p.x, slot, tid, p.threads), y = load(p.y, slot, tid, p.threads);
#if ECC_PACKED_POLY_STATE
            P131 dp = load(denominators, slot, tid, p.threads), ii;
            const int j = 3 + ((dp.v[4] >> 3) & 7);
            dp.v[4] &= 7;
            const P131 normalY = fromPolynomial131(y);
            P131 ep = toPolynomial131(add131(normalY, sigma131(normalY, j)));
#elif ECC_PACKED_CACHE_DENOM
            P131 d = load(denominators, slot, tid, p.threads), ii;
            const int j = 3 + ((d.v[4] >> 3) & 7);
            d.v[4] &= 7;
            P131 e = add131(y, sigma131(y, j));
#else
            const int j = js[slot];
            P131 d = add131(x, sigma131(x, j)), e = add131(y, sigma131(y, j)), ii;
#endif
            if (slot) {
#if ECC_PACKED_POLY_CHAIN
#if ECC_PACKED_PAIR_PRODUCTS
                PolynomialPair pair=mulPolynomialPair131(inv,
                    load(p.pchain,slot-1,tid,p.threads),
#if ECC_PACKED_POLY_STATE
                    dp);
#else
                    load(polyDenominators,slot,tid,p.threads));
#endif
                ii=pair.first; inv=pair.second;
#else
                ii = mulPolynomial131(inv, load(p.pchain, slot - 1, tid, p.threads));
#if ECC_PACKED_POLY_STATE
                inv = mulPolynomial131(inv, dp);
#else
                inv = mulPolynomial131(inv, load(polyDenominators, slot, tid, p.threads));
#endif
#endif
#else
                ii = mul131(inv, load(p.pchain, slot - 1, tid, p.threads));
                inv = mul131(inv, d);
#endif
            } else ii = inv;
#if ECC_PACKED_POLY_STATE
            P131 lambdaPoly = mulPolynomial131(ep, ii);
            P131 nx = add131(add131(squarePolynomial131(lambdaPoly), lambdaPoly), dp);
            P131 product = mulPolynomial131(lambdaPoly, add131(x, nx));
            P131 ny = add131(add131(product, nx), y);
#else
#if ECC_PACKED_POLY_CHAIN
            P131 lambdaPoly = mulPolynomial131(toPolynomial131(e), ii);
            P131 lambda = fromPolynomial131(lambdaPoly);
#else
            P131 lambda = mul131(e, ii);
#endif
            P131 nx = add131(add131(sqr131(lambda), lambda), d);
#if ECC_PACKED_POLY_CHAIN
            P131 product = mulPolynomial131(lambdaPoly, toPolynomial131(add131(x, nx)));
            P131 ny = add131(add131(fromPolynomial131(product), nx), y);
#else
            P131 ny = add131(add131(mul131(lambda, add131(x, nx)), nx), y);
#endif
#endif
            store(p.x, slot, tid, p.threads, nx);
            store(p.y, slot, tid, p.threads, ny);
        }
    }
}
} // namespace eccPacked131
