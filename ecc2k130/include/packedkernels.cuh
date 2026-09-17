// One scalar walk per slot, retaining the existing iteration and DP record.
#pragma once
#include "kernel.h"
#include "packed131.h"
#ifndef ECC_PACKED_BLOCK_INVERSE
#define ECC_PACKED_BLOCK_INVERSE 0
#endif
#if ECC_PACKED_BLOCK_INVERSE != 0 && ECC_PACKED_BLOCK_INVERSE != 1
#error "ECC_PACKED_BLOCK_INVERSE must be 0 or 1"
#endif
#if ECC_PACKED_LOGICAL_PAIR_INVERSE && !ECC_PACKED_BLOCK_INVERSE
#error "LOGICAL_PAIR_INVERSE requires BLOCK_INVERSE"
#endif
#if ECC_PACKED_BLOCK_INVERSE
#include "packedblockinverse131.cuh"
#endif
#ifndef ECC_PACKED_COMPACT_STATE
#define ECC_PACKED_COMPACT_STATE 0
#endif
#if ECC_PACKED_COMPACT_STATE != 0 && ECC_PACKED_COMPACT_STATE != 1
#error "ECC_PACKED_COMPACT_STATE must be 0 or 1"
#endif
#if ECC_PACKED_COMPACT_STATE
#include "packedcompactstate.cuh"
#endif

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
#if ECC_PACKED_BLOCK_INVERSE && !ECC_PACKED_POLY_CHAIN
#error "ECC_PACKED_BLOCK_INVERSE requires polynomial chains"
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
#if ECC_PACKED_WEIGHTED_PREFIX && (!ECC_PACKED_POLY_STATE || !ECC_PACKED_POLY_CHAIN || !ECC_PACKED_CACHE_DENOM || !ECC_PACKED_PAIR_PRODUCTS)
#error "ECC_PACKED_WEIGHTED_PREFIX requires polynomial state, polynomial chains, denominator cache and paired products"
#endif
#if ECC_PACKED_WEIGHTED_PREFIX == 2 && !(ECC_PACKED_PERM_SIGMA & 1)
#error "ECC_PACKED_WEIGHTED_PREFIX=2 requires the walk permutation network"
#endif
#ifndef ECC_PACKED_STATE_TILE
#define ECC_PACKED_STATE_TILE 0
#endif
#if ECC_PACKED_STATE_TILE != 0 && ECC_PACKED_STATE_TILE != 256
#error "ECC_PACKED_STATE_TILE must be 0 or 256"
#endif
#if ECC_PACKED_STATE_TILE && (!ECC_PACKED_POLY_STATE || !ECC_PACKED_CACHE_DENOM || !ECC_PACKED_POLY_CHAIN)
#error "ECC_PACKED_STATE_TILE requires polynomial state, denominator cache and polynomial chains"
#endif
#if ECC_PACKED_STATE_TILE && ECC_THREADS != 256
#error "ECC_PACKED_STATE_TILE requires ECC_THREADS=256"
#endif
#if ECC_PACKED_COMPACT_STATE && (ECC_PACKED_STATE_TILE != 256 || !ECC_PACKED_POLY_STATE || !ECC_PACKED_CACHE_DENOM || !ECC_PACKED_POLY_CHAIN)
#error "ECC_PACKED_COMPACT_STATE requires TILE256, polynomial state, denominator cache and polynomial chains"
#endif
#if ECC_PACKED_FUSED_SIGMA && (!ECC_PACKED_POLY_STATE || ECC_PACKED_WEIGHTED_PREFIX != 2)
#error "ECC_PACKED_FUSED_SIGMA requires polynomial state and weighted-prefix mode 2"
#endif
#ifndef ECC_PACKED_BATCH_SPLIT
#define ECC_PACKED_BATCH_SPLIT 1
#endif
#if ECC_PACKED_BATCH_SPLIT != 1 && ECC_PACKED_BATCH_SPLIT != 2
#error "ECC_PACKED_BATCH_SPLIT must be 1 or 2"
#endif
#if ECC_PACKED_BATCH_SPLIT == 2 && ((ECC_BATCH < 16 || ECC_BATCH > 32 || (ECC_BATCH & 1)) || ECC_THREADS != 256 || !ECC_PACKED_BLOCK_INVERSE || ECC_PACKED_WEIGHTED_PREFIX != 2)
#error "Split batches require an even batch16..32, threads256, block inversion and weighted-prefix mode2"
#endif
#ifndef ECC_PACKED_LAST_SLOT_CACHE
#define ECC_PACKED_LAST_SLOT_CACHE 0
#endif
#if ECC_PACKED_LAST_SLOT_CACHE < 0 || ECC_PACKED_LAST_SLOT_CACHE > 2
#error "ECC_PACKED_LAST_SLOT_CACHE must be 0, 1 or 2"
#endif
#if ECC_PACKED_LAST_SLOT_CACHE && ((ECC_BATCH < 16 || ECC_BATCH > 32 || (ECC_BATCH & 1)) || ECC_PACKED_BATCH_SPLIT != 2 || !ECC_PACKED_BLOCK_INVERSE || ECC_PACKED_WEIGHTED_PREFIX != 2)
#error "Last-slot cache requires an even batch16..32, split2, block inversion and weighted-prefix mode2"
#endif
#if ECC_PACKED_LAST_SLOT_CACHE
static constexpr int lastLocalSlot131 = ECC_BATCH / ECC_PACKED_BATCH_SPLIT - 1;
#endif
#ifndef ECC_PACKED_SHARED_X_SLOTS
#define ECC_PACKED_SHARED_X_SLOTS 0
#endif
#if ECC_PACKED_SHARED_X_SLOTS != 0 && ECC_PACKED_SHARED_X_SLOTS != 2 && ECC_PACKED_SHARED_X_SLOTS != 4
#error "SHARED_X_SLOTS must be 0, 2 or 4"
#endif
#if ECC_PACKED_SHARED_X_SLOTS && ((ECC_BATCH < 16 || ECC_BATCH > 32 || (ECC_BATCH & 1)) || ECC_THREADS != 256 || ECC_PACKED_BATCH_SPLIT != 2 || !ECC_PACKED_BLOCK_INVERSE || !ECC_PACKED_COMPACT_STATE || !ECC_PACKED_POLY_STATE || ECC_PACKED_WEIGHTED_PREFIX != 2)
#error "SHARED_X_SLOTS requires an even B16..32 split2 compact polynomial block-inverse layout"
#endif
#if ECC_PACKED_SHARED_X_SLOTS
static constexpr int firstSharedXLocalSlot131=ECC_BATCH/ECC_PACKED_BATCH_SPLIT-ECC_PACKED_SHARED_X_SLOTS;
__device__ __forceinline__ P131 loadSharedX131(const uint4*low,const unsigned char*tail,int index) {
    const uint4 v=low[index];return P131{{v.x,v.y,v.z,v.w,unsigned(tail[index])}};
}
__device__ __forceinline__ void storeSharedX131(uint4*low,unsigned char*tail,int index,P131 value) {
    low[index]=uint4{value.v[0],value.v[1],value.v[2],value.v[3]};
    tail[index]=static_cast<unsigned char>(value.v[4]);
}
#endif
static constexpr int walkWorkersPerBlock131 = ECC_THREADS / ECC_PACKED_BATCH_SPLIT;
#if ECC_PACKED_STATE_TILE
ECC_HD size_t physicalStateThreads(size_t threads) {
    return ((threads + 255) / 256) * 256;
}
ECC_HD size_t stateWordIndex(int slot, int word, int tid) {
    return ((size_t(tid) / 256 * ECC_BATCH * 5 + size_t(slot) * 5 + word) * 256)
           + size_t(tid) % 256;
}
#endif
static __constant__ P131 orbitX[128], orbitY[128], targetX, targetY;

__device__ __forceinline__ P131 load(const unsigned *p, int slot, int tid, int threads) {
#if ECC_PACKED_COMPACT_STATE
    return compactLoad131(p, slot, tid);
#else
    P131 a;
#pragma unroll
#if ECC_PACKED_STATE_TILE
    for (int i = 0; i < 5; ++i) a.v[i] = p[stateWordIndex(slot, i, tid)];
#else
    for (int i = 0; i < 5; ++i) a.v[i] = p[(size_t(slot) * 5 + i) * threads + tid];
#endif
    return a;
#endif
}
__device__ __forceinline__ void store(unsigned *p, int slot, int tid, int threads, P131 a) {
#if ECC_PACKED_COMPACT_STATE
    compactStore131(p, slot, tid, a);
#else
#pragma unroll
#if ECC_PACKED_STATE_TILE
    for (int i = 0; i < 5; ++i) p[stateWordIndex(slot, i, tid)] = a.v[i];
#else
    for (int i = 0; i < 5; ++i) p[(size_t(slot) * 5 + i) * threads + tid] = a.v[i];
#endif
#endif
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

#if !ECC_PACKED_XONLY_23
static __global__ void ECC_BOUNDS walk(WalkParams<unsigned> p, unsigned *denominators) {
#if ECC_PACKED_BATCH_SPLIT == 2
    // Paired physical threads own alternating slots of one logical worker.
    const int tid = blockIdx.x * walkWorkersPerBlock131 + threadIdx.x % walkWorkersPerBlock131;
    const int firstSlot = threadIdx.x / walkWorkersPerBlock131;
#else
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int firstSlot = 0;
#endif
#if ECC_PACKED_SHARED_SIGMA && !ECC_PACKED_FUSED_SIGMA
    // All block threads participate, including inactive partial-tile workers.
    initSigmaWalkShared131();
#endif
#if ECC_PACKED_FUSED_SIGMA
    initFusedDelta131();
#endif
#if ECC_PACKED_BLOCK_INVERSE
    const bool active = tid < p.threads;
    __shared__ uint32_t inverseTree[blockInverseWords131 + (ECC_PACKED_BATCH_SPLIT == 2 ? blockInverseFlagWords131 : 0)];
#else
    if (tid >= p.threads) return;
#endif
#if ECC_PACKED_POLY_CHAIN && !ECC_PACKED_POLY_STATE
    unsigned *polyDenominators=denominators+size_t(p.threads)*ECC_BATCH*5;
#endif
#if !ECC_PACKED_CACHE_DENOM
    unsigned char js[ECC_BATCH];
#endif
#if ECC_PACKED_SHARED_X_SLOTS
    __shared__ uint4 sharedXLow[ECC_PACKED_SHARED_X_SLOTS*ECC_THREADS];
    __shared__ unsigned char sharedXTail[ECC_PACKED_SHARED_X_SLOTS*ECC_THREADS];
    if (active && p.steps > 0) {
#pragma unroll
        for (int cacheSlot=0;cacheSlot<ECC_PACKED_SHARED_X_SLOTS;++cacheSlot) {
            const int slot=(firstSharedXLocalSlot131+cacheSlot)*ECC_PACKED_BATCH_SPLIT+firstSlot;
            const int index=cacheSlot*ECC_THREADS+threadIdx.x;
            storeSharedX131(sharedXLow,sharedXTail,index,load(p.x,slot,tid,p.threads));
        }
    }
#endif
    P131 prod, inv;
#pragma unroll 1
    for (int step = 0; step < p.steps; ++step) {
        const unsigned long long now = p.iterBase + step;
        const bool guard = p.maxIters && now % ECC_GUARD_PERIOD == 0;
#if ECC_PACKED_LAST_SLOT_CACHE
        P131 lastDenominator{};
#if ECC_PACKED_LAST_SLOT_CACHE == 2
        P131 lastWeightedPrefix{};
#endif
#endif
#if ECC_PACKED_BLOCK_INVERSE
        prod = P131{{1,0,0,0,0}};
        if (active) {
#endif
#pragma unroll 1
        for (int localSlot = 0; localSlot < ECC_BATCH / ECC_PACKED_BATCH_SPLIT; ++localSlot) {
            const int slot = localSlot * ECC_PACKED_BATCH_SPLIT + firstSlot;
#if ECC_PACKED_SHARED_X_SLOTS
            P131 x;
            if (localSlot >= firstSharedXLocalSlot131)
                x=loadSharedX131(sharedXLow,sharedXTail,(localSlot-firstSharedXLocalSlot131)*ECC_THREADS+threadIdx.x);
            else x=load(p.x,slot,tid,p.threads);
#else
            P131 x = load(p.x, slot, tid, p.threads);
#endif
#if ECC_PACKED_FUSED_SIGMA
            const P131 polyX=x;
#endif
#if ECC_PACKED_POLY_STATE
            x = fromPolynomial131(x);
#endif
            const int hw = weight(x);
            const size_t id = size_t(slot) * p.threads + tid;
            if (!p.dead[id]) {
                if (hw <= p.dpWeight) {
                    if ((p.seed[id] & 0xffffull) == 0xffffull) atomicAdd(p.dpCount + 2, 1u);
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
                    if ((p.seed[id] & 0xffffull) == 0xffffull) atomicAdd(p.dpCount + 2, 1u);
                    // Overdue walks need a restart, not a false DP report.
                    p.dead[id] = 1;
                    atomicAdd(p.dpCount + 1, 1u);
                }
            }
            const int j = 3 + ((hw >> 1) & 7);
#if !ECC_PACKED_CACHE_DENOM
            js[slot] = j;
#endif
#if ECC_PACKED_WEIGHTED_PREFIX
#if ECC_PACKED_FUSED_SIGMA
            const auto deltas=fusedDeltaPair131(polyX,load(p.y,slot,tid,p.threads),j-3);
            P131 dp=deltas.first, ep=deltas.second;
#else
            const P131 normalY = fromPolynomial131(load(p.y, slot, tid, p.threads));
#if ECC_PACKED_WEIGHTED_PREFIX == 2
#if ECC_PACKED_SHARED_SIGMA
            const SigmaWalkPair131 sigmas = sigmaWalkNetworkPairShared131(x, normalY, j - 3);
#else
            const SigmaWalkPair131 sigmas = sigmaWalkNetworkPair131(x, normalY, j - 3);
#endif
            P131 d = add131(x, sigmas.first);
            P131 ep = toPolynomial131(add131(normalY, sigmas.second));
#else
            P131 d = add131(x, sigma131(x, j));
            P131 ep = toPolynomial131(add131(normalY, sigma131(normalY, j)));
#endif
            P131 dp = toPolynomial131(d);
#endif
            // W_i = E_i * product_{k<i}(D_k). Every prefix slot is used.
            if (localSlot) {
                PolynomialPair pair = mulPolynomialPair131(prod, ep, dp);
#if ECC_PACKED_LAST_SLOT_CACHE == 2
                if (localSlot == lastLocalSlot131) lastWeightedPrefix = pair.first;
                else
#endif
                store(p.pchain, slot, tid, p.threads, pair.first);
                prod = pair.second;
            } else {
                prod = dp;
                store(p.pchain, slot, tid, p.threads, ep);
            }
#else
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
#endif
#if ECC_PACKED_CACHE_DENOM
            // The high 29 bits are unused by the field. Keep the jump index
            // beside the denominator, eliminating the separate local array.
#if ECC_PACKED_POLY_STATE
            dp.v[4] |= unsigned(j - 3) << 3;
#if ECC_PACKED_LAST_SLOT_CACHE
            if (localSlot == lastLocalSlot131) lastDenominator = dp;
            else
#endif
            store(denominators, slot, tid, p.threads, dp);
#else
            d.v[4] |= unsigned(j - 3) << 3;
            store(denominators, slot, tid, p.threads, d);
#endif
#endif
        }
#if ECC_PACKED_BLOCK_INVERSE
        } // active forward state access
        inv = blockInverse131<(ECC_PACKED_BATCH_SPLIT == 2)>(prod, inverseTree);
#elif ECC_PACKED_POLY_CHAIN
        inv = toPolynomial131(inv131(fromPolynomial131(prod)));
#else
        inv = inv131(prod);
#endif
#if ECC_PACKED_BLOCK_INVERSE
        if (active) {
#endif
#pragma unroll 1
        for (int localSlot = ECC_BATCH / ECC_PACKED_BATCH_SPLIT - 1; localSlot >= 0; --localSlot) {
            const int slot = localSlot * ECC_PACKED_BATCH_SPLIT + firstSlot;
#if ECC_PACKED_SHARED_X_SLOTS
            P131 x;
            if (localSlot >= firstSharedXLocalSlot131)
                x=loadSharedX131(sharedXLow,sharedXTail,(localSlot-firstSharedXLocalSlot131)*ECC_THREADS+threadIdx.x);
            else x=load(p.x,slot,tid,p.threads);
            P131 y=load(p.y,slot,tid,p.threads);
#else
            P131 x = load(p.x, slot, tid, p.threads), y = load(p.y, slot, tid, p.threads);
#endif
#if ECC_PACKED_WEIGHTED_PREFIX
#if ECC_PACKED_LAST_SLOT_CACHE
            P131 dp = localSlot == lastLocalSlot131 ? lastDenominator : load(denominators, slot, tid, p.threads);
#else
            P131 dp = load(denominators, slot, tid, p.threads);
#endif
            dp.v[4] &= 7;
            P131 lambdaPoly;
            if (localSlot) {
                PolynomialPair pair = mulPolynomialPair131(inv,
#if ECC_PACKED_LAST_SLOT_CACHE == 2
                    localSlot == lastLocalSlot131 ? lastWeightedPrefix : load(p.pchain, slot, tid, p.threads), dp);
#else
                    load(p.pchain, slot, tid, p.threads), dp);
#endif
                lambdaPoly = pair.first;
                inv = pair.second;
            } else {
                lambdaPoly = mulPolynomial131(inv, load(p.pchain, firstSlot, tid, p.threads));
            }
#else
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
            if (localSlot) {
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
#endif
#endif
#if ECC_PACKED_POLY_STATE
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
#if ECC_PACKED_SHARED_X_SLOTS
            if (localSlot >= firstSharedXLocalSlot131)
                storeSharedX131(sharedXLow,sharedXTail,(localSlot-firstSharedXLocalSlot131)*ECC_THREADS+threadIdx.x,nx);
            else store(p.x,slot,tid,p.threads,nx);
#else
            store(p.x, slot, tid, p.threads, nx);
#endif
            store(p.y, slot, tid, p.threads, ny);
        }
#if ECC_PACKED_BLOCK_INVERSE
        } // active backward state access
#endif
    }
#if ECC_PACKED_SHARED_X_SLOTS
    if (active && p.steps > 0) {
#pragma unroll
        for (int cacheSlot=0;cacheSlot<ECC_PACKED_SHARED_X_SLOTS;++cacheSlot) {
            const int slot=(firstSharedXLocalSlot131+cacheSlot)*ECC_PACKED_BATCH_SPLIT+firstSlot;
            const int index=cacheSlot*ECC_THREADS+threadIdx.x;
            store(p.x,slot,tid,p.threads,loadSharedX131(sharedXLow,sharedXTail,index));
        }
    }
#endif
}
#else
#include "packedxonly23.cuh"
#endif
} // namespace eccPacked131
