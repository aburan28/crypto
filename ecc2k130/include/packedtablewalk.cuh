// Device side of the table walk (tablewalk.h) for the packed GF(2^131) kernel.
//
// Per block, dynamic shared memory holds the walk's tables:
//
//   table[k][h]  sigma^k(T_h) in the polynomial basis, 10 words (x then y)
//   maskLt[k]    normal-basis coordinates e with L(e) < k, 5 words
//   fromRow[p]   row p of the polynomial->normal transform, 5 words: the
//                normal-basis coordinate p of y is the parity of y_p & fromRow[p]
//   inv[w]       w^-1 mod 131
//
// and the eight L bit-planes, identical for every lane, live in constant
// memory.  The host fills one flat buffer in this order (twOffsets) from the
// reference TableWalk and the kernel copies it to shared memory on entry.
#pragma once
#include "tablewalk.h"
#include "packed131.h"

namespace eccPacked131 {

static const int TW_H = ECC_TABLE_BRANCHES;
static const int TW_TABLE_WORDS = 131 * TW_H * 10;
static const int TW_MASK_OFF = TW_TABLE_WORDS;
static const int TW_ROW_OFF = TW_MASK_OFF + 131 * 5;
static const int TW_INV_OFF = TW_ROW_OFF + 131 * 5;
static const int TW_WORDS = TW_INV_OFF + 132;
static const size_t TW_SHARED_BYTES = size_t(TW_WORDS) * sizeof(uint32_t);

#ifdef __CUDACC__
static __constant__ uint32_t twPlane[8][5];

__device__ __forceinline__ void twLoadShared(uint32_t *shared, const uint32_t *global) {
    for (int i = threadIdx.x; i < TW_WORDS; i += blockDim.x) shared[i] = global[i];
    __syncthreads();
}

// Frobenius phase k(x) = (sum_e L(e) x_e) * HW(x)^-1 mod 131 of a normal-basis x.
__device__ __forceinline__ int twPhase(const P131 &x, int hw, const uint32_t *inv) {
    unsigned w = 0;
#pragma unroll
    for (int b = 0; b < 8; ++b) {
        unsigned c = __popc(x.v[0] & twPlane[b][0]);
        c += __popc(x.v[1] & twPlane[b][1]);
        c += __popc(x.v[2] & twPlane[b][2]);
        c += __popc(x.v[3] & twPlane[b][3]);
        c += __popc(x.v[4] & twPlane[b][4]);
        w += c << b;
    }
    return int(((w % 131u) * inv[hw]) % 131u);
}

// Index of the support element whose L is last before k in cyclic order:
// among set bits with L < k if any, else among all set bits, the largest L.
__device__ __forceinline__ int twPivot(const P131 &x, int k, const uint32_t *maskLt) {
    const uint32_t *m = maskLt + k * 5;
    uint32_t s[5];
#pragma unroll
    for (int i = 0; i < 5; ++i) s[i] = x.v[i] & m[i];
    const bool any = (s[0] | s[1] | s[2] | s[3] | s[4]) != 0;
#pragma unroll
    for (int i = 0; i < 5; ++i) s[i] = any ? s[i] : x.v[i];
#pragma unroll
    for (int b = 7; b >= 0; --b) {
        uint32_t t[5];
#pragma unroll
        for (int i = 0; i < 5; ++i) t[i] = s[i] & twPlane[b][i];
        const bool anyT = (t[0] | t[1] | t[2] | t[3] | t[4]) != 0;
#pragma unroll
        for (int i = 0; i < 5; ++i) s[i] = anyT ? t[i] : s[i];
    }
    int p = 0;
#pragma unroll
    for (int i = 4; i >= 0; --i) p = s[i] ? 32 * i + __ffs(s[i]) - 1 : p;
    return p;
}

// Coordinate p of the normal-basis image of a polynomial-basis y.
__device__ __forceinline__ int twCoordinate(const P131 &yp, int p, const uint32_t *fromRow) {
    const uint32_t *r = fromRow + p * 5;
    unsigned c = __popc(yp.v[0] & r[0]) + __popc(yp.v[1] & r[1]) + __popc(yp.v[2] & r[2]) +
                 __popc(yp.v[3] & r[3]) + __popc(yp.v[4] & r[4]);
    return int(c & 1u);
}

// The whole selection for one point: (h, k, eps) after the cycle rule, with
// the history advanced.  x is in the normal basis, yp in the polynomial basis.
__device__ __forceinline__ unsigned twSelect(const P131 &x, const P131 &yp, int hw,
                                             unsigned long long *hist, const uint32_t *shared) {
    const int k = twPhase(x, hw, shared + TW_INV_OFF);
    const int eps = twCoordinate(yp, twPivot(x, k, shared + TW_MASK_OFF), shared + TW_ROW_OFF);
    int h = (hw >> 1) & (TW_H - 1);
    unsigned tag = eccTag(h, k, eps);
    const unsigned long long old = *hist;
    while (eccTagFruitless(tag, old)) {
        h = (h + 1) & (TW_H - 1);
        tag = eccTag(h, k, eps);
    }
    *hist = eccHistPush(old, tag);
    return tag;
}

// d = x + x_T and e = y + y_T (+ x_T when the table point is negated), in the
// polynomial basis, for the selected tag.
__device__ __forceinline__ void twAddend(unsigned tag, const P131 &xp, const P131 &yp,
                                         const uint32_t *shared, P131 *d, P131 *e) {
    const uint32_t *t = shared + (eccTagK(tag) * TW_H + eccTagH(tag)) * 10;
    const uint32_t negMask = 0u - unsigned(eccTagEps(tag));
#pragma unroll
    for (int i = 0; i < 5; ++i) {
        const uint32_t tx = t[i];
        d->v[i] = xp.v[i] ^ tx;
        e->v[i] = yp.v[i] ^ t[5 + i] ^ (tx & negMask);
    }
}
#endif  // __CUDACC__

// Host: fill the flat constant buffer from the reference walk.
template <class TW>
inline void twFillConsts(const TW &walk, uint32_t *out, uint32_t planes[8][5]) {
    auto pack = [](const unsigned long long *v, uint32_t *w) {
        for (int i = 0; i < 5; ++i) w[i] = uint32_t(v[i / 2] >> (32 * (i & 1)));
    };
    for (int k = 0; k < 131; ++k)
        for (int h = 0; h < TW_H; ++h) {
            P131 x, y;
            pack(walk.table[h][k].x.v, x.v);
            pack(walk.table[h][k].y.v, y.v);
            x = toPolynomial131(x);
            y = toPolynomial131(y);
            uint32_t *t = out + (k * TW_H + h) * 10;
            for (int i = 0; i < 5; ++i) { t[i] = x.v[i]; t[5 + i] = y.v[i]; }
        }
    for (int k = 0; k < 131; ++k) pack(walk.consts.maskLt[k], out + TW_MASK_OFF + k * 5);
    for (int b = 0; b < 8; ++b) pack(walk.consts.plane[b], planes[b]);
    for (int p = 0; p < 131 * 5; ++p) out[TW_ROW_OFF + p] = 0;
    for (int j = 0; j < 131; ++j) {
        P131 e = {{0, 0, 0, 0, 0}};
        e.v[j >> 5] = 1u << (j & 31);
        const P131 n = fromPolynomial131(e);
        for (int p = 0; p < 131; ++p)
            if ((n.v[p >> 5] >> (p & 31)) & 1u) out[TW_ROW_OFF + p * 5 + (j >> 5)] |= 1u << (j & 31);
    }
    for (int w = 0; w < 132; ++w) out[TW_INV_OFF + w] = uint32_t(walk.consts.inv[w]);
}

} // namespace eccPacked131
