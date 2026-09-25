// The table walk, ITERATION-FUNCTION.md §4:
//
//     R' = R + (-1)^eps(R) * sigma^k(R) ( T_h(R) ),      T_h = a_h P + b_h Q
//
//     h(R)   = (HW(x_n) / 2) mod H                       normal-basis weight
//     k(R)   = (sum_e L(e) x_e) * HW(x_n)^-1  mod m      Frobenius phase
//     eps(R) = bit p(R) of y_n,  p = argmax_{e in supp x_n} (L(e) - k) mod m
//
// L(e) is the discrete logarithm base 2 of the folded coordinate index e in
// (Z/n)*/±1, so squaring shifts every L by one and k(sigma R) = k(R) + 1;
// eps is sigma-invariant and flips under negation.  Hence f(sigma R) = sigma
// f(R) and f(-R) = -f(R): the walk is a function on the 2m-element classes
// without any canonical representative, exactly as the sigma^j + 1 walk is.
//
// Unlike that walk this one is additive, so it has fruitless cycles: runs of
// steps whose addends sum to zero.  Pairs that cancel are the obvious ones,
// but Frobenius satisfies s^2 + s + 2 = 0 on these curves, so four steps in
// one branch adding sigma^(k+2) T, sigma^(k+1) T, sigma^k T, sigma^k T also
// return a walk to where it started, as do ten more families of up to six
// steps (WALK-CONSTANT.md section 5; benchmarks/walk-constant/fruitless_patterns.py).
// Every one of them is an identity in Z[tau], tau^2 + tau + 2 = 0, per branch.
// The rule maps a step to phi = +-c_h r^k in Z/2^16, r the odd root of
// x^2 + x + 2 there: a ring homomorphism from Z[tau], so any run of steps that
// returns formally maps to zero, and a run that does not maps to zero with
// probability 2^-16.  A step is refused when it would close such a run of 2, 4
// or 6 steps (5 tags of history at H = 8, 4 at H = 16, where 6-step runs
// remain); a run of odd length cannot close, since every phi is odd.  A
// refused step advances h.  The decision depends only on the last few steps,
// so two trails that merge re-synchronise within a few steps: merging, and
// with it the rho collision structure, survives.  A spurious refusal is just
// another deterministic step, so the 2^-16 costs nothing but that.
//
// This header carries what host and device share: the tag encoding, the cycle
// rule and the coordinate tables for one field size.  The reference step and
// the table itself are in TableWalk<Cfg> below (host only).
#pragma once
#include <stdint.h>

#ifndef ECC_WALK_TABLE
#define ECC_WALK_TABLE 0
#endif
#ifndef ECC_TABLE_BRANCHES
#define ECC_TABLE_BRANCHES 8
#endif
#if ECC_TABLE_BRANCHES != 8 && ECC_TABLE_BRANCHES != 16
#error "ECC_TABLE_BRANCHES must be 8 or 16"
#endif
#include "bitslice.h"

// Step tag: h in the low ECC_TAG_HBITS bits, k in the next 8, eps above.
// The history word holds the last ECC_HIST_DEPTH tags, newest in the low slot.
// A lane that has not stepped yet holds all ones in every slot, whose k (255)
// is no phase; phi reads it as zero.
#if ECC_TABLE_BRANCHES == 8
#define ECC_TAG_HBITS 3
#define ECC_HIST_SLOT 12
#else
#define ECC_TAG_HBITS 4
#define ECC_HIST_SLOT 16
#endif
#define ECC_HIST_DEPTH (64 / ECC_HIST_SLOT)
#define ECC_CYCLE_WINDOWS ((ECC_HIST_DEPTH + 1) / 2)   // closing lengths 2, 4 (, 6)
#define ECC_TAG_EPS (1u << (ECC_TAG_HBITS + 8))
#define ECC_TAG_MASK ((1u << ECC_HIST_SLOT) - 1u)
#define ECC_TAG_NONE ECC_TAG_MASK
#define ECC_HIST_EMPTY 0xFFFFFFFFFFFFFFFFull

ECC_HD unsigned eccTag(int h, int k, int eps) {
    return unsigned(h) | (unsigned(k) << ECC_TAG_HBITS) | (unsigned(eps) << (ECC_TAG_HBITS + 8));
}
ECC_HD int eccTagH(unsigned t) { return int(t & ((1u << ECC_TAG_HBITS) - 1u)); }
ECC_HD int eccTagK(unsigned t) { return int((t >> ECC_TAG_HBITS) & 255u); }
ECC_HD int eccTagEps(unsigned t) { return int((t >> (ECC_TAG_HBITS + 8)) & 1u); }
ECC_HD unsigned long long eccHistPush(unsigned long long hist, unsigned t) {
    return (hist << ECC_HIST_SLOT) | (unsigned long long)(t & ECC_TAG_MASK);
}

// The odd root of x^2 + x + 2 modulo 2^16, lifted bit by bit (the derivative
// 2x + 1 is odd), and its powers.
constexpr uint32_t eccCycleRoot() {
    uint32_t r = 1;
    for (int j = 1; j < 16; ++j)
        if ((r * r + r + 2u) & (1u << j)) r += 1u << j;
    return r & 0xFFFFu;
}
constexpr uint32_t eccCyclePow(int e) {
    uint32_t v = 1;
    for (int i = 0; i < e; ++i) v = (v * eccCycleRoot()) & 0xFFFFu;
    return v;
}
static_assert(((eccCycleRoot() * eccCycleRoot() + eccCycleRoot() + 2u) & 0xFFFFu) == 0, "cycle root");

// phi of one step in the frame cut at phase b: exponents run (k - b) mod m, so
// r^k picks up r^m below the cut, and the common factor r^-b is dropped.  The
// frame puts the new step's own phase opposite the cut, so every run that
// could close through it is uncut.  rpow[k] = r^k for k < m and rpow[m] = 0.
ECC_HD uint32_t eccCyclePhi(unsigned t, int b, const uint16_t *rpow, int m, uint32_t rm) {
    int k = eccTagK(t);
    k = k < m ? k : m;
    uint32_t v = uint32_t(rpow[k]) * (uint32_t(2 * eccTagH(t) + 1) * 0x9E3779B9u);
    if (k < b) v *= rm;
    return eccTagEps(t) ? 0u - v : v;
}
// What the rule needs from the history for a step of phase k: the cut, and phi
// summed over the last 1, 3 (and 5) steps.  Retries change only h, so this is
// computed once per step.
struct EccCycleWindow {
    int b;
    uint32_t s[ECC_CYCLE_WINDOWS];
};
ECC_HD EccCycleWindow eccCycleWindow(int k, unsigned long long hist, const uint16_t *rpow, int m, uint32_t rm) {
    EccCycleWindow w;
    w.b = k + (m + 1) / 2;
    if (w.b >= m) w.b -= m;
    uint32_t acc = 0;
#pragma unroll
    for (int i = 0; i < 2 * ECC_CYCLE_WINDOWS - 1; ++i) {
        acc += eccCyclePhi(unsigned(hist >> (ECC_HIST_SLOT * i)) & ECC_TAG_MASK, w.b, rpow, m, rm);
        if ((i & 1) == 0) w.s[i >> 1] = acc;
    }
    return w;
}
// A step with tag t is fruitless when it closes a run of 2, 4 (or 6) steps.
ECC_HD bool eccTagFruitless(unsigned t, const EccCycleWindow &w, const uint16_t *rpow, int m, uint32_t rm) {
    const uint32_t v = eccCyclePhi(t, w.b, rpow, m, rm);
    bool closes = false;
#pragma unroll
    for (int i = 0; i < ECC_CYCLE_WINDOWS; ++i) closes |= ((v + w.s[i]) & 0xFFFFu) == 0;
    return closes;
}
// The tag after the rule: advance the branch while the step would be
// fruitless, at most H times.
ECC_HD unsigned eccResolveTag(unsigned t, unsigned long long hist, const uint16_t *rpow, int m, uint32_t rm) {
    const EccCycleWindow w = eccCycleWindow(eccTagK(t), hist, rpow, m, rm);
    for (int i = 0; i < ECC_TABLE_BRANCHES && eccTagFruitless(t, w, rpow, m, rm); ++i)
        t = eccTag((eccTagH(t) + 1) & (ECC_TABLE_BRANCHES - 1), eccTagK(t), eccTagEps(t));
    return t;
}

// Coordinate tables for GF(2^M) in the permuted type-II ONB, coordinate i in
// 1..M stored at bit i-1.  Host-side construction; the device gets copies.
template <int M>
struct TableWalkConsts {
    static const int N = 2 * M + 1;
    static const int NL = (M + 63) / 64;
    int L[M + 1];                  // L[i], i = 1..M; L[fold(2^t)] = t
    int inv[M + 1];                // inverse of w mod M, inv[0] = 0
    unsigned long long maskLt[M][NL];        // {i : L(i) < k}
    unsigned long long plane[8][NL];         // {i : bit b of L(i) set}
    uint16_t rpow[M + 1];                    // cycle rule: r^k, and rpow[M] = 0
    uint32_t rm;                             // r^M

    static int fold(long long e) {
        e %= N;
        if (e < 0) e += N;
        return e > M ? int(N - e) : int(e);
    }
    static void setBit(unsigned long long *v, int i) { v[(i - 1) >> 6] |= 1ull << ((i - 1) & 63); }
    void build() {
        for (int i = 0; i <= M; ++i) L[i] = -1;
        long long e = 1;
        for (int t = 0; t < M; ++t) {
            L[fold(e)] = t;
            e = (2 * e) % N;
        }
        for (int i = 0; i <= M; ++i) {
            inv[i] = 0;
            if (i == 0) continue;
            for (int j = 1; j < M; ++j)
                if ((i * j) % M == 1) { inv[i] = j; break; }
        }
        for (int k = 0; k < M; ++k)
            for (int l = 0; l < NL; ++l) maskLt[k][l] = 0;
        for (int b = 0; b < 8; ++b)
            for (int l = 0; l < NL; ++l) plane[b][l] = 0;
        for (int k = 0; k < M; ++k) rpow[k] = uint16_t(eccCyclePow(k));
        rpow[M] = 0;
        rm = eccCyclePow(M);
        for (int i = 1; i <= M; ++i) {
            for (int k = L[i] + 1; k < M; ++k) setBit(maskLt[k], i);
            for (int b = 0; b < 8; ++b)
                if ((L[i] >> b) & 1) setBit(plane[b], i);
        }
    }
    bool consistent() const {
        bool seen[M + 1] = {};
        for (int i = 1; i <= M; ++i) {
            if (L[i] < 0 || L[i] >= M || seen[L[i]]) return false;
            seen[L[i]] = true;
        }
        return true;
    }

    // Frobenius phase of a normal-basis coordinate vector of weight hw.
    int phase(const unsigned long long *xn, int hw) const {
        int w = 0;
        for (int b = 0; b < 8; ++b) {
            int c = 0;
            for (int l = 0; l < NL; ++l) c += __builtin_popcountll(xn[l] & plane[b][l]);
            w += c << b;
        }
        return int((long long)(w % M) * inv[hw % M] % M);
    }
    // The support element whose L is last before the phase in cyclic order:
    // set bits with L < k if there are any, otherwise all set bits, and the
    // largest L among them.  Returns a one-hot vector in `p`.
    void pivot(const unsigned long long *xn, int k, unsigned long long *p) const {
        bool any = false;
        for (int l = 0; l < NL; ++l) { p[l] = xn[l] & maskLt[k][l]; any |= p[l] != 0; }
        if (!any) for (int l = 0; l < NL; ++l) p[l] = xn[l];
        for (int b = 7; b >= 0; --b) {
            unsigned long long t[NL];
            bool anyT = false;
            for (int l = 0; l < NL; ++l) { t[l] = p[l] & plane[b][l]; anyT |= t[l] != 0; }
            if (anyT) for (int l = 0; l < NL; ++l) p[l] = t[l];
        }
    }
    int negationBit(const unsigned long long *xn, const unsigned long long *yn, int k) const {
        unsigned long long p[NL];
        pivot(xn, k, p);
        int c = 0;
        for (int l = 0; l < NL; ++l) c += __builtin_popcountll(yn[l] & p[l]);
        return c & 1;
    }
};
