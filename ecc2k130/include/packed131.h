// Packed GF(2^131) in the same permuted type-II normal basis as eccF131.
// Each element occupies five uint32_t words. Squaring is a bit permutation;
// multiplication uses gamma_i*gamma_j = gamma_(i+j) + gamma_(i-j).
#pragma once
#ifndef ECC_PACKED_WEIGHTED_PREFIX
#define ECC_PACKED_WEIGHTED_PREFIX 0
#endif
#if ECC_PACKED_WEIGHTED_PREFIX < 0 || ECC_PACKED_WEIGHTED_PREFIX > 2
#error "ECC_PACKED_WEIGHTED_PREFIX must be 0, 1 or 2"
#endif
#ifndef ECC_PACKED_SHARED_SIGMA
#define ECC_PACKED_SHARED_SIGMA 0
#endif
#if ECC_PACKED_SHARED_SIGMA != 0 && ECC_PACKED_SHARED_SIGMA != 1
#error "ECC_PACKED_SHARED_SIGMA must be 0 or 1"
#endif
#if ECC_PACKED_SHARED_SIGMA && ECC_PACKED_WEIGHTED_PREFIX != 2
#error "ECC_PACKED_SHARED_SIGMA requires weighted-prefix mode 2"
#endif
#ifndef ECC_PACKED_CLMAD
#define ECC_PACKED_CLMAD 0
#endif
#if ECC_PACKED_CLMAD != 0 && ECC_PACKED_CLMAD != 1
#error "ECC_PACKED_CLMAD must be 0 or 1"
#endif
#if ECC_PACKED_CLMAD && defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ < 13 || (__CUDACC_VER_MAJOR__ == 13 && __CUDACC_VER_MINOR__ < 3))
#error "ECC_PACKED_CLMAD requires CUDA 13.3 or newer (PTX 9.3)"
#endif
/* Fat ARCHES="75 89 120" CLMAD=1: emit clmad on sm_80+ slices and keep the
   software product on Turing (sm_75, g4dn). An #error here would make that
   mixed-fleet client unbuildable. */
#if ECC_PACKED_CLMAD && defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 800
#define ECC_USE_CLMAD_INSN 1
#else
#define ECC_USE_CLMAD_INSN 0
#endif
#ifndef ECC_PACKED_TOP_CLMAD
#define ECC_PACKED_TOP_CLMAD 0
#endif
#if ECC_PACKED_TOP_CLMAD != 0 && ECC_PACKED_TOP_CLMAD != 1
#error "ECC_PACKED_TOP_CLMAD must be 0 or 1"
#endif
#if ECC_PACKED_TOP_CLMAD && !ECC_PACKED_CLMAD
#error "ECC_PACKED_TOP_CLMAD requires ECC_PACKED_CLMAD"
#endif
#ifndef ECC_PACKED_CLMAD_SQUARE
#define ECC_PACKED_CLMAD_SQUARE 0
#endif
#if ECC_PACKED_CLMAD_SQUARE != 0 && ECC_PACKED_CLMAD_SQUARE != 1
#error "ECC_PACKED_CLMAD_SQUARE must be 0 or 1"
#endif
#if ECC_PACKED_CLMAD_SQUARE && !ECC_PACKED_CLMAD
#error "ECC_PACKED_CLMAD_SQUARE requires ECC_PACKED_CLMAD"
#endif
#ifndef ECC_PACKED_KARAT3
#define ECC_PACKED_KARAT3 0
#endif
#if ECC_PACKED_KARAT3 != 0 && ECC_PACKED_KARAT3 != 1
#error "ECC_PACKED_KARAT3 must be 0 or 1"
#endif
#if ECC_PACKED_KARAT3 && !ECC_PACKED_CLMAD
#error "ECC_PACKED_KARAT3 requires ECC_PACKED_CLMAD"
#endif
#include "bitslice.h"
namespace eccPacked131 {
#ifndef ECC_PACKED_ADD_COMBINE
#define ECC_PACKED_ADD_COMBINE 0
#endif
#if ECC_PACKED_ADD_COMBINE != 0 && ECC_PACKED_ADD_COMBINE != 1
#error "ECC_PACKED_ADD_COMBINE must be 0 or 1"
#endif
#if ECC_PACKED_ADD_COMBINE
// The four masked class products of clmul32 occupy disjoint bit positions, so
// their union is also their integer sum. Spelling the union as PTX adds keeps
// the compiler from folding it back into LOP3 ORs (nvvm canonicalises adds of
// provably disjoint values to ors), so ptxas emits IADD3/IMAD.IADD instead:
// benchmarks/hardware-limits measured that adds can co-issue with the logic
// pipe that the rest of this multiplier saturates. Two combines per word, per
// clmul32. The host path adds too, so a non-disjoint operand fails the tests.
ECC_HD uint32_t disjointUnion4(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
#ifdef __CUDA_ARCH__
    uint32_t ab, cd, r;
    asm("add.u32 %0, %1, %2;" : "=r"(ab) : "r"(a), "r"(b));
    asm("add.u32 %0, %1, %2;" : "=r"(cd) : "r"(c), "r"(d));
    asm("add.u32 %0, %1, %2;" : "=r"(r) : "r"(ab), "r"(cd));
    return r;
#else
    return a + b + c + d;
#endif
}
#endif
// Integer-mask carryless primitives adapted from gpu/ecc2k/f2m.cuh.
ECC_HD uint32_t clmul32(uint32_t x, uint32_t y, uint32_t *hi) {
#if !ECC_PACKED_ADD_COMBINE
    const uint64_t m0 = 0x1111111111111111ull;
    const uint64_t m1 = 0x2222222222222222ull;
    const uint64_t m2 = 0x4444444444444444ull;
    const uint64_t m3 = 0x8888888888888888ull;
#endif

    uint64_t x0 = x & 0x11111111u, x1 = x & 0x22222222u;
    uint64_t x2 = x & 0x44444444u, x3 = x & 0x88888888u;
    uint64_t y0 = y & 0x11111111u, y1 = y & 0x22222222u;
    uint64_t y2 = y & 0x44444444u, y3 = y & 0x88888888u;

    uint64_t z0 = (x0 * y0) ^ (x1 * y3) ^ (x2 * y2) ^ (x3 * y1);
    uint64_t z1 = (x0 * y1) ^ (x1 * y0) ^ (x2 * y3) ^ (x3 * y2);
    uint64_t z2 = (x0 * y2) ^ (x1 * y1) ^ (x2 * y0) ^ (x3 * y3);
    uint64_t z3 = (x0 * y3) ^ (x1 * y2) ^ (x2 * y1) ^ (x3 * y0);

#if ECC_PACKED_ADD_COMBINE
    *hi = disjointUnion4(uint32_t(z0 >> 32) & 0x11111111u, uint32_t(z1 >> 32) & 0x22222222u,
                         uint32_t(z2 >> 32) & 0x44444444u, uint32_t(z3 >> 32) & 0x88888888u);
    return disjointUnion4(uint32_t(z0) & 0x11111111u, uint32_t(z1) & 0x22222222u,
                          uint32_t(z2) & 0x44444444u, uint32_t(z3) & 0x88888888u);
#else
    uint64_t z = (z0 & m0) | (z1 & m1) | (z2 & m2) | (z3 & m3);
    *hi = (uint32_t)(z >> 32);
    return (uint32_t)z;
#endif
}

/* 2 x 2 words -> 4 words. The software path uses three clmul32 leaves;
   CLMAD retains both halves of the native 64-bit carryless product. */
ECC_HD void clmul64(uint32_t r[4], const uint32_t a[2], const uint32_t b[2]) {
#if ECC_USE_CLMAD_INSN
    const uint64_t aa = uint64_t(a[0]) | (uint64_t(a[1]) << 32);
    const uint64_t bb = uint64_t(b[0]) | (uint64_t(b[1]) << 32);
    uint64_t lo, hi;
    asm("clmad.lo.u64 %0, %1, %2, 0;" : "=l"(lo) : "l"(aa), "l"(bb));
    asm("clmad.hi.u64 %0, %1, %2, 0;" : "=l"(hi) : "l"(aa), "l"(bb));
    r[0] = uint32_t(lo); r[1] = uint32_t(lo >> 32);
    r[2] = uint32_t(hi); r[3] = uint32_t(hi >> 32);
#else
    uint32_t lo0, hi0, lo1, hi1, lom, him;
    lo0 = clmul32(a[0], b[0], &hi0);
    lo1 = clmul32(a[1], b[1], &hi1);
    lom = clmul32(a[0] ^ a[1], b[0] ^ b[1], &him);
    /* middle = (a0+a1)(b0+b1) - a0b0 - a1b1, XOR in characteristic 2 */
    lom ^= lo0 ^ lo1;
    him ^= hi0 ^ hi1;
    r[0] = lo0;
    r[1] = hi0 ^ lom;
    r[2] = lo1 ^ him;
    r[3] = hi1;
#endif
}

/* 4 x 4 words -> 8 words, Karatsuba again: 3 clmul64 = 9 clmul32. */
ECC_HD void clmul128(uint32_t r[8], const uint32_t a[4], const uint32_t b[4]) {
    uint32_t lo[4], hi[4], mid[4], as[2], bs[2];
    clmul64(lo, a, b);
    clmul64(hi, a + 2, b + 2);
    as[0] = a[0] ^ a[2]; as[1] = a[1] ^ a[3];
    bs[0] = b[0] ^ b[2]; bs[1] = b[1] ^ b[3];
    clmul64(mid, as, bs);
#pragma unroll
    for (int i = 0; i < 4; i++) mid[i] ^= lo[i] ^ hi[i];
    r[0] = lo[0]; r[1] = lo[1];
    r[2] = lo[2] ^ mid[0];
    r[3] = lo[3] ^ mid[1];
    r[4] = hi[0] ^ mid[2];
    r[5] = hi[1] ^ mid[3];
    r[6] = hi[2];
    r[7] = hi[3];
}

struct P131 { uint32_t v[5]; };

#if ECC_PACKED_TOP_CLMAD
/* The low 64 bits of x*y, plus an addend.  On the device this is one clmad
   whose third operand is free; the host emulates it through the software
   clmul64 so that testpacked.cpp checks the device formulation bit for bit. */
ECC_HD uint64_t clmadLo64(uint64_t x, uint64_t y, uint64_t add) {
#if defined(__CUDA_ARCH__)
    uint64_t r;
    asm("clmad.lo.u64 %0, %1, %2, %3;" : "=l"(r) : "l"(x), "l"(y), "l"(add));
    return r;
#else
    uint32_t r[4], a[2] = {uint32_t(x), uint32_t(x >> 32)}, b[2] = {uint32_t(y), uint32_t(y >> 32)};
    clmul64(r, a, b);
    return (r[0] | (uint64_t(r[1]) << 32)) ^ add;
#endif
}
/* x*y for x and y both below 8: three masked shifts, no multiply.  Used for
   the two bits a 3x64 product pushes past its limb (y is then below 4) and
   for a4*b4 itself. */
ECC_HD uint32_t clmulTiny(uint32_t x, uint32_t y) {
    return (x & (0u - (y & 1u))) ^ ((x << 1) & (0u - ((y >> 1) & 1u))) ^ ((x << 2) & (0u - ((y >> 2) & 1u)));
}
/* c[4..8] ^= a4*B_lo + A_lo*b4 + a4*b4*x^256, the terms the 3-bit top word
   adds to the 128x128 clmul128 product.

   Measured in SASS at the shipping preset, this correction is 65 of
   product131's 77 instructions; the whole four-word clmad product is 12.
   Written as masked shifts it is 6 masks, 12 masked ANDs and 10 funnel-shift
   accumulations per product, and ptxas already has it near that floor.

   So spend the other pipe.  a4*B_lo is two 3x64-bit products whose low 64
   bits each fit one clmad.lo, and clmad's addend is free, so the accumulate
   into c costs nothing: c[4..5] = a4*B0 + b4*A0 + c[4..5] is two instructions
   on the carryless unit and zero on the ALU.  What does not fit is the 2 bits
   each 3x64 product pushes past bit 63; those are (a4 * (B0 >> 62)) >> 2, a
   3x2-bit product.  Measured with the same compiler: product131 goes from 77
   ALU and 6 clmads to 52 and 10, so the residue is 40 against the 65 it
   replaces, and mulPolynomial131 from 164 to 128 because ptxas then fuses more
   of the reducer.  Per update at batch 16 the routines sass_cost.py weights
   fall 1621 -> 1534 ALU slots (-5.4%) for 29.1 -> 41.2 clmads (+41%).

   The price is paid on a unit the profile puts at 51.4%: scaled by 41.2/29.1
   that is 73% at today's rate and about 76% if the ALU cut converts, both
   under the 87% the ALU pipe already sustains, so the pipe model predicts the
   whole saving, about +4%.  What it cannot see is a product issuing its ten
   clmads back to back into a unit at 76% average.  A knob until a card
   decides, off by default; TOP-CLMAD.md has the tables. */
ECC_HD void topCrossClmad131(const P131 &a, const P131 &b, uint32_t *c) {
    const uint64_t A0 = a.v[0] | (uint64_t(a.v[1]) << 32), A1 = a.v[2] | (uint64_t(a.v[3]) << 32);
    const uint64_t B0 = b.v[0] | (uint64_t(b.v[1]) << 32), B1 = b.v[2] | (uint64_t(b.v[3]) << 32);
    const uint32_t a4 = a.v[4] & 7u, b4 = b.v[4] & 7u;
    // bits 64 and 65 of each 3x64 product, from the limb's top two bits
    const uint32_t h0 = (clmulTiny(a4, uint32_t(B0 >> 62)) ^ clmulTiny(b4, uint32_t(A0 >> 62))) >> 2;
    const uint32_t h1 = (clmulTiny(a4, uint32_t(B1 >> 62)) ^ clmulTiny(b4, uint32_t(A1 >> 62))) >> 2;
    uint64_t C2 = c[4] | (uint64_t(c[5]) << 32), C3 = c[6] | (uint64_t(c[7]) << 32);
    C2 = clmadLo64(a4, B0, clmadLo64(b4, A0, C2));
    C3 = clmadLo64(a4, B1, clmadLo64(b4, A1, C3 ^ h0));
    c[4] = uint32_t(C2); c[5] = uint32_t(C2 >> 32);
    c[6] = uint32_t(C3); c[7] = uint32_t(C3 >> 32);
    c[8] ^= h1 ^ clmulTiny(a4, b4);
}
#endif
ECC_HD uint32_t reverse32(uint32_t x) {
#ifdef __CUDA_ARCH__
    return __brev(x);
#else
    x=((x&0x55555555u)<<1)|((x>>1)&0x55555555u);
    x=((x&0x33333333u)<<2)|((x>>2)&0x33333333u);
    x=((x&0x0f0f0f0fu)<<4)|((x>>4)&0x0f0f0f0fu);
    x=((x&0x00ff00ffu)<<8)|((x>>8)&0x00ff00ffu);
    return (x<<16)|(x>>16);
#endif
}
ECC_HD P131 reverse131(const P131 &a) {
    P131 r;
#pragma unroll
    for(int i=0;i<5;i++) r.v[i]=(reverse32(a.v[4-i])>>29) | (i<4?reverse32(a.v[3-i])<<3:0);
    return r;
}
// Three limbs (64,64,3 bits): reuse the low diagonal products when
// reconstructing the top cross terms. Five full 64-bit products replace
// the old three products plus the bit-by-bit 128x3 tail.
ECC_HD void product131Karat3(const P131 &a, const P131 &b, uint32_t *c) {
    uint32_t lo[4], hi[4], mid[4], aa[2], bb[2];
    clmul64(lo, a.v, b.v);
    clmul64(hi, a.v + 2, b.v + 2);
    const uint32_t a2 = a.v[4], b2 = b.v[4];
    const uint32_t top = (b2 & (0u - (a2 & 1u)))
        ^ ((b2 & (0u - ((a2 >> 1) & 1u))) << 1)
        ^ ((b2 & (0u - ((a2 >> 2) & 1u))) << 2);
#pragma unroll
    for (int i = 0; i < 4; ++i) { c[i] = lo[i]; c[i + 4] = hi[i]; }
    c[8] = top;
    aa[0] = a.v[0] ^ a.v[2]; aa[1] = a.v[1] ^ a.v[3];
    bb[0] = b.v[0] ^ b.v[2]; bb[1] = b.v[1] ^ b.v[3];
    clmul64(mid, aa, bb);
#pragma unroll
    for (int i = 0; i < 4; ++i) c[i + 2] ^= mid[i] ^ lo[i] ^ hi[i];
    aa[0] = a.v[0] ^ a2; aa[1] = a.v[1];
    bb[0] = b.v[0] ^ b2; bb[1] = b.v[1];
    clmul64(mid, aa, bb);
#pragma unroll
    for (int i = 0; i < 4; ++i) c[i + 4] ^= mid[i] ^ lo[i] ^ (i == 0 ? top : 0u);
    aa[0] = a.v[2] ^ a2; aa[1] = a.v[3];
    bb[0] = b.v[2] ^ b2; bb[1] = b.v[3];
    clmul64(mid, aa, bb);
    c[6] ^= mid[0] ^ hi[0] ^ top;
    c[7] ^= mid[1] ^ hi[1];
    c[8] ^= mid[2] ^ hi[2];
    // The omitted tenth word is zero: a 64x3 cross term has degree <=65.
}
ECC_HD void product131(const P131 &a,const P131 &b,uint32_t *c) {
#if ECC_PACKED_KARAT3
    product131Karat3(a, b, c);
#else
    clmul128(c,a.v,b.v); c[8]=0;
#if ECC_PACKED_TOP_CLMAD
    topCrossClmad131(a,b,c);
#else
#pragma unroll
    for(int k=0;k<3;k++) {
        uint32_t ma=0u-((a.v[4]>>k)&1u), mb=0u-((b.v[4]>>k)&1u);
#pragma unroll
        for(int i=0;i<4;i++) {
            uint32_t t=(a.v[i]&mb)^(b.v[i]&ma);
            c[4+i]^=t<<k;
            if(k) c[5+i]^=t>>(32-k);
        }
        c[8]^=(b.v[4]&ma)<<k;
    }
#endif
#endif
}
ECC_HD P131 add131(const P131 &a,const P131 &b) {
    P131 r;
#pragma unroll
    for(int i=0;i<5;i++) r.v[i]=a.v[i]^b.v[i];
    return r;
}
// The generated linear transforms convert to/from the polynomial basis used
// by codegen/build.py, allowing a single product while retaining ONB storage.
#include "packedtransform131.h"
#ifndef ECC_PACKED_DIRECT_REDUCE
#define ECC_PACKED_DIRECT_REDUCE 0
#endif
#if ECC_PACKED_DIRECT_REDUCE != 0 && ECC_PACKED_DIRECT_REDUCE != 1
#error "ECC_PACKED_DIRECT_REDUCE must be 0 or 1"
#endif
#if ECC_PACKED_DIRECT_REDUCE
#include "packeddirectreduce131.h"
#else
#include "packedpolyreduce131.h"
#endif
#ifndef ECC_PACKED_GENERATED_PRODUCT
#define ECC_PACKED_GENERATED_PRODUCT 0
#endif
#if ECC_PACKED_GENERATED_PRODUCT != 0 && ECC_PACKED_GENERATED_PRODUCT != 1
#error "ECC_PACKED_GENERATED_PRODUCT must be 0 or 1"
#endif
#if ECC_PACKED_GENERATED_PRODUCT && !ECC_PACKED_DIRECT_REDUCE
#error "ECC_PACKED_GENERATED_PRODUCT requires ECC_PACKED_DIRECT_REDUCE"
#endif
#if ECC_PACKED_GENERATED_PRODUCT
#include "packedgeneratedproduct131.h"
#endif
ECC_HD P131 fromPolynomial131(const P131 &a) {
    const uint32_t h[9]={a.v[0],a.v[1],a.v[2],a.v[3],a.v[4],0,0,0,0};
    return fromPolynomialProduct131(h);
}
static ECC_BIG P131 mulPolynomial131(P131 a, P131 b) {
// Native carryless products supersede the generated software multiplier.
#if ECC_PACKED_GENERATED_PRODUCT && !ECC_PACKED_CLMAD
    return generatedProduct131(a,b);
#else
    uint32_t h[9]; product131(a,b,h);
    return reducePolynomial131(h);
#endif
}
struct PolynomialPair { P131 first,second; };
static ECC_BIG PolynomialPair mulPolynomialPair131(P131 a,P131 b,P131 c) {
#if ECC_PACKED_GENERATED_PRODUCT && !ECC_PACKED_CLMAD
    P131 first=generatedProduct131(a,b);
    return PolynomialPair{first,generatedProduct131(a,c)};
#else
    uint32_t h[9];
    product131(a,b,h);
    P131 first=reducePolynomial131(h);
    product131(a,c,h);
    return PolynomialPair{first,reducePolynomial131(h)};
#endif
}
#ifndef ECC_PACKED_SINGLE_PRODUCT
#define ECC_PACKED_SINGLE_PRODUCT 0
#endif
#if ECC_PACKED_SINGLE_PRODUCT != 0 && ECC_PACKED_SINGLE_PRODUCT != 1
#error "ECC_PACKED_SINGLE_PRODUCT must be 0 or 1"
#endif
#ifndef ECC_PACKED_BY_VALUE
#define ECC_PACKED_BY_VALUE 0
#endif
#if ECC_PACKED_BY_VALUE != 0 && ECC_PACKED_BY_VALUE != 1
#error "ECC_PACKED_BY_VALUE must be 0 or 1"
#endif
#if ECC_PACKED_BY_VALUE
// Passing these five-word aggregates by value lets the device ABI use
// registers instead of materializing the caller's operands in local memory.
using MulArg = P131;
#else
using MulArg = const P131 &;
#endif
static ECC_BIG P131 mul131(MulArg a, MulArg b) {
#if ECC_PACKED_SINGLE_PRODUCT
    const P131 pa = toPolynomial131(a), pb = toPolynomial131(b);
    uint32_t h[9];
    product131(pa,pb,h);
    return fromPolynomialProduct131(h);
#else
    // gamma_i gamma_j = gamma_(i+j) + gamma_(i-j), gamma_0=0,
    // gamma_k=gamma_(263-k): the original two-product multiplier.
    uint32_t c[9],d[9];
    P131 rb=reverse131(b),r;
    product131(a,b,c); product131(a,rb,d);
#pragma unroll
    for(int i=0;i<5;i++) {
        uint32_t low=(c[i]<<1)|(i?c[i-1]>>31:0);
        uint32_t high=(reverse32(c[8-i])>>27)|(reverse32(c[7-i])<<5);
        uint32_t corr=(d[4+i]>>3)|(i<4?d[5+i]<<29:0);
        uint32_t rev=(reverse32(d[4-i])>>30)|(i<4?reverse32(d[3-i])<<2:0);
        r.v[i]=low^high^corr^rev;
    }
    r.v[4]&=7;
    return r;
#endif
}
#ifndef ECC_PACKED_ALU_SQUARE
#define ECC_PACKED_ALU_SQUARE 0
#endif
#if ECC_PACKED_ALU_SQUARE != 0 && ECC_PACKED_ALU_SQUARE != 1
#error "ECC_PACKED_ALU_SQUARE must be 0 or 1"
#endif
// The logic-op spread, kept callable beside the clmad one: with the table walk
// the ALU pipe is no longer the only saturated one (ITERATION-FUNCTION.md
// section 6), so the per-update polynomial squaring can be moved back here
// with ECC_PACKED_ALU_SQUARE=1, trading five clmad for about sixty logic ops.
ECC_HD uint64_t spread32alu(uint32_t x){
 uint64_t r=x;
 r=(r|(r<<16))&0x0000ffff0000ffffull;
 r=(r|(r<<8))&0x00ff00ff00ff00ffull;
 r=(r|(r<<4))&0x0f0f0f0f0f0f0f0full;
 r=(r|(r<<2))&0x3333333333333333ull;
 r=(r|(r<<1))&0x5555555555555555ull;
 return r;
}
ECC_HD uint64_t spread32p(uint32_t x){
#if ECC_USE_CLMAD_INSN
 /* A carryless square is a bit spread.  x*x = sum_(i,j) x_i x_j t^(i+j) and
    every i != j term appears twice, so in characteristic two only the doubled
    positions survive, which is exactly what the stages below build.  Degree
    2*31 fits the low half, so clmad.hi is not needed.

    One clmad is dearer than the ~22 logic ops it replaces -- the price
    benchmarks/clmad-price measures is 37.7 -- and is still the right trade
    because the two pipes are not equally loaded.  Nsight Compute puts the ALU
    pipe at 87.3% here and the FP64 pipe that carries clmad at 51.4%, so the
    walk is short of ALU and has carryless capacity to spend.  The same
    direction holds on Ada: benchmarks/preblackwell measured CLMAD=1 at
    1.811× the software product on an L40S and 1.881× on an L4. */
 uint64_t r;
 asm("clmad.lo.u64 %0, %1, %1, 0;" : "=l"(r) : "l"((uint64_t)x));
 return r;
#else
 uint64_t r=x;
 r=(r|(r<<16))&0x0000ffff0000ffffull;
 r=(r|(r<<8))&0x00ff00ff00ff00ffull;
 r=(r|(r<<4))&0x0f0f0f0f0f0f0f0full;
 r=(r|(r<<2))&0x3333333333333333ull;
 r=(r|(r<<1))&0x5555555555555555ull;
 return r;
#endif
}
// Polynomial coefficients square into the even positions of a degree-260
// product. This is distinct from sqr131's normal-basis permutation.
ECC_HD P131 squarePolynomial131(P131 a) {
    uint32_t h[9];
#pragma unroll
    for (int i = 0; i < 4; ++i) {
        const uint64_t w = ECC_PACKED_ALU_SQUARE ? spread32alu(a.v[i]) : spread32p(a.v[i]);
        h[2 * i] = uint32_t(w);
        h[2 * i + 1] = uint32_t(w >> 32);
    }
    h[8] = uint32_t(ECC_PACKED_ALU_SQUARE ? spread32alu(a.v[4]) : spread32p(a.v[4]));
    return reducePolynomial131(h);
}
ECC_HD P131 sqr131(const P131 &a){
 P131 rev=reverse131(a),r;
 uint64_t lo=(spread32p(a.v[0])<<1)^spread32p(rev.v[0]);
 uint64_t hi=(spread32p(a.v[1])<<1)^spread32p(rev.v[1]);
 r.v[0]=uint32_t(lo);r.v[1]=uint32_t(lo>>32);
 r.v[2]=uint32_t(hi);r.v[3]=uint32_t(hi>>32);
 r.v[4]=(rev.v[2]&1u)|((a.v[2]&1u)<<1)|((rev.v[2]&2u)<<1);
 return r;
}
#ifndef ECC_PACKED_PERM_SIGMA
#define ECC_PACKED_PERM_SIGMA 0
#endif
#if ECC_PACKED_PERM_SIGMA < 0 || ECC_PACKED_PERM_SIGMA > 3
#error "ECC_PACKED_PERM_SIGMA must be a bit mask from 0 to 3"
#endif
#if ECC_PACKED_SHARED_SIGMA && !(ECC_PACKED_PERM_SIGMA & 1)
#error "ECC_PACKED_SHARED_SIGMA requires the walk permutation network"
#endif
#if ECC_PACKED_PERM_SIGMA
#include "packedsigma131.h"
#endif
ECC_HD P131 sigma131(P131 a,int k){
#if ECC_PACKED_PERM_SIGMA & 1
 if(k>=3 && k<=10) return sigmaWalkNetwork131(a,k-3);
#endif
#if ECC_PACKED_PERM_SIGMA & 2
 if(k==16 || k==32 || k==65) return sigmaInvNetwork131(a,k==16?0:k==32?1:2);
#endif
#pragma unroll 1
 for(int i=0;i<k;i++)a=sqr131(a);
 return a;
}
#ifndef ECC_PACKED_UNROLL_INV
#define ECC_PACKED_UNROLL_INV 0
#endif
ECC_HD P131 inv131(P131 a){
#if ECC_PACKED_UNROLL_INV
 // The same Itoh–Tsujii chain with explicit powers: beta_2,4,8,16,32,64,65,130.
 P131 acc=mul131(sqr131(a),a);
 acc=mul131(sqr131(sqr131(acc)),acc);
 acc=mul131(sigma131(acc,4),acc);
 acc=mul131(sigma131(acc,8),acc);
 acc=mul131(sigma131(acc,16),acc);
 acc=mul131(sigma131(acc,32),acc);
 acc=mul131(sqr131(acc),a);
 acc=mul131(sigma131(acc,65),acc);
 return sqr131(acc);
#else
 P131 acc=a;int k=1;
#pragma unroll 1
 for(int bit=6;bit>=0;--bit){
  acc=mul131(sigma131(acc,k),acc);k*=2;
  if((130>>bit)&1){acc=mul131(sqr131(acc),a);k++;}
 }
 return sqr131(acc);
#endif
}

} // namespace eccPacked131
