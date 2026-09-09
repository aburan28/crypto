// Packed GF(2^131) in the same permuted type-II normal basis as eccF131.
// Each element occupies five uint32_t words. Squaring is a bit permutation;
// multiplication uses gamma_i*gamma_j = gamma_(i+j) + gamma_(i-j).
#pragma once
#include "bitslice.h"
namespace eccPacked131 {
// Integer-mask carryless primitives adapted from gpu/ecc2k/f2m.cuh.
ECC_HD uint32_t clmul32(uint32_t x, uint32_t y, uint32_t *hi) {
    const uint64_t m0 = 0x1111111111111111ull;
    const uint64_t m1 = 0x2222222222222222ull;
    const uint64_t m2 = 0x4444444444444444ull;
    const uint64_t m3 = 0x8888888888888888ull;

    uint64_t x0 = x & 0x11111111u, x1 = x & 0x22222222u;
    uint64_t x2 = x & 0x44444444u, x3 = x & 0x88888888u;
    uint64_t y0 = y & 0x11111111u, y1 = y & 0x22222222u;
    uint64_t y2 = y & 0x44444444u, y3 = y & 0x88888888u;

    uint64_t z0 = (x0 * y0) ^ (x1 * y3) ^ (x2 * y2) ^ (x3 * y1);
    uint64_t z1 = (x0 * y1) ^ (x1 * y0) ^ (x2 * y3) ^ (x3 * y2);
    uint64_t z2 = (x0 * y2) ^ (x1 * y1) ^ (x2 * y0) ^ (x3 * y3);
    uint64_t z3 = (x0 * y3) ^ (x1 * y2) ^ (x2 * y1) ^ (x3 * y0);

    uint64_t z = (z0 & m0) | (z1 & m1) | (z2 & m2) | (z3 & m3);
    *hi = (uint32_t)(z >> 32);
    return (uint32_t)z;
}

/* 2 x 2 words -> 4 words, Karatsuba (3 clmul32 instead of 4). */
ECC_HD void clmul64(uint32_t r[4], const uint32_t a[2], const uint32_t b[2]) {
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
ECC_HD void product131(const P131 &a,const P131 &b,uint32_t *c) {
    clmul128(c,a.v,b.v); c[8]=0;
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
#ifndef ECC_PACKED_SINGLE_PRODUCT
#define ECC_PACKED_SINGLE_PRODUCT 0
#endif
#if ECC_PACKED_SINGLE_PRODUCT != 0 && ECC_PACKED_SINGLE_PRODUCT != 1
#error "ECC_PACKED_SINGLE_PRODUCT must be 0 or 1"
#endif
static ECC_BIG P131 mul131(const P131 &a,const P131 &b) {
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
ECC_HD uint64_t spread32p(uint32_t x){
 uint64_t r=x;
 r=(r|(r<<16))&0x0000ffff0000ffffull;
 r=(r|(r<<8))&0x00ff00ff00ff00ffull;
 r=(r|(r<<4))&0x0f0f0f0f0f0f0f0full;
 r=(r|(r<<2))&0x3333333333333333ull;
 r=(r|(r<<1))&0x5555555555555555ull;
 return r;
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
ECC_HD P131 sigma131(P131 a,int k){
#pragma unroll 1
 for(int i=0;i<k;i++)a=sqr131(a);
 return a;
}
ECC_HD P131 inv131(P131 a){
 P131 acc=a;int k=1;
#pragma unroll 1
 for(int bit=6;bit>=0;--bit){
  acc=mul131(sigma131(acc,k),acc);k*=2;
  if((130>>bit)&1){acc=mul131(sqr131(acc),a);k++;}
 }
 return sqr131(acc);
}

} // namespace eccPacked131
