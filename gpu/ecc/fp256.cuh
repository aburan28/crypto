/* fp256.cuh -- 256-bit prime-field arithmetic for CUDA (and host emulation).
 *
 * Representation: 8 x 32-bit little-endian limbs.  32-bit limbs are the
 * right choice on every NVIDIA architecture through Blackwell: the integer
 * datapath is 32-bit (IMAD / IMAD.WIDE.U32), 64-bit multiplies are emulated
 * from 32-bit ones, and the extended-precision carry chain
 * (mad.lo.cc / madc.hi.cc / addc) is exposed at the 32-bit granularity only.
 *
 * Two reduction strategies, selected per modulus descriptor:
 *
 *   FAST = true   secp256k1 only.  Elements are canonical integers < p and a
 *                 512-bit product is folded with 2^256 = 2^32 + 977 (mod p).
 *                 ~40% fewer multiply-adds than Montgomery.
 *
 *   FAST = false  Any odd modulus < 2^256.  Elements are in Montgomery form
 *                 (x * 2^256 mod p) and multiplication is CIOS Montgomery.
 *                 Used for generic curves, toy curves, and the scalar ring
 *                 mod n on the host.
 *
 * Every function is `__host__ __device__` so the *same* code runs in the
 * CPU test harness (compiled by g++/clang++ without CUDA) and in kernels.
 * The portable path uses 64-bit accumulators; nvcc turns
 * `(uint64_t)a*b + c + d` into IMAD.WIDE.U32 + carry handling and the result
 * is within ~15% of hand-written PTX.  Define FP_PTX=1 to use the inline-PTX
 * carry chains instead (device only; the host still uses the portable path).
 *
 * Nothing here is constant-time with respect to *control flow* except where
 * noted -- this is a cryptanalysis / throughput library, not a signing
 * library.  Data-dependent branches are marked.
 */
#ifndef GPU_ECC_FP256_CUH
#define GPU_ECC_FP256_CUH

#include <stdint.h>
#include <string.h>

#ifndef GPU_ECC_CURVE_HEADER
#define GPU_ECC_CURVE_HEADER "curve_secp256k1.h"
#endif
#include GPU_ECC_CURVE_HEADER

/* FP_HD marks the small primitives that must be inlined into the caller's
 * register allocation (add, sub, mul, limb shuffling).  FP_BIG marks the
 * heavyweight routines -- Fermat inversion, scalar multiplication, walk
 * seeding -- which must NOT be force-inlined: each one expands to a few
 * hundred field multiplications, and inlining them at every call site turns
 * a two-second compile into a ten-minute one and blows out the instruction
 * cache.  Letting the compiler emit them as real device functions costs a
 * call and gains everything else. */
#ifdef __CUDACC__
#define FP_HD __host__ __device__ __forceinline__
#define FP_BIG __host__ __device__ __noinline__
#define FP_DEV __device__ __forceinline__
#else
#define FP_HD inline
#define FP_BIG inline
#define FP_DEV inline
#endif

#ifndef FP_PTX
#define FP_PTX 0
#endif
/* Use the secp256k1 special reduction for the base field when the curve is
 * secp256k1.  Override with -DFP_FAST=0 to force Montgomery (used by the
 * tests to exercise both paths on the same vectors). */
#ifndef FP_FAST
#define FP_FAST CURVE_IS_SECP256K1
#endif

#define FP_LIMBS 8

struct fp256 {
    uint32_t v[FP_LIMBS];
};

/* ------------------------------------------------------------------------
 * Modulus descriptors.  Limbs are returned through functions that build a
 * local const array: after `#pragma unroll` every index is a compile-time
 * constant and the values become immediates -- no constant-memory loads,
 * and the same source works for host and device compilation.
 * ---------------------------------------------------------------------- */
#define FP_DEFINE_LIMB_FN(name, INIT) \
    static FP_HD uint32_t name(int i) { const uint32_t t[FP_LIMBS] = INIT; return t[i]; }

struct ModP {
    FP_DEFINE_LIMB_FN(limb, FP_P_LIMBS)
    FP_DEFINE_LIMB_FN(r1, FP_R1_LIMBS)
    FP_DEFINE_LIMB_FN(r2, FP_R2_LIMBS)
    FP_DEFINE_LIMB_FN(half, FP_HALF_LIMBS)
    FP_DEFINE_LIMB_FN(pm2, FP_PM2_LIMBS)
    static FP_HD uint32_t nprime() { return FP_NPRIME; }
    static FP_HD int bits() { return CURVE_BITS; }
};

struct ModN {
    FP_DEFINE_LIMB_FN(limb, CURVE_N_LIMBS)
    FP_DEFINE_LIMB_FN(r1, N_R1_LIMBS)
    FP_DEFINE_LIMB_FN(r2, N_R2_LIMBS)
    FP_DEFINE_LIMB_FN(half, N_HALF_LIMBS)
    FP_DEFINE_LIMB_FN(pm2, N_PM2_LIMBS)
    static FP_HD uint32_t nprime() { return N_NPRIME; }
    static FP_HD int bits() { return N_BITS; }
};

/* ------------------------------------------------------------------------
 * Multi-precision primitives (portable).
 * ---------------------------------------------------------------------- */

/* r = a + b, returns carry. */
FP_HD uint32_t mp_add(uint32_t r[8], const uint32_t a[8], const uint32_t b[8]) {
#if FP_PTX && defined(__CUDA_ARCH__)
    uint32_t c;
    asm("add.cc.u32  %0, %9,  %17;\n\t"
        "addc.cc.u32 %1, %10, %18;\n\t"
        "addc.cc.u32 %2, %11, %19;\n\t"
        "addc.cc.u32 %3, %12, %20;\n\t"
        "addc.cc.u32 %4, %13, %21;\n\t"
        "addc.cc.u32 %5, %14, %22;\n\t"
        "addc.cc.u32 %6, %15, %23;\n\t"
        "addc.cc.u32 %7, %16, %24;\n\t"
        "addc.u32    %8, 0, 0;"
        : "=&r"(r[0]), "=&r"(r[1]), "=&r"(r[2]), "=&r"(r[3]),
          "=&r"(r[4]), "=&r"(r[5]), "=&r"(r[6]), "=&r"(r[7]), "=&r"(c)
        : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
          "r"(a[4]), "r"(a[5]), "r"(a[6]), "r"(a[7]),
          "r"(b[0]), "r"(b[1]), "r"(b[2]), "r"(b[3]),
          "r"(b[4]), "r"(b[5]), "r"(b[6]), "r"(b[7]));
    return c;
#else
    uint64_t c = 0;
#pragma unroll
    for (int j = 0; j < 8; j++) {
        c += (uint64_t)a[j] + b[j];
        r[j] = (uint32_t)c;
        c >>= 32;
    }
    return (uint32_t)c;
#endif
}

/* r = a - b, returns borrow (1 if a < b). */
FP_HD uint32_t mp_sub(uint32_t r[8], const uint32_t a[8], const uint32_t b[8]) {
#if FP_PTX && defined(__CUDA_ARCH__)
    uint32_t bw;
    asm("sub.cc.u32  %0, %9,  %17;\n\t"
        "subc.cc.u32 %1, %10, %18;\n\t"
        "subc.cc.u32 %2, %11, %19;\n\t"
        "subc.cc.u32 %3, %12, %20;\n\t"
        "subc.cc.u32 %4, %13, %21;\n\t"
        "subc.cc.u32 %5, %14, %22;\n\t"
        "subc.cc.u32 %6, %15, %23;\n\t"
        "subc.cc.u32 %7, %16, %24;\n\t"
        "subc.u32    %8, 0, 0;"
        : "=&r"(r[0]), "=&r"(r[1]), "=&r"(r[2]), "=&r"(r[3]),
          "=&r"(r[4]), "=&r"(r[5]), "=&r"(r[6]), "=&r"(r[7]), "=&r"(bw)
        : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
          "r"(a[4]), "r"(a[5]), "r"(a[6]), "r"(a[7]),
          "r"(b[0]), "r"(b[1]), "r"(b[2]), "r"(b[3]),
          "r"(b[4]), "r"(b[5]), "r"(b[6]), "r"(b[7]));
    /* subc.u32 d, 0, 0 yields 0 - 0 - borrow = 0xFFFFFFFF when a borrow is
     * pending, so normalise to 0/1. */
    return bw & 1u;
#else
    uint64_t bw = 0;
#pragma unroll
    for (int j = 0; j < 8; j++) {
        uint64_t d = (uint64_t)a[j] - b[j] - bw;
        r[j] = (uint32_t)d;
        bw = d >> 63;
    }
    return (uint32_t)bw;
#endif
}

/* Constant-time select: r = flag ? a : r  (flag in {0,1}). */
FP_HD void mp_cmov(uint32_t r[8], const uint32_t a[8], uint32_t flag) {
    uint32_t mask = 0u - (flag & 1u);
#pragma unroll
    for (int j = 0; j < 8; j++) r[j] = (r[j] & ~mask) | (a[j] & mask);
}

FP_HD int mp_is_zero(const uint32_t a[8]) {
    uint32_t acc = 0;
#pragma unroll
    for (int j = 0; j < 8; j++) acc |= a[j];
    return acc == 0;
}

FP_HD int mp_eq(const uint32_t a[8], const uint32_t b[8]) {
    uint32_t acc = 0;
#pragma unroll
    for (int j = 0; j < 8; j++) acc |= a[j] ^ b[j];
    return acc == 0;
}

/* t[0..9] += a[0..7] * b.  The row primitive of CIOS Montgomery
 * multiplication, where the accumulator really can be nine limbs wide at
 * row entry.  The caller guarantees the sum fits in the 10-limb window
 * (see bounds in mont_mul). */
FP_HD void mp_mac_row(uint32_t t[10], const uint32_t a[8], uint32_t b) {
#if FP_PTX && defined(__CUDA_ARCH__)
    /* Two carry chains: the low halves of the products, then the high
     * halves shifted up one limb.  Each chain is a single asm statement
     * because the PTX carry flag is not preserved across statements. */
    asm("mad.lo.cc.u32  %0, %10, %18, %0;\n\t"
        "madc.lo.cc.u32 %1, %11, %18, %1;\n\t"
        "madc.lo.cc.u32 %2, %12, %18, %2;\n\t"
        "madc.lo.cc.u32 %3, %13, %18, %3;\n\t"
        "madc.lo.cc.u32 %4, %14, %18, %4;\n\t"
        "madc.lo.cc.u32 %5, %15, %18, %5;\n\t"
        "madc.lo.cc.u32 %6, %16, %18, %6;\n\t"
        "madc.lo.cc.u32 %7, %17, %18, %7;\n\t"
        "addc.cc.u32    %8, %8, 0;\n\t"
        "addc.u32       %9, %9, 0;"
        : "+r"(t[0]), "+r"(t[1]), "+r"(t[2]), "+r"(t[3]), "+r"(t[4]),
          "+r"(t[5]), "+r"(t[6]), "+r"(t[7]), "+r"(t[8]), "+r"(t[9])
        : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
          "r"(a[4]), "r"(a[5]), "r"(a[6]), "r"(a[7]), "r"(b));
    asm("mad.hi.cc.u32  %1, %10, %18, %1;\n\t"
        "madc.hi.cc.u32 %2, %11, %18, %2;\n\t"
        "madc.hi.cc.u32 %3, %12, %18, %3;\n\t"
        "madc.hi.cc.u32 %4, %13, %18, %4;\n\t"
        "madc.hi.cc.u32 %5, %14, %18, %5;\n\t"
        "madc.hi.cc.u32 %6, %15, %18, %6;\n\t"
        "madc.hi.cc.u32 %7, %16, %18, %7;\n\t"
        "madc.hi.cc.u32 %8, %17, %18, %8;\n\t"
        "addc.u32       %9, %9, 0;"
        : "+r"(t[0]), "+r"(t[1]), "+r"(t[2]), "+r"(t[3]), "+r"(t[4]),
          "+r"(t[5]), "+r"(t[6]), "+r"(t[7]), "+r"(t[8]), "+r"(t[9])
        : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
          "r"(a[4]), "r"(a[5]), "r"(a[6]), "r"(a[7]), "r"(b));
#else
    uint64_t c = 0;
#pragma unroll
    for (int j = 0; j < 8; j++) {
        /* max: (2^32-1) + (2^32-1)^2 + (2^32-1) = 2^64 - 1, no overflow */
        c += (uint64_t)t[j] + (uint64_t)a[j] * b;
        t[j] = (uint32_t)c;
        c >>= 32;
    }
    c += t[8];
    t[8] = (uint32_t)c;
    t[9] += (uint32_t)(c >> 32);
#endif
}

/* t[0..8] = a[0..7] * b.  The first row of a schoolbook product has nothing
 * to accumulate onto, so its low half is plain multiplies and it needs no
 * zeroed accumulator to start from. */
FP_HD void mp_mul_row0(uint32_t t[9], const uint32_t a[8], uint32_t b) {
#if FP_PTX && defined(__CUDA_ARCH__)
    asm("mul.lo.u32 %0, %8,  %16;\n\t"
        "mul.lo.u32 %1, %9,  %16;\n\t"
        "mul.lo.u32 %2, %10, %16;\n\t"
        "mul.lo.u32 %3, %11, %16;\n\t"
        "mul.lo.u32 %4, %12, %16;\n\t"
        "mul.lo.u32 %5, %13, %16;\n\t"
        "mul.lo.u32 %6, %14, %16;\n\t"
        "mul.lo.u32 %7, %15, %16;"
        : "=&r"(t[0]), "=&r"(t[1]), "=&r"(t[2]), "=&r"(t[3]),
          "=&r"(t[4]), "=&r"(t[5]), "=&r"(t[6]), "=&r"(t[7])
        : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
          "r"(a[4]), "r"(a[5]), "r"(a[6]), "r"(a[7]), "r"(b));
    asm("mad.hi.cc.u32  %0, %8,  %16, %0;\n\t"
        "madc.hi.cc.u32 %1, %9,  %16, %1;\n\t"
        "madc.hi.cc.u32 %2, %10, %16, %2;\n\t"
        "madc.hi.cc.u32 %3, %11, %16, %3;\n\t"
        "madc.hi.cc.u32 %4, %12, %16, %4;\n\t"
        "madc.hi.cc.u32 %5, %13, %16, %5;\n\t"
        "madc.hi.cc.u32 %6, %14, %16, %6;\n\t"
        "madc.hi.u32    %7, %15, %16, 0;"
        : "+r"(t[1]), "+r"(t[2]), "+r"(t[3]), "+r"(t[4]),
          "+r"(t[5]), "+r"(t[6]), "+r"(t[7]), "=&r"(t[8])
        : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
          "r"(a[4]), "r"(a[5]), "r"(a[6]), "r"(a[7]), "r"(b));
#else
    uint64_t c = 0;
#pragma unroll
    for (int j = 0; j < 8; j++) {
        c += (uint64_t)a[j] * b;
        t[j] = (uint32_t)c;
        c >>= 32;
    }
    t[8] = (uint32_t)c;
#endif
}

/* t[0..7] += a[0..7] * b, with the carry out *written* to t[8] rather than
 * accumulated into it.
 *
 * The row never carries past limb 8: the largest value it can produce is
 * (2^256 - 1) + (2^256 - 1)(2^32 - 1) = (2^256 - 1) * 2^32 < 2^288. So a
 * schoolbook product needs no pre-zeroed accumulator above the row window
 * either, which is what separates this from `mp_mac_row`. */
FP_HD void mp_mac_row9(uint32_t t[9], const uint32_t a[8], uint32_t b) {
#if FP_PTX && defined(__CUDA_ARCH__)
    asm("mad.lo.cc.u32  %0, %9,  %17, %0;\n\t"
        "madc.lo.cc.u32 %1, %10, %17, %1;\n\t"
        "madc.lo.cc.u32 %2, %11, %17, %2;\n\t"
        "madc.lo.cc.u32 %3, %12, %17, %3;\n\t"
        "madc.lo.cc.u32 %4, %13, %17, %4;\n\t"
        "madc.lo.cc.u32 %5, %14, %17, %5;\n\t"
        "madc.lo.cc.u32 %6, %15, %17, %6;\n\t"
        "madc.lo.cc.u32 %7, %16, %17, %7;\n\t"
        "addc.u32       %8, 0, 0;"
        : "+r"(t[0]), "+r"(t[1]), "+r"(t[2]), "+r"(t[3]),
          "+r"(t[4]), "+r"(t[5]), "+r"(t[6]), "+r"(t[7]), "=&r"(t[8])
        : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
          "r"(a[4]), "r"(a[5]), "r"(a[6]), "r"(a[7]), "r"(b));
    asm("mad.hi.cc.u32  %0, %8,  %16, %0;\n\t"
        "madc.hi.cc.u32 %1, %9,  %16, %1;\n\t"
        "madc.hi.cc.u32 %2, %10, %16, %2;\n\t"
        "madc.hi.cc.u32 %3, %11, %16, %3;\n\t"
        "madc.hi.cc.u32 %4, %12, %16, %4;\n\t"
        "madc.hi.cc.u32 %5, %13, %16, %5;\n\t"
        "madc.hi.cc.u32 %6, %14, %16, %6;\n\t"
        "madc.hi.u32    %7, %15, %16, %7;"
        : "+r"(t[1]), "+r"(t[2]), "+r"(t[3]), "+r"(t[4]),
          "+r"(t[5]), "+r"(t[6]), "+r"(t[7]), "+r"(t[8])
        : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
          "r"(a[4]), "r"(a[5]), "r"(a[6]), "r"(a[7]), "r"(b));
#else
    uint64_t c = 0;
#pragma unroll
    for (int j = 0; j < 8; j++) {
        c += (uint64_t)t[j] + (uint64_t)a[j] * b;
        t[j] = (uint32_t)c;
        c >>= 32;
    }
    t[8] = (uint32_t)c;
#endif
}

/* t[1..9] += a[0..7], i.e. t += a << 32.
 *
 * The secp256k1 fold needs `hi << 32` added to a running 10-word value, and
 * that shifted add is pure carry propagation -- exactly the shape the
 * portable path is worst at, since it has to rebuild the carry out of a
 * 64-bit accumulator each limb. */
FP_HD void mp_add_shift32(uint32_t t[10], const uint32_t a[8]) {
#if FP_PTX && defined(__CUDA_ARCH__)
    asm("add.cc.u32  %0, %0, %9;\n\t"
        "addc.cc.u32 %1, %1, %10;\n\t"
        "addc.cc.u32 %2, %2, %11;\n\t"
        "addc.cc.u32 %3, %3, %12;\n\t"
        "addc.cc.u32 %4, %4, %13;\n\t"
        "addc.cc.u32 %5, %5, %14;\n\t"
        "addc.cc.u32 %6, %6, %15;\n\t"
        "addc.cc.u32 %7, %7, %16;\n\t"
        "addc.u32    %8, %8, 0;"
        : "+r"(t[1]), "+r"(t[2]), "+r"(t[3]), "+r"(t[4]),
          "+r"(t[5]), "+r"(t[6]), "+r"(t[7]), "+r"(t[8]), "+r"(t[9])
        : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]),
          "r"(a[4]), "r"(a[5]), "r"(a[6]), "r"(a[7]));
#else
    uint64_t c = 0;
#pragma unroll
    for (int j = 0; j < 8; j++) {
        c += (uint64_t)t[j + 1] + a[j];
        t[j + 1] = (uint32_t)c;
        c >>= 32;
    }
    t[9] += (uint32_t)c;
#endif
}

/* t[0..7] += (a0, a1, a2) placed at limbs 0, 1, 2.  Returns the carry out
 * of limb 7.
 *
 * The later secp256k1 folds all have this shape: a value below 2^66 added
 * to the low half, then carry rippling to the top.  Writing them this way
 * rather than interleaving the multiply into the carry loop turns each
 * fold into one straight carry chain. */
FP_HD uint32_t mp_add_small(uint32_t t[8], uint32_t a0, uint32_t a1, uint32_t a2) {
#if FP_PTX && defined(__CUDA_ARCH__)
    uint32_t c;
    asm("add.cc.u32  %0, %0, %9;\n\t"
        "addc.cc.u32 %1, %1, %10;\n\t"
        "addc.cc.u32 %2, %2, %11;\n\t"
        "addc.cc.u32 %3, %3, 0;\n\t"
        "addc.cc.u32 %4, %4, 0;\n\t"
        "addc.cc.u32 %5, %5, 0;\n\t"
        "addc.cc.u32 %6, %6, 0;\n\t"
        "addc.cc.u32 %7, %7, 0;\n\t"
        "addc.u32    %8, 0, 0;"
        : "+r"(t[0]), "+r"(t[1]), "+r"(t[2]), "+r"(t[3]),
          "+r"(t[4]), "+r"(t[5]), "+r"(t[6]), "+r"(t[7]), "=&r"(c)
        : "r"(a0), "r"(a1), "r"(a2));
    return c;
#else
    uint64_t c = (uint64_t)t[0] + a0;
    t[0] = (uint32_t)c; c >>= 32;
    c += (uint64_t)t[1] + a1;
    t[1] = (uint32_t)c; c >>= 32;
    c += (uint64_t)t[2] + a2;
    t[2] = (uint32_t)c; c >>= 32;
#pragma unroll
    for (int j = 3; j < 8; j++) {
        c += t[j];
        t[j] = (uint32_t)c; c >>= 32;
    }
    return (uint32_t)c;
#endif
}

/* r[0..15] = a * b (schoolbook, row by row). */
FP_HD void mp_mul_full(uint32_t r[16], const uint32_t a[8], const uint32_t b[8]) {
    uint32_t t[16];
    /* Row i reads t[i..i+7] and writes t[i..i+8], so the limb each row
     * writes last is the one the next row's carry-out overwrites, and the
     * accumulator is defined limb by limb as the rows advance.  Nothing has
     * to be zeroed first, and the product is exactly 16 limbs because the
     * last row cannot carry past t[15]. */
    mp_mul_row0(t, a, b[0]);
#pragma unroll
    for (int i = 1; i < 8; i++) mp_mac_row9(t + i, a, b[i]);
#pragma unroll
    for (int j = 0; j < 16; j++) r[j] = t[j];
}

/* ------------------------------------------------------------------------
 * Field operations, templated on the modulus descriptor.
 * ---------------------------------------------------------------------- */
template <class M, bool FAST>
struct Field {
    typedef fp256 elt;

    static FP_HD elt zero() { elt r; for (int j = 0; j < 8; j++) r.v[j] = 0; return r; }
    static FP_HD elt one() {
        elt r;
        if (FAST) { r = zero(); r.v[0] = 1; }
        else { for (int j = 0; j < 8; j++) r.v[j] = M::r1(j); }
        return r;
    }
    static FP_HD elt modulus() { elt r; for (int j = 0; j < 8; j++) r.v[j] = M::limb(j); return r; }

    static FP_HD int is_zero(const elt &a) { return mp_is_zero(a.v); }
    static FP_HD int eq(const elt &a, const elt &b) { return mp_eq(a.v, b.v); }
    static FP_HD void cmov(elt &r, const elt &a, uint32_t flag) { mp_cmov(r.v, a.v, flag); }

    static FP_HD elt add(const elt &a, const elt &b) {
        elt s, t;
        uint32_t c = mp_add(s.v, a.v, b.v);
        uint32_t bw = mp_sub(t.v, s.v, modulus().v);
        /* take s - p when the sum overflowed 2^256 or is >= p */
        mp_cmov(s.v, t.v, c | (bw ^ 1u));
        return s;
    }

    static FP_HD elt sub(const elt &a, const elt &b) {
        elt d, t;
        uint32_t bw = mp_sub(d.v, a.v, b.v);
        mp_add(t.v, d.v, modulus().v);
        mp_cmov(d.v, t.v, bw);
        return d;
    }

    static FP_HD elt neg(const elt &a) { return sub(zero(), a); }

    static FP_HD elt dbl(const elt &a) { return add(a, a); }

    /* a > (p-1)/2, i.e. a is the "larger" of {a, -a}.  Used by the
     * negation map.  Note this compares the *internal* representation. */
    static FP_HD int gt_half(const elt &a) {
        elt h, t;
        for (int j = 0; j < 8; j++) h.v[j] = M::half(j);
        return (int)mp_sub(t.v, h.v, a.v);   /* borrow <=> half < a */
    }

    /* --- Montgomery CIOS: r = a * b * 2^-256 mod p ------------------- */
    static FP_HD elt mont_mul(const elt &a, const elt &b) {
        uint32_t t[10];
        uint32_t p[8];
#pragma unroll
        for (int j = 0; j < 10; j++) t[j] = 0;
#pragma unroll
        for (int j = 0; j < 8; j++) p[j] = M::limb(j);
        const uint32_t np = M::nprime();
#pragma unroll
        for (int i = 0; i < 8; i++) {
            /* invariant: t < 2p < 2^257 at loop entry (t[8] <= 1, t[9] = 0);
             * t + a*b_i < 2^289 and t + a*b_i + m*p < 2^290 so the 10-limb
             * window suffices. */
            mp_mac_row(t, a.v, b.v[i]);
            uint32_t m = t[0] * np;
            mp_mac_row(t, p, m);            /* now t[0] == 0 */
#pragma unroll
            for (int j = 0; j < 9; j++) t[j] = t[j + 1];
            t[9] = 0;
        }
        /* t[0..8] < 2p: one conditional subtraction */
        elt r, s;
#pragma unroll
        for (int j = 0; j < 8; j++) r.v[j] = t[j];
        uint32_t bw = mp_sub(s.v, r.v, p);
        mp_cmov(r.v, s.v, (t[8] != 0) | (bw ^ 1u));
        return r;
    }

    /* --- secp256k1 special reduction: 2^256 = 2^32 + 977 (mod p) ------ */
    static FP_HD elt secp_reduce(const uint32_t r[16]) {
        uint32_t t[10];
        /* fold 1: T = lo + hi*977 + hi<<32   (< 2^289)
         *
         * Both halves go through the same carry-chain primitives the
         * multiplication rows use.  This fold is the bulk of the reduction
         * and the reduction is ~40% of a field multiply, so leaving it on
         * the portable path wasted much of what FP_PTX buys elsewhere.
         *
         * Bounds: lo + hi*977 < 2^266, so the mac writes t[8] < 2^10 and
         * leaves t[9] to be cleared here; adding hi << 32 takes
         * the total under 2^289, so t[9] ends at 0 or 1. */
#pragma unroll
        for (int j = 0; j < 8; j++) t[j] = r[j];
        mp_mac_row9(t, r + 8, 977u);      /* T = lo + hi * 977, writes t[8] */
        t[9] = 0;
        mp_add_shift32(t, r + 8);         /* T += hi << 32 */

        /* fold 2: hi2 = T >> 256 < 2^33, and hi2 * (2^32 + 977) < 2^65 + 2^43
         * is three limbs, so the whole fold is one carry chain over the low
         * half rather than a multiply interleaved with the carry loop. */
        uint64_t hi2 = (uint64_t)t[8] | ((uint64_t)t[9] << 32);
        uint64_t m = hi2 * 977u;                   /* < 2^43 */
        uint64_t mid = (m >> 32) + (uint32_t)hi2;  /* < 2^32 + 2^11 */
        uint32_t a0 = (uint32_t)m;
        uint32_t a1 = (uint32_t)mid;
        uint32_t a2 = (uint32_t)(mid >> 32) + (uint32_t)(hi2 >> 32);  /* <= 2 */
        uint32_t carry = mp_add_small(t, a0, a1, a2);

        /* fold 3: the carry is 0 or 1, and when it is 1 the low half is
         * below 2^66, so adding 2^32 + 977 cannot carry out again. */
        mp_add_small(t, carry * 977u, carry, 0u);

        /* value < 2^256 < 2p: final conditional subtraction */
        elt out, s;
#pragma unroll
        for (int j = 0; j < 8; j++) out.v[j] = t[j];
        uint32_t bw = mp_sub(s.v, out.v, modulus().v);
        mp_cmov(out.v, s.v, bw ^ 1u);
        return out;
    }

    static FP_HD elt mul(const elt &a, const elt &b) {
        if (FAST) {
            uint32_t r[16];
            mp_mul_full(r, a.v, b.v);
            return secp_reduce(r);
        }
        return mont_mul(a, b);
    }

    static FP_HD elt sqr(const elt &a) { return mul(a, a); }

    /* Conversions between canonical integers and the internal form. */
    static FP_HD elt from_canonical(const elt &x) {
        if (FAST) return x;
        elt r2;
        for (int j = 0; j < 8; j++) r2.v[j] = M::r2(j);
        return mont_mul(x, r2);
    }
    static FP_HD elt to_canonical(const elt &x) {
        if (FAST) return x;
        elt one_plain = zero();
        one_plain.v[0] = 1;
        return mont_mul(x, one_plain);
    }
    static FP_HD elt from_u32(uint32_t k) {
        elt r = zero();
        r.v[0] = k;
        return from_canonical(r);
    }
    static FP_HD elt from_limbs(const uint32_t l[8]) {
        elt r;
        for (int j = 0; j < 8; j++) r.v[j] = l[j];
        return from_canonical(r);
    }

    /* Squaring chain helper: a^(2^n). */
    static FP_HD elt sqr_n(elt a, int n) {
#pragma unroll 1
        for (int i = 0; i < n; i++) a = sqr(a);
        return a;
    }

    /* Addition chain for a^(p-2) on secp256k1, where p - 2 = 2^256 - 2^32 - 979.
     * Builds a^(2^k - 1) for k = 2, 3, 6, 9, 11, 22, 44, 88, 176, 220, 223 and
     * assembles the tail with a sliding window (the chain used by
     * libsecp256k1).  255 squarings + 15 multiplications = 270 operations
     * against 316 for a 4-bit window, and -- the reason it matters on a GPU --
     * it needs 11 live field elements instead of a 16-entry table, cutting
     * 512 bytes off every thread's stack frame. */
    static FP_BIG elt inv_addchain_secp256k1(const elt &a) {
        elt x2 = mul(sqr(a), a);                       /* 2^2 - 1 */
        elt x3 = mul(sqr(x2), a);                      /* 2^3 - 1 */
        elt x6 = mul(sqr_n(x3, 3), x3);                /* 2^6 - 1 */
        elt x9 = mul(sqr_n(x6, 3), x3);                /* 2^9 - 1 */
        elt x11 = mul(sqr_n(x9, 2), x2);               /* 2^11 - 1 */
        elt x22 = mul(sqr_n(x11, 11), x11);            /* 2^22 - 1 */
        elt x44 = mul(sqr_n(x22, 22), x22);            /* 2^44 - 1 */
        elt x88 = mul(sqr_n(x44, 44), x44);            /* 2^88 - 1 */
        elt x176 = mul(sqr_n(x88, 88), x88);           /* 2^176 - 1 */
        elt x220 = mul(sqr_n(x176, 44), x44);          /* 2^220 - 1 */
        elt x223 = mul(sqr_n(x220, 3), x3);            /* 2^223 - 1 */
        elt t = sqr_n(x223, 23);
        t = mul(t, x22);
        t = sqr_n(t, 5);
        t = mul(t, a);
        t = sqr_n(t, 3);
        t = mul(t, x2);
        t = sqr_n(t, 2);
        return mul(t, a);
    }

    /* Fermat inversion a^(p-2), 4-bit fixed window.  Branch structure is
     * uniform across the warp (depends only on the constant exponent), so
     * no divergence; the window table lives in local memory (L1).
     * inv(0) = 0. */
    static FP_BIG elt inv_window(const elt &a) {
        elt tbl[16];
        tbl[0] = one();
        tbl[1] = a;
#pragma unroll
        for (int i = 2; i < 16; i++) tbl[i] = mul(tbl[i - 1], a);
        elt acc = one();
        const int top = (M::bits() + 3) / 4 - 1;
#pragma unroll 1
        for (int w = top; w >= 0; w--) {
            if (w != top) {
                acc = sqr(acc); acc = sqr(acc); acc = sqr(acc); acc = sqr(acc);
            }
            uint32_t nib = (M::pm2(w >> 3) >> (4 * (w & 7))) & 15u;
            if (nib) acc = mul(acc, tbl[nib]);
        }
        return acc;
    }

    /* Inversion entry point: the short addition chain where we have one,
     * the generic window otherwise. */
    static FP_HD elt inv(const elt &a) {
        if (FAST) return inv_addchain_secp256k1(a);
        return inv_window(a);
    }

    /* Montgomery's simultaneous-inversion trick for n elements:
     * one inversion + 3(n-1) multiplications.  `scratch` needs n slots.
     * Zero inputs are not allowed (caller filters them). */
    static FP_BIG void batch_inv(elt *x, int n, elt *scratch) {
        scratch[0] = x[0];
        for (int i = 1; i < n; i++) scratch[i] = mul(scratch[i - 1], x[i]);
        elt t = inv(scratch[n - 1]);
        for (int i = n - 1; i > 0; i--) {
            elt xi = mul(t, scratch[i - 1]);
            t = mul(t, x[i]);
            x[i] = xi;
        }
        x[0] = t;
    }
};

typedef Field<ModP, (FP_FAST != 0)> Fp;   /* base field */
typedef Field<ModN, false> Fn;            /* scalar ring mod group order */

#endif /* GPU_ECC_FP256_CUH */
