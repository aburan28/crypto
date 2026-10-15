/* fp64.cuh -- Montgomery arithmetic mod an odd prime p < 2^62 on one
 * 64-bit word, for CUDA devices and for host emulation of the same code.
 *
 * This is a line-for-line mirror of `Fp64` in
 * `src/cryptanalysis/prime_fast.rs` (R = 2^64, values held as a*R mod p,
 * reduction `t -> (t + m*p) >> 64` with m = lo(t) * (-p^-1 mod 2^64)),
 * so a Montgomery-form value produced here is bit-identical to the one
 * the Rust engine produces for the same canonical input.  That is what
 * lets the Rust harness compare the GPU oracle against `prime_fast`
 * exactly rather than modulo re-encoding.
 *
 * Every function is `__host__ __device__`: the CPU test harness and the
 * emulation backend compile *this file* with g++/clang++, the kernels
 * compile it with nvcc.  The only device/host split is the 64x64->128
 * product: `__umul64hi` on the device, `unsigned __int128` on the host.
 *
 * Blackwell (sm_120) notes
 * ------------------------
 * - The integer datapath is 32-bit.  A 64x64->64 low product is three
 *   IMAD.WIDE.U32, `__umul64hi` is four plus carry handling, so one
 *   Montgomery multiplication here is roughly 12-16 IMAD-class
 *   instructions plus a handful of IADD3/ISETP.  Reductions are branch
 *   free (the final conditional subtraction is a select), so a warp never
 *   diverges inside the field layer.
 * - A thread holds a few dozen 64-bit values at most (see rootfind.cuh),
 *   which keeps the kernel well under the 255-register ceiling and lets
 *   ptxas place >= 8 warps per scheduler on sm_120.  Register pressure
 *   comes from the batched-inversion scratch, not from the field.
 * - p < 2^62 means a + b < 2^63 and t + m*p < 2^128 never overflow, the
 *   same bound `prime_fast::MAX_BITS` enforces on the CPU.
 */
#ifndef GPU_ROOTFIND_FP64_CUH
#define GPU_ROOTFIND_FP64_CUH

#include <stdint.h>

#ifdef __CUDACC__
#define RF_HD __host__ __device__ __forceinline__
#define RF_BIG __host__ __device__ __noinline__
#else
#define RF_HD inline
#define RF_BIG inline
#endif

/* Field context.  Plain-old-data so it can be passed to a kernel by value
 * (it lands in the constant parameter bank, which is broadcast to every
 * thread for free). */
struct Fp64Ctx {
    uint64_t p;      /* odd prime, 3 < p < 2^62 */
    uint64_t n_inv;  /* -p^-1 mod 2^64 */
    uint64_t one;    /* R mod p  (Montgomery 1) */
    uint64_t r2;     /* R^2 mod p (canonical -> Montgomery) */
    uint64_t pm2;    /* p - 2, the Fermat exponent */
};

RF_HD uint64_t rf_mulhi64(uint64_t a, uint64_t b) {
#ifdef __CUDA_ARCH__
    return __umul64hi(a, b);
#else
    return (uint64_t)(((unsigned __int128)a * b) >> 64);
#endif
}

/* Montgomery reduction of the 128-bit value (hi:lo) < p * 2^64:
 * returns (hi:lo) * R^-1 mod p.  Written without a 128-bit type so the
 * device path is exactly the arithmetic the hardware offers:
 *   m  = lo * n_inv                 (mod 2^64)
 *   u  = hi + mulhi(m, p) + carry(lo + lo(m*p))
 * and lo + lo(m*p) == 0 mod 2^64 by construction, so the carry is simply
 * (lo != 0).  Result is < 2p, one conditional subtraction finishes. */
RF_HD uint64_t fp64_reduce(const Fp64Ctx &f, uint64_t hi, uint64_t lo) {
    uint64_t m = lo * f.n_inv;
    uint64_t u = hi + rf_mulhi64(m, f.p) + (lo != 0 ? 1u : 0u);
    return u >= f.p ? u - f.p : u;
}

RF_HD uint64_t fp64_mul(const Fp64Ctx &f, uint64_t a, uint64_t b) {
    return fp64_reduce(f, rf_mulhi64(a, b), a * b);
}

RF_HD uint64_t fp64_sqr(const Fp64Ctx &f, uint64_t a) { return fp64_mul(f, a, a); }

RF_HD uint64_t fp64_add(const Fp64Ctx &f, uint64_t a, uint64_t b) {
    uint64_t s = a + b;
    return s >= f.p ? s - f.p : s;
}

RF_HD uint64_t fp64_sub(const Fp64Ctx &f, uint64_t a, uint64_t b) {
    return a >= b ? a - b : a + f.p - b;
}

RF_HD uint64_t fp64_neg(const Fp64Ctx &f, uint64_t a) { return a == 0 ? 0 : f.p - a; }

/* canonical -> Montgomery (input must already be < p) */
RF_HD uint64_t fp64_to_mont(const Fp64Ctx &f, uint64_t a) { return fp64_mul(f, a, f.r2); }

/* Montgomery -> canonical */
RF_HD uint64_t fp64_from_mont(const Fp64Ctx &f, uint64_t a) { return fp64_reduce(f, 0, a); }

/* Fermat inversion a^(p-2).  The exponent is the public modulus, so the
 * branch pattern is identical for every thread of the warp: no
 * divergence, ~62 squarings + popcount(p-2) multiplications.  This is
 * deliberately *not* the binary extended Euclid the Rust side uses --
 * EEA has data-dependent loop trip counts and would serialise a warp.
 * Both compute the same field element, and the Montgomery representation
 * of a residue is unique, so the result matches `Fp64::inv` bit for bit.
 * inv(0) = 0. */
RF_BIG uint64_t fp64_inv(const Fp64Ctx &f, uint64_t a) {
    uint64_t acc = f.one;
    uint64_t base = a;
    uint64_t e = f.pm2;
#pragma unroll 1
    while (e) {
        if (e & 1) acc = fp64_mul(f, acc, base);
        base = fp64_sqr(f, base);
        e >>= 1;
    }
    return acc;
}

/* Montgomery's simultaneous inversion of x[0..n) in place using scratch
 * of n words: one fp64_inv plus 3(n-1) multiplications.  Inputs must be
 * non-zero (the caller substitutes `one` for degenerate entries and
 * flags them). */
RF_HD void fp64_batch_inv(const Fp64Ctx &f, uint64_t *x, int n, uint64_t *scratch) {
    uint64_t acc = f.one;
    for (int i = 0; i < n; i++) {
        scratch[i] = acc;
        acc = fp64_mul(f, acc, x[i]);
    }
    uint64_t inv_acc = fp64_inv(f, acc);
    for (int i = n - 1; i >= 0; i--) {
        uint64_t xi = x[i];
        x[i] = fp64_mul(f, inv_acc, scratch[i]);
        inv_acc = fp64_mul(f, inv_acc, xi);
    }
}

/* Host-side context construction (mirrors `Fp64::new`). */
static inline Fp64Ctx fp64_ctx_new(uint64_t p) {
    Fp64Ctx f;
    f.p = p;
    uint64_t inv = 1;
    for (int i = 0; i < 6; i++) inv *= 2u - p * inv;   /* Newton: p^-1 mod 2^64 */
    f.n_inv = (uint64_t)0 - inv;
    unsigned __int128 R = ((unsigned __int128)1) << 64;
    f.one = (uint64_t)(R % p);
    f.r2 = (uint64_t)(((unsigned __int128)f.one * f.one) % p);
    f.pm2 = p - 2;
    return f;
}

#endif /* GPU_ROOTFIND_FP64_CUH */
