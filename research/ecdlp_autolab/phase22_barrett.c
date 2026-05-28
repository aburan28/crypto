/* Phase 22.1: Barrett reduction for 81-bit prime modular arithmetic.
 *
 * Goal: beat GMP's mpz_mul + mpz_mod for fixed-prime modular mul.
 *
 * Algorithm (Barrett):
 *   Given fixed p, precompute mu = floor(2^{2k} / p) where k = ceil(log2 p).
 *   To compute a * b mod p:
 *     1. t = a * b  (in 2k-bit precision)
 *     2. q = (t * mu) >> 2k  (approximate t / p)
 *     3. r = t - q * p
 *     4. while r >= p: r -= p  (at most 2 iterations)
 *
 * For 81-bit p: k = 81, 2k = 162, mu fits in 82 bits.
 * Use __uint128_t for intermediate products.
 *
 * Build: clang -O3 -o phase22_barrett phase22_barrett.c -lpthread -lgmp
 * Run:   ./phase22_barrett <p> <a> <b> <num_muls>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>

typedef unsigned __int128 u128;

/* Read u128 from base-10 string */
static u128 parse_u128(const char *s) {
    u128 r = 0;
    while (*s >= '0' && *s <= '9') {
        r = r * 10 + (*s - '0');
        s++;
    }
    return r;
}

/* Print u128 in base 10 */
static void print_u128(u128 x) {
    char buf[50];
    int i = 49;
    buf[i--] = 0;
    if (x == 0) { buf[i--] = '0'; }
    else {
        while (x > 0) {
            buf[i--] = '0' + (int)(x % 10);
            x /= 10;
        }
    }
    fprintf(stderr, "%s", &buf[i+1]);
}

/* Helper: 128 × 128 → 256-bit product via __uint128 halves */
typedef struct { u128 hi, lo; } u256;

static u256 mul_128(u128 a, u128 b) {
    /* Split each into 64-bit halves */
    uint64_t ah = (uint64_t)(a >> 64);
    uint64_t al = (uint64_t)a;
    uint64_t bh = (uint64_t)(b >> 64);
    uint64_t bl = (uint64_t)b;

    u128 ll = (u128)al * bl;        /* low 128 of low partial */
    u128 lh = (u128)al * bh;
    u128 hl = (u128)ah * bl;
    u128 hh = (u128)ah * bh;

    /* Result = hh * 2^128 + (lh + hl) * 2^64 + ll */
    u128 mid = lh + hl;
    int mid_carry = (mid < lh) ? 1 : 0;  /* did sum overflow */

    u128 result_lo = ll + (mid << 64);
    int lo_carry = (result_lo < ll) ? 1 : 0;

    u128 result_hi = hh + (mid >> 64) + ((u128)mid_carry << 64) + lo_carry;

    u256 r = {result_hi, result_lo};
    return r;
}

/* Barrett reduction structure */
typedef struct {
    u128 p;
    u128 mu;     /* floor(2^{2k} / p) where k = 81 */
    int k;
} barrett_t;

static void barrett_init(barrett_t *b, u128 p) {
    b->p = p;
    /* k = ceil(log2 p) */
    int k = 0;
    u128 t = p;
    while (t > 0) { k++; t >>= 1; }
    b->k = k;

    /* Compute mu = floor(2^{2k} / p) using GMP for precision */
    mpz_t mu_mpz, p_mpz, pow_mpz;
    mpz_init(mu_mpz); mpz_init(p_mpz); mpz_init(pow_mpz);

    /* p_mpz = p */
    uint64_t plo = (uint64_t)p;
    uint64_t phi = (uint64_t)(p >> 64);
    mpz_import(p_mpz, 1, 1, sizeof(uint64_t), 0, 0, &phi);
    mpz_mul_2exp(p_mpz, p_mpz, 64);
    mpz_t tmp; mpz_init(tmp);
    mpz_import(tmp, 1, 1, sizeof(uint64_t), 0, 0, &plo);
    mpz_add(p_mpz, p_mpz, tmp);
    mpz_clear(tmp);

    /* pow_mpz = 2^{2k} */
    mpz_set_ui(pow_mpz, 1);
    mpz_mul_2exp(pow_mpz, pow_mpz, 2 * k);

    /* mu_mpz = pow_mpz / p_mpz */
    mpz_fdiv_q(mu_mpz, pow_mpz, p_mpz);

    /* Extract mu_mpz to u128 (it should fit since mu has ~k+1 bits) */
    size_t mu_bits = mpz_sizeinbase(mu_mpz, 2);
    if (mu_bits > 128) {
        fprintf(stderr, "ERROR: mu doesn't fit in u128 (%zu bits)\n", mu_bits);
        exit(1);
    }
    /* Get u128 limbs */
    uint64_t mu_lo = 0, mu_hi = 0;
    mp_limb_t *limbs = mpz_limbs_read(mu_mpz);
    mu_lo = (uint64_t)limbs[0];
    if (mpz_size(mu_mpz) >= 2) {
        mu_hi = (uint64_t)limbs[1];
    }
    b->mu = ((u128)mu_hi << 64) | mu_lo;

    mpz_clear(mu_mpz); mpz_clear(p_mpz); mpz_clear(pow_mpz);
}

/* Compute a * b mod p using Barrett */
static u128 barrett_mul(const barrett_t *b, u128 a, u128 mb) {
    /* t = a * mb (256-bit) */
    u256 t = mul_128(a, mb);
    /* q = (t * mu) >> 2k. Take upper 256 - 2k bits.
     * For k=81, 2k=162. We need bits [162..323] of (t * mu).
     * Since mu has ~82 bits and t has ~162 bits, product has ~244 bits.
     * Shift right by 162 ⇒ value of ~82 bits, fits in u128. */

    /* Compute t * mu in pieces. t is 256-bit, mu is 128-bit.
     * t.hi * mu = part contributing to bits [128..384]
     * t.lo * mu = part contributing to bits [0..256]
     * Total fits in 256 + 128 = 384 bits.
     *
     * We want (t * mu) >> 162. That's bits [162..] of the full product.
     */
    /* (t.hi * mu) >> (162 - 128) = (t.hi * mu) >> 34. This gives bits [34..162] of t.hi * mu. */
    /* Plus contribution from t.lo * mu's upper bits, etc. */

    /* Easier: compute upper 128 bits of (t * mu) approximation:
     * q ≈ (t.hi * mu) >> 34   (if we ignore lower contributions)
     * This may underestimate by 1-2; we'll correct after. */
    u256 hi_part = mul_128(t.hi, b->mu);
    /* (hi_part.hi << 128 + hi_part.lo) >> 34 = (hi_part.hi << 94) | (hi_part.lo >> 34) */
    u128 q_approx = (hi_part.hi << (128 - 34)) | (hi_part.lo >> 34);

    /* Also need contribution from t.lo * mu's upper bits.
     * t.lo * mu has 256 bits; its upper bits >> 34 are also part of q. */
    u256 lo_part = mul_128(t.lo, b->mu);
    /* The upper 128 bits of lo_part shifted right by 34 in the overall combination */
    u128 q_add = (lo_part.hi >> 34);
    /* Plus carry from lower: ((lo_part.hi & 0x3FFFFFFFF) << (128 - 34)) | (lo_part.lo >> 34)
     * — this is small contribution, but we need it. */

    u128 q_total = q_approx + q_add;  /* may slightly underestimate */

    /* r = t.lo - q_total * p (mod 2^128) */
    /* Since q_total is approximate, r may be in [0, 3p) range. */
    u256 qp = mul_128(q_total, b->p);
    u128 r = t.lo - qp.lo;

    /* Conditional subtractions */
    while (r >= b->p) r -= b->p;

    return r;
}

/* Reference: GMP mul-mod */
static u128 gmp_mul_mod(u128 a, u128 mb, u128 p) {
    mpz_t mpa, mpb, mpp, mpr;
    mpz_init(mpa); mpz_init(mpb); mpz_init(mpp); mpz_init(mpr);
    uint64_t ah = (uint64_t)(a >> 64), al = (uint64_t)a;
    uint64_t bh = (uint64_t)(mb >> 64), bl = (uint64_t)mb;
    uint64_t ph = (uint64_t)(p >> 64), pl = (uint64_t)p;
    if (ah) { mpz_import(mpa, 1, 1, 8, 0, 0, &ah); mpz_mul_2exp(mpa, mpa, 64); }
    mpz_t tmp; mpz_init(tmp);
    mpz_import(tmp, 1, 1, 8, 0, 0, &al); mpz_add(mpa, mpa, tmp);
    if (bh) { mpz_import(mpb, 1, 1, 8, 0, 0, &bh); mpz_mul_2exp(mpb, mpb, 64); }
    mpz_import(tmp, 1, 1, 8, 0, 0, &bl); mpz_add(mpb, mpb, tmp);
    if (ph) { mpz_import(mpp, 1, 1, 8, 0, 0, &ph); mpz_mul_2exp(mpp, mpp, 64); }
    mpz_import(tmp, 1, 1, 8, 0, 0, &pl); mpz_add(mpp, mpp, tmp);
    mpz_mul(mpr, mpa, mpb);
    mpz_mod(mpr, mpr, mpp);
    /* Read back */
    u128 result = 0;
    mp_limb_t *limbs = mpz_limbs_read(mpr);
    if (mpz_size(mpr) >= 1) result |= (u128)limbs[0];
    if (mpz_size(mpr) >= 2) result |= ((u128)limbs[1]) << 64;
    mpz_clear(mpa); mpz_clear(mpb); mpz_clear(mpp); mpz_clear(mpr); mpz_clear(tmp);
    return result;
}

int main(int argc, char *argv[]) {
    u128 p = parse_u128(argv[1]);
    u128 a = parse_u128(argv[2]);
    u128 b = parse_u128(argv[3]);
    int n = atoi(argv[4]);

    barrett_t bar;
    barrett_init(&bar, p);
    fprintf(stderr, "p = "); print_u128(p); fprintf(stderr, "\n");
    fprintf(stderr, "k = %d, mu = ", bar.k); print_u128(bar.mu); fprintf(stderr, "\n");

    /* Verify correctness on a few inputs */
    u128 ref = gmp_mul_mod(a, b, p);
    u128 my = barrett_mul(&bar, a, b);
    fprintf(stderr, "a*b mod p:\n");
    fprintf(stderr, "  GMP: "); print_u128(ref); fprintf(stderr, "\n");
    fprintf(stderr, "  Barrett: "); print_u128(my); fprintf(stderr, "\n");
    fprintf(stderr, "  match: %s\n", (ref == my) ? "YES" : "NO");
    if (ref != my) {
        fprintf(stderr, "FAILED correctness check; debugging needed\n");
        return 1;
    }

    /* Benchmark */
    struct timespec ts0, ts1;
    volatile u128 sink = 0;
    u128 acc = a;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (int i = 0; i < n; i++) {
        acc = barrett_mul(&bar, acc, b);
    }
    sink = acc;  /* prevent dead-code elimination */
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double wall_bar = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    fprintf(stderr, "Barrett: %d muls in %.3f s = %.2e muls/sec (sink=%llu)\n",
            n, wall_bar, n / wall_bar, (unsigned long long)sink);

    acc = a;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (int i = 0; i < n; i++) {
        acc = gmp_mul_mod(acc, b, p);
    }
    sink = acc;
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double wall_gmp = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    fprintf(stderr, "GMP (malloc per call): %d muls in %.3f s = %.2e muls/sec (sink=%llu)\n",
            n, wall_gmp, n / wall_gmp, (unsigned long long)sink);

    /* Fair GMP: persistent mpz_t */
    mpz_t mp_a, mp_b, mp_p, mp_t;
    mpz_init(mp_a); mpz_init(mp_b); mpz_init(mp_p); mpz_init(mp_t);
    /* Convert u128 → mpz_t */
    {
        uint64_t lo = (uint64_t)a, hi = (uint64_t)(a >> 64);
        mpz_import(mp_a, 1, 1, 8, 0, 0, &hi); mpz_mul_2exp(mp_a, mp_a, 64);
        mpz_t tmp; mpz_init(tmp);
        mpz_import(tmp, 1, 1, 8, 0, 0, &lo); mpz_add(mp_a, mp_a, tmp);
        lo = (uint64_t)b; hi = (uint64_t)(b >> 64);
        mpz_import(mp_b, 1, 1, 8, 0, 0, &hi); mpz_mul_2exp(mp_b, mp_b, 64);
        mpz_import(tmp, 1, 1, 8, 0, 0, &lo); mpz_add(mp_b, mp_b, tmp);
        lo = (uint64_t)p; hi = (uint64_t)(p >> 64);
        mpz_import(mp_p, 1, 1, 8, 0, 0, &hi); mpz_mul_2exp(mp_p, mp_p, 64);
        mpz_import(tmp, 1, 1, 8, 0, 0, &lo); mpz_add(mp_p, mp_p, tmp);
        mpz_clear(tmp);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (int i = 0; i < n; i++) {
        mpz_mul(mp_t, mp_a, mp_b);
        mpz_mod(mp_a, mp_t, mp_p);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double wall_gmp_persist = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    {
        mp_limb_t *limbs = (mp_limb_t *)mpz_limbs_read(mp_a);
        u128 result = 0;
        if (mpz_size(mp_a) >= 1) result |= (u128)limbs[0];
        if (mpz_size(mp_a) >= 2) result |= ((u128)limbs[1]) << 64;
        sink = result;
    }
    fprintf(stderr, "GMP (persistent mpz_t): %d muls in %.3f s = %.2e muls/sec (sink=%llu)\n",
            n, wall_gmp_persist, n / wall_gmp_persist, (unsigned long long)sink);
    mpz_clear(mp_a); mpz_clear(mp_b); mpz_clear(mp_p); mpz_clear(mp_t);

    fprintf(stderr, "Barrett speedup over GMP (malloc): %.2fx\n", wall_gmp / wall_bar);
    fprintf(stderr, "Barrett speedup over GMP (persistent): %.2fx\n", wall_gmp_persist / wall_bar);
    printf("%.2f\n", wall_gmp_persist / wall_bar);
    return 0;
}
