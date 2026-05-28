/* Phase 22.6: Pollard rho random walk with Barrett-reduction modular mul.
 *
 * Builds on phase21_rho_v3.c (affine + pthreads + GMP) by replacing
 * mpz_mul + mpz_mod with Barrett reduction using __uint128_t.
 *
 * Modular inverse is kept via GMP (used once per affine point add),
 * because Fermat's via Barrett would be ~120 Barrett muls per inverse,
 * which is similar cost to GMP's optimized binary-GCD inverse.
 *
 * Build: clang -O3 -o phase22_rho_barrett phase22_rho_barrett.c -lgmp -lpthread
 *
 * Hypothesis: end-to-end >=5x over phase21_rho_v3.c.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <pthread.h>
#include <gmp.h>

typedef unsigned __int128 u128;

#define NUM_PARTITIONS 16

typedef struct {
    u128 p;
    u128 a;
    u128 b;
    /* Barrett structure */
    u128 mu;
    int k;
} curve_t;

typedef struct {
    u128 x, y;
    int is_infinity;
} point_t;

/* 256-bit product */
typedef struct { u128 hi, lo; } u256;
static inline u256 mul_128(u128 a, u128 b) {
    uint64_t ah = (uint64_t)(a >> 64), al = (uint64_t)a;
    uint64_t bh = (uint64_t)(b >> 64), bl = (uint64_t)b;
    u128 ll = (u128)al * bl;
    u128 lh = (u128)al * bh;
    u128 hl = (u128)ah * bl;
    u128 hh = (u128)ah * bh;
    u128 mid = lh + hl;
    int mid_carry = (mid < lh) ? 1 : 0;
    u128 result_lo = ll + (mid << 64);
    int lo_carry = (result_lo < ll) ? 1 : 0;
    u128 result_hi = hh + (mid >> 64) + ((u128)mid_carry << 64) + lo_carry;
    return (u256){result_hi, result_lo};
}

/* Barrett: a * mb mod p (where p, mu, k are precomputed in curve_t) */
static inline u128 mod_mul(u128 a, u128 mb, const curve_t *c) {
    u256 t = mul_128(a, mb);
    u256 hi_part = mul_128(t.hi, c->mu);
    int shift = 2 * c->k - 128;  /* for k=81, shift=34 */
    u128 q_approx = (hi_part.hi << (128 - shift)) | (hi_part.lo >> shift);
    u256 lo_part = mul_128(t.lo, c->mu);
    u128 q_add = (lo_part.hi >> shift);
    u128 q_total = q_approx + q_add;
    u256 qp = mul_128(q_total, c->p);
    u128 r = t.lo - qp.lo;
    while (r >= c->p) r -= c->p;
    return r;
}

static inline u128 mod_add(u128 a, u128 b, u128 p) {
    u128 r = a + b;
    if (r < a || r >= p) r -= p;
    return r;
}

static inline u128 mod_sub(u128 a, u128 b, u128 p) {
    if (a >= b) return a - b;
    return p - (b - a);
}

/* Thread-local persistent mpz_t for modular inversion (kept as comparison) */
typedef struct {
    mpz_t ma, mp, mr, tmp;
    int initialized;
} inv_ctx_t;

static __thread inv_ctx_t inv_ctx = {0};

static void init_inv_ctx(u128 p) {
    mpz_init(inv_ctx.ma); mpz_init(inv_ctx.mp); mpz_init(inv_ctx.mr); mpz_init(inv_ctx.tmp);
    uint64_t plo = (uint64_t)p, phi = (uint64_t)(p >> 64);
    if (phi) { mpz_import(inv_ctx.mp, 1, 1, 8, 0, 0, &phi); mpz_mul_2exp(inv_ctx.mp, inv_ctx.mp, 64); }
    mpz_import(inv_ctx.tmp, 1, 1, 8, 0, 0, &plo);
    mpz_add(inv_ctx.mp, inv_ctx.mp, inv_ctx.tmp);
    inv_ctx.initialized = 1;
}

/* Modular inverse via GMP — persistent mpz_t */
static u128 mod_inv_gmp(u128 a, u128 p) {
    if (!inv_ctx.initialized) init_inv_ctx(p);
    uint64_t alo = (uint64_t)a, ahi = (uint64_t)(a >> 64);
    mpz_set_ui(inv_ctx.ma, 0);
    if (ahi) {
        mpz_set_ui(inv_ctx.ma, (unsigned long)ahi);
        mpz_mul_2exp(inv_ctx.ma, inv_ctx.ma, 64);
    }
    mpz_add_ui(inv_ctx.ma, inv_ctx.ma, (unsigned long)alo);
    mpz_invert(inv_ctx.mr, inv_ctx.ma, inv_ctx.mp);
    u128 result = 0;
    const mp_limb_t *limbs = mpz_limbs_read(inv_ctx.mr);
    if (mpz_size(inv_ctx.mr) >= 1) result |= (u128)limbs[0];
    if (mpz_size(inv_ctx.mr) >= 2) result |= ((u128)limbs[1]) << 64;
    return result;
}

/* Modular inverse via Fermat's: a^(p-2) mod p, using Barrett mul */
static u128 mod_inv_fermat(u128 a, const curve_t *c) {
    /* exp = p - 2 */
    u128 exp = c->p - 2;
    u128 result = 1;
    u128 base = a;
    while (exp > 0) {
        if (exp & 1) result = mod_mul(result, base, c);
        base = mod_mul(base, base, c);
        exp >>= 1;
    }
    return result;
}

/* Choose inverse implementation: 0 = GMP, 1 = Fermat */
#ifndef USE_FERMAT_INV
#define USE_FERMAT_INV 0
#endif

static inline u128 mod_inv(u128 a, const curve_t *c) {
#if USE_FERMAT_INV
    return mod_inv_fermat(a, c);
#else
    return mod_inv_gmp(a, c->p);
#endif
}

static void point_double(point_t *r, const point_t *P, const curve_t *c) {
    if (P->is_infinity || P->y == 0) { r->is_infinity = 1; return; }
    /* lambda_num = 3*x^2 + a */
    u128 x_sq = mod_mul(P->x, P->x, c);
    u128 t = mod_add(x_sq, x_sq, c->p);
    t = mod_add(t, x_sq, c->p);
    u128 lambda_num = mod_add(t, c->a, c->p);
    /* lambda_den = 2y */
    u128 lambda_den = mod_add(P->y, P->y, c->p);
    u128 lambda_den_inv = mod_inv(lambda_den, c);
    u128 lambda = mod_mul(lambda_num, lambda_den_inv, c);
    /* x3 = lambda^2 - 2x */
    u128 lambda_sq = mod_mul(lambda, lambda, c);
    u128 x3 = mod_sub(lambda_sq, P->x, c->p);
    x3 = mod_sub(x3, P->x, c->p);
    /* y3 = lambda*(x - x3) - y */
    u128 dx = mod_sub(P->x, x3, c->p);
    u128 y3 = mod_mul(lambda, dx, c);
    y3 = mod_sub(y3, P->y, c->p);
    r->x = x3; r->y = y3; r->is_infinity = 0;
}

static void point_add(point_t *r, const point_t *P, const point_t *Q, const curve_t *c) {
    if (P->is_infinity) { *r = *Q; return; }
    if (Q->is_infinity) { *r = *P; return; }
    if (P->x == Q->x) {
        if (mod_add(P->y, Q->y, c->p) == 0) { r->is_infinity = 1; return; }
        point_double(r, P, c); return;
    }
    u128 lambda_num = mod_sub(Q->y, P->y, c->p);
    u128 lambda_den = mod_sub(Q->x, P->x, c->p);
    u128 lambda_den_inv = mod_inv(lambda_den, c);
    u128 lambda = mod_mul(lambda_num, lambda_den_inv, c);
    u128 lambda_sq = mod_mul(lambda, lambda, c);
    u128 x3 = mod_sub(lambda_sq, P->x, c->p);
    x3 = mod_sub(x3, Q->x, c->p);
    u128 dx = mod_sub(P->x, x3, c->p);
    u128 y3 = mod_mul(lambda, dx, c);
    y3 = mod_sub(y3, P->y, c->p);
    r->x = x3; r->y = y3; r->is_infinity = 0;
}

static int point_partition(const point_t *P) {
    if (P->is_infinity) return 0;
    return (int)((uint64_t)P->x & (NUM_PARTITIONS - 1));
}

static void apply_neg(point_t *P, const curve_t *c) {
    if (P->is_infinity) return;
    u128 half = c->p >> 1;
    if (P->y > half) P->y = c->p - P->y;
}

typedef struct {
    int tid, num_steps;
    const curve_t *c;
    const point_t *M;
} thread_t;

static void *worker(void *arg) {
    thread_t *T = (thread_t *)arg;
    point_t cur = T->M[T->tid % NUM_PARTITIONS];
    point_t next;
    for (int i = 0; i < T->num_steps; i++) {
        int part = point_partition(&cur);
        point_add(&next, &cur, &T->M[part], T->c);
        apply_neg(&next, T->c);
        cur = next;
    }
    return NULL;
}

static u128 parse_u128(const char *s) {
    u128 r = 0;
    while (*s >= '0' && *s <= '9') { r = r * 10 + (*s - '0'); s++; }
    return r;
}

/* Initialize Barrett structure */
static void init_barrett(curve_t *c) {
    int k = 0; u128 t = c->p;
    while (t > 0) { k++; t >>= 1; }
    c->k = k;
    /* mu = floor(2^{2k} / p) using GMP */
    mpz_t mu, p_mpz, pow;
    mpz_init(mu); mpz_init(p_mpz); mpz_init(pow);
    uint64_t plo = (uint64_t)c->p, phi = (uint64_t)(c->p >> 64);
    mpz_t tmp; mpz_init(tmp);
    if (phi) { mpz_import(p_mpz, 1, 1, 8, 0, 0, &phi); mpz_mul_2exp(p_mpz, p_mpz, 64); }
    mpz_import(tmp, 1, 1, 8, 0, 0, &plo); mpz_add(p_mpz, p_mpz, tmp);
    mpz_set_ui(pow, 1); mpz_mul_2exp(pow, pow, 2 * k);
    mpz_fdiv_q(mu, pow, p_mpz);
    const mp_limb_t *limbs = mpz_limbs_read(mu);
    u128 mu_val = 0;
    if (mpz_size(mu) >= 1) mu_val |= (u128)limbs[0];
    if (mpz_size(mu) >= 2) mu_val |= ((u128)limbs[1]) << 64;
    c->mu = mu_val;
    mpz_clear(mu); mpz_clear(p_mpz); mpz_clear(pow); mpz_clear(tmp);
}

int main(int argc, char *argv[]) {
    if (argc < 7) {
        fprintf(stderr, "Usage: %s p a b xP yP num_steps [num_threads]\n", argv[0]);
        return 1;
    }
    curve_t c;
    c.p = parse_u128(argv[1]);
    c.a = parse_u128(argv[2]);
    c.b = parse_u128(argv[3]);
    init_barrett(&c);

    point_t M[NUM_PARTITIONS];
    M[0].x = parse_u128(argv[4]);
    M[0].y = parse_u128(argv[5]);
    M[0].is_infinity = 0;
    for (int i = 1; i < NUM_PARTITIONS; i++) {
        point_double(&M[i], &M[i-1], &c);
    }
    int num_steps = atoi(argv[6]);
    int n_threads = (argc >= 8) ? atoi(argv[7]) : 1;

    thread_t T[16];
    pthread_t tids[16];
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (int i = 0; i < n_threads; i++) {
        T[i].tid = i; T[i].num_steps = num_steps; T[i].c = &c; T[i].M = M;
        pthread_create(&tids[i], NULL, worker, &T[i]);
    }
    for (int i = 0; i < n_threads; i++) pthread_join(tids[i], NULL);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double wall = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    long total = (long)n_threads * num_steps;
    fprintf(stderr, "Barrett-rho: %d threads x %d = %ld ops in %.3f s = %.2e ops/sec\n",
            n_threads, num_steps, total, wall, total / wall);
    printf("%.2e\n", total / wall);
    return 0;
}
