/* Phase 22.9: Pollard rho with Montgomery batch inversion across N parallel walks.
 *
 * Each thread runs N parallel walks in lockstep. Their inverses are batched:
 * for N inversions, do 1 inverse + (3N-1) multiplications instead of N inverses.
 *
 * Per-step amortized inverse cost drops from ~74ns to ~74/N + ~25ns ≈ 50ns for N=4.
 *
 * Build: clang -O3 -o phase22_rho_batch_inv phase22_rho_batch_inv.c -lgmp -lpthread
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
#define WALKS_PER_THREAD 4  /* N parallel walks per thread */

typedef struct {
    u128 p;
    u128 a, b;
    u128 mu;
    int k;
} curve_t;

typedef struct {
    u128 x, y;
    int is_infinity;
} point_t;

/* 256-bit product helper */
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

/* Barrett: a * mb mod p */
static inline u128 mod_mul(u128 a, u128 mb, const curve_t *c) {
    u256 t = mul_128(a, mb);
    u256 hi_part = mul_128(t.hi, c->mu);
    int shift = 2 * c->k - 128;
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

/* GMP inverse using persistent thread-local mpz_t */
typedef struct {
    mpz_t ma, mp, mr;
    int initialized;
} inv_ctx_t;

static __thread inv_ctx_t inv_ctx = {.initialized=0};

static void init_inv_ctx(u128 p) {
    mpz_init(inv_ctx.ma); mpz_init(inv_ctx.mp); mpz_init(inv_ctx.mr);
    uint64_t plo = (uint64_t)p, phi = (uint64_t)(p >> 64);
    mpz_set_ui(inv_ctx.mp, 0);
    if (phi) { mpz_set_ui(inv_ctx.mp, (unsigned long)phi); mpz_mul_2exp(inv_ctx.mp, inv_ctx.mp, 64); }
    mpz_add_ui(inv_ctx.mp, inv_ctx.mp, (unsigned long)plo);
    inv_ctx.initialized = 1;
}

static u128 mod_inv(u128 a, u128 p) {
    if (!inv_ctx.initialized) init_inv_ctx(p);
    uint64_t alo = (uint64_t)a, ahi = (uint64_t)(a >> 64);
    mpz_set_ui(inv_ctx.ma, 0);
    if (ahi) { mpz_set_ui(inv_ctx.ma, (unsigned long)ahi); mpz_mul_2exp(inv_ctx.ma, inv_ctx.ma, 64); }
    mpz_add_ui(inv_ctx.ma, inv_ctx.ma, (unsigned long)alo);
    mpz_invert(inv_ctx.mr, inv_ctx.ma, inv_ctx.mp);
    u128 result = 0;
    const mp_limb_t *limbs = mpz_limbs_read(inv_ctx.mr);
    if (mpz_size(inv_ctx.mr) >= 1) result |= (u128)limbs[0];
    if (mpz_size(inv_ctx.mr) >= 2) result |= ((u128)limbs[1]) << 64;
    return result;
}

/* Batch inversion: invert N values a[0..N-1] in-place using 1 inverse + 3N-1 muls */
static void batch_invert(u128 a[], int N, const curve_t *c) {
    static __thread u128 prod[WALKS_PER_THREAD];
    /* prod[i] = a[0] * a[1] * ... * a[i] */
    prod[0] = a[0];
    for (int i = 1; i < N; i++) {
        prod[i] = mod_mul(prod[i-1], a[i], c);
    }
    /* Invert the total product */
    u128 inv_total = mod_inv(prod[N-1], c->p);
    /* Work backwards */
    for (int i = N - 1; i > 0; i--) {
        u128 inv_i = mod_mul(inv_total, prod[i-1], c);
        inv_total = mod_mul(inv_total, a[i], c);
        a[i] = inv_i;
    }
    a[0] = inv_total;
}

static int point_partition(const point_t *P) {
    if (P->is_infinity) return 0;
    return (int)((uint64_t)P->x & (NUM_PARTITIONS - 1));
}

static u128 parse_u128(const char *s) {
    u128 r = 0;
    while (*s >= '0' && *s <= '9') { r = r * 10 + (*s - '0'); s++; }
    return r;
}

static void init_barrett(curve_t *c) {
    int k = 0; u128 t = c->p;
    while (t > 0) { k++; t >>= 1; }
    c->k = k;
    mpz_t mu, p_mpz, pow, tmp;
    mpz_init(mu); mpz_init(p_mpz); mpz_init(pow); mpz_init(tmp);
    uint64_t plo = (uint64_t)c->p, phi = (uint64_t)(c->p >> 64);
    mpz_set_ui(p_mpz, 0);
    if (phi) { mpz_set_ui(p_mpz, (unsigned long)phi); mpz_mul_2exp(p_mpz, p_mpz, 64); }
    mpz_add_ui(p_mpz, p_mpz, (unsigned long)plo);
    mpz_set_ui(pow, 1); mpz_mul_2exp(pow, pow, 2 * k);
    mpz_fdiv_q(mu, pow, p_mpz);
    const mp_limb_t *limbs = mpz_limbs_read(mu);
    u128 mu_val = 0;
    if (mpz_size(mu) >= 1) mu_val |= (u128)limbs[0];
    if (mpz_size(mu) >= 2) mu_val |= ((u128)limbs[1]) << 64;
    c->mu = mu_val;
    mpz_clear(mu); mpz_clear(p_mpz); mpz_clear(pow); mpz_clear(tmp);
}

/* Apply negation map to canonical form */
static inline void apply_neg(point_t *P, const curve_t *c) {
    if (P->is_infinity) return;
    u128 half = c->p >> 1;
    if (P->y > half) P->y = c->p - P->y;
}

typedef struct {
    int tid, num_beats;  /* beat = WALKS_PER_THREAD steps */
    const curve_t *c;
    const point_t *M;
} thread_t;

static void *worker(void *arg) {
    thread_t *T = (thread_t *)arg;
    const curve_t *c = T->c;
    const point_t *M = T->M;

    /* Initialize N parallel walks with different starting points */
    point_t cur[WALKS_PER_THREAD];
    for (int w = 0; w < WALKS_PER_THREAD; w++) {
        cur[w] = M[(T->tid * WALKS_PER_THREAD + w) % NUM_PARTITIONS];
    }

    u128 dx[WALKS_PER_THREAD];     /* x-differences */
    u128 dy[WALKS_PER_THREAD];     /* y-differences */
    int parts[WALKS_PER_THREAD];

    for (int beat = 0; beat < T->num_beats; beat++) {
        /* For each walk: determine partition and compute dx, dy */
        for (int w = 0; w < WALKS_PER_THREAD; w++) {
            parts[w] = point_partition(&cur[w]);
            const point_t *Mp = &M[parts[w]];
            /* Same-x check (rare); fallback to GMP inverse just for safety */
            if (cur[w].x == Mp->x) {
                /* Degenerate; restart this walk */
                cur[w] = M[(beat + w + 1) % NUM_PARTITIONS];
                dx[w] = mod_sub(M[parts[w]].x, cur[w].x, c->p);
                dy[w] = mod_sub(M[parts[w]].y, cur[w].y, c->p);
            } else {
                dx[w] = mod_sub(Mp->x, cur[w].x, c->p);
                dy[w] = mod_sub(Mp->y, cur[w].y, c->p);
            }
        }
        /* Batch invert dx[0..N-1] */
        batch_invert(dx, WALKS_PER_THREAD, c);
        /* dx[w] is now (Mp.x - cur.x)^(-1); compute slope and new point */
        for (int w = 0; w < WALKS_PER_THREAD; w++) {
            const point_t *Mp = &M[parts[w]];
            u128 lambda = mod_mul(dy[w], dx[w], c);
            u128 lambda_sq = mod_mul(lambda, lambda, c);
            u128 new_x = mod_sub(lambda_sq, cur[w].x, c->p);
            new_x = mod_sub(new_x, Mp->x, c->p);
            u128 ddx = mod_sub(cur[w].x, new_x, c->p);
            u128 new_y = mod_mul(lambda, ddx, c);
            new_y = mod_sub(new_y, cur[w].y, c->p);
            cur[w].x = new_x;
            cur[w].y = new_y;
            cur[w].is_infinity = 0;
            apply_neg(&cur[w], c);
        }
    }
    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc < 7) {
        fprintf(stderr, "Usage: %s p a b xP yP num_beats [num_threads]\n", argv[0]);
        return 1;
    }
    curve_t c;
    c.p = parse_u128(argv[1]);
    c.a = parse_u128(argv[2]);
    c.b = parse_u128(argv[3]);
    init_barrett(&c);

    /* Build M[0..NUM_PARTITIONS-1] via repeated doubling */
    point_t M[NUM_PARTITIONS];
    M[0].x = parse_u128(argv[4]);
    M[0].y = parse_u128(argv[5]);
    M[0].is_infinity = 0;
    for (int i = 1; i < NUM_PARTITIONS; i++) {
        /* M[i] = 2 * M[i-1] (affine point doubling using mod_inv) */
        u128 lambda_num = mod_mul(M[i-1].x, M[i-1].x, &c);
        lambda_num = mod_mul(lambda_num, 3, &c);
        lambda_num = mod_add(lambda_num, c.a, c.p);
        u128 two_y = mod_add(M[i-1].y, M[i-1].y, c.p);
        u128 lambda = mod_mul(lambda_num, mod_inv(two_y, c.p), &c);
        u128 lambda_sq = mod_mul(lambda, lambda, &c);
        u128 nx = mod_sub(lambda_sq, M[i-1].x, c.p);
        nx = mod_sub(nx, M[i-1].x, c.p);
        u128 dx = mod_sub(M[i-1].x, nx, c.p);
        u128 ny = mod_mul(lambda, dx, &c);
        ny = mod_sub(ny, M[i-1].y, c.p);
        M[i].x = nx; M[i].y = ny; M[i].is_infinity = 0;
    }
    int num_beats = atoi(argv[6]);
    int n_threads = (argc >= 8) ? atoi(argv[7]) : 1;

    thread_t T[16];
    pthread_t tids[16];
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (int i = 0; i < n_threads; i++) {
        T[i].tid = i; T[i].num_beats = num_beats; T[i].c = &c; T[i].M = M;
        pthread_create(&tids[i], NULL, worker, &T[i]);
    }
    for (int i = 0; i < n_threads; i++) pthread_join(tids[i], NULL);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double wall = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    /* Each beat = WALKS_PER_THREAD steps */
    long total_steps = (long)n_threads * num_beats * WALKS_PER_THREAD;
    fprintf(stderr, "Batch-inv-rho (N=%d): %d threads x %d beats x %d walks = %ld steps in %.3f s = %.2e steps/sec\n",
            WALKS_PER_THREAD, n_threads, num_beats, WALKS_PER_THREAD, total_steps, wall, total_steps / wall);
    printf("%.2e\n", total_steps / wall);
    return 0;
}
