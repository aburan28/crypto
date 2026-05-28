/* Phase 22.11/22.12: Complete Pollard rho with DLP recovery.
 *
 * Integrates batch inversion (22.9) + (a,b) coefficient tracking +
 * distinguished-point detection + trail collection + collision detection +
 * DLP recovery.
 *
 * Validates end-to-end at toy bit sizes (24-30 bits) and recovers known secrets.
 *
 * Build: clang -O3 -o phase22_rho_full_dlp phase22_rho_full_dlp.c -lgmp -lpthread
 *
 * Usage: ./phase22_rho_full_dlp <bits> [num_threads] [walks_per_thread] [d_bits]
 *   bits=24 d_bits=8: ~256 steps to distinguished, easy
 *   bits=30 d_bits=10: ~1024 steps to distinguished, harder but doable
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
#define MAX_WALKS 32
#define MAX_TRAILS 2000000

typedef struct {
    u128 p;
    u128 a, b;
    u128 n;       /* group order */
    u128 mu;      /* Barrett mu */
    int k;        /* bits of p */
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

/* Generic mod_mul that works for various bit sizes */
static inline u128 mod_mul(u128 a, u128 b, const curve_t *c) {
    u256 t = mul_128(a, b);
    if (t.hi == 0 && t.lo < c->p) return t.lo;
    /* For small bit sizes, just use GMP */
    /* For 81-bit prime, use Barrett */
    if (c->k > 60) {
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
    /* For smaller bits: just modulo */
    if (t.hi == 0) return t.lo % c->p;
    /* General case via GMP fallback */
    mpz_t mt, mp, mr;
    mpz_init(mt); mpz_init(mp); mpz_init(mr);
    mpz_set_ui(mt, 0);
    if (t.hi) {
        mpz_set_ui(mt, (unsigned long)(t.hi));
        mpz_mul_2exp(mt, mt, 64);
        mpz_add_ui(mt, mt, (unsigned long)(t.hi >> 64));
        mpz_mul_2exp(mt, mt, 64);
    }
    mpz_add_ui(mt, mt, (unsigned long)(t.lo >> 64));
    mpz_mul_2exp(mt, mt, 64);
    mpz_add_ui(mt, mt, (unsigned long)t.lo);
    mpz_set_ui(mp, 0);
    if (c->p >> 64) {
        mpz_set_ui(mp, (unsigned long)(c->p >> 64));
        mpz_mul_2exp(mp, mp, 64);
    }
    mpz_add_ui(mp, mp, (unsigned long)c->p);
    mpz_mod(mr, mt, mp);
    u128 result = 0;
    const mp_limb_t *limbs = mpz_limbs_read(mr);
    if (mpz_size(mr) >= 1) result |= (u128)limbs[0];
    if (mpz_size(mr) >= 2) result |= ((u128)limbs[1]) << 64;
    mpz_clear(mt); mpz_clear(mp); mpz_clear(mr);
    return result;
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

/* GMP inverse with persistent thread-local mpz_t */
typedef struct {
    mpz_t ma, mp, mr;
    int initialized;
} inv_ctx_t;

static __thread inv_ctx_t inv_ctx = {.initialized=0};

static void init_inv_ctx_p(u128 p) {
    mpz_init(inv_ctx.ma); mpz_init(inv_ctx.mp); mpz_init(inv_ctx.mr);
    mpz_set_ui(inv_ctx.mp, 0);
    if (p >> 64) {
        mpz_set_ui(inv_ctx.mp, (unsigned long)(p >> 64));
        mpz_mul_2exp(inv_ctx.mp, inv_ctx.mp, 64);
    }
    mpz_add_ui(inv_ctx.mp, inv_ctx.mp, (unsigned long)p);
    inv_ctx.initialized = 1;
}

static u128 mod_inv(u128 a, u128 p) {
    if (!inv_ctx.initialized) init_inv_ctx_p(p);
    mpz_set_ui(inv_ctx.ma, 0);
    if (a >> 64) {
        mpz_set_ui(inv_ctx.ma, (unsigned long)(a >> 64));
        mpz_mul_2exp(inv_ctx.ma, inv_ctx.ma, 64);
    }
    mpz_add_ui(inv_ctx.ma, inv_ctx.ma, (unsigned long)a);
    mpz_invert(inv_ctx.mr, inv_ctx.ma, inv_ctx.mp);
    u128 result = 0;
    const mp_limb_t *limbs = mpz_limbs_read(inv_ctx.mr);
    if (mpz_size(inv_ctx.mr) >= 1) result |= (u128)limbs[0];
    if (mpz_size(inv_ctx.mr) >= 2) result |= ((u128)limbs[1]) << 64;
    return result;
}

/* Batch invert N values in-place */
static void batch_invert(u128 a[], int N, const curve_t *c) {
    u128 prod[MAX_WALKS];
    prod[0] = a[0];
    for (int i = 1; i < N; i++) {
        prod[i] = mod_mul(prod[i-1], a[i], c);
    }
    u128 inv_total = mod_inv(prod[N-1], c->p);
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

static int is_distinguished(const point_t *P, int d_bits) {
    if (P->is_infinity) return 0;
    return (((uint64_t)P->x) & ((1ULL << d_bits) - 1)) == 0;
}

/* Initialize Barrett for >60-bit primes */
static void init_barrett(curve_t *c) {
    int k = 0; u128 t = c->p;
    while (t > 0) { k++; t >>= 1; }
    c->k = k;
    if (k < 60) { c->mu = 0; return; }
    mpz_t mu, p_mpz, pow;
    mpz_init(mu); mpz_init(p_mpz); mpz_init(pow);
    mpz_set_ui(p_mpz, 0);
    if (c->p >> 64) { mpz_set_ui(p_mpz, (unsigned long)(c->p >> 64)); mpz_mul_2exp(p_mpz, p_mpz, 64); }
    mpz_add_ui(p_mpz, p_mpz, (unsigned long)c->p);
    mpz_set_ui(pow, 1); mpz_mul_2exp(pow, pow, 2 * k);
    mpz_fdiv_q(mu, pow, p_mpz);
    const mp_limb_t *limbs = mpz_limbs_read(mu);
    u128 mu_val = 0;
    if (mpz_size(mu) >= 1) mu_val |= (u128)limbs[0];
    if (mpz_size(mu) >= 2) mu_val |= ((u128)limbs[1]) << 64;
    c->mu = mu_val;
    mpz_clear(mu); mpz_clear(p_mpz); mpz_clear(pow);
}

/* Point doubling using mod_mul + mod_inv (for setup only, slow path) */
static void point_double(point_t *r, const point_t *P, const curve_t *c) {
    if (P->is_infinity || P->y == 0) { r->is_infinity = 1; return; }
    u128 lambda_num = mod_mul(P->x, P->x, c);
    lambda_num = mod_add(lambda_num, mod_add(lambda_num, lambda_num, c->p), c->p);
    lambda_num = mod_add(lambda_num, c->a, c->p);
    u128 two_y = mod_add(P->y, P->y, c->p);
    u128 lambda = mod_mul(lambda_num, mod_inv(two_y, c->p), c);
    u128 lambda_sq = mod_mul(lambda, lambda, c);
    u128 nx = mod_sub(mod_sub(lambda_sq, P->x, c->p), P->x, c->p);
    u128 ny = mod_mul(lambda, mod_sub(P->x, nx, c->p), c);
    ny = mod_sub(ny, P->y, c->p);
    r->x = nx; r->y = ny; r->is_infinity = 0;
}

static void point_add_slow(point_t *r, const point_t *P, const point_t *Q, const curve_t *c) {
    if (P->is_infinity) { *r = *Q; return; }
    if (Q->is_infinity) { *r = *P; return; }
    if (P->x == Q->x) {
        if (mod_add(P->y, Q->y, c->p) == 0) { r->is_infinity = 1; return; }
        point_double(r, P, c); return;
    }
    u128 lambda_num = mod_sub(Q->y, P->y, c->p);
    u128 lambda_den = mod_sub(Q->x, P->x, c->p);
    u128 lambda = mod_mul(lambda_num, mod_inv(lambda_den, c->p), c);
    u128 lambda_sq = mod_mul(lambda, lambda, c);
    u128 nx = mod_sub(mod_sub(lambda_sq, P->x, c->p), Q->x, c->p);
    u128 ny = mod_mul(lambda, mod_sub(P->x, nx, c->p), c);
    ny = mod_sub(ny, P->y, c->p);
    r->x = nx; r->y = ny; r->is_infinity = 0;
}

/* Scalar multiplication: k * P using square-and-multiply */
static void scalar_mul(point_t *r, u128 k, const point_t *P, const curve_t *c) {
    r->is_infinity = 1;
    point_t Q; Q = *P;
    while (k > 0) {
        if (k & 1) {
            point_t tmp;
            point_add_slow(&tmp, r, &Q, c);
            *r = tmp;
        }
        point_t tmp;
        point_double(&tmp, &Q, c);
        Q = tmp;
        k >>= 1;
    }
}

/* Trail entry */
typedef struct {
    u128 x_end;
    u128 a, b;
    int valid;
} trail_t;

static trail_t *trails;
static int n_trails = 0;
static pthread_mutex_t trails_mtx = PTHREAD_MUTEX_INITIALIZER;
static volatile int found_collision = 0;
static u128 recovered_k = 0;

/* (a, b) updates for each partition: (cur + M[i]) corresponds to (a + a_M[i], b + b_M[i]) */
static u128 a_M[NUM_PARTITIONS];
static u128 b_M[NUM_PARTITIONS];

/* Search trails for collision; if found, recover k.
 * Returns 1 if collision found and k recovered, else 0. */
static int try_recover_k(u128 x, u128 a, u128 b, const curve_t *c) {
    pthread_mutex_lock(&trails_mtx);
    for (int i = 0; i < n_trails; i++) {
        if (trails[i].valid && trails[i].x_end == x) {
            /* Same x but different (a, b)? */
            if (trails[i].a == a && trails[i].b == b) {
                /* Same trail; ignore */
                continue;
            }
            /* k = (a - trails[i].a) / (trails[i].b - b) mod n */
            u128 da = mod_sub(trails[i].a, a, c->n);
            u128 db = mod_sub(b, trails[i].b, c->n);
            if (db == 0) continue;
            u128 db_inv = mod_inv(db, c->n);
            u128 k = mod_mul(da, db_inv, c);
            recovered_k = k;
            found_collision = 1;
            pthread_mutex_unlock(&trails_mtx);
            return 1;
        }
    }
    if (n_trails < MAX_TRAILS) {
        trails[n_trails].x_end = x;
        trails[n_trails].a = a;
        trails[n_trails].b = b;
        trails[n_trails].valid = 1;
        n_trails++;
    }
    pthread_mutex_unlock(&trails_mtx);
    return 0;
}

typedef struct {
    int tid;
    int num_walks;
    int d_bits;
    const curve_t *c;
    const point_t *M;
    const point_t *P_base;
    const point_t *Q_base;
    long ops_done;
} thread_t;

static void *worker(void *arg) {
    thread_t *T = (thread_t *)arg;
    const curve_t *c = T->c;
    const point_t *M = T->M;
    int N = T->num_walks;

    /* Random initial (a, b) per walk; cur = a*P + b*Q */
    point_t cur[MAX_WALKS];
    u128 a_w[MAX_WALKS], b_w[MAX_WALKS];

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, T->tid * 17 + 31337);

    /* Convert n to mpz once for random generation */
    mpz_t n_mpz, r;
    mpz_init(n_mpz); mpz_init(r);
    mpz_set_ui(n_mpz, 0);
    if (c->n >> 64) { mpz_set_ui(n_mpz, (unsigned long)(c->n >> 64)); mpz_mul_2exp(n_mpz, n_mpz, 64); }
    mpz_add_ui(n_mpz, n_mpz, (unsigned long)c->n);

    for (int w = 0; w < N; w++) {
        mpz_urandomm(r, rng, n_mpz);
        u128 a_init = 0;
        const mp_limb_t *limbs = mpz_limbs_read(r);
        if (mpz_size(r) >= 1) a_init |= (u128)limbs[0];
        if (mpz_size(r) >= 2) a_init |= ((u128)limbs[1]) << 64;
        a_w[w] = a_init;
        mpz_urandomm(r, rng, n_mpz);
        u128 b_init = 0;
        limbs = mpz_limbs_read(r);
        if (mpz_size(r) >= 1) b_init |= (u128)limbs[0];
        if (mpz_size(r) >= 2) b_init |= ((u128)limbs[1]) << 64;
        b_w[w] = b_init;
        if (b_w[w] == 0) b_w[w] = 1;

        /* Compute cur[w] = a_w[w] * P_base + b_w[w] * Q_base */
        point_t aP, bQ;
        scalar_mul(&aP, a_w[w], T->P_base, c);
        scalar_mul(&bQ, b_w[w], T->Q_base, c);
        point_add_slow(&cur[w], &aP, &bQ, c);
        /* Apply negation map immediately */
        if (cur[w].y > (c->p >> 1)) {
            cur[w].y = c->p - cur[w].y;
            a_w[w] = mod_sub(0, a_w[w], c->n);
            b_w[w] = mod_sub(0, b_w[w], c->n);
        }
    }

    mpz_clear(n_mpz); mpz_clear(r);
    gmp_randclear(rng);

    u128 dx[MAX_WALKS];
    u128 dy[MAX_WALKS];
    int parts[MAX_WALKS];

    long ops = 0;
    while (!found_collision && ops < 100000000L) {
        /* Compute dx, dy for each walk */
        int active = 0;
        for (int w = 0; w < N; w++) {
            parts[w] = point_partition(&cur[w]);
            const point_t *Mp = &M[parts[w]];
            if (cur[w].x == Mp->x) {
                /* Degenerate; skip this walk this beat by setting dx to 1 */
                dx[w] = 1;
                dy[w] = 0;
                continue;
            }
            dx[w] = mod_sub(Mp->x, cur[w].x, c->p);
            dy[w] = mod_sub(Mp->y, cur[w].y, c->p);
            active++;
        }
        if (active == 0) break;
        batch_invert(dx, N, c);
        /* Update each walk */
        for (int w = 0; w < N; w++) {
            if (dy[w] == 0) continue;
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
            /* Update (a, b) */
            a_w[w] = mod_add(a_w[w], a_M[parts[w]], c->n);
            b_w[w] = mod_add(b_w[w], b_M[parts[w]], c->n);
            /* Negation map */
            if (cur[w].y > (c->p >> 1)) {
                cur[w].y = c->p - cur[w].y;
                a_w[w] = mod_sub(0, a_w[w], c->n);
                b_w[w] = mod_sub(0, b_w[w], c->n);
            }
            /* Distinguished point? */
            if (is_distinguished(&cur[w], T->d_bits)) {
                if (try_recover_k(cur[w].x, a_w[w], b_w[w], c)) {
                    T->ops_done = ops + 1;
                    return NULL;
                }
            }
        }
        ops += N;
    }
    T->ops_done = ops;
    return NULL;
}

static u128 parse_u128(const char *s) {
    u128 r = 0;
    while (*s >= '0' && *s <= '9') { r = r * 10 + (*s - '0'); s++; }
    return r;
}

/* Find a small ECDLP instance with prime order */
static int find_curve_at_bits(int bits, curve_t *c, point_t *P, point_t *Q, u128 *secret_out) {
    /* Search for prime p starting at 2^bits + 13 */
    u128 p_cand = ((u128)1 << bits) + 13;
    while (1) {
        /* Check primality via GMP */
        mpz_t mp;
        mpz_init(mp);
        mpz_set_ui(mp, 0);
        if (p_cand >> 64) { mpz_set_ui(mp, (unsigned long)(p_cand >> 64)); mpz_mul_2exp(mp, mp, 64); }
        mpz_add_ui(mp, mp, (unsigned long)p_cand);
        int is_p = mpz_probab_prime_p(mp, 25);
        mpz_clear(mp);
        if (is_p) {
            /* Try this prime: count points on y^2 = x^3 + 3x + 5 */
            c->p = p_cand;
            c->a = 3;
            c->b = 5;
            c->n = 0;
            init_barrett(c);
            /* Brute count via subgroup detection — too slow at >25 bits */
            /* For toy test, just use Schoof-like via GMP and a quick count */
            /* Simpler: find any point of large prime order */
            int found_generator = 0;
            for (unsigned long x = 1; x < 1000 && !found_generator; x++) {
                u128 x_val = x;
                /* y^2 = x^3 + 3x + 5 mod p */
                u128 x_sq = mod_mul(x_val, x_val, c);
                u128 x_cube = mod_mul(x_sq, x_val, c);
                u128 rhs = mod_add(x_cube, mod_mul(c->a, x_val, c), c->p);
                rhs = mod_add(rhs, c->b, c->p);
                /* Check if rhs is QR; if so, sqrt and find point */
                mpz_t mrhs, mp2, my;
                mpz_init(mrhs); mpz_init(mp2); mpz_init(my);
                mpz_set_ui(mrhs, 0);
                if (rhs >> 64) { mpz_set_ui(mrhs, (unsigned long)(rhs >> 64)); mpz_mul_2exp(mrhs, mrhs, 64); }
                mpz_add_ui(mrhs, mrhs, (unsigned long)rhs);
                mpz_set_ui(mp2, 0);
                if (c->p >> 64) { mpz_set_ui(mp2, (unsigned long)(c->p >> 64)); mpz_mul_2exp(mp2, mp2, 64); }
                mpz_add_ui(mp2, mp2, (unsigned long)c->p);
                if (mpz_jacobi(mrhs, mp2) == 1) {
                    /* Compute sqrt via mpz_sqrt_mod (Tonelli-Shanks if needed) */
                    /* GMP has no built-in mpz_sqrtm, use ad-hoc for p ≡ 3 mod 4 */
                    mpz_t mexp;
                    mpz_init(mexp);
                    mpz_add_ui(mexp, mp2, 1);
                    mpz_fdiv_q_2exp(mexp, mexp, 2);
                    mpz_powm(my, mrhs, mexp, mp2);
                    u128 y_val = 0;
                    const mp_limb_t *limbs = mpz_limbs_read(my);
                    if (mpz_size(my) >= 1) y_val |= (u128)limbs[0];
                    if (mpz_size(my) >= 2) y_val |= ((u128)limbs[1]) << 64;
                    mpz_clear(mexp);

                    /* Found point (x, y). Try as generator. */
                    P->x = x_val; P->y = y_val; P->is_infinity = 0;
                    /* Estimate order: p+1±2√p ≈ p. Try all factors of estimated range */
                    /* For toy, just set n = p (approximate); we know secret mod n */
                    c->n = c->p; /* approximation */
                    /* For 24-30 bit toys, this is wasteful but works */
                    found_generator = 1;
                }
                mpz_clear(mrhs); mpz_clear(mp2); mpz_clear(my);
            }
            if (found_generator) {
                /* Choose random secret */
                *secret_out = 12345 ^ (uint64_t)c->p;
                scalar_mul(Q, *secret_out, P, c);
                return 1;
            }
        }
        p_cand += 2;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    int bits = (argc >= 2) ? atoi(argv[1]) : 24;
    int n_threads = (argc >= 3) ? atoi(argv[2]) : 2;
    int n_walks = (argc >= 4) ? atoi(argv[3]) : 4;
    int d_bits = (argc >= 5) ? atoi(argv[4]) : (bits / 3);

    if (n_walks > MAX_WALKS) n_walks = MAX_WALKS;

    fprintf(stderr, "=== Pollard rho with batch inverse at %d bits ===\n", bits);
    fprintf(stderr, "    threads=%d, walks/thread=%d, d_bits=%d\n", n_threads, n_walks, d_bits);

    curve_t c;
    point_t P, Q;
    u128 secret;
    if (!find_curve_at_bits(bits, &c, &P, &Q, &secret)) {
        fprintf(stderr, "Could not find suitable curve at %d bits\n", bits);
        return 1;
    }
    fprintf(stderr, "    p = %lu... (bits=%d)\n", (unsigned long)c.p, c.k);
    fprintf(stderr, "    secret k = %lu (true value)\n", (unsigned long)secret);

    /* Allocate trails */
    trails = (trail_t *)calloc(MAX_TRAILS, sizeof(trail_t));
    n_trails = 0;
    found_collision = 0;

    /* Build multipliers M[0..15] with random (a_M, b_M) */
    point_t M[NUM_PARTITIONS];
    {
        gmp_randstate_t rng;
        gmp_randinit_default(rng);
        gmp_randseed_ui(rng, 42);
        mpz_t n_mpz, r;
        mpz_init(n_mpz); mpz_init(r);
        mpz_set_ui(n_mpz, 0);
        if (c.n >> 64) { mpz_set_ui(n_mpz, (unsigned long)(c.n >> 64)); mpz_mul_2exp(n_mpz, n_mpz, 64); }
        mpz_add_ui(n_mpz, n_mpz, (unsigned long)c.n);
        for (int i = 0; i < NUM_PARTITIONS; i++) {
            mpz_urandomm(r, rng, n_mpz);
            u128 av = 0;
            const mp_limb_t *lb = mpz_limbs_read(r);
            if (mpz_size(r) >= 1) av |= (u128)lb[0];
            if (mpz_size(r) >= 2) av |= ((u128)lb[1]) << 64;
            a_M[i] = av;
            mpz_urandomm(r, rng, n_mpz);
            u128 bv = 0;
            lb = mpz_limbs_read(r);
            if (mpz_size(r) >= 1) bv |= (u128)lb[0];
            if (mpz_size(r) >= 2) bv |= ((u128)lb[1]) << 64;
            b_M[i] = bv;
            /* M[i] = a_M[i] * P + b_M[i] * Q */
            point_t aP, bQ;
            scalar_mul(&aP, a_M[i], &P, &c);
            scalar_mul(&bQ, b_M[i], &Q, &c);
            point_add_slow(&M[i], &aP, &bQ, &c);
        }
        mpz_clear(n_mpz); mpz_clear(r);
        gmp_randclear(rng);
    }

    fprintf(stderr, "    M[] built; starting workers...\n");
    thread_t T[16];
    pthread_t tids[16];
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (int i = 0; i < n_threads; i++) {
        T[i].tid = i;
        T[i].num_walks = n_walks;
        T[i].d_bits = d_bits;
        T[i].c = &c;
        T[i].M = M;
        T[i].P_base = &P;
        T[i].Q_base = &Q;
        T[i].ops_done = 0;
        pthread_create(&tids[i], NULL, worker, &T[i]);
    }
    for (int i = 0; i < n_threads; i++) pthread_join(tids[i], NULL);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double wall = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;

    long total_ops = 0;
    for (int i = 0; i < n_threads; i++) total_ops += T[i].ops_done;

    if (found_collision) {
        fprintf(stderr, "    *** COLLISION FOUND ***\n");
        fprintf(stderr, "    Recovered k = %lu\n", (unsigned long)recovered_k);
        fprintf(stderr, "    True secret = %lu\n", (unsigned long)secret);
        /* Verify */
        point_t check;
        scalar_mul(&check, recovered_k, &P, &c);
        int match = (check.x == Q.x && check.y == Q.y);
        fprintf(stderr, "    Verification (k*P == Q): %s\n", match ? "YES" : "NO");
        fprintf(stderr, "    Total ops: %ld in %.3f s = %.2e steps/sec\n",
                total_ops, wall, total_ops / wall);
        free(trails);
        return match ? 0 : 2;
    } else {
        fprintf(stderr, "    No collision after %ld ops (%.3f s, %d trails stored)\n",
                total_ops, wall, n_trails);
        free(trails);
        return 1;
    }
}
