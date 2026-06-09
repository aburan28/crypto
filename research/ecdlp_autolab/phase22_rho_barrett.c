/* Phase 22.2: Pollard rho walk with Barrett reduction + batched Montgomery
 * inversion.
 *
 * This fuses the two engineering wins that the v7 roadmap flagged but never
 * integrated into the affine walk:
 *
 *   1. Barrett reduction (phase22_barrett.c proved 9.5x on mulmod in
 *      isolation, but it was never wired into the rho point arithmetic).
 *
 *   2. Batched modular inversion (Montgomery's trick).  The v3 affine walk
 *      pays one GMP mpz_invert per step, which dominates the per-add cost
 *      (~80 muls via Fermat).  By running W independent walks per thread and
 *      batching their denominators, W inversions collapse to 1 inversion +
 *      ~3W multiplies.  Amortized, each walk's inversion cost -> ~3 muls.
 *
 * This is the CPU embodiment of the GPU strategy: many independent walks,
 * special-form (Barrett) reduction, and a single shared inversion per batch.
 *
 * Field arithmetic is fixed to the 81-bit benchmark prime via __uint128_t,
 * so there is no GMP in the hot loop.  GMP is used only at startup for the
 * one-time Barrett mu precomputation and the correctness self-checks.
 *
 * Build: gcc -O3 -o phase22_rho_barrett phase22_rho_barrett.c -lgmp -lpthread
 * Run:   ./phase22_rho_barrett <p> <a> <b> <xP> <yP> <num_steps> \
 *                              [num_threads] [walks_per_thread]
 *
 * Example (67.a1 benchmark, 4 threads, 512 walks/thread, 200k steps):
 *   ./phase22_rho_barrett 1208925819614629469615699 12 5 \
 *     942402533238732514353214 807227645196494903606326 200000 4 512
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
#define MAX_THREADS 64

/* ----- u128 parse / print helpers (shared with phase22_barrett.c) ------- */

static u128 parse_u128(const char *s) {
    u128 r = 0;
    while (*s >= '0' && *s <= '9') { r = r * 10 + (u128)(*s - '0'); s++; }
    return r;
}

static void print_u128(u128 x) {
    char buf[50];
    int i = 49;
    buf[i--] = 0;
    if (x == 0) buf[i--] = '0';
    else while (x > 0) { buf[i--] = '0' + (int)(x % 10); x /= 10; }
    fprintf(stderr, "%s", &buf[i + 1]);
}

/* ----- 128x128 -> 256-bit product (shared with phase22_barrett.c) ------- */

typedef struct { u128 hi, lo; } u256;

static u256 mul_128(u128 a, u128 b) {
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
    u256 r = { result_hi, result_lo };
    return r;
}

/* ----- Barrett reduction state ------------------------------------------ */

typedef struct {
    u128 p;
    u128 mu;   /* floor(2^{2k} / p) */
    int k;     /* ceil(log2 p) */
    u128 half; /* floor(p/2), for the negation map */
} field_t;

static void field_init(field_t *f, u128 p) {
    f->p = p;
    f->half = p >> 1;
    int k = 0; u128 t = p;
    while (t > 0) { k++; t >>= 1; }
    f->k = k;

    mpz_t mu_mpz, p_mpz, pow_mpz, tmp;
    mpz_init(mu_mpz); mpz_init(p_mpz); mpz_init(pow_mpz); mpz_init(tmp);
    uint64_t plo = (uint64_t)p, phi = (uint64_t)(p >> 64);
    mpz_import(p_mpz, 1, 1, sizeof(uint64_t), 0, 0, &phi);
    mpz_mul_2exp(p_mpz, p_mpz, 64);
    mpz_import(tmp, 1, 1, sizeof(uint64_t), 0, 0, &plo);
    mpz_add(p_mpz, p_mpz, tmp);
    mpz_set_ui(pow_mpz, 1);
    mpz_mul_2exp(pow_mpz, pow_mpz, 2 * k);
    mpz_fdiv_q(mu_mpz, pow_mpz, p_mpz);
    if (mpz_sizeinbase(mu_mpz, 2) > 128) {
        fprintf(stderr, "ERROR: mu doesn't fit in u128\n");
        exit(1);
    }
    uint64_t mu_lo = 0, mu_hi = 0;
    const mp_limb_t *limbs = mpz_limbs_read(mu_mpz);
    mu_lo = (uint64_t)limbs[0];
    if (mpz_size(mu_mpz) >= 2) mu_hi = (uint64_t)limbs[1];
    f->mu = ((u128)mu_hi << 64) | mu_lo;
    mpz_clear(mu_mpz); mpz_clear(p_mpz); mpz_clear(pow_mpz); mpz_clear(tmp);
}

/* a * b mod p, Barrett.  Requires a, b < p (81-bit). 2k = 162 for k=81. */
static inline u128 fmul(const field_t *f, u128 a, u128 b) {
    u256 t = mul_128(a, b);
    int sh = 2 * f->k - 128;                  /* = 34 for k=81 */
    u256 hi_part = mul_128(t.hi, f->mu);
    u128 q_approx = (hi_part.hi << (128 - sh)) | (hi_part.lo >> sh);
    u256 lo_part = mul_128(t.lo, f->mu);
    u128 q_total = q_approx + (lo_part.hi >> sh);
    u256 qp = mul_128(q_total, f->p);
    u128 r = t.lo - qp.lo;
    while (r >= f->p) r -= f->p;
    return r;
}

static inline u128 fadd(const field_t *f, u128 a, u128 b) {
    u128 s = a + b;
    if (s >= f->p) s -= f->p;
    return s;
}

static inline u128 fsub(const field_t *f, u128 a, u128 b) {
    return (a >= b) ? (a - b) : (a + f->p - b);
}

/* x^e mod p via Barrett square-and-multiply. */
static u128 fpow(const field_t *f, u128 x, u128 e) {
    u128 r = 1 % f->p;
    while (e) {
        if (e & 1) r = fmul(f, r, x);
        x = fmul(f, x, x);
        e >>= 1;
    }
    return r;
}

/* Modular inverse via Fermat: x^{p-2} mod p.  Used once per batch. */
static inline u128 finv(const field_t *f, u128 x) {
    return fpow(f, x, f->p - 2);
}

/* Batched inversion (Montgomery's trick).  d[0..n) all nonzero -> out[0..n).
 * Cost: 1 inversion + 3(n-1) muls.  acc is caller-provided scratch[n]. */
static void batch_invert(const field_t *f, const u128 *d, u128 *out,
                         u128 *acc, int n) {
    acc[0] = d[0];
    for (int j = 1; j < n; j++) acc[j] = fmul(f, acc[j - 1], d[j]);
    u128 inv = finv(f, acc[n - 1]);
    for (int j = n - 1; j >= 1; j--) {
        out[j] = fmul(f, inv, acc[j - 1]);
        inv = fmul(f, inv, d[j]);
    }
    out[0] = inv;
}

/* ----- Affine points ---------------------------------------------------- */

typedef struct { u128 x, y; int inf; } apoint;

/* Single affine addition with one Fermat inversion (bootstrap only). */
static apoint affine_add(const field_t *f, u128 a_coef, apoint P, apoint Q) {
    if (P.inf) return Q;
    if (Q.inf) return P;
    u128 lambda, num, den;
    if (P.x == Q.x) {
        if (fadd(f, P.y, Q.y) == 0) { apoint R = {0, 0, 1}; return R; }
        num = fadd(f, fmul(f, 3 % f->p, fmul(f, P.x, P.x)), a_coef);
        den = fadd(f, P.y, P.y);
    } else {
        num = fsub(f, Q.y, P.y);
        den = fsub(f, Q.x, P.x);
    }
    lambda = fmul(f, num, finv(f, den));
    u128 x3 = fsub(f, fsub(f, fmul(f, lambda, lambda), P.x), Q.x);
    u128 y3 = fsub(f, fmul(f, lambda, fsub(f, P.x, x3)), P.y);
    apoint R = { x3, y3, 0 };
    return R;
}

/* ----- Worker ----------------------------------------------------------- */

typedef struct {
    int tid;
    long num_steps;
    int walks;
    const field_t *f;
    u128 a_coef;
    const apoint *M;     /* [NUM_PARTITIONS] */
    const apoint *start; /* [walks] distinct seeds for this thread */
    long ops_done;
} thread_t;

/* kind codes for the two-pass step */
enum { K_ADD = 0, K_DOUBLE = 1, K_NEUTRAL = 2 };

static void *worker(void *arg) {
    thread_t *T = (thread_t *)arg;
    const field_t *f = T->f;
    int W = T->walks;
    u128 a_coef = T->a_coef;
    const apoint *M = T->M;

    apoint *cur   = malloc(sizeof(apoint) * W);
    u128 *den     = malloc(sizeof(u128) * W);
    u128 *num     = malloc(sizeof(u128) * W);
    u128 *invd    = malloc(sizeof(u128) * W);
    u128 *acc     = malloc(sizeof(u128) * W);
    int  *kind    = malloc(sizeof(int) * W);
    apoint *nxt   = malloc(sizeof(apoint) * W);  /* neutral-lane results */

    for (int j = 0; j < W; j++) cur[j] = T->start[j];

    for (long s = 0; s < T->num_steps; s++) {
        /* Pass 1: classify each lane and gather denominators. */
        for (int j = 0; j < W; j++) {
            apoint c = cur[j];
            if (c.inf) {
                nxt[j] = M[0];          /* infinity + M[0] */
                kind[j] = K_NEUTRAL;
                den[j] = 1;
                continue;
            }
            int part = (int)((uint64_t)c.x & (NUM_PARTITIONS - 1));
            apoint Q = M[part];
            if (c.x == Q.x) {
                if (fadd(f, c.y, Q.y) == 0) {
                    apoint inf = {0, 0, 1};
                    nxt[j] = inf;
                    kind[j] = K_NEUTRAL;
                    den[j] = 1;
                } else {
                    kind[j] = K_DOUBLE;
                    num[j] = fadd(f, fmul(f, 3 % f->p, fmul(f, c.x, c.x)), a_coef);
                    den[j] = fadd(f, c.y, c.y);
                }
            } else {
                kind[j] = K_ADD;
                num[j] = fsub(f, Q.y, c.y);
                den[j] = fsub(f, Q.x, c.x);
            }
        }

        /* One inversion for the whole batch. */
        batch_invert(f, den, invd, acc, W);

        /* Pass 2: finish each lane and fold in the negation map. */
        for (int j = 0; j < W; j++) {
            if (kind[j] == K_NEUTRAL) {
                apoint r = nxt[j];
                if (!r.inf && r.y > f->half) r.y = f->p - r.y;
                cur[j] = r;
                continue;
            }
            apoint c = cur[j];
            int part = (int)((uint64_t)c.x & (NUM_PARTITIONS - 1));
            apoint Q = M[part];
            u128 lambda = fmul(f, num[j], invd[j]);
            u128 x3, y3;
            if (kind[j] == K_ADD)
                x3 = fsub(f, fsub(f, fmul(f, lambda, lambda), c.x), Q.x);
            else
                x3 = fsub(f, fsub(f, fmul(f, lambda, lambda), c.x), c.x);
            y3 = fsub(f, fmul(f, lambda, fsub(f, c.x, x3)), c.y);
            if (y3 > f->half) y3 = f->p - y3;  /* negation map */
            cur[j].x = x3; cur[j].y = y3; cur[j].inf = 0;
        }
    }

    T->ops_done = T->num_steps * (long)W;
    free(cur); free(den); free(num); free(invd); free(acc); free(kind); free(nxt);
    return NULL;
}

/* ----- Correctness self-checks ------------------------------------------ */

static void u128_to_mpz(mpz_t out, u128 v) {
    uint64_t lo = (uint64_t)v, hi = (uint64_t)(v >> 64);
    mpz_t tmp; mpz_init(tmp);
    mpz_import(out, 1, 1, 8, 0, 0, &hi);
    mpz_mul_2exp(out, out, 64);
    mpz_import(tmp, 1, 1, 8, 0, 0, &lo);
    mpz_add(out, out, tmp);
    mpz_clear(tmp);
}

static u128 mpz_to_u128(const mpz_t v) {
    u128 r = 0;
    const mp_limb_t *limbs = mpz_limbs_read(v);
    if (mpz_size(v) >= 1) r |= (u128)limbs[0];
    if (mpz_size(v) >= 2) r |= ((u128)limbs[1]) << 64;
    return r;
}

/* Check batched inversion against GMP on random inputs. */
static int check_batch_invert(const field_t *f) {
    int n = 64;
    u128 *d = malloc(sizeof(u128) * n);
    u128 *out = malloc(sizeof(u128) * n);
    u128 *acc = malloc(sizeof(u128) * n);
    srand(12345);
    for (int j = 0; j < n; j++) {
        u128 v = ((u128)rand() << 96) ^ ((u128)rand() << 64) ^
                 ((u128)rand() << 32) ^ (u128)rand();
        v %= f->p;
        if (v == 0) v = 1;
        d[j] = v;
    }
    batch_invert(f, d, out, acc, n);
    mpz_t pm, dm, inv; mpz_init(pm); mpz_init(dm); mpz_init(inv);
    u128_to_mpz(pm, f->p);
    int ok = 1;
    for (int j = 0; j < n; j++) {
        u128_to_mpz(dm, d[j]);
        mpz_invert(inv, dm, pm);
        if (mpz_to_u128(inv) != out[j]) { ok = 0; break; }
    }
    mpz_clear(pm); mpz_clear(dm); mpz_clear(inv);
    free(d); free(out); free(acc);
    return ok;
}

/* GMP reference walk: one lane, N steps, returns final affine point. */
static void gmp_ref_walk(const field_t *f, u128 a_coef, const apoint *M,
                         apoint start, long steps, mpz_t out_x, mpz_t out_y,
                         int *out_inf) {
    mpz_t p, a, cx, cy, qx, qy, t1, t2, lam, x3, y3, half;
    mpz_init(p); mpz_init(a); mpz_init(cx); mpz_init(cy);
    mpz_init(qx); mpz_init(qy); mpz_init(t1); mpz_init(t2);
    mpz_init(lam); mpz_init(x3); mpz_init(y3); mpz_init(half);
    u128_to_mpz(p, f->p); u128_to_mpz(a, a_coef); u128_to_mpz(half, f->half);
    u128_to_mpz(cx, start.x); u128_to_mpz(cy, start.y);
    int inf = start.inf;
    for (long s = 0; s < steps; s++) {
        int part = inf ? 0 : (int)((uint64_t)mpz_to_u128(cx) & (NUM_PARTITIONS - 1));
        u128_to_mpz(qx, M[part].x); u128_to_mpz(qy, M[part].y);
        if (inf) { mpz_set(cx, qx); mpz_set(cy, qy); inf = M[part].inf; }
        else {
            if (mpz_cmp(cx, qx) == 0) {
                mpz_add(t1, cy, qy); mpz_mod(t1, t1, p);
                if (mpz_sgn(t1) == 0) { inf = 1; goto neg; }
                mpz_mul(t1, cx, cx); mpz_mul_ui(t1, t1, 3);
                mpz_add(t1, t1, a); mpz_mod(t1, t1, p);
                mpz_add(t2, cy, cy); mpz_invert(t2, t2, p);
            } else {
                mpz_sub(t1, qy, cy);
                mpz_sub(t2, qx, cx); mpz_invert(t2, t2, p);
            }
            mpz_mul(lam, t1, t2); mpz_mod(lam, lam, p);
            mpz_mul(x3, lam, lam); mpz_sub(x3, x3, cx); mpz_sub(x3, x3, qx);
            mpz_mod(x3, x3, p);
            mpz_sub(y3, cx, x3); mpz_mul(y3, y3, lam); mpz_sub(y3, y3, cy);
            mpz_mod(y3, y3, p);
            mpz_set(cx, x3); mpz_set(cy, y3); inf = 0;
        }
neg:
        if (!inf && mpz_cmp(cy, half) > 0) mpz_sub(cy, p, cy);
    }
    mpz_set(out_x, cx); mpz_set(out_y, cy); *out_inf = inf;
    mpz_clear(p); mpz_clear(a); mpz_clear(cx); mpz_clear(cy);
    mpz_clear(qx); mpz_clear(qy); mpz_clear(t1); mpz_clear(t2);
    mpz_clear(lam); mpz_clear(x3); mpz_clear(y3); mpz_clear(half);
}

/* Run a single Barrett lane for N steps (mirrors worker, W=1). */
static apoint barrett_walk_one(const field_t *f, u128 a_coef, const apoint *M,
                               apoint start, long steps) {
    apoint c = start;
    for (long s = 0; s < steps; s++) {
        if (c.inf) { c = M[0]; if (!c.inf && c.y > f->half) c.y = f->p - c.y; continue; }
        int part = (int)((uint64_t)c.x & (NUM_PARTITIONS - 1));
        apoint Q = M[part];
        u128 num, den;
        int dbl = 0;
        if (c.x == Q.x) {
            if (fadd(f, c.y, Q.y) == 0) { apoint inf = {0,0,1}; c = inf; continue; }
            num = fadd(f, fmul(f, 3 % f->p, fmul(f, c.x, c.x)), a_coef);
            den = fadd(f, c.y, c.y); dbl = 1;
        } else {
            num = fsub(f, Q.y, c.y);
            den = fsub(f, Q.x, c.x);
        }
        u128 lambda = fmul(f, num, finv(f, den));
        u128 x3 = dbl
            ? fsub(f, fsub(f, fmul(f, lambda, lambda), c.x), c.x)
            : fsub(f, fsub(f, fmul(f, lambda, lambda), c.x), Q.x);
        u128 y3 = fsub(f, fmul(f, lambda, fsub(f, c.x, x3)), c.y);
        if (y3 > f->half) y3 = f->p - y3;
        c.x = x3; c.y = y3; c.inf = 0;
    }
    return c;
}

static int check_walk(const field_t *f, u128 a_coef, const apoint *M, apoint P) {
    long steps = 2000;
    apoint got = barrett_walk_one(f, a_coef, M, P, steps);
    mpz_t rx, ry; mpz_init(rx); mpz_init(ry); int rinf;
    gmp_ref_walk(f, a_coef, M, P, steps, rx, ry, &rinf);
    int ok = (got.inf == rinf) &&
             (rinf || (mpz_to_u128(rx) == got.x && mpz_to_u128(ry) == got.y));
    mpz_clear(rx); mpz_clear(ry);
    return ok;
}

/* ----- main ------------------------------------------------------------- */

int main(int argc, char *argv[]) {
    if (argc < 7) {
        fprintf(stderr,
            "Usage: %s p a b xP yP num_steps [num_threads] [walks_per_thread]\n",
            argv[0]);
        return 1;
    }
    u128 p = parse_u128(argv[1]);
    u128 a_coef = parse_u128(argv[2]);
    /* argv[3] = b: not needed for affine add/double, kept for CLI parity. */
    u128 xP = parse_u128(argv[4]);
    u128 yP = parse_u128(argv[5]);
    long num_steps = atol(argv[6]);
    int n_threads = (argc >= 8) ? atoi(argv[7]) : 4;
    int walks = (argc >= 9) ? atoi(argv[8]) : 256;
    if (n_threads < 1) n_threads = 1;
    if (n_threads > MAX_THREADS) n_threads = MAX_THREADS;
    if (walks < 1) walks = 1;

    field_t f;
    field_init(&f, p);

    /* Build multipliers M[0..15] = P, 2P, 4P, ... (matches v3 bootstrap). */
    apoint M[NUM_PARTITIONS];
    M[0].x = xP; M[0].y = yP; M[0].inf = 0;
    for (int i = 1; i < NUM_PARTITIONS; i++)
        M[i] = affine_add(&f, a_coef, M[i - 1], M[i - 1]);

    fprintf(stderr, "p = "); print_u128(p);
    fprintf(stderr, " (%d-bit), mu = ", f.k); print_u128(f.mu); fprintf(stderr, "\n");

    /* Self-checks before benchmarking. */
    if (!check_batch_invert(&f)) {
        fprintf(stderr, "FAIL: batched inversion != GMP mpz_invert\n");
        return 1;
    }
    if (!check_walk(&f, a_coef, M, M[0])) {
        fprintf(stderr, "FAIL: Barrett walk trajectory != GMP reference walk\n");
        return 1;
    }
    fprintf(stderr, "self-check OK: batched inverse and walk match GMP\n");

    /* Distinct seeds: thread t, lane j starts at (t*walks + j + 1) * P,
     * built by a running affine addition so every walk is independent. */
    int total_seeds = n_threads * walks;
    apoint *seeds = malloc(sizeof(apoint) * total_seeds);
    apoint running = M[0];
    seeds[0] = running;
    for (int i = 1; i < total_seeds; i++) {
        running = affine_add(&f, a_coef, running, M[0]);
        seeds[i] = running;
    }

    pthread_t tids[MAX_THREADS];
    thread_t T[MAX_THREADS];
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (int t = 0; t < n_threads; t++) {
        T[t].tid = t;
        T[t].num_steps = num_steps;
        T[t].walks = walks;
        T[t].f = &f;
        T[t].a_coef = a_coef;
        T[t].M = M;
        T[t].start = &seeds[t * walks];
        T[t].ops_done = 0;
        pthread_create(&tids[t], NULL, worker, &T[t]);
    }
    long total = 0;
    for (int t = 0; t < n_threads; t++) {
        pthread_join(tids[t], NULL);
        total += T[t].ops_done;
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double wall = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    double rate = total / wall;

    fprintf(stderr,
        "%d threads x %d walks x %ld steps = %ld point-ops in %.3f s = %.3e ops/sec\n",
        n_threads, walks, num_steps, total, wall, rate);
    printf("%.3e\n", rate);

    free(seeds);
    return 0;
}
