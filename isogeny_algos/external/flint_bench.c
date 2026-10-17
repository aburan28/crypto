/* FLINT side of the external comparison: the workloads, inputs and digests of
 * examples/ext_compare.rs (see that file for the definitions), one JSON line per case.
 *
 *   cc -O2 -o flint_bench flint_bench.c -lflint -lgmp
 *   ./flint_bench poly <p61|p127|p255|p511>
 *   ./flint_bench powx <p61|p127|p255|p511>
 *
 * p61 uses nmod_poly (FLINT's word-size type), the others fmpz_mod_poly. Single-threaded
 * (FLINT's default). X^p mod g includes the pre-inverse of g, as our side includes building
 * its reducer. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/nmod_poly.h>

static const char *P61 = "2305843009213693951";
static const char *P127 = "170141183460469231731687303715884105727";
static const char *P255 = "57896044618658097711785492504343953926634992332820282019728792003956564819949";
static const char *P511 = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";
static const long LENS[] = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384};
static const long DEGS[] = {32, 64, 90, 128, 256, 512, 1024};

static double now_us(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec * 1e6 + t.tv_nsec / 1e3;
}
static int cmpd(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
static long reps_for(long n) { return n <= 256 ? 50 : (n <= 2048 ? 9 : 5); }

/* ---- multi-word primes: fmpz_mod_poly ---- */
#define GOLD 11400714819323198485UL /* G of examples/ext_compare.rs */

/* coefficient i = G^(i + s) mod p (full-size residues) */
static void seq_fmpz(fmpz_mod_poly_t a, long n, ulong s, const fmpz_mod_ctx_t ctx) {
    fmpz_t g, c;
    fmpz_init(g);
    fmpz_init(c);
    fmpz_set_ui(g, GOLD);
    fmpz_mod_set_fmpz(g, g, ctx);
    fmpz_mod_pow_ui(c, g, s, ctx);
    fmpz_mod_poly_zero(a, ctx);
    fmpz_mod_poly_fit_length(a, n, ctx);
    for (long i = 0; i < n; i++) {
        fmpz_mod_poly_set_coeff_fmpz(a, i, c, ctx);
        fmpz_mod_mul(c, c, g, ctx);
    }
    fmpz_clear(g);
    fmpz_clear(c);
}
static void digest_fmpz(const fmpz_mod_poly_t a, const fmpz_mod_ctx_t ctx, char **out) {
    fmpz_t r, v;
    fmpz_init_set_ui(r, 123456789);
    fmpz_init(v);
    fmpz_mod_poly_evaluate_fmpz(v, a, r, ctx);
    *out = fmpz_get_str(NULL, 10, v);
    fmpz_clear(r);
    fmpz_clear(v);
}

static void poly_fmpz(const char *tag, const char *pdec, int max_n) {
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_str(p, pdec, 10);
    fmpz_mod_ctx_t ctx;
    fmpz_mod_ctx_init(ctx, p);
    for (size_t k = 0; k < sizeof LENS / sizeof *LENS; k++) {
        long n = LENS[k];
        if (n > max_n) continue;
        fmpz_mod_poly_t a, b, r;
        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(r, ctx);
        seq_fmpz(a, n, 16, ctx);
        seq_fmpz(b, n, 1000003, ctx);
        long reps = reps_for(n);
        double *ts = malloc(sizeof(double) * reps);
        for (int op = 0; op < 2; op++) {
            for (long i = 0; i < reps; i++) {
                double t0 = now_us();
                if (op == 0) fmpz_mod_poly_mul(r, a, b, ctx);
                else fmpz_mod_poly_sqr(r, a, ctx);
                ts[i] = now_us() - t0;
            }
            qsort(ts, reps, sizeof(double), cmpd);
            char *d;
            digest_fmpz(r, ctx, &d);
            printf("{\"tool\":\"flint-%s\",\"op\":\"%s\",\"p\":\"%s\",\"n\":%ld,\"reps\":%ld,"
                   "\"best_us\":%.2f,\"median_us\":%.2f,\"digest\":\"%s\"}\n",
                   FLINT_VERSION, op == 0 ? "mul" : "sqr", tag, n, reps, ts[0], ts[reps / 2], d);
            flint_free(d);
        }
        free(ts);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(r, ctx);
    }
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(p);
}

static void powx_fmpz(const char *tag, const char *pdec) {
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_str(p, pdec, 10);
    fmpz_mod_ctx_t ctx;
    fmpz_mod_ctx_init(ctx, p);
    for (size_t k = 0; k < sizeof DEGS / sizeof *DEGS; k++) {
        long d = DEGS[k];
        fmpz_mod_poly_t g, ginv, r;
        fmpz_mod_poly_init(g, ctx);
        fmpz_mod_poly_init(ginv, ctx);
        fmpz_mod_poly_init(r, ctx);
        seq_fmpz(g, d, 2000003, ctx);
        fmpz_mod_poly_set_coeff_ui(g, d, 1, ctx);
        long reps = d <= 128 ? 9 : 3;
        double *ts = malloc(sizeof(double) * reps);
        for (long i = 0; i < reps; i++) {
            double t0 = now_us();
            fmpz_mod_poly_reverse(ginv, g, g->length, ctx);
            fmpz_mod_poly_inv_series(ginv, ginv, g->length, ctx);
            fmpz_mod_poly_powmod_x_fmpz_preinv(r, p, g, ginv, ctx);
            ts[i] = now_us() - t0;
        }
        qsort(ts, reps, sizeof(double), cmpd);
        char *dg;
        digest_fmpz(r, ctx, &dg);
        printf("{\"tool\":\"flint-%s\",\"op\":\"powx\",\"p\":\"%s\",\"n\":%ld,\"reps\":%ld,"
               "\"best_us\":%.2f,\"median_us\":%.2f,\"digest\":\"%s\"}\n",
               FLINT_VERSION, tag, d, reps, ts[0], ts[reps / 2], dg);
        flint_free(dg);
        free(ts);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_mod_poly_clear(ginv, ctx);
        fmpz_mod_poly_clear(r, ctx);
    }
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(p);
}

/* ---- word-size prime: nmod_poly ---- */
static void seq_nmod(nmod_poly_t a, long n, ulong s) {
    ulong p = a->mod.n, g = GOLD % p, c = n_powmod2(g, s, p);
    nmod_poly_zero(a);
    for (long i = 0; i < n; i++) {
        nmod_poly_set_coeff_ui(a, i, c);
        c = n_mulmod2_preinv(c, g, p, a->mod.ninv);
    }
}
static ulong digest_nmod(const nmod_poly_t a) { return nmod_poly_evaluate_nmod(a, 123456789 % a->mod.n); }

static void poly_nmod(const char *tag, ulong p) {
    for (size_t k = 0; k < sizeof LENS / sizeof *LENS; k++) {
        long n = LENS[k];
        nmod_poly_t a, b, r;
        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(r, p);
        seq_nmod(a, n, 16);
        seq_nmod(b, n, 1000003);
        long reps = reps_for(n);
        double *ts = malloc(sizeof(double) * reps);
        for (int op = 0; op < 2; op++) {
            for (long i = 0; i < reps; i++) {
                double t0 = now_us();
                if (op == 0) nmod_poly_mul(r, a, b);
                else nmod_poly_mul(r, a, a);
                ts[i] = now_us() - t0;
            }
            qsort(ts, reps, sizeof(double), cmpd);
            printf("{\"tool\":\"flint-%s\",\"op\":\"%s\",\"p\":\"%s\",\"n\":%ld,\"reps\":%ld,"
                   "\"best_us\":%.2f,\"median_us\":%.2f,\"digest\":\"%lu\"}\n",
                   FLINT_VERSION, op == 0 ? "mul" : "sqr", tag, n, reps, ts[0], ts[reps / 2],
                   digest_nmod(r));
        }
        free(ts);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(r);
    }
}

static void powx_nmod(const char *tag, ulong p) {
    for (size_t k = 0; k < sizeof DEGS / sizeof *DEGS; k++) {
        long d = DEGS[k];
        nmod_poly_t g, ginv, r;
        nmod_poly_init(g, p);
        nmod_poly_init(ginv, p);
        nmod_poly_init(r, p);
        seq_nmod(g, d, 2000003);
        nmod_poly_set_coeff_ui(g, d, 1);
        long reps = d <= 128 ? 9 : 3;
        double *ts = malloc(sizeof(double) * reps);
        for (long i = 0; i < reps; i++) {
            double t0 = now_us();
            nmod_poly_reverse(ginv, g, g->length);
            nmod_poly_inv_series(ginv, ginv, g->length);
            nmod_poly_powmod_x_ui_preinv(r, p, g, ginv);
            ts[i] = now_us() - t0;
        }
        qsort(ts, reps, sizeof(double), cmpd);
        printf("{\"tool\":\"flint-%s\",\"op\":\"powx\",\"p\":\"%s\",\"n\":%ld,\"reps\":%ld,"
               "\"best_us\":%.2f,\"median_us\":%.2f,\"digest\":\"%lu\"}\n",
               FLINT_VERSION, tag, d, reps, ts[0], ts[reps / 2], digest_nmod(r));
        free(ts);
        nmod_poly_clear(g);
        nmod_poly_clear(ginv);
        nmod_poly_clear(r);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <poly|powx> <p61|p127|p255|p511>\n", argv[0]);
        return 2;
    }
    const char *op = argv[1], *tag = argv[2];
    if (!strcmp(tag, "p61")) {
        ulong p = strtoull(P61, NULL, 10);
        if (!strcmp(op, "poly")) poly_nmod(tag, p);
        else powx_nmod(tag, p);
    } else {
        const char *pdec = !strcmp(tag, "p127") ? P127 : !strcmp(tag, "p255") ? P255 : P511;
        int max_n = !strcmp(tag, "p511") ? 8192 : 16384;
        if (!strcmp(op, "poly")) poly_fmpz(tag, pdec, max_n);
        else powx_fmpz(tag, pdec);
    }
    flint_cleanup();
    return 0;
}
