// Host calibration and simulator controls for RESEARCH_CACHEGRIND_N41_20260929.md.
// Single-threaded, 4 KiB pages (no madvise), one cache line per node.
//
//   mem_calib freq                          effective core clock from a dependent add chain
//   mem_calib chase W k steps reps warm     k independent pointer chains over W bytes total;
//                                           one JSON line, ns per group step (k loads)
//   mem_calib seq W passes reps             sequential 8-byte reads over W bytes; ns per line
//   mem_calib branch                        cycles per mispredicted branch (random vs constant)
//   mem_calib antagonist W k seconds        run k chains over W bytes for a wall-clock duration
//
// Build: gcc -O2 -o mem_calib mem_calib.c
#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct node { struct node *next; uint64_t pad[7]; } node;

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;
static uint64_t rng(void) {
    uint64_t x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    return rng_state = x;
}

static int cmp_d(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

// k disjoint random cycles (Sattolo) of n/k nodes each; heads[c] starts chain c.
static node *build(size_t bytes, int k, node **heads) {
    size_t n = bytes / sizeof(node);
    size_t per = n / k;
    node *nodes = aligned_alloc(64, n * sizeof(node));
    uint32_t *perm = malloc(per * sizeof(uint32_t));
    if (!nodes || !perm) { fprintf(stderr, "alloc failed\n"); exit(2); }
    memset(nodes, 0, n * sizeof(node));
    for (int c = 0; c < k; c++) {
        node *base = nodes + (size_t)c * per;
        for (size_t i = 0; i < per; i++) perm[i] = (uint32_t)i;
        for (size_t i = per - 1; i > 0; i--) {
            size_t j = rng() % i; // Sattolo: j < i gives a single n-cycle
            uint32_t t = perm[i]; perm[i] = perm[j]; perm[j] = t;
        }
        for (size_t i = 0; i < per; i++) base[perm[i]].next = &base[perm[(i + 1) % per]];
        heads[c] = &base[perm[0]];
    }
    free(perm);
    return nodes;
}

#define CHASE_CASE(K)                                                        \
    case K: {                                                                \
        node *p[K];                                                          \
        for (int c = 0; c < K; c++) p[c] = heads[c];                         \
        for (uint64_t s = 0; s < steps; s++) {                               \
            _Pragma("GCC unroll 16")                                         \
            for (int c = 0; c < K; c++) p[c] = p[c]->next;                   \
        }                                                                    \
        for (int c = 0; c < K; c++) { heads[c] = p[c]; sink ^= (uintptr_t)p[c]; } \
        break;                                                               \
    }

static uintptr_t sink;
static void run_chains(node **heads, int k, uint64_t steps) {
    switch (k) {
        CHASE_CASE(1) CHASE_CASE(2) CHASE_CASE(4) CHASE_CASE(6) CHASE_CASE(8)
        CHASE_CASE(10) CHASE_CASE(12) CHASE_CASE(16)
        default: fprintf(stderr, "unsupported k=%d\n", k); exit(2);
    }
}

static int cmd_chase(size_t bytes, int k, uint64_t steps, int reps, uint64_t warm) {
    node *heads[16];
    node *nodes = build(bytes, k, heads);
    if (warm) run_chains(heads, k, warm);
    double t[64];
    for (int r = 0; r < reps; r++) {
        double t0 = now();
        run_chains(heads, k, steps);
        t[r] = (now() - t0) * 1e9 / (double)steps;
    }
    qsort(t, reps, sizeof(double), cmp_d);
    printf("{\"mode\":\"chase\",\"bytes\":%zu,\"k\":%d,\"steps\":%llu,\"reps\":%d,\"warm\":%llu,"
           "\"ns_per_step_min\":%.4f,\"ns_per_step_median\":%.4f,\"ns_per_load_median\":%.4f,\"sink\":%llu}\n",
           bytes, k, (unsigned long long)steps, reps, (unsigned long long)warm, t[0], t[reps / 2],
           t[reps / 2] / k, (unsigned long long)sink);
    free(nodes);
    return 0;
}

static int cmd_seq(size_t bytes, int passes, int reps) {
    size_t words = bytes / 8;
    uint64_t *a = aligned_alloc(64, bytes);
    for (size_t i = 0; i < words; i++) a[i] = i * 0x9E3779B97F4A7C15ULL;
    double t[64];
    uint64_t acc = 0;
    for (int r = 0; r < reps; r++) {
        double t0 = now();
        for (int p = 0; p < passes; p++)
            for (size_t i = 0; i < words; i += 8) {
                uint64_t s = a[i] + a[i + 1] + a[i + 2] + a[i + 3] + a[i + 4] + a[i + 5] + a[i + 6] + a[i + 7];
                acc += s;
            }
        t[r] = (now() - t0) * 1e9 / ((double)passes * (double)(words / 8));
    }
    qsort(t, reps, sizeof(double), cmp_d);
    printf("{\"mode\":\"seq\",\"bytes\":%zu,\"passes\":%d,\"reps\":%d,\"ns_per_line_min\":%.4f,"
           "\"ns_per_line_median\":%.4f,\"sink\":%llu}\n",
           bytes, passes, reps, t[0], t[reps / 2], (unsigned long long)acc);
    free(a);
    return 0;
}

static int cmd_freq(void) {
    double best = 1e30;
    uint64_t x = 0;
    const uint64_t iters = 200000000ULL; // x8 adds per iteration
    for (int r = 0; r < 7; r++) {
        double t0 = now();
        for (uint64_t i = 0; i < iters; i++)
            __asm__ volatile("add $1,%0\n\tadd $1,%0\n\tadd $1,%0\n\tadd $1,%0\n\t"
                             "add $1,%0\n\tadd $1,%0\n\tadd $1,%0\n\tadd $1,%0" : "+r"(x));
        double dt = now() - t0;
        if (dt < best) best = dt;
    }
    printf("{\"mode\":\"freq\",\"f_hz\":%.0f,\"sink\":%llu}\n", (double)iters * 8.0 / best,
           (unsigned long long)x);
    return 0;
}

// One branch per element: random bytes (50% taken) vs all-ones (always taken).
static double branch_pass(const uint8_t *v, size_t n, int passes) {
    uint64_t s = 1;
    double t0 = now();
    for (int p = 0; p < passes; p++)
        for (size_t i = 0; i < n; i++) {
            uint64_t b = v[i];
            __asm__ volatile("test $1,%1\n\tjz 1f\n\tadd $3,%0\n\tjmp 2f\n1:\tsub $1,%0\n2:"
                             : "+r"(s) : "r"(b) : "cc");
        }
    double dt = now() - t0;
    sink ^= s;
    return dt;
}

static int cmd_branch(void) {
    const size_t n = 1 << 16; // 64 KiB of bytes; L2-resident
    const int passes = 6000;
    uint8_t *rnd = malloc(n), *cst = malloc(n);
    for (size_t i = 0; i < n; i++) { rnd[i] = rng() & 1; cst[i] = 1; }
    double tr[7], tc[7];
    for (int r = 0; r < 7; r++) { tc[r] = branch_pass(cst, n, passes); tr[r] = branch_pass(rnd, n, passes); }
    qsort(tr, 7, sizeof(double), cmp_d); qsort(tc, 7, sizeof(double), cmp_d);
    double iters = (double)n * passes;
    // 50% of the random branches are mispredicted by an ideal-random predictor
    printf("{\"mode\":\"branch\",\"iters\":%.0f,\"ns_per_iter_random\":%.4f,\"ns_per_iter_const\":%.4f,"
           "\"ns_per_mispredict\":%.4f,\"sink\":%llu}\n",
           iters, tr[3] * 1e9 / iters, tc[3] * 1e9 / iters,
           (tr[3] - tc[3]) * 1e9 / (0.5 * iters), (unsigned long long)sink);
    return 0;
}

static int cmd_antagonist(size_t bytes, int k, double seconds) {
    node *heads[16];
    node *nodes = build(bytes, k, heads);
    double t0 = now();
    while (now() - t0 < seconds) run_chains(heads, k, 1000000);
    printf("{\"mode\":\"antagonist\",\"bytes\":%zu,\"k\":%d,\"seconds\":%.1f,\"sink\":%llu}\n", bytes, k,
           seconds, (unsigned long long)sink);
    free(nodes);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "freq")) return cmd_freq();
    if (argc >= 7 && !strcmp(argv[1], "chase"))
        return cmd_chase(strtoull(argv[2], 0, 10), atoi(argv[3]), strtoull(argv[4], 0, 10),
                         atoi(argv[5]), strtoull(argv[6], 0, 10));
    if (argc >= 5 && !strcmp(argv[1], "seq"))
        return cmd_seq(strtoull(argv[2], 0, 10), atoi(argv[3]), atoi(argv[4]));
    if (argc >= 2 && !strcmp(argv[1], "branch")) return cmd_branch();
    if (argc >= 5 && !strcmp(argv[1], "antagonist"))
        return cmd_antagonist(strtoull(argv[2], 0, 10), atoi(argv[3]), atof(argv[4]));
    fprintf(stderr, "usage: see header\n");
    return 2;
}
