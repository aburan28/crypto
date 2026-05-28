/* Profile GMP modular inverse cost at 81-bit */
#include <stdio.h>
#include <time.h>
#include <gmp.h>

int main() {
    mpz_t p, a, r;
    mpz_init(p); mpz_init(a); mpz_init(r);
    mpz_set_str(p, "1208925819614629469615699", 10);
    mpz_set_ui(a, 12345);

    int N = 1000000;
    struct timespec ts0, ts1;

    /* Warmup */
    for (int i = 0; i < 1000; i++) mpz_invert(r, a, p);

    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (int i = 0; i < N; i++) {
        mpz_invert(r, a, p);
        /* Slightly modify a to prevent over-caching */
        if (i % 100 == 0) mpz_add_ui(a, a, 1);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double dt = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    printf("GMP mpz_invert: %d calls in %.3f s = %.1f ns/call\n",
           N, dt, dt * 1e9 / N);

    /* mul comparison */
    mpz_t b;
    mpz_init(b);
    mpz_set_ui(b, 999999);
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (int i = 0; i < N; i++) {
        mpz_mul(r, a, b);
        mpz_mod(r, r, p);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    dt = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;
    printf("GMP mpz_mul+mod: %d calls in %.3f s = %.1f ns/call\n",
           N, dt, dt * 1e9 / N);

    mpz_clear(p); mpz_clear(a); mpz_clear(b); mpz_clear(r);
    return 0;
}
