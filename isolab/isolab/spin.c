/* isolab-spin: the A/A calibration kernel.
 *
 * A fixed amount of work with a serial dependency chain (xorshift64* mixed
 * with an integer multiply), so the instruction count is identical on every
 * run and the wall time measures the host, not the program.  Prints a
 * checksum, which must be identical between repeats, and writes metrics.json
 * into $ISOLAB_OUTPUT_DIR when that is set.
 *
 *   cc -O2 -static -o isolab-spin spin.c      (static so it runs in any image)
 *   isolab-spin [iterations]                  default 2e9, about 1-3 s per core
 */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char **argv) {
    uint64_t n = argc > 1 ? strtoull(argv[1], NULL, 10) : 2000000000ULL;
    uint64_t x = 0x9E3779B97F4A7C15ULL, acc = 0;
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (uint64_t i = 0; i < n; i++) {
        x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
        acc += x * 0x2545F4914F6CDD1DULL;
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double secs = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    printf("iterations=%" PRIu64 " checksum=%016" PRIx64 " seconds=%.6f\n", n, acc, secs);
    const char *out = getenv("ISOLAB_OUTPUT_DIR");
    if (out && *out) {
        char path[4096];
        snprintf(path, sizeof path, "%s/metrics.json", out);
        FILE *f = fopen(path, "w");
        if (f) {
            fprintf(f, "{\"iterations\": %" PRIu64 ", \"checksum\": \"%016" PRIx64 "\", \"seconds\": %.6f}\n",
                    n, acc, secs);
            fclose(f);
        }
    }
    return 0;
}
