/* bench2.cu -- device driver for the pairs-and-solve sweep.
 *
 *   ./bench2 selftest      device sweep vs the verified host sweep
 *   ./bench2 throughput    pairs per second, swept over l
 *   ./bench2 occupancy     launch configuration report
 *
 * One thread owns one `X₁` and walks every `X₂ >= X₁`.  That is the
 * decomposition `RESEARCH_SEMAEV_DECOMPOSITION.md` describes as
 * "embarrassingly parallel over `X₁`, and nothing in it shares state":
 * no shared memory, no atomics except the single result slot, and no
 * communication between lanes at all.
 *
 * The load is triangular -- thread 0 walks `2^l` pairs and the last
 * thread walks one -- so a static one-row-per-thread mapping wastes
 * about half the lanes.  `sweep_kernel` therefore takes a grid-stride
 * over rows, which lets a smaller grid re-balance; `sweep_kernel_flat`
 * is the alternative that maps threads to *pairs* instead and is what
 * to compare against on real hardware.
 *
 * NOTE: this file has never been run.  The environment it was written
 * in had no NVIDIA device.  Treat every number it would print as
 * unmeasured until it has been.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "params.h"
#include "decomp.cuh"
#include "vectors.h"

#ifndef SEM_THREADS
#  define SEM_THREADS 256
#endif

#define CUDA_OK(call)                                                        \
    do {                                                                     \
        cudaError_t e_ = (call);                                             \
        if (e_ != cudaSuccess) {                                             \
            printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_),       \
                   __FILE__, __LINE__);                                      \
            exit(1);                                                         \
        }                                                                    \
    } while (0)

/* Result slot.  `found` is claimed with an atomic so the first writer
 * wins and the witness is never half-written by two lanes. */
struct SweepResult {
    unsigned int found;
    uint64_t w[3];
};

/* One thread per `X₁`, grid-stride so the triangular load can be
 * re-balanced by choosing a grid smaller than `2^l`. */
__global__ void sweep_kernel(uint64_t xr, int l, uint64_t irr_low, int n,
                             const uint64_t* lv, int lv_len, SweepResult* out) {
    const Gf2n f = {n, irr_low};
    const TargetPowers t = target_powers(xr, f);
    const uint64_t span = 1ull << l;
    for (uint64_t x1 = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x; x1 < span;
         x1 += (uint64_t)gridDim.x * blockDim.x) {
        if (out->found) return;          /* someone already answered */
        uint64_t w[3];
        if (decompose_row(x1, l, t, lv, lv_len, f, w)) {
            if (atomicCAS(&out->found, 0u, 1u) == 0u) {
                out->w[0] = w[0];
                out->w[1] = w[1];
                out->w[2] = w[2];
            }
            return;
        }
    }
}

/* One thread per *pair*, which removes the triangular imbalance at the
 * cost of recomputing the quartic setup per pair.  Which wins is a
 * hardware question; this is here to be measured against the above. */
__global__ void sweep_kernel_flat(uint64_t xr, int l, uint64_t irr_low, int n,
                                  const uint64_t* lv, int lv_len,
                                  SweepResult* out) {
    const Gf2n f = {n, irr_low};
    const TargetPowers t = target_powers(xr, f);
    const uint64_t span = 1ull << l;
    const uint64_t npairs = span * (span + 1) / 2;
    for (uint64_t idx = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
         idx < npairs; idx += (uint64_t)gridDim.x * blockDim.x) {
        if (out->found) return;
        /* Unrank idx into (x1, x2) with x1 <= x2.  Linear scan over
         * rows is fine: l is small and this is once per thread. */
        uint64_t x1 = 0, rem = idx;
        while (rem >= span - x1) { rem -= span - x1; x1++; }
        uint64_t x2 = x1 + rem;

        SemPoly q = quartic_with(x1, x2, t, f);
        int d = poly_deg(&q);
        uint64_t w[3] = {x1, x2, 0};
        int hit = 0;
        if (d < 0) {
            hit = 1;
        } else {
            uint64_t inv = gf_inv(q.c[d], f);
            for (int i = 0; i <= d; i++) q.c[i] = gf_mul(q.c[i], inv, f);
            SemPoly g = roots_in_subspace(&q, d, lv, lv_len, f);
            int dg = poly_deg(&g);
            if (dg == 1) {
                w[2] = gf_mul(g.c[0], gf_inv(g.c[1], f), f);
                hit = 1;
            } else if (dg > 1) {
                for (uint64_t tt = 0; tt < span && !hit; tt++) {
                    uint64_t v = 0;
                    for (int j = SEM_MAX_DEG; j >= 0; j--)
                        v = gf_mul(v, tt, f) ^ g.c[j];
                    if (v == 0) { w[2] = tt; hit = 1; }
                }
            }
        }
        if (hit && atomicCAS(&out->found, 0u, 1u) == 0u) {
            out->w[0] = w[0]; out->w[1] = w[1]; out->w[2] = w[2];
            return;
        }
    }
}

static const Gf2n HF = {SEM_N, SEM_IRR};

/* Host sweep, using the functions `make test` verified against sref.py.
 * Returns 1 and fills `w` if a decomposition exists. */
static int host_sweep(uint64_t xr, const uint64_t* lv, uint64_t w[3]) {
    TargetPowers t = target_powers(xr, HF);
    for (uint64_t x1 = 0; x1 < (1ull << SEM_L); x1++) {
        if (decompose_row(x1, SEM_L, t, lv, SEM_L + 1, HF, w)) return 1;
    }
    return 0;
}

static int witness_ok(uint64_t xr, const uint64_t w[3]) {
    TargetPowers t = target_powers(xr, HF);
    SemPoly q = quartic_with(w[0], w[1], t, HF);
    uint64_t v = 0;
    for (int j = SEM_MAX_DEG; j >= 0; j--) v = gf_mul(v, w[2], HF) ^ q.c[j];
    const uint64_t span = 1ull << SEM_L;
    return v == 0 && w[0] < span && w[1] < span && w[2] < span;
}

static int selftest() {
    printf("selftest: n=%d l=%d threads=%d\n", SEM_N, SEM_L, SEM_THREADS);
    std::vector<uint64_t> lv(SEM_L + 1);
    subspace_poly(SEM_L, lv.data(), HF);

    uint64_t* d_lv = nullptr;
    SweepResult* d_out = nullptr;
    CUDA_OK(cudaMalloc(&d_lv, lv.size() * sizeof(uint64_t)));
    CUDA_OK(cudaMalloc(&d_out, sizeof(SweepResult)));
    CUDA_OK(cudaMemcpy(d_lv, lv.data(), lv.size() * sizeof(uint64_t),
                       cudaMemcpyHostToDevice));

    int bad = 0, yes = 0, no = 0;
    for (int i = 0; i < TARGETS_N; i++) {
        uint64_t xr = TARGETS[i][0];
        int want = (int)TARGETS[i][1];   /* from exhaustive triples */
        uint64_t hw[3] = {0, 0, 0};
        int host = host_sweep(xr, lv.data(), hw);
        if (host != want) {
            printf("  target %llu: HOST sweep disagrees with the oracle\n",
                   (unsigned long long)xr);
            bad++;
        }

        for (int which = 0; which < 2; which++) {
            SweepResult init{0u, {0, 0, 0}};
            CUDA_OK(cudaMemcpy(d_out, &init, sizeof(init), cudaMemcpyHostToDevice));
            int grid = (int)(((1ull << SEM_L) + SEM_THREADS - 1) / SEM_THREADS);
            if (grid < 1) grid = 1;
            if (which == 0)
                sweep_kernel<<<grid, SEM_THREADS>>>(xr, SEM_L, SEM_IRR, SEM_N,
                                                    d_lv, SEM_L + 1, d_out);
            else
                sweep_kernel_flat<<<grid * 4, SEM_THREADS>>>(
                    xr, SEM_L, SEM_IRR, SEM_N, d_lv, SEM_L + 1, d_out);
            CUDA_OK(cudaGetLastError());
            CUDA_OK(cudaDeviceSynchronize());
            SweepResult got{};
            CUDA_OK(cudaMemcpy(&got, d_out, sizeof(got), cudaMemcpyDeviceToHost));
            const char* nm = which ? "flat" : "row ";
            if ((int)got.found != want) {
                printf("  target %llu: %s kernel says %u, exhaustive says %d\n",
                       (unsigned long long)xr, nm, got.found, want);
                bad++;
            } else if (got.found && !witness_ok(xr, got.w)) {
                printf("  target %llu: %s kernel witness does not satisfy f3\n",
                       (unsigned long long)xr, nm);
                bad++;
            }
        }
        if (want) yes++; else no++;
    }
    /* A run that only ever saw one verdict would pass by always saying
     * the same thing. */
    if (!yes || !no) {
        printf("  only one verdict occurred (%d yes, %d no): vacuous\n", yes, no);
        bad++;
    }
    cudaFree(d_lv);
    cudaFree(d_out);
    printf("%s (%d targets, %d yes, %d no, %d failures)\n",
           bad ? "SELFTEST FAILED" : "selftest ok", TARGETS_N, yes, no, bad);
    return bad ? 1 : 0;
}

static void throughput() {
    printf("throughput: n=%d l=%d threads=%d\n", SEM_N, SEM_L, SEM_THREADS);
    printf("  %6s %12s %16s %16s\n", "kernel", "ms", "pairs", "pairs/s");
    std::vector<uint64_t> lv(SEM_L + 1);
    subspace_poly(SEM_L, lv.data(), HF);
    uint64_t* d_lv = nullptr;
    SweepResult* d_out = nullptr;
    CUDA_OK(cudaMalloc(&d_lv, lv.size() * sizeof(uint64_t)));
    CUDA_OK(cudaMalloc(&d_out, sizeof(SweepResult)));
    CUDA_OK(cudaMemcpy(d_lv, lv.data(), lv.size() * sizeof(uint64_t),
                       cudaMemcpyHostToDevice));

    /* A target chosen NOT to decompose, so every pair is walked: the
     * note measures the same way, because rejection is ~85% of an
     * attack's work and a lucky early hit measures nothing. */
    uint64_t xr = 0;
    for (int i = 0; i < TARGETS_N; i++) {
        if (!TARGETS[i][1]) { xr = TARGETS[i][0]; break; }
    }
    if (!xr) { printf("  no undecomposable target in the vectors\n"); return; }

    const uint64_t span = 1ull << SEM_L;
    const uint64_t npairs = span * (span + 1) / 2;
    for (int which = 0; which < 2; which++) {
        SweepResult init{0u, {0, 0, 0}};
        int grid = (int)((span + SEM_THREADS - 1) / SEM_THREADS);
        if (grid < 1) grid = 1;
        cudaEvent_t t0, t1;
        CUDA_OK(cudaEventCreate(&t0));
        CUDA_OK(cudaEventCreate(&t1));
        for (int rep = 0; rep < 2; rep++) {
            /* Reset the result slot each rep.  Both kernels early-exit
             * once `found` is set, so a warm-up that answered would
             * leave the timed rep measuring an immediate return.  The
             * target above is chosen not to decompose, which already
             * prevents that, but the reset means the measurement does
             * not silently depend on that choice. */
            CUDA_OK(cudaMemcpy(d_out, &init, sizeof(init), cudaMemcpyHostToDevice));
            CUDA_OK(cudaDeviceSynchronize());
            if (rep == 1) CUDA_OK(cudaEventRecord(t0));
            if (which == 0)
                sweep_kernel<<<grid, SEM_THREADS>>>(xr, SEM_L, SEM_IRR, SEM_N,
                                                    d_lv, SEM_L + 1, d_out);
            else
                sweep_kernel_flat<<<grid * 4, SEM_THREADS>>>(
                    xr, SEM_L, SEM_IRR, SEM_N, d_lv, SEM_L + 1, d_out);
            if (rep == 1) CUDA_OK(cudaEventRecord(t1));
            CUDA_OK(cudaDeviceSynchronize());
        }
        float ms = 0;
        CUDA_OK(cudaEventElapsedTime(&ms, t0, t1));
        printf("  %6s %12.3f %16llu %16.0f\n", which ? "flat" : "row", ms,
               (unsigned long long)npairs, npairs / (ms / 1000.0));
        cudaEventDestroy(t0);
        cudaEventDestroy(t1);
    }
    cudaFree(d_lv);
    cudaFree(d_out);
    printf("\nThe CPU reference is about 1 us per pair at l = 8 (roughly 190\n");
    printf("field multiplications), measured in RESEARCH_SEMAEV_DECOMPOSITION.md.\n");
    printf("Quote the ratio, not the pairs/s: the note's own conclusion is that a\n");
    printf("faster oracle does not move the attack's exponent.\n");
}

static void occupancy() {
    int dev = 0;
    cudaDeviceProp prop{};
    CUDA_OK(cudaGetDevice(&dev));
    CUDA_OK(cudaGetDeviceProperties(&prop, dev));
    printf("device: %s, sm_%d%d, %d SMs\n", prop.name, prop.major, prop.minor,
           prop.multiProcessorCount);
    int blocks = 0;
    CUDA_OK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&blocks, sweep_kernel,
                                                          SEM_THREADS, 0));
    printf("blocks per SM at %d threads: %d (%d resident rows)\n", SEM_THREADS,
           blocks, blocks * prop.multiProcessorCount * SEM_THREADS);
    printf("rows in this sweep: %llu\n", (unsigned long long)(1ull << SEM_L));
    printf("no shared memory, no atomics except the single result slot.\n");
}

int main(int argc, char** argv) {
    const char* cmd = argc > 1 ? argv[1] : "selftest";
    if (!strcmp(cmd, "selftest")) return selftest();
    if (!strcmp(cmd, "throughput")) { throughput(); return 0; }
    if (!strcmp(cmd, "occupancy")) { occupancy(); return 0; }
    printf("usage: %s [selftest | throughput | occupancy]\n", argv[0]);
    return 2;
}
