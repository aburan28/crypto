/* bench.cu -- driver for the elliptic-curve GPU kernels.
 *
 *   ./bench selftest        run every kernel against the host implementation
 *   ./bench field           microbenchmark mul / sqr / inv
 *   ./bench mul             batch scalar multiplication throughput
 *   ./bench rho [opts]      Pollard-rho walk throughput (and, on a small
 *                           curve, an actual DLP solve)
 *   ./bench bsgs [opts]     baby-step giant-step: build the table on the
 *                           device, then solve a planted interval log
 *
 * rho options:  --walks N --iters N --w W --rbits R --dp BITS --neg 0|1
 *               --variant reg|lowmem|ref  --solve
 * bsgs options: --walks N (threads) --iters N --w W --neg 0|1 --wbits B
 *               --variant reg|ref
 *
 * Everything the kernels compute is checked against the same host code that
 * test_cpu.cpp verifies against the Python reference, so `selftest` is a
 * true end-to-end check of the device path.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

#include "kernels.cuh"
#include "kernels_bsgs.cuh"
#include "rho_host.hpp"
#include "bsgs_host.hpp"

#define CU(call) do { \
    cudaError_t e_ = (call); \
    if (e_ != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s at %s:%d: %s\n", #call, __FILE__, __LINE__, \
                cudaGetErrorString(e_)); \
        exit(1); \
    } } while (0)

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

struct DevProps {
    cudaDeviceProp p;
    int sms, maxThreadsPerSM, smemPerSM, smemPerBlockOptin;
};

static DevProps g_dev;

static void print_device() {
    int dev = 0;
    CU(cudaGetDevice(&dev));
    CU(cudaGetDeviceProperties(&g_dev.p, dev));
    g_dev.sms = g_dev.p.multiProcessorCount;
    g_dev.maxThreadsPerSM = g_dev.p.maxThreadsPerMultiProcessor;
    g_dev.smemPerSM = (int)g_dev.p.sharedMemPerMultiprocessor;
    CU(cudaDeviceGetAttribute(&g_dev.smemPerBlockOptin,
                              cudaDevAttrMaxSharedMemoryPerBlockOptin, dev));
    printf("device: %s  sm_%d%d  %d SMs  %.1f GHz  %.1f GB  L2 %.1f MB\n",
           g_dev.p.name, g_dev.p.major, g_dev.p.minor, g_dev.sms,
           g_dev.p.clockRate / 1e6, g_dev.p.totalGlobalMem / 1e9,
           g_dev.p.l2CacheSize / 1e6);
    printf("        %d threads/SM, %d KB shared/SM, %d KB shared/block (opt-in)\n",
           g_dev.maxThreadsPerSM, g_dev.smemPerSM >> 10, g_dev.smemPerBlockOptin >> 10);
    if (g_dev.p.major >= 10)
        printf("        Blackwell class: see OPTIMIZATION_BLACKWELL.md\n");
}

/* ---------------------------------------------------------------- */
struct Timer {
    cudaEvent_t a, b;
    Timer() { CU(cudaEventCreate(&a)); CU(cudaEventCreate(&b)); }
    ~Timer() { cudaEventDestroy(a); cudaEventDestroy(b); }
    void start() { CU(cudaEventRecord(a)); }
    double stop() {
        CU(cudaEventRecord(b));
        CU(cudaEventSynchronize(b));
        float ms = 0;
        CU(cudaEventElapsedTime(&ms, a, b));
        return ms / 1e3;
    }
};

static void selftest_bsgs();

/* ---------------------------------------------------------------- */
static void selftest() {
    printf("[selftest] curve %s, %s reduction\n", CURVE_NAME, FP_FAST ? "special" : "Montgomery");
    const uint32_t n = 4096;

    /* --- scalar multiplication against the host --- */
    std::vector<affine_pt> h_in(n), h_out(n), h_ref(n);
    std::vector<uint32_t> h_k(n * 8);
    affine_pt G = Curve::generator();
    uint64_t seed = 0x1234;
    for (uint32_t i = 0; i < n; i++) {
        uint32_t s[8];
        rho_scalar_from_seed(seed, s);
        h_in[i] = Curve::to_affine(Curve::scalar_mul(G, s, 0));
        rho_scalar_from_seed(seed, s);
        for (int l = 0; l < 8; l++) h_k[i * 8 + l] = s[l];
    }
    for (uint32_t i = 0; i < n; i++)
        h_ref[i] = Curve::to_affine(Curve::scalar_mul(h_in[i], &h_k[i * 8], 0));

    affine_pt *d_in, *d_out;
    uint32_t *d_k;
    CU(cudaMalloc(&d_in, n * sizeof(affine_pt)));
    CU(cudaMalloc(&d_out, n * sizeof(affine_pt)));
    CU(cudaMalloc(&d_k, n * 8 * sizeof(uint32_t)));
    CU(cudaMemcpy(d_in, h_in.data(), n * sizeof(affine_pt), cudaMemcpyHostToDevice));
    CU(cudaMemcpy(d_k, h_k.data(), n * 8 * sizeof(uint32_t), cudaMemcpyHostToDevice));

    for (int ct = 0; ct < 2; ct++) {
        k_scalar_mul<<<(n + 127) / 128, 128>>>(d_out, d_in, d_k, n, ct);
        CU(cudaGetLastError());
        CU(cudaMemcpy(h_out.data(), d_out, n * sizeof(affine_pt), cudaMemcpyDeviceToHost));
        uint32_t bad = 0;
        for (uint32_t i = 0; i < n; i++) if (!Curve::affine_eq(h_out[i], h_ref[i])) bad++;
        CHECK(bad == 0, "scalar_mul(ct=%d): %u/%u mismatches", ct, bad, n);
        if (!bad) printf("  scalar_mul(ct=%d): %u points match the host\n", ct, n);
    }

    /* --- rho: kernel state must track the host stepper exactly --- */
    RhoHost h;
    h.prm.r_bits = 8; h.prm.dp_mask = 0x1F; h.prm.neg_map = 1;
    h.prm.max_steps = 1u << 20; h.prm.table_seed = 99;
    h.P = Curve::generator();
    uint32_t sk[8] = {0x9e3779b9u, 0x85ebca6bu, 0, 0, 0, 0, 0, 0};
    h.Q = Curve::to_affine(Curve::scalar_mul(h.P, sk, 0));
    h.build_table();

    const uint32_t T = 256, W = 8, nw = T * W, dpcap = 1 << 14;
    rho_ctx hc{}, dc{};
    std::vector<uint32_t> hX(8 * nw), hY(8 * nw), hH(2 * nw), hE(nw), hS(nw), hR(nw);
    std::vector<rho_dp> hdp(dpcap);
    uint32_t hcount = 0;
    hc.X = hX.data(); hc.Y = hY.data(); hc.H = hH.data(); hc.esc = hE.data();
    hc.steps = hS.data(); hc.restarts = hR.data();
    hc.nthreads = T; hc.walks_per_thread = W;
    hc.table = h.table.data(); hc.P = h.P; hc.Q = h.Q; hc.prm = h.prm;
    hc.dp_out = hdp.data(); hc.dp_count = &hcount; hc.dp_cap = dpcap;
    hc.cycle_counter = nullptr;

    dc = hc;
    CU(cudaMalloc(&dc.X, 8 * nw * sizeof(uint32_t)));
    CU(cudaMalloc(&dc.Y, 8 * nw * sizeof(uint32_t)));
    CU(cudaMalloc(&dc.H, 2 * nw * sizeof(uint32_t)));
    CU(cudaMalloc(&dc.esc, nw * sizeof(uint32_t)));
    CU(cudaMalloc(&dc.steps, nw * sizeof(uint32_t)));
    CU(cudaMalloc(&dc.restarts, nw * sizeof(uint32_t)));
    affine_pt *d_table;
    CU(cudaMalloc(&d_table, h.table.size() * sizeof(affine_pt)));
    CU(cudaMemcpy(d_table, h.table.data(), h.table.size() * sizeof(affine_pt),
                  cudaMemcpyHostToDevice));
    dc.table = d_table;
    rho_dp *d_dp;
    uint32_t *d_cnt;
    CU(cudaMalloc(&d_dp, dpcap * sizeof(rho_dp)));
    CU(cudaMalloc(&d_cnt, sizeof(uint32_t)));
    CU(cudaMemset(d_cnt, 0, sizeof(uint32_t)));
    dc.dp_out = d_dp; dc.dp_count = d_cnt;
    dc.cycle_counter = nullptr;

    size_t smem = rho_smem_bytes(h.prm.r_bits);
    k_rho_init<<<(T + 127) / 128, 128>>>(dc);
    CU(cudaGetLastError());
    for (uint32_t t = 0; t < T; t++) rho_init_thread(hc, t);

    const uint32_t iters = 64;
    k_rho_walk<W><<<(T + RHO_BLOCK - 1) / RHO_BLOCK, RHO_BLOCK, smem>>>(dc, iters);
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());
    for (uint32_t it = 0; it < iters; it++)
        for (uint32_t t = 0; t < T; t++) rho_step_batch<W>(hc, t);

    std::vector<uint32_t> gX(8 * nw), gY(8 * nw), gS(nw), gR(nw);
    CU(cudaMemcpy(gX.data(), dc.X, gX.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(gY.data(), dc.Y, gY.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(gS.data(), dc.steps, gS.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(gR.data(), dc.restarts, gR.size() * 4, cudaMemcpyDeviceToHost));
    CHECK(gX == hX && gY == hY, "rho walk state after %u iterations", iters);
    CHECK(gS == hS && gR == hR, "rho step/restart counters");
    uint32_t gcount = 0;
    CU(cudaMemcpy(&gcount, d_cnt, 4, cudaMemcpyDeviceToHost));
    CHECK(gcount == hcount, "dp count: device %u vs host %u", gcount, hcount);
    if (!failures) printf("  rho: %u walks x %u steps identical to the host, %u DPs\n",
                          nw, iters, gcount);

    /* every device DP must replay on the host to the reported x */
    std::vector<rho_dp> gdp(gcount < dpcap ? gcount : dpcap);
    if (!gdp.empty()) {
        CU(cudaMemcpy(gdp.data(), d_dp, gdp.size() * sizeof(rho_dp), cudaMemcpyDeviceToHost));
        uint32_t bad = 0;
        for (size_t i = 0; i < gdp.size() && i < 64; i++) {
            fp256 a, b; affine_pt e;
            if (!h.replay(gdp[i].walk, gdp[i].restart, gdp[i].steps, a, b, e) ||
                !mp_eq(e.x.v, gdp[i].x)) bad++;
        }
        CHECK(bad == 0, "%u device DPs failed host replay", bad);
        if (!bad) printf("  rho: device DPs replay correctly on the host\n");
    }

    /* the lowmem variant must produce the same state as the register one */
    {
        CU(cudaMemset(d_cnt, 0, sizeof(uint32_t)));
        k_rho_init<<<(T + 127) / 128, 128>>>(dc);
        k_rho_walk_lowmem<W><<<(T + RHO_BLOCK - 1) / RHO_BLOCK, RHO_BLOCK, smem>>>(dc, iters);
        CU(cudaGetLastError());
        CU(cudaDeviceSynchronize());
        std::vector<uint32_t> lX(8 * nw), lY(8 * nw);
        CU(cudaMemcpy(lX.data(), dc.X, lX.size() * 4, cudaMemcpyDeviceToHost));
        CU(cudaMemcpy(lY.data(), dc.Y, lY.size() * 4, cudaMemcpyDeviceToHost));
        CHECK(lX == hX && lY == hY, "lowmem variant state");
        if (lX == hX && lY == hY) printf("  rho lowmem variant matches\n");
    }

    cudaFree(d_in); cudaFree(d_out); cudaFree(d_k);
    cudaFree(dc.X); cudaFree(dc.Y); cudaFree(dc.H); cudaFree(dc.esc);
    cudaFree(dc.steps); cudaFree(dc.restarts);
    cudaFree(d_table); cudaFree(d_dp); cudaFree(d_cnt);
    selftest_bsgs();
    printf(failures ? "SELFTEST FAILED\n" : "selftest OK\n");
}

/* ---------------------------------------------------------------- */
static void bench_field() {
    const int blocks = g_dev.sms * 8, threads = 256;
    const uint32_t iters = 2000;
    fp256 *d_in, *d_out;
    CU(cudaMalloc(&d_in, 64 * sizeof(fp256)));
    CU(cudaMalloc(&d_out, (size_t)blocks * threads * sizeof(fp256)));
    std::vector<fp256> h_in(64);
    uint64_t seed = 7;
    for (int i = 0; i < 64; i++) {
        uint32_t s[8];
        rho_scalar_from_seed(seed, s);
        h_in[i] = Fp::from_limbs(s);
    }
    CU(cudaMemcpy(d_in, h_in.data(), 64 * sizeof(fp256), cudaMemcpyHostToDevice));
    double total = (double)blocks * threads;
    Timer tm;

    k_bench_mul<<<blocks, threads>>>(d_out, d_in, 10, 1);
    CU(cudaDeviceSynchronize());

    tm.start();
    k_bench_mul<<<blocks, threads>>>(d_out, d_in, iters, 1);
    double t_mul = tm.stop();
    printf("  mul (throughput): %.2f Gop/s\n", total * iters * 4 / t_mul / 1e9);

    tm.start();
    k_bench_mul<<<blocks, threads>>>(d_out, d_in, iters, 0);
    double t_dep = tm.stop();
    printf("  mul (dependent):  %.2f Gop/s\n", total * iters / t_dep / 1e9);

    tm.start();
    k_bench_sqr<<<blocks, threads>>>(d_out, d_in, iters);
    double t_sqr = tm.stop();
    printf("  sqr:              %.2f Gop/s\n", total * iters / t_sqr / 1e9);

    tm.start();
    k_bench_inv<<<blocks, threads>>>(d_out, d_in, iters / 100);
    double t_inv = tm.stop();
    printf("  inv:              %.2f Mop/s  (%.0f mul-equivalents)\n",
           total * (iters / 100) / t_inv / 1e6,
           (t_inv / (iters / 100)) / (t_dep / iters));
    cudaFree(d_in); cudaFree(d_out);
}

static void bench_mul() {
    const uint32_t n = 1u << 20;
    affine_pt *d_in, *d_out;
    uint32_t *d_k;
    CU(cudaMalloc(&d_in, n * sizeof(affine_pt)));
    CU(cudaMalloc(&d_out, n * sizeof(affine_pt)));
    CU(cudaMalloc(&d_k, n * 8 * sizeof(uint32_t)));
    std::vector<affine_pt> h_in(n);
    std::vector<uint32_t> h_k(n * 8);
    affine_pt G = Curve::generator();
    uint64_t seed = 3;
    affine_pt cur = G;
    for (uint32_t i = 0; i < n; i++) {
        h_in[i] = cur;
        cur = Curve::to_affine(Curve::madd(Curve::to_jac(cur), G));
        uint32_t s[8];
        rho_scalar_from_seed(seed, s);
        for (int l = 0; l < 8; l++) h_k[i * 8 + l] = s[l];
    }
    CU(cudaMemcpy(d_in, h_in.data(), n * sizeof(affine_pt), cudaMemcpyHostToDevice));
    CU(cudaMemcpy(d_k, h_k.data(), n * 8 * sizeof(uint32_t), cudaMemcpyHostToDevice));
    Timer tm;
    k_scalar_mul<<<(n + 127) / 128, 128>>>(d_out, d_in, d_k, n, 0);
    CU(cudaDeviceSynchronize());
    tm.start();
    k_scalar_mul<<<(n + 127) / 128, 128>>>(d_out, d_in, d_k, n, 0);
    double t = tm.stop();
    printf("  variable-base scalar mul: %.2f M/s (%u points in %.3f s)\n", n / t / 1e6, n, t);
    cudaFree(d_in); cudaFree(d_out); cudaFree(d_k);
}

/* ---------------------------------------------------------------- */
struct RhoOpts {
    uint32_t threads = 0, iters = 256, w = 8, rbits = 8, dpbits = 20, neg = 1;
    std::string variant = "lowmem";
    bool solve = false;
};

template <int W>
static void launch_rho(const RhoOpts &o, rho_ctx &dc, uint32_t blocks, size_t smem) {
    if (o.variant == "reg")
        k_rho_walk<W><<<blocks, RHO_BLOCK, smem>>>(dc, o.iters);
    else if (o.variant == "ref")
        k_rho_walk_ref<<<blocks, RHO_BLOCK, smem>>>(dc, o.iters);
    else
        k_rho_walk_lowmem<W><<<blocks, RHO_BLOCK, smem>>>(dc, o.iters);
}

static void bench_rho(RhoOpts o) {
    if (!o.threads) o.threads = (uint32_t)g_dev.sms * 2048 / o.w;
    o.threads = (o.threads + RHO_BLOCK - 1) / RHO_BLOCK * RHO_BLOCK;
    RhoHost h;
    h.prm.r_bits = o.rbits;
    h.prm.dp_mask = (o.dpbits >= 24) ? 0xFFFFFFu : ((1u << o.dpbits) - 1u);
    h.prm.neg_map = o.neg;
    h.prm.max_steps = 100u << 20;
    h.prm.table_seed = 2024;
    h.P = Curve::generator();
    /* A 48-bit secret.  On a toy curve --solve finds it; on secp256k1 the
     * search space is 2^128 and it never will -- the run still measures
     * throughput, which is the point of the benchmark. */
    uint32_t sk[8] = {0xdeadbeefu, 0x1234u, 0, 0, 0, 0, 0, 0};
    h.Q = Curve::to_affine(Curve::scalar_mul(h.P, sk, 0));
    if (o.solve && ModN::bits() > 96)
        printf("  note: --solve on a %d-bit group will not terminate; "
               "measuring throughput only\n", ModN::bits());
    h.build_table();

    uint32_t T = o.threads, W = o.w, nw = T * W;
    const uint32_t dpcap = 1u << 20;
    rho_ctx dc{};
    dc.nthreads = T; dc.walks_per_thread = W; dc.P = h.P; dc.Q = h.Q; dc.prm = h.prm;
    dc.dp_cap = dpcap; dc.cycle_counter = nullptr;
    CU(cudaMalloc(&dc.X, 8 * (size_t)nw * 4));
    CU(cudaMalloc(&dc.Y, 8 * (size_t)nw * 4));
    CU(cudaMalloc(&dc.H, 2 * (size_t)nw * 4));
    CU(cudaMalloc(&dc.esc, (size_t)nw * 4));
    CU(cudaMalloc(&dc.steps, (size_t)nw * 4));
    CU(cudaMalloc(&dc.restarts, (size_t)nw * 4));
    affine_pt *d_table;
    CU(cudaMalloc(&d_table, h.table.size() * sizeof(affine_pt)));
    CU(cudaMemcpy(d_table, h.table.data(), h.table.size() * sizeof(affine_pt),
                  cudaMemcpyHostToDevice));
    dc.table = d_table;
    CU(cudaMalloc(&dc.dp_out, (size_t)dpcap * sizeof(rho_dp)));
    CU(cudaMalloc(&dc.dp_count, 4));
    CU(cudaMemset(dc.dp_count, 0, 4));
    CU(cudaMalloc(&dc.cycle_counter, 8));
    CU(cudaMemset(dc.cycle_counter, 0, 8));

    size_t smem = rho_smem_bytes(o.rbits);
    uint32_t blocks = (T + RHO_BLOCK - 1) / RHO_BLOCK;
    printf("  %u walks (%u threads x %u), r=2^%u, dp=2^%u, neg=%u, variant=%s\n",
           nw, T, W, o.rbits, o.dpbits, o.neg, o.variant.c_str());
    printf("  %u blocks x %d threads, %zu B shared\n", blocks, RHO_BLOCK, smem);

    k_rho_init<<<(T + 127) / 128, 128>>>(dc);
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());

    Timer tm;
    /* warm up */
    if (W == 8) launch_rho<8>(o, dc, blocks, smem);
    else if (W == 16) launch_rho<16>(o, dc, blocks, smem);
    else if (W == 32) launch_rho<32>(o, dc, blocks, smem);
    else { fprintf(stderr, "unsupported --w %u (use 8, 16 or 32)\n", W); exit(1); }
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());

    tm.start();
    if (W == 8) launch_rho<8>(o, dc, blocks, smem);
    else if (W == 16) launch_rho<16>(o, dc, blocks, smem);
    else launch_rho<32>(o, dc, blocks, smem);
    double t = tm.stop();
    CU(cudaGetLastError());

    double steps = (double)nw * o.iters;
    unsigned long long cyc = 0;
    CU(cudaMemcpy(&cyc, dc.cycle_counter, 8, cudaMemcpyDeviceToHost));
    printf("  %.3f Gstep/s  (%.0f M steps in %.3f s)\n", steps / t / 1e9, steps / 1e6, t);
    printf("  %.1f ns/step/walk, %llu cycle escapes (%.2f%% of steps)\n",
           t / steps * 1e9, cyc, 100.0 * cyc / steps);

    if (o.solve) {
        RhoHost solver = h;
        std::vector<rho_dp> host_dp(dpcap);
        uint32_t k[8], consumed = 0;
        bool ok = false;
        double elapsed = t;
        for (int round = 0; round < 100000 && !ok; round++) {
            uint32_t cnt = 0;
            CU(cudaMemcpy(&cnt, dc.dp_count, 4, cudaMemcpyDeviceToHost));
            if (cnt > dpcap) cnt = dpcap;
            if (cnt > consumed) {
                CU(cudaMemcpy(host_dp.data() + consumed, dc.dp_out + consumed,
                              (cnt - consumed) * sizeof(rho_dp), cudaMemcpyDeviceToHost));
                while (consumed < cnt) if (solver.add_dp(host_dp[consumed++], k)) { ok = true; break; }
            }
            if (ok) break;
            tm.start();
            if (W == 8) launch_rho<8>(o, dc, blocks, smem);
            else if (W == 16) launch_rho<16>(o, dc, blocks, smem);
            else launch_rho<32>(o, dc, blocks, smem);
            elapsed += tm.stop();
            steps += (double)nw * o.iters;
        }
        if (ok) {
            printf("  SOLVED: k = 0x");
            for (int l = 7; l >= 0; l--) printf("%08x", k[l]);
            printf("\n  after %.0f M steps in %.1f s\n", steps / 1e6, elapsed);
        } else {
            printf("  not solved within the step budget\n");
        }
    }

    cudaFree(dc.X); cudaFree(dc.Y); cudaFree(dc.H); cudaFree(dc.esc);
    cudaFree(dc.steps); cudaFree(dc.restarts); cudaFree(d_table);
    cudaFree(dc.dp_out); cudaFree(dc.dp_count); cudaFree(dc.cycle_counter);
}


/* ---------------------------------------------------------------- *
 * Baby-step giant-step
 * ---------------------------------------------------------------- */
struct DevBsgs {
    bsgs_ctx c{};
    uint64_t *table = nullptr;
    uint32_t nchains = 0;

    void alloc_chains(uint32_t T, uint32_t W, uint32_t cand_cap) {
        nchains = T * W;
        CU(cudaMalloc(&c.X, 8 * (size_t)nchains * 4));
        CU(cudaMalloc(&c.Y, 8 * (size_t)nchains * 4));
        CU(cudaMalloc(&c.inf, (size_t)nchains * 4));
        CU(cudaMalloc(&c.pos, (size_t)nchains * 8));
        CU(cudaMalloc(&c.overflow, 4));
        CU(cudaMalloc(&c.cand, (size_t)cand_cap * sizeof(bsgs_cand)));
        CU(cudaMalloc(&c.cand_count, 4));
        c.cand_cap = cand_cap;
        reset_counters();
    }
    void alloc_table(uint32_t bits) {
        size_t bytes = (size_t)8 << bits;
        CU(cudaMalloc(&table, bytes));
        CU(cudaMemset(table, 0xFF, bytes));
        c.table = table;
        c.table_bits = bits;
    }
    void reset_counters() {
        CU(cudaMemset(c.overflow, 0, 4));
        CU(cudaMemset(c.cand_count, 0, 4));
    }
    uint32_t cand_count() {
        uint32_t n = 0;
        CU(cudaMemcpy(&n, c.cand_count, 4, cudaMemcpyDeviceToHost));
        return n;
    }
    uint32_t overflow() {
        uint32_t n = 0;
        CU(cudaMemcpy(&n, c.overflow, 4, cudaMemcpyDeviceToHost));
        return n;
    }
    unsigned long long table_count() {
        unsigned long long *d, h = 0;
        CU(cudaMalloc(&d, 8));
        CU(cudaMemset(d, 0, 8));
        k_bsgs_count<<<g_dev.sms * 4, 256>>>(table, 1ull << c.table_bits, d);
        CU(cudaGetLastError());
        CU(cudaMemcpy(&h, d, 8, cudaMemcpyDeviceToHost));
        cudaFree(d);
        return h;
    }
    void release() {
        cudaFree(c.X); cudaFree(c.Y); cudaFree(c.inf); cudaFree(c.pos);
        cudaFree(c.overflow); cudaFree(c.cand); cudaFree(c.cand_count);
        if (table) cudaFree(table);
    }
};

template <int W>
static void bsgs_launch_seed(DevBsgs &d, uint32_t T) {
    k_bsgs_seed<W><<<(T + BSGS_BLOCK - 1) / BSGS_BLOCK, BSGS_BLOCK>>>(d.c);
    CU(cudaGetLastError());
}

template <int W>
static void bsgs_launch_run(DevBsgs &d, uint32_t T, uint32_t iters, bool ref) {
    uint32_t blocks = (T + BSGS_BLOCK - 1) / BSGS_BLOCK;
    if (ref) k_bsgs_run_ref<<<blocks, BSGS_BLOCK>>>(d.c, iters);
    else k_bsgs_run<W><<<blocks, BSGS_BLOCK>>>(d.c, iters);
    CU(cudaGetLastError());
}

static void bsgs_seed_w(DevBsgs &d, uint32_t T, uint32_t W) {
    if (W == 4) bsgs_launch_seed<4>(d, T);
    else if (W == 8) bsgs_launch_seed<8>(d, T);
    else if (W == 16) bsgs_launch_seed<16>(d, T);
    else { fprintf(stderr, "unsupported --w %u (use 4, 8 or 16)\n", W); exit(1); }
}

static void bsgs_run_w(DevBsgs &d, uint32_t T, uint32_t W, uint32_t iters, bool ref) {
    if (W == 4) bsgs_launch_run<4>(d, T, iters, ref);
    else if (W == 8) bsgs_launch_run<8>(d, T, iters, ref);
    else bsgs_launch_run<16>(d, T, iters, ref);
}

/* Device state -> host vectors, for comparison with the CPU driver. */
static void bsgs_fetch(DevBsgs &d, std::vector<uint32_t> &X, std::vector<uint32_t> &Y,
                       std::vector<uint32_t> &inf, std::vector<uint64_t> &pos) {
    X.resize(8 * (size_t)d.nchains); Y.resize(8 * (size_t)d.nchains);
    inf.resize(d.nchains); pos.resize(d.nchains);
    CU(cudaMemcpy(X.data(), d.c.X, X.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(Y.data(), d.c.Y, Y.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(inf.data(), d.c.inf, inf.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(pos.data(), d.c.pos, pos.size() * 8, cudaMemcpyDeviceToHost));
}

static std::vector<uint64_t> bsgs_sorted_entries(const std::vector<uint64_t> &t) {
    std::vector<uint64_t> v;
    for (uint64_t e : t) if (e != BSGS_EMPTY) v.push_back(e);
    std::sort(v.begin(), v.end());
    return v;
}

static std::vector<bsgs_cand> bsgs_sorted_cands(std::vector<bsgs_cand> v) {
    std::sort(v.begin(), v.end(), [](const bsgs_cand &a, const bsgs_cand &b) {
        return a.i != b.i ? a.i < b.i : a.j < b.j;
    });
    return v;
}

/* The device table, chain state and candidates must match the host driver
 * that test_bsgs.cpp verifies.  The table is compared as a set: insertion
 * order across threads is not deterministic under atomicCAS, so the slot
 * an entry lands in may differ, but the entries may not. */
static void selftest_bsgs() {
    const uint32_t T = 256, W = 8, iters = 5;
    printf("[selftest bsgs] %u chains\n", T * W);
    affine_pt G = Curve::generator();
    uint32_t x0[8] = {0x12345678u, 0x9abcdef0u, 0x13579bdfu, 0, 0, 0, 0, 0};
    BsgsPlan p = bsgs_plan(1ull << 30, 1, T, W, x0);
    BsgsHost h;
    h.setup(G, p);

    /* host table */
    std::vector<uint64_t> htable = bsgs_new_table(p.table_bits);
    bsgs_ctx hc{};
    BsgsChains hb;
    hb.bind(hc, T, W, htable, 1024);
    BsgsStats hs;
    bsgs_cpu_build_table<W>(h, hc, iters, hs);

    /* device table */
    DevBsgs d;
    d.alloc_chains(T, W, 1024);
    d.alloc_table(p.table_bits);
    h.fill_baby_ctx(d.c);
    bsgs_seed_w(d, T, W);
    for (uint32_t r = 0; r < bsgs_rounds(p.Lb, iters); r++) bsgs_run_w(d, T, W, iters, false);
    CU(cudaDeviceSynchronize());
    std::vector<uint64_t> dtable((size_t)1 << p.table_bits);
    CU(cudaMemcpy(dtable.data(), d.table, dtable.size() * 8, cudaMemcpyDeviceToHost));
    CHECK(bsgs_sorted_entries(dtable) == bsgs_sorted_entries(htable), "bsgs baby table differs from the host");
    CHECK(d.overflow() == 0, "bsgs table overflow on the device");
    {
        std::vector<uint32_t> X, Y, inf; std::vector<uint64_t> pos;
        bsgs_fetch(d, X, Y, inf, pos);
        CHECK(X == hb.X && Y == hb.Y && inf == hb.inf && pos == hb.pos, "bsgs baby chain state");
    }
    if (!failures) printf("  bsgs: device baby table (%llu entries) and chain state match the host\n",
                          (unsigned long long)bsgs_sorted_entries(dtable).size());

    /* giant phase: one round on both, compare state and candidates; then
     * finish on the device and verify the planted secret. */
    uint32_t xs[8]; memcpy(xs, x0, sizeof(xs));
    uint64_t r = 0x2F1E3D4Cull;                     /* < 2^30 */
    uint64_t lo = (uint64_t)xs[0] + r; xs[0] = (uint32_t)lo; xs[1] += (uint32_t)(lo >> 32);
    affine_pt Q = Curve::to_affine(Curve::scalar_mul(G, xs, 0));
    h.set_target(Q);

    bsgs_ctx gc{};
    BsgsChains gb;
    gb.bind(gc, T, W, htable, 1024);
    h.fill_giant_ctx(gc);
    BsgsStats gs;
    bsgs_cpu_seed<W>(gc, gs);
    bsgs_cpu_round<W>(gc, iters, 0);

    d.reset_counters();
    h.fill_giant_ctx(d.c);
    bsgs_seed_w(d, T, W);
    bsgs_run_w(d, T, W, iters, false);
    CU(cudaDeviceSynchronize());
    {
        std::vector<uint32_t> X, Y, inf; std::vector<uint64_t> pos;
        bsgs_fetch(d, X, Y, inf, pos);
        CHECK(X == gb.X && Y == gb.Y && inf == gb.inf && pos == gb.pos, "bsgs giant chain state after one round");
        uint32_t dn = d.cand_count();
        std::vector<bsgs_cand> dc(dn < 1024 ? dn : 1024);
        if (!dc.empty()) CU(cudaMemcpy(dc.data(), d.c.cand, dc.size() * sizeof(bsgs_cand), cudaMemcpyDeviceToHost));
        std::vector<bsgs_cand> hcands(gb.cand.begin(), gb.cand.begin() + (gb.cand_count < 1024 ? gb.cand_count : 1024));
        auto a = bsgs_sorted_cands(dc), b = bsgs_sorted_cands(hcands);
        bool same = a.size() == b.size();
        for (size_t k = 0; same && k < a.size(); k++) same = a[k].i == b[k].i && a[k].j == b[k].j;
        CHECK(same, "bsgs candidates after one round: device %zu vs host %zu", a.size(), b.size());
    }
    /* the reference kernel must agree with the batched one */
    {
        d.reset_counters();
        bsgs_seed_w(d, T, W);
        bsgs_run_w(d, T, W, iters, true);
        CU(cudaDeviceSynchronize());
        std::vector<uint32_t> X, Y, inf; std::vector<uint64_t> pos;
        bsgs_fetch(d, X, Y, inf, pos);
        CHECK(X == gb.X && Y == gb.Y && inf == gb.inf && pos == gb.pos, "bsgs reference kernel state");
    }
    /* finish: run rounds until a candidate verifies */
    d.reset_counters();
    bsgs_seed_w(d, T, W);
    uint32_t k[8], consumed = 0;
    bool ok = false;
    for (uint32_t rd = 0; rd < bsgs_rounds(p.Lg, iters) && !ok; rd++) {
        bsgs_run_w(d, T, W, iters, false);
        CU(cudaDeviceSynchronize());
        uint32_t cnt = d.cand_count();
        if (cnt > 1024) cnt = 1024;
        if (cnt > consumed) {
            std::vector<bsgs_cand> dc(cnt - consumed);
            CU(cudaMemcpy(dc.data(), d.c.cand + consumed, dc.size() * sizeof(bsgs_cand), cudaMemcpyDeviceToHost));
            for (auto &cd : dc) { consumed++; if (h.verify(cd, k)) { ok = true; break; } }
        }
    }
    CHECK(ok, "bsgs: planted interval log not recovered on the device");
    if (ok) CHECK(Fn::eq(Fn::from_limbs(xs), Fn::from_limbs(k)), "bsgs: recovered k != planted secret");
    if (ok) printf("  bsgs: device giant phase recovered the planted secret\n");
    d.release();
}

struct BsgsOpts {
    uint32_t threads = 0, iters = 64, w = 8, neg = 1, wbits = 40;
    std::string variant = "reg";
};

static void bench_bsgs(BsgsOpts o) {
    if (!o.threads) o.threads = (uint32_t)g_dev.sms * 2048 / o.w;
    o.threads = (o.threads + BSGS_BLOCK - 1) / BSGS_BLOCK * BSGS_BLOCK;
    const uint32_t T = o.threads, W = o.w;
    const bool ref = (o.variant == "ref");
    if (o.wbits > 63) { fprintf(stderr, "--wbits must be <= 63\n"); exit(1); }
    uint64_t width = 1ull << o.wbits;
    if (ModN::bits() <= 63) {
        uint64_t n = 0;
        for (int l = 1; l >= 0; l--) n = (n << 32) | ModN::limb(l);
        if (width > n) { width = n; printf("  width clamped to the group order n = %llu\n", (unsigned long long)n); }
    }

    /* interval start and a planted secret inside it */
    uint64_t seed = 0xB5B5 + o.wbits;
    uint32_t x0[8];
    rho_scalar_from_seed(seed, x0);
    x0[7] &= 0x3FFFFFFFu;                 /* comfortably below n on a 256-bit curve */
    if (ModN::bits() <= 63) { x0[0] = 0; x0[1] = 0; for (int l = 2; l < 8; l++) x0[l] = 0; }
    uint64_t r = rho_splitmix64(seed) % width;
    uint32_t xs[8]; memcpy(xs, x0, sizeof(xs));
    {
        uint64_t carry = r;
        for (int l = 0; l < 8 && carry; l++) {
            uint64_t s = (uint64_t)xs[l] + (carry & 0xFFFFFFFFu);
            xs[l] = (uint32_t)s;
            carry = (carry >> 32) + (s >> 32);
        }
    }
    affine_pt G = Curve::generator();
    affine_pt Q = Curve::to_affine(Curve::scalar_mul(G, xs, 0));

    BsgsPlan p = bsgs_plan(width, o.neg, T, W, x0);
    BsgsHost h;
    h.setup(G, p);
    double sqrt_w = sqrt((double)width);
    printf("  interval 2^%u, neg=%u: m = %llu baby, stride M = %llu, %llu giants; table 2^%u slots "
           "(%.1f MB, load %.2f)\n", o.wbits, o.neg, (unsigned long long)p.m,
           (unsigned long long)p.M, (unsigned long long)p.giant_count, p.table_bits,
           (double)(8ull << p.table_bits) / 1e6, (double)p.m / (double)(1ull << p.table_bits));
    printf("  %u chains (%u threads x %u), Lb = %llu, Lg = %llu, %u steps per launch, variant=%s\n",
           T * W, T, W, (unsigned long long)p.Lb, (unsigned long long)p.Lg, o.iters, o.variant.c_str());
    if ((8ull << p.table_bits) > g_dev.p.totalGlobalMem) {
        printf("  table does not fit device memory; lower --wbits\n");
        return;
    }

    DevBsgs d;
    d.alloc_chains(T, W, 1u << 16);
    d.alloc_table(p.table_bits);
    Timer tm;
    BsgsStats st;

    /* baby phase */
    h.fill_baby_ctx(d.c);
    tm.start();
    bsgs_seed_w(d, T, W);
    double t_seed = tm.stop();
    uint32_t rounds_b = bsgs_rounds(p.Lb, o.iters);
    tm.start();
    for (uint32_t rd = 0; rd < rounds_b; rd++) bsgs_run_w(d, T, W, o.iters, ref);
    double t_baby = tm.stop();
    st.baby_steps = (unsigned long long)rounds_b * o.iters * T * W;
    bsgs_account_seed(d.c, st);
    unsigned long long entries = d.table_count();
    printf("  baby:  %llu entries in %.3f s (+%.3f s seed): %.3f Gstep/s, %.1f ns/step, overflow %u\n",
           entries, t_baby, t_seed, st.baby_steps / t_baby / 1e9, t_baby / st.baby_steps * 1e9,
           d.overflow());
    if (entries != p.m - 1) printf("  WARNING: expected %llu entries\n", (unsigned long long)(p.m - 1));

    /* giant phase, stopping at the first verified candidate */
    h.set_target(Q);
    d.reset_counters();
    h.fill_giant_ctx(d.c);
    tm.start();
    bsgs_seed_w(d, T, W);
    double t_gseed = tm.stop();
    bsgs_account_seed(d.c, st);
    uint32_t k[8], consumed = 0;
    bool ok = false;
    double t_giant = 0;
    uint32_t rounds_g = bsgs_rounds(p.Lg, o.iters);
    for (uint32_t rd = 0; rd < rounds_g && !ok; rd++) {
        tm.start();
        bsgs_run_w(d, T, W, o.iters, ref);
        t_giant += tm.stop();
        st.giant_steps += (unsigned long long)o.iters * T * W;
        st.rounds++;
        uint32_t cnt = d.cand_count();
        if (cnt > d.c.cand_cap) cnt = d.c.cand_cap;
        if (cnt > consumed) {
            std::vector<bsgs_cand> dc(cnt - consumed);
            CU(cudaMemcpy(dc.data(), d.c.cand + consumed, dc.size() * sizeof(bsgs_cand), cudaMemcpyDeviceToHost));
            for (auto &cd : dc) {
                consumed++; st.candidates++;
                if (h.verify(cd, k)) { ok = true; break; }
                st.false_candidates++;
            }
        }
    }
    printf("  giant: %llu steps in %llu launches, %.3f s (+%.3f s seed): %.3f Gstep/s, %.1f ns/step; "
           "%u candidates, %u false\n", st.giant_steps, st.rounds, t_giant, t_gseed,
           st.giant_steps / t_giant / 1e9, t_giant / st.giant_steps * 1e9, st.candidates, st.false_candidates);
    if (ok) {
        bool right = Fn::eq(Fn::from_limbs(xs), Fn::from_limbs(k));
        printf("  %s: k = 0x", right ? "SOLVED" : "WRONG ANSWER");
        for (int l = 7; l >= 0; l--) printf("%08x", k[l]);
        printf("\n");
        if (!right) failures++;
    } else {
        printf("  not solved: the giant phase ran to the end without a verified hit\n");
        failures++;
    }
    unsigned long long ops = st.baby_steps + st.giant_steps + st.seed_ops + st.candidates;
    printf("  S = %.3f  (ops / sqrt(width): baby %.3f + giant %.3f + seed %.4f), "
           "wall %.3f s\n", ops / sqrt_w, st.baby_steps / sqrt_w, st.giant_steps / sqrt_w,
           st.seed_ops / sqrt_w, t_seed + t_baby + t_gseed + t_giant);
    d.release();
}

int main(int argc, char **argv) {
    print_device();
    std::string cmd = argc > 1 ? argv[1] : "selftest";
    RhoOpts o;
    BsgsOpts b;
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { return (i + 1 < argc) ? (uint32_t)strtoul(argv[++i], nullptr, 0) : 0u; };
        if (a == "--walks") o.threads = b.threads = next();
        else if (a == "--iters") o.iters = b.iters = next();
        else if (a == "--w") o.w = b.w = next();
        else if (a == "--rbits") o.rbits = next();
        else if (a == "--dp") o.dpbits = next();
        else if (a == "--neg") o.neg = b.neg = next();
        else if (a == "--wbits") b.wbits = next();
        else if (a == "--variant" && i + 1 < argc) o.variant = b.variant = argv[++i];
        else if (a == "--solve") o.solve = true;
        else { fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    if (cmd == "selftest") selftest();
    else if (cmd == "field") bench_field();
    else if (cmd == "mul") bench_mul();
    else if (cmd == "rho") bench_rho(o);
    else if (cmd == "bsgs") bench_bsgs(b);
    else { fprintf(stderr, "usage: %s {selftest|field|mul|rho|bsgs}\n", argv[0]); return 2; }
    return failures ? 1 : 0;
}
