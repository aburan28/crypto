/* bench2k.cu -- driver for the Koblitz-curve GPU kernels.
 *
 *   ./bench2k selftest      every kernel against the host implementation
 *   ./bench2k field         microbenchmark mul / sqr / inv / class weight
 *   ./bench2k rho [opts]    Frobenius-class walk throughput
 *
 * rho options: --walks N --iters N --w W --dp T --variant reg|lowmem|ref
 *
 * The host code it checks against is the same code test_cpu2k.cpp verifies
 * against the Python oracle, so `selftest` closes the loop on the device.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

#include "kernels2k.cuh"
#include "rho2k_host.hpp"

#define CU(call) do { \
    cudaError_t e_ = (call); \
    if (e_ != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s at %s:%d: %s\n", #call, __FILE__, __LINE__, \
                cudaGetErrorString(e_)); \
        exit(1); \
    } } while (0)

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

static int g_sms = 1;

static const uint32_t *cb_flat() { return &f2m_cb_table[0][0][0]; }
static const size_t cb_words = (size_t)F2M_CB_WINDOWS * 16 * F2M_CB_STRIDE;

static void print_device() {
    int dev = 0;
    cudaDeviceProp p;
    CU(cudaGetDevice(&dev));
    CU(cudaGetDeviceProperties(&p, dev));
    g_sms = p.multiProcessorCount;
    int smem_optin = 0;
    CU(cudaDeviceGetAttribute(&smem_optin, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev));
    printf("device: %s  sm_%d%d  %d SMs  %.1f GHz\n", p.name, p.major, p.minor,
           g_sms, p.clockRate / 1e6);
    printf("        %d threads/SM, %d KB shared/SM, %d KB shared/block (opt-in)\n",
           p.maxThreadsPerMultiProcessor, (int)(p.sharedMemPerMultiprocessor >> 10),
           smem_optin >> 10);
    printf("curve:  %s, F_2^%d, r = %d bits, class size 2m = %d "
           "(rho speedup %.1fx)\n",
           CURVE2K_NAME, F2M_M, CURVE2K_R_BITS, 2 * F2M_M, sqrt(2.0 * F2M_M));
    printf("        class-weight table %zu B shared\n", r2k_smem_bytes());
}

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

/* ---------------------------------------------------------------- */
static void selftest() {
    printf("[selftest] %s\n", CURVE2K_NAME);
    const uint32_t n = 2048;

    std::vector<pt2k> h_in(n), h_out(n), h_ref(n), h_tau(n);
    std::vector<uint32_t> h_k(n * SC_WORDS);
    pt2k G = Koblitz::generator();
    uint64_t seed = 0xABCD;
    for (uint32_t i = 0; i < n; i++) {
        uint32_t s[SC_WORDS];
        r2k_scalar_from_seed(seed, s);
        h_in[i] = Koblitz::mul(G, s);
        r2k_scalar_from_seed(seed, s);
        for (int l = 0; l < SC_WORDS; l++) h_k[i * SC_WORDS + l] = s[l];
    }
    for (uint32_t i = 0; i < n; i++) {
        h_ref[i] = Koblitz::mul(h_in[i], &h_k[i * SC_WORDS]);
        h_tau[i] = Koblitz::frob(h_in[i], 5);
    }

    pt2k *d_in, *d_out;
    uint32_t *d_k;
    CU(cudaMalloc(&d_in, n * sizeof(pt2k)));
    CU(cudaMalloc(&d_out, n * sizeof(pt2k)));
    CU(cudaMalloc(&d_k, n * SC_WORDS * sizeof(uint32_t)));
    CU(cudaMemcpy(d_in, h_in.data(), n * sizeof(pt2k), cudaMemcpyHostToDevice));
    CU(cudaMemcpy(d_k, h_k.data(), n * SC_WORDS * 4, cudaMemcpyHostToDevice));

    k2k_scalar_mul<<<(n + 127) / 128, 128>>>(d_out, d_in, d_k, n);
    CU(cudaGetLastError());
    CU(cudaMemcpy(h_out.data(), d_out, n * sizeof(pt2k), cudaMemcpyDeviceToHost));
    uint32_t bad = 0;
    for (uint32_t i = 0; i < n; i++) if (!Koblitz::eq(h_out[i], h_ref[i])) bad++;
    CHECK(bad == 0, "scalar_mul: %u/%u mismatches", bad, n);
    if (!bad) printf("  scalar_mul: %u points match the host\n", n);

    k2k_frob<<<(n + 127) / 128, 128>>>(d_out, d_in, 5, n);
    CU(cudaGetLastError());
    CU(cudaMemcpy(h_out.data(), d_out, n * sizeof(pt2k), cudaMemcpyDeviceToHost));
    bad = 0;
    for (uint32_t i = 0; i < n; i++) if (!Koblitz::eq(h_out[i], h_tau[i])) bad++;
    CHECK(bad == 0, "tau^5: %u/%u mismatches", bad, n);
    if (!bad) printf("  tau^5: %u points match the host\n", n);

    /* canonical class representatives */
    {
        f2e *d_c;
        CU(cudaMalloc(&d_c, n * sizeof(f2e)));
        k2k_canonical<<<(n + 127) / 128, 128>>>(d_c, d_in, n);
        CU(cudaGetLastError());
        std::vector<f2e> h_c(n);
        CU(cudaMemcpy(h_c.data(), d_c, n * sizeof(f2e), cudaMemcpyDeviceToHost));
        bad = 0;
        for (uint32_t i = 0; i < n; i++)
            if (!F2::eq(h_c[i], r2k_canonical_x(h_in[i]))) bad++;
        CHECK(bad == 0, "canonical: %u/%u mismatches", bad, n);
        if (!bad) printf("  canonical class representative: %u match\n", n);
        cudaFree(d_c);
    }

    /* the walk must track the host stepper exactly */
    Rho2kHost h;
    h.prm.nj = 8; h.prm.jmin = 3;
    h.prm.dp_threshold = (uint32_t)(F2M_M / 2 - 3);
    h.prm.max_steps = 1u << 20;
    h.cb = cb_flat();
    h.P = G;
    uint32_t sk[SC_WORDS] = {0x9e3779b9u, 0x85ebca6bu, 0, 0};
    h.Q = Koblitz::mul(h.P, sk);
    h.build();

    const uint32_t T = 256, W = 8, nw = T * W, dpcap = 1 << 14;
    rho2k_ctx hc{}, dc{};
    std::vector<uint32_t> hX(F2M_WORDS * nw), hY(F2M_WORDS * nw), hS(nw), hR(nw);
    std::vector<rho2k_dp> hdp(dpcap);
    uint32_t hcount = 0;
    hc.X = hX.data(); hc.Y = hY.data(); hc.steps = hS.data(); hc.restarts = hR.data();
    hc.nthreads = T; hc.walks_per_thread = W;
    hc.cb = cb_flat(); hc.P = h.P; hc.Q = h.Q; hc.prm = h.prm;
    hc.dp_out = hdp.data(); hc.dp_count = &hcount; hc.dp_cap = dpcap;

    dc = hc;
    CU(cudaMalloc(&dc.X, F2M_WORDS * (size_t)nw * 4));
    CU(cudaMalloc(&dc.Y, F2M_WORDS * (size_t)nw * 4));
    CU(cudaMalloc(&dc.steps, (size_t)nw * 4));
    CU(cudaMalloc(&dc.restarts, (size_t)nw * 4));
    uint32_t *d_cb;
    CU(cudaMalloc(&d_cb, cb_words * 4));
    CU(cudaMemcpy(d_cb, cb_flat(), cb_words * 4, cudaMemcpyHostToDevice));
    dc.cb = d_cb;
    rho2k_dp *d_dp;
    uint32_t *d_cnt;
    CU(cudaMalloc(&d_dp, dpcap * sizeof(rho2k_dp)));
    CU(cudaMalloc(&d_cnt, 4));
    CU(cudaMemset(d_cnt, 0, 4));
    dc.dp_out = d_dp; dc.dp_count = d_cnt;

    size_t smem = r2k_smem_bytes();
    k2k_rho_init<<<(T + 127) / 128, 128>>>(dc);
    CU(cudaGetLastError());
    for (uint32_t t = 0; t < T; t++) r2k_init_thread(hc, t);

    const uint32_t iters = 64;
    k2k_rho_walk<W><<<(T + R2K_BLOCK - 1) / R2K_BLOCK, R2K_BLOCK, smem>>>(dc, iters);
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());
    for (uint32_t it = 0; it < iters; it++)
        for (uint32_t t = 0; t < T; t++) r2k_step_batch<W>(hc, t);

    std::vector<uint32_t> gX(F2M_WORDS * nw), gY(F2M_WORDS * nw), gS(nw), gR(nw);
    CU(cudaMemcpy(gX.data(), dc.X, gX.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(gY.data(), dc.Y, gY.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(gS.data(), dc.steps, gS.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(gR.data(), dc.restarts, gR.size() * 4, cudaMemcpyDeviceToHost));
    CHECK(gX == hX && gY == hY, "rho walk state after %u iterations", iters);
    CHECK(gS == hS && gR == hR, "rho step/restart counters");
    uint32_t gcount = 0;
    CU(cudaMemcpy(&gcount, d_cnt, 4, cudaMemcpyDeviceToHost));
    CHECK(gcount == hcount, "dp count: device %u vs host %u", gcount, hcount);
    if (!failures)
        printf("  rho: %u class walks x %u steps identical to the host, %u DPs\n",
               nw, iters, gcount);

    /* device DPs must replay on the host to the same class */
    std::vector<rho2k_dp> gdp(gcount < dpcap ? gcount : dpcap);
    if (!gdp.empty()) {
        CU(cudaMemcpy(gdp.data(), d_dp, gdp.size() * sizeof(rho2k_dp),
                      cudaMemcpyDeviceToHost));
        uint32_t badr = 0;
        for (size_t i = 0; i < gdp.size() && i < 32; i++) {
            sc_t a, b;
            pt2k end;
            if (!h.replay(gdp[i].walk, gdp[i].restart, gdp[i].steps, a, b, end)) { badr++; continue; }
            f2e c = r2k_canonical_x(end);
            if (memcmp(c.v, gdp[i].x, sizeof(c.v)) != 0) badr++;
        }
        CHECK(badr == 0, "%u device DPs failed host replay", badr);
        if (!badr) printf("  rho: device DPs replay to the same class on the host\n");
    }

    /* lowmem variant must match */
    {
        CU(cudaMemset(d_cnt, 0, 4));
        k2k_rho_init<<<(T + 127) / 128, 128>>>(dc);
        k2k_rho_walk_lowmem<W><<<(T + R2K_BLOCK - 1) / R2K_BLOCK, R2K_BLOCK, smem>>>(dc, iters);
        CU(cudaGetLastError());
        CU(cudaDeviceSynchronize());
        std::vector<uint32_t> lX(F2M_WORDS * nw), lY(F2M_WORDS * nw);
        CU(cudaMemcpy(lX.data(), dc.X, lX.size() * 4, cudaMemcpyDeviceToHost));
        CU(cudaMemcpy(lY.data(), dc.Y, lY.size() * 4, cudaMemcpyDeviceToHost));
        CHECK(lX == hX && lY == hY, "lowmem variant state");
        if (lX == hX && lY == hY) printf("  rho lowmem variant matches\n");
    }

    cudaFree(d_in); cudaFree(d_out); cudaFree(d_k);
    cudaFree(dc.X); cudaFree(dc.Y); cudaFree(dc.steps); cudaFree(dc.restarts);
    cudaFree(d_cb); cudaFree(d_dp); cudaFree(d_cnt);
    printf(failures ? "SELFTEST FAILED\n" : "selftest OK\n");
}

/* ---------------------------------------------------------------- */
static void bench_field() {
    const int blocks = g_sms * 8, threads = 256;
    const uint32_t iters = 2000;
    f2e *d_in, *d_out;
    uint32_t *d_cb, *d_u32;
    CU(cudaMalloc(&d_in, 64 * sizeof(f2e)));
    CU(cudaMalloc(&d_out, (size_t)blocks * threads * sizeof(f2e)));
    CU(cudaMalloc(&d_u32, (size_t)blocks * threads * 4));
    CU(cudaMalloc(&d_cb, cb_words * 4));
    CU(cudaMemcpy(d_cb, cb_flat(), cb_words * 4, cudaMemcpyHostToDevice));
    std::vector<f2e> h_in(64);
    uint64_t seed = 11;
    for (int i = 0; i < 64; i++) {
        uint32_t s[SC_WORDS];
        r2k_scalar_from_seed(seed, s);
        h_in[i] = F2::from_limbs(s);
    }
    CU(cudaMemcpy(d_in, h_in.data(), 64 * sizeof(f2e), cudaMemcpyHostToDevice));
    double total = (double)blocks * threads;
    Timer tm;

    k2k_bench_mul<<<blocks, threads>>>(d_out, d_in, 10, 1);
    CU(cudaDeviceSynchronize());

    tm.start();
    k2k_bench_mul<<<blocks, threads>>>(d_out, d_in, iters, 1);
    double t_mul = tm.stop();
    printf("  mul (throughput): %.2f Gop/s\n", total * iters * 4 / t_mul / 1e9);

    tm.start();
    k2k_bench_mul<<<blocks, threads>>>(d_out, d_in, iters, 0);
    double t_dep = tm.stop();
    printf("  mul (dependent):  %.2f Gop/s\n", total * iters / t_dep / 1e9);

    tm.start();
    k2k_bench_sqr<<<blocks, threads>>>(d_out, d_in, iters);
    double t_sqr = tm.stop();
    printf("  sqr:              %.2f Gop/s  (%.2f mul-equivalents)\n",
           total * iters / t_sqr / 1e9, t_sqr / t_dep);

    tm.start();
    k2k_bench_inv<<<blocks, threads>>>(d_out, d_in, iters / 50);
    double t_inv = tm.stop();
    printf("  inv:              %.2f Mop/s  (%.0f mul-equivalents)\n",
           total * (iters / 50) / t_inv / 1e6, (t_inv / (iters / 50)) / (t_dep / iters));

    tm.start();
    k2k_bench_weight<<<blocks, threads, r2k_smem_bytes()>>>(d_u32, d_in, d_cb, iters);
    double t_w = tm.stop();
    printf("  class weight:     %.2f Gop/s  (%.2f mul-equivalents)\n",
           total * iters / t_w / 1e9, t_w / t_dep);

    cudaFree(d_in); cudaFree(d_out); cudaFree(d_cb); cudaFree(d_u32);
}

/* ---------------------------------------------------------------- */
struct RhoOpts {
    uint32_t threads = 0, iters = 256, w = 8, dp = 0;
    std::string variant = "lowmem";
};

template <int W>
static void launch(const RhoOpts &o, rho2k_ctx &dc, uint32_t blocks, size_t smem) {
    if (o.variant == "reg")
        k2k_rho_walk<W><<<blocks, R2K_BLOCK, smem>>>(dc, o.iters);
    else if (o.variant == "ref")
        k2k_rho_walk_ref<<<blocks, R2K_BLOCK, smem>>>(dc, o.iters);
    else
        k2k_rho_walk_lowmem<W><<<blocks, R2K_BLOCK, smem>>>(dc, o.iters);
}

static void bench_rho(RhoOpts o) {
    if (!o.threads) o.threads = (uint32_t)g_sms * 2048 / o.w;
    o.threads = (o.threads + R2K_BLOCK - 1) / R2K_BLOCK * R2K_BLOCK;
    if (!o.dp) o.dp = (uint32_t)(F2M_M / 2 - 3.5 * sqrt(F2M_M / 4.0));

    Rho2kHost h;
    h.prm.nj = 8; h.prm.jmin = 3;
    h.prm.dp_threshold = o.dp;
    h.prm.max_steps = 100u << 20;
    h.cb = cb_flat();
    h.P = Koblitz::generator();
    uint32_t sk[SC_WORDS] = {0xdeadbeefu, 0x1234u, 0, 0};
    h.Q = Koblitz::mul(h.P, sk);
    h.build();

    uint32_t T = o.threads, W = o.w, nw = T * W;
    const uint32_t dpcap = 1u << 20;
    rho2k_ctx dc{};
    dc.nthreads = T; dc.walks_per_thread = W;
    dc.P = h.P; dc.Q = h.Q; dc.prm = h.prm; dc.dp_cap = dpcap;
    CU(cudaMalloc(&dc.X, F2M_WORDS * (size_t)nw * 4));
    CU(cudaMalloc(&dc.Y, F2M_WORDS * (size_t)nw * 4));
    CU(cudaMalloc(&dc.steps, (size_t)nw * 4));
    CU(cudaMalloc(&dc.restarts, (size_t)nw * 4));
    uint32_t *d_cb;
    CU(cudaMalloc(&d_cb, cb_words * 4));
    CU(cudaMemcpy(d_cb, cb_flat(), cb_words * 4, cudaMemcpyHostToDevice));
    dc.cb = d_cb;
    CU(cudaMalloc(&dc.dp_out, (size_t)dpcap * sizeof(rho2k_dp)));
    CU(cudaMalloc(&dc.dp_count, 4));
    CU(cudaMemset(dc.dp_count, 0, 4));

    size_t smem = r2k_smem_bytes();
    uint32_t blocks = (T + R2K_BLOCK - 1) / R2K_BLOCK;
    printf("  %u class walks (%u threads x %u), dp g <= %u, variant=%s\n",
           nw, T, W, o.dp, o.variant.c_str());
    printf("  %u blocks x %d threads, %zu B shared\n", blocks, R2K_BLOCK, smem);

    k2k_rho_init<<<(T + 127) / 128, 128>>>(dc);
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());

    Timer tm;
    if (W == 8) launch<8>(o, dc, blocks, smem);
    else if (W == 16) launch<16>(o, dc, blocks, smem);
    else if (W == 32) launch<32>(o, dc, blocks, smem);
    else { fprintf(stderr, "unsupported --w %u (use 8, 16 or 32)\n", W); exit(1); }
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());

    tm.start();
    if (W == 8) launch<8>(o, dc, blocks, smem);
    else if (W == 16) launch<16>(o, dc, blocks, smem);
    else launch<32>(o, dc, blocks, smem);
    double t = tm.stop();
    CU(cudaGetLastError());

    double steps = (double)nw * o.iters;
    printf("  %.3f Gstep/s  (%.0f M steps in %.3f s)\n", steps / t / 1e9, steps / 1e6, t);
    printf("  %.1f ns/step/walk\n", t / steps * 1e9);
    /* what that means for the challenge */
    double need = sqrt(3.14159265 * ldexp(1.0, CURVE2K_R_BITS) / (4.0 * F2M_M));
    printf("  ECC2K-%d needs ~%.3e class-walk steps: %.1f device-years at this rate\n",
           CURVE2K_R_BITS, need, need / (steps / t) / (365.25 * 24 * 3600));

    cudaFree(dc.X); cudaFree(dc.Y); cudaFree(dc.steps); cudaFree(dc.restarts);
    cudaFree(d_cb); cudaFree(dc.dp_out); cudaFree(dc.dp_count);
}

int main(int argc, char **argv) {
    print_device();
    std::string cmd = argc > 1 ? argv[1] : "selftest";
    RhoOpts o;
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { return (i + 1 < argc) ? (uint32_t)strtoul(argv[++i], nullptr, 0) : 0u; };
        if (a == "--walks") o.threads = next();
        else if (a == "--iters") o.iters = next();
        else if (a == "--w") o.w = next();
        else if (a == "--dp") o.dp = next();
        else if (a == "--variant" && i + 1 < argc) o.variant = argv[++i];
        else { fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    if (cmd == "selftest") selftest();
    else if (cmd == "field") bench_field();
    else if (cmd == "rho") bench_rho(o);
    else { fprintf(stderr, "usage: %s {selftest|field|rho}\n", argv[0]); return 2; }
    return failures ? 1 : 0;
}
