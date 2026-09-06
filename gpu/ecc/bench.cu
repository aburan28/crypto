/* bench.cu -- driver for the elliptic-curve GPU kernels.
 *
 *   ./bench selftest        run every kernel against the host implementation
 *   ./bench field           microbenchmark mul / sqr / inv
 *   ./bench mul             batch scalar multiplication throughput
 *   ./bench rho [opts]      Pollard-rho walk throughput (and, on a small
 *                           curve, an actual DLP solve)
 *
 * rho options:  --walks N --iters N --w W --rbits R --dp BITS --neg 0|1
 *               --variant reg|lowmem|ref  --solve
 *
 * Everything the kernels compute is checked against the same host code that
 * test_cpu.cpp verifies against the Python reference, so `selftest` is a
 * true end-to-end check of the device path.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "kernels.cuh"
#include "rho_host.hpp"

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
        else if (a == "--rbits") o.rbits = next();
        else if (a == "--dp") o.dpbits = next();
        else if (a == "--neg") o.neg = next();
        else if (a == "--variant" && i + 1 < argc) o.variant = argv[++i];
        else if (a == "--solve") o.solve = true;
        else { fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    if (cmd == "selftest") selftest();
    else if (cmd == "field") bench_field();
    else if (cmd == "mul") bench_mul();
    else if (cmd == "rho") bench_rho(o);
    else { fprintf(stderr, "usage: %s {selftest|field|mul|rho}\n", argv[0]); return 2; }
    return failures ? 1 : 0;
}
