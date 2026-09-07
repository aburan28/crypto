/* kangaroo.cu -- GPU solver and benchmark for the interval ECDLP.
 *
 *   ./bench selftest                 kernels against the host implementation
 *   ./bench list                     the puzzle registry and its cost model
 *   ./bench bench [--bits N]         walk throughput
 *   ./bench solve --bits N --pubkey <hex>
 *   ./bench solve --puzzle N         target taken from puzzles.txt
 *   ./bench solve --bits N --self    generate a key and solve it (a real test)
 *
 * Common options: --walks N --w W --dp BITS --jumps BITS --iters N
 *                 --variant reg|lowmem|ref
 *
 * The host code this checks against is the same code test_kangaroo.cpp
 * verifies, so `selftest` closes the loop on the device.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

#include "kernels_kangaroo.cuh"
#include "kangaroo_host.hpp"

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

static void print_device() {
    int dev = 0;
    cudaDeviceProp p;
    CU(cudaGetDevice(&dev));
    CU(cudaGetDeviceProperties(&p, dev));
    g_sms = p.multiProcessorCount;
    printf("device: %s  sm_%d%d  %d SMs  %.1f GHz  %d threads/SM  %d KB shared/SM\n",
           p.name, p.major, p.minor, g_sms, p.clockRate / 1e6,
           p.maxThreadsPerMultiProcessor, (int)(p.sharedMemPerMultiprocessor >> 10));
}

/* ---------------------------------------------------------------- *
 * device-side context
 * ---------------------------------------------------------------- */
struct DeviceRun {
    kg_ctx dc{};
    uint32_t nkang = 0;
    kg_jump *d_jumps = nullptr;
    uint32_t dpcap = 1u << 20;
    size_t smem = 0;
    uint32_t blocks = 0;

    void alloc(KangarooHost &h, uint32_t T, uint32_t W) {
        nkang = T * W;
        dc.nthreads = T;
        dc.kang_per_thread = W;
        dc.Qshift = h.Qshift;
        dc.prm = h.prm;
        dc.dp_cap = dpcap;
        CU(cudaMalloc(&dc.X, 8 * (size_t)nkang * 4));
        CU(cudaMalloc(&dc.Y, 8 * (size_t)nkang * 4));
        CU(cudaMalloc(&dc.D, 8 * (size_t)nkang * 4));
        CU(cudaMalloc(&dc.steps, (size_t)nkang * 4));
        CU(cudaMalloc(&dc.restarts, (size_t)nkang * 4));
        CU(cudaMalloc(&d_jumps, h.jumps.size() * sizeof(kg_jump)));
        CU(cudaMemcpy(d_jumps, h.jumps.data(), h.jumps.size() * sizeof(kg_jump),
                      cudaMemcpyHostToDevice));
        dc.jumps = d_jumps;
        CU(cudaMalloc(&dc.dp_out, (size_t)dpcap * sizeof(kg_dp)));
        CU(cudaMalloc(&dc.dp_count, 4));
        CU(cudaMemset(dc.dp_count, 0, 4));
        smem = kg_smem_bytes(h.prm.njump_bits);
        blocks = (T + KG_BLOCK - 1) / KG_BLOCK;
    }

    void free_all() {
        cudaFree(dc.X); cudaFree(dc.Y); cudaFree(dc.D);
        cudaFree(dc.steps); cudaFree(dc.restarts);
        cudaFree(d_jumps); cudaFree(dc.dp_out); cudaFree(dc.dp_count);
    }
};

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

struct Opts {
    uint32_t threads = 0, iters = 256, w = 8, dp = 0, jumps = 6;
    int bits = 40, puzzle = 0;
    bool self = false;
    std::string variant = "lowmem", pubkey;
};

template <int W>
static void launch(const Opts &o, kg_ctx &dc, uint32_t blocks, size_t smem, uint32_t iters) {
    if (o.variant == "reg") k_kang_walk<W><<<blocks, KG_BLOCK, smem>>>(dc, iters);
    else if (o.variant == "ref") k_kang_walk_ref<<<blocks, KG_BLOCK, smem>>>(dc, iters);
    else k_kang_walk_lowmem<W><<<blocks, KG_BLOCK, smem>>>(dc, iters);
}

static void launch_any(const Opts &o, kg_ctx &dc, uint32_t blocks, size_t smem,
                       uint32_t iters, uint32_t W) {
    if (W == 8) launch<8>(o, dc, blocks, smem, iters);
    else if (W == 16) launch<16>(o, dc, blocks, smem, iters);
    else if (W == 32) launch<32>(o, dc, blocks, smem, iters);
    else { fprintf(stderr, "unsupported --w %u (use 8, 16 or 32)\n", W); exit(1); }
}

/* ---------------------------------------------------------------- */
static void selftest() {
    printf("[selftest] kangaroo kernels vs the host\n");
    KangarooHost h;
    h.prm.njump_bits = 6;
    h.prm.dp_mask = (1u << 8) - 1u;
    h.prm.max_steps = 1u << 20;
    h.prm.seed = 17;
    h.prm.reseed_on_dp = 1;
    u256 secret = u256_pow2(39);
    secret.v[0] ^= 0x1234567u;
    affine_pt Q = Curve::to_affine(Curve::scalar_mul(Curve::generator(), secret.v, 0));
    h.setup(Q, u256_pow2(39), 39);

    const uint32_t T = 256, W = 8, nk = T * W, cap = 1 << 14;
    kg_ctx hc{};
    std::vector<uint32_t> hX(8 * nk), hY(8 * nk), hD(8 * nk), hS(nk), hR(nk);
    std::vector<kg_dp> hdp(cap);
    uint32_t hcount = 0;
    hc.X = hX.data(); hc.Y = hY.data(); hc.D = hD.data();
    hc.steps = hS.data(); hc.restarts = hR.data();
    hc.nthreads = T; hc.kang_per_thread = W;
    hc.jumps = h.jumps.data(); hc.Qshift = h.Qshift; hc.prm = h.prm;
    hc.dp_out = hdp.data(); hc.dp_count = &hcount; hc.dp_cap = cap;
    for (uint32_t t = 0; t < T; t++) kg_init_thread(hc, t);

    DeviceRun dr;
    dr.dpcap = cap;
    dr.alloc(h, T, W);
    k_kang_init<<<(T + 127) / 128, 128>>>(dr.dc);
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());

    const uint32_t iters = 64;
    Opts o;
    o.variant = "reg";
    launch<W>(o, dr.dc, dr.blocks, dr.smem, iters);
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());
    for (uint32_t it = 0; it < iters; it++)
        for (uint32_t t = 0; t < T; t++) kg_step_batch<W>(hc, t);

    std::vector<uint32_t> gX(8 * nk), gY(8 * nk), gD(8 * nk), gS(nk), gR(nk);
    CU(cudaMemcpy(gX.data(), dr.dc.X, gX.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(gY.data(), dr.dc.Y, gY.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(gD.data(), dr.dc.D, gD.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(gS.data(), dr.dc.steps, gS.size() * 4, cudaMemcpyDeviceToHost));
    CU(cudaMemcpy(gR.data(), dr.dc.restarts, gR.size() * 4, cudaMemcpyDeviceToHost));
    CHECK(gX == hX && gY == hY, "kangaroo positions after %u iterations", iters);
    CHECK(gD == hD, "kangaroo distances");
    CHECK(gS == hS && gR == hR, "step and restart counters");
    uint32_t gcount = 0;
    CU(cudaMemcpy(&gcount, dr.dc.dp_count, 4, cudaMemcpyDeviceToHost));
    CHECK(gcount == hcount, "dp count: device %u vs host %u", gcount, hcount);
    if (!failures)
        printf("  %u kangaroos x %u steps identical to the host, %u DPs\n",
               nk, iters, gcount);

    /* the walk invariant, checked on the device's own output */
    {
        uint32_t bad = 0;
        for (uint32_t i = 0; i < nk && i < 256; i++) {
            u256 d;
            for (int l = 0; l < 8; l++) d.v[l] = gD[l * nk + i];
            affine_pt want = Curve::to_affine(
                Curve::scalar_mul(Curve::generator(), d.v, 0));
            if (kg_herd_of(i) == KG_HERD_WILD)
                want = Curve::to_affine(Curve::madd(Curve::to_jac(want), h.Qshift));
            fp256 wx = want.x;
            for (int l = 0; l < 8; l++)
                if (wx.v[l] != gX[l * nk + i]) { bad++; break; }
        }
        CHECK(bad == 0, "%u device kangaroos violate distance*G == position", bad);
        if (!bad) printf("  distance/position invariant holds on device state\n");
    }

    /* the low-memory variant must agree */
    {
        CU(cudaMemset(dr.dc.dp_count, 0, 4));
        k_kang_init<<<(T + 127) / 128, 128>>>(dr.dc);
        Opts o2;
        o2.variant = "lowmem";
        launch<W>(o2, dr.dc, dr.blocks, dr.smem, iters);
        CU(cudaGetLastError());
        CU(cudaDeviceSynchronize());
        std::vector<uint32_t> lX(8 * nk), lD(8 * nk);
        CU(cudaMemcpy(lX.data(), dr.dc.X, lX.size() * 4, cudaMemcpyDeviceToHost));
        CU(cudaMemcpy(lD.data(), dr.dc.D, lD.size() * 4, cudaMemcpyDeviceToHost));
        CHECK(lX == hX && lD == hD, "lowmem variant state");
        if (lX == hX && lD == hD) printf("  lowmem variant matches\n");
    }

    dr.free_all();
    printf(failures ? "SELFTEST FAILED\n" : "selftest OK\n");
}

/* ---------------------------------------------------------------- */
static void cmd_list() {
    PuzzleRegistry reg;
    if (!reg.load("puzzles.txt")) {
        printf("puzzles.txt not found (run from gpu/btcpuzzle)\n");
        return;
    }
    for (const auto &p : reg.problems) printf("REGISTRY PROBLEM: %s\n", p.c_str());
    printf("%-4s %-8s %-10s %s\n", "n", "interval", "kangaroo", "status");
    for (const auto &e : reg.entries) {
        double steps = 2.0 * ldexp(1.0, (e.n - 1) / 2.0);
        printf("%-4d 2^%-6d 2^%-8.1f %s%s\n", e.n, e.n - 1, log2(steps),
               e.has_pub() ? "public key present" : "no public key: not attackable by kangaroo",
               e.known_key.empty() ? "" : ", key known (regression target)");
    }
    printf("\nkangaroo cost is ~2*sqrt(W) group operations for an interval of\n"
           "width W = 2^(n-1); measured constant on this implementation is about\n"
           "3 (see README).  Without a public key the search is over the whole\n"
           "interval instead, 2^(n-1) rather than 2^((n-1)/2).\n");
}

/* ---------------------------------------------------------------- */
static void run_solver(Opts &o, KangarooHost &h, const u256 *known) {
    if (!o.threads) o.threads = (uint32_t)g_sms * 2048 / o.w;
    o.threads = (o.threads + KG_BLOCK - 1) / KG_BLOCK * KG_BLOCK;
    if (o.threads < KG_BLOCK) o.threads = KG_BLOCK;

    DeviceRun dr;
    dr.alloc(h, o.threads, o.w);
    printf("  %u kangaroos (%u threads x %u), dp %u bits, %u jumps, %zu B shared\n",
           dr.nkang, o.threads, o.w, 32 - __builtin_clz(h.prm.dp_mask + 1) - 1,
           1u << h.prm.njump_bits, dr.smem);
    printf("  expected ~%.3e group operations\n", h.expected_steps());

    k_kang_init<<<(o.threads + 127) / 128, 128>>>(dr.dc);
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());

    std::vector<kg_dp> host_dp(dr.dpcap);
    uint32_t consumed = 0;
    double steps = 0, elapsed = 0;
    Timer tm;
    u256 key;
    bool solved = false;

    for (int round = 0; round < 1000000 && !solved; round++) {
        tm.start();
        launch_any(o, dr.dc, dr.blocks, dr.smem, o.iters, o.w);
        elapsed += tm.stop();
        CU(cudaGetLastError());
        steps += (double)dr.nkang * o.iters;

        uint32_t cnt = 0;
        CU(cudaMemcpy(&cnt, dr.dc.dp_count, 4, cudaMemcpyDeviceToHost));
        if (cnt > dr.dpcap) {
            printf("  distinguished-point buffer overflowed; raise --dp\n");
            break;
        }
        if (cnt > consumed) {
            CU(cudaMemcpy(host_dp.data() + consumed, dr.dc.dp_out + consumed,
                          (cnt - consumed) * sizeof(kg_dp), cudaMemcpyDeviceToHost));
            while (consumed < cnt)
                if (h.add_dp(host_dp[consumed++], key)) { solved = true; break; }
        }
        if (round % 64 == 0 && round)
            printf("  %.3e steps, %.2f Gstep/s, %zu distinguished points\n",
                   steps, steps / elapsed / 1e9, h.seen.size());
    }

    printf("  %.3e steps in %.1f s (%.3f Gstep/s), %llu same-herd collisions\n",
           steps, elapsed, steps / elapsed / 1e9, h.same_herd_collisions);
    if (solved) {
        printf("  SOLVED: k = %s\n", u256_hex(key).c_str());
        printf("  took %.2f * sqrt(W)\n", steps / ldexp(1.0, (int)h.prm.w_bits / 2.0));
        if (known) CHECK(u256_cmp(key, *known) == 0, "recovered key differs from the planted one");
    } else {
        printf("  not solved within the round budget\n");
    }
    dr.free_all();
}

static void cmd_solve(Opts &o) {
    KangarooHost h;
    h.prm.njump_bits = o.jumps;
    h.prm.max_steps = 1u << 24;
    h.prm.seed = 1;
    h.prm.reseed_on_dp = 1;

    affine_pt Q;
    u256 known = u256_zero();
    bool have_known = false;

    if (o.puzzle) {
        PuzzleRegistry reg;
        if (!reg.load("puzzles.txt")) { printf("puzzles.txt not found\n"); return; }
        const Puzzle *p = reg.find(o.puzzle);
        if (!p) { printf("puzzle %d is not in puzzles.txt\n", o.puzzle); return; }
        if (!p->has_pub()) {
            printf("puzzle %d has no public key in the registry, so kangaroo cannot\n"
                   "attack it; only the address is known and the search would be\n"
                   "over the whole 2^%d interval.\n", o.puzzle, o.puzzle - 1);
            return;
        }
        if (!pubkey_from_hex(p->pubkey.c_str(), Q)) {
            printf("puzzle %d has an invalid public key\n", o.puzzle);
            return;
        }
        o.bits = p->n;
        if (!p->known_key.empty() && u256_from_hex(p->known_key.c_str(), known))
            have_known = true;
        printf("[solve] puzzle %d\n", p->n);
    } else if (o.self) {
        known = u256_pow2(o.bits - 1);
        uint64_t s = 0xC0FFEEull * o.bits + 12345;
        u256 off = u256_zero();
        for (int i = 0; i < 4; i++) {
            uint64_t z = kg_splitmix64(s);
            off.v[2 * i] = (uint32_t)z;
            off.v[2 * i + 1] = (uint32_t)(z >> 32);
        }
        for (int l = 0; l < 8; l++) {
            int lo = 32 * l;
            if (lo >= o.bits - 1) off.v[l] = 0;
            else if (lo + 32 > o.bits - 1) off.v[l] &= (1u << (o.bits - 1 - lo)) - 1u;
        }
        known = u256_add(known, off);
        have_known = true;
        Q = Curve::to_affine(Curve::scalar_mul(Curve::generator(), known.v, 0));
        printf("[solve] generated a %d-bit key and will look for it\n", o.bits);
    } else if (!o.pubkey.empty()) {
        if (!pubkey_from_hex(o.pubkey.c_str(), Q)) {
            printf("public key is not a point on secp256k1\n");
            return;
        }
        printf("[solve] %d-bit interval, supplied public key\n", o.bits);
    } else {
        printf("need --puzzle N, --pubkey <hex>, or --self\n");
        return;
    }

    if (o.bits < 8 || o.bits > 160) { printf("--bits out of range\n"); return; }
    uint32_t wbits = (uint32_t)(o.bits - 1);
    /* Pick the distinguished-point rate so the per-kangaroo tail stays well
     * under the birthday term: total tail is nkang * 2^dp_bits. */
    if (!o.dp) {
        double root = ldexp(1.0, wbits / 2.0);
        uint32_t nk = ((o.threads ? o.threads : (uint32_t)g_sms * 2048 / o.w)) * o.w;
        double target = root / (8.0 * nk);
        int b = 0;
        while (b < 24 && ldexp(1.0, b + 1) < target) b++;
        o.dp = (uint32_t)b;
    }
    h.prm.dp_mask = (1u << o.dp) - 1u;
    h.setup(Q, u256_pow2(o.bits - 1), wbits);
    run_solver(o, h, have_known ? &known : nullptr);
}

static void cmd_bench(Opts &o) {
    KangarooHost h;
    h.prm.njump_bits = o.jumps;
    h.prm.dp_mask = (1u << (o.dp ? o.dp : 20)) - 1u;
    h.prm.max_steps = 1u << 24;
    h.prm.seed = 3;
    h.prm.reseed_on_dp = 1;
    u256 secret = u256_pow2(o.bits - 1);
    secret.v[0] ^= 0xbeefu;
    affine_pt Q = Curve::to_affine(Curve::scalar_mul(Curve::generator(), secret.v, 0));
    h.setup(Q, u256_pow2(o.bits - 1), (uint32_t)(o.bits - 1));

    if (!o.threads) o.threads = (uint32_t)g_sms * 2048 / o.w;
    o.threads = (o.threads + KG_BLOCK - 1) / KG_BLOCK * KG_BLOCK;
    DeviceRun dr;
    dr.alloc(h, o.threads, o.w);
    printf("  %u kangaroos, variant=%s, %zu B shared\n",
           dr.nkang, o.variant.c_str(), dr.smem);

    k_kang_init<<<(o.threads + 127) / 128, 128>>>(dr.dc);
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());

    Timer tm;
    launch_any(o, dr.dc, dr.blocks, dr.smem, o.iters, o.w);
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());
    tm.start();
    launch_any(o, dr.dc, dr.blocks, dr.smem, o.iters, o.w);
    double t = tm.stop();
    CU(cudaGetLastError());
    double steps = (double)dr.nkang * o.iters;
    printf("  %.3f Gstep/s  (%.1f ns per step per kangaroo)\n",
           steps / t / 1e9, t / steps * 1e9);
    printf("  at this rate a %d-bit interval takes ~%.2e s = %.2f device-days\n",
           o.bits, h.expected_steps() / (steps / t),
           h.expected_steps() / (steps / t) / 86400.0);
    dr.free_all();
}

int main(int argc, char **argv) {
    print_device();
    std::string cmd = argc > 1 ? argv[1] : "selftest";
    Opts o;
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        auto nx = [&]() { return (i + 1 < argc) ? (uint32_t)strtoul(argv[++i], nullptr, 0) : 0u; };
        if (a == "--walks") o.threads = nx();
        else if (a == "--iters") o.iters = nx();
        else if (a == "--w") o.w = nx();
        else if (a == "--dp") o.dp = nx();
        else if (a == "--jumps") o.jumps = nx();
        else if (a == "--bits") o.bits = (int)nx();
        else if (a == "--puzzle") o.puzzle = (int)nx();
        else if (a == "--self") o.self = true;
        else if (a == "--pubkey" && i + 1 < argc) o.pubkey = argv[++i];
        else if (a == "--variant" && i + 1 < argc) o.variant = argv[++i];
        else { fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    if (cmd == "selftest") selftest();
    else if (cmd == "list") cmd_list();
    else if (cmd == "solve") cmd_solve(o);
    else if (cmd == "bench") cmd_bench(o);
    else { fprintf(stderr, "usage: %s {selftest|list|bench|solve}\n", argv[0]); return 2; }
    return failures ? 1 : 0;
}
