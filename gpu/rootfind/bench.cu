/* bench.cu -- device driver for the root-finding oracle.
 *
 *   ./bench selftest [bits] [fb] [targets]   every (lookup, chunk) CUDA
 *                                            configuration vs host emulation
 *   ./bench throughput [bits] [fb] [targets] targets/s and (target,i) pairs/s
 *   ./bench p256 [fb] [targets]              same two checks on P-256
 *
 * Build: make bench [ARCH=sm_120]        (needs nvcc; CUDA 12.8+ for sm_120)
 *
 * UNTESTED ON HARDWARE: written on a machine without an NVIDIA GPU.  The
 * host side of every check is the emulation backend that
 * test_cpu.cpp and examples/gpu_rootfind_validate.rs verify; `selftest`
 * is the step that closes the gap between that and the device.  Exit
 * code is non-zero if any configuration disagrees with the emulation.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "rootfind_host.hpp"

typedef unsigned __int128 u128;
typedef std::set<std::tuple<uint32_t, uint32_t, uint32_t, int, int>> RelSet;

static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t p) { return (uint64_t)(((u128)a * b) % p); }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t p) {
    uint64_t r = 1;
    while (e) { if (e & 1) r = mulmod(r, a, p); a = mulmod(a, a, p); e >>= 1; }
    return r;
}
static bool is_prime(uint64_t n) {
    if (n < 2) return false;
    for (uint64_t q : {2ull, 3ull, 5ull, 7ull, 11ull, 13ull}) if (n % q == 0) return n == q;
    uint64_t d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (uint64_t a : {2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull}) {
        uint64_t x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { comp = false; break; } }
        if (comp) return false;
    }
    return true;
}

/* a = -3 curve over the largest p = 3 mod 4 below 2^bits; small-x factor
 * base; targets = positive controls F_i +- F_j followed by random points. */
struct Instance {
    uint64_t p;
    std::vector<uint64_t> fbx, fby, tx, ty;
};

static Instance make_instance(int bits, uint32_t n_fb, uint32_t n_targets, uint64_t seed) {
    Instance I;
    uint64_t p = (1ull << bits) - 1;
    while (!(p % 4 == 3 && is_prime(p))) p--;
    I.p = p;
    Fp64Ctx f = fp64_ctx_new(p);
    uint64_t a = fp64_to_mont(f, p - 3), b = fp64_to_mont(f, 7);
    auto lift = [&](uint64_t xc, uint64_t &xm, uint64_t &ym) {
        xm = fp64_to_mont(f, xc);
        uint64_t rhs = fp64_add(f, fp64_add(f, fp64_mul(f, fp64_sqr(f, xm), xm), fp64_mul(f, a, xm)), b);
        uint64_t rc = fp64_from_mont(f, rhs);
        uint64_t yc = powmod(rc, (p + 1) / 4, p);
        if (mulmod(yc, yc, p) != rc) return false;
        ym = fp64_to_mont(f, yc);
        if (ym > (p - 1) / 2) ym = fp64_neg(f, ym);
        return true;
    };
    std::vector<uint64_t> mx, my;
    for (uint64_t x = 1; mx.size() < n_fb; x++) {
        uint64_t xm, ym;
        if (lift(x, xm, ym)) { mx.push_back(xm); my.push_back(ym); I.fbx.push_back(x); I.fby.push_back(fp64_from_mont(f, ym)); }
    }
    std::mt19937_64 rng(seed);
    uint32_t n_ctl = n_targets / 4;
    for (uint32_t k = 0; k < n_ctl; k++) {
        uint32_t i = rng() % n_fb, j = rng() % n_fb;
        if (i == j) { k--; continue; }
        uint64_t yj = (rng() & 1) ? fp64_neg(f, my[j]) : my[j];
        uint64_t lam = fp64_mul(f, fp64_sub(f, yj, my[i]), fp64_inv(f, fp64_sub(f, mx[j], mx[i])));
        uint64_t x3 = fp64_sub(f, fp64_sub(f, fp64_sqr(f, lam), mx[i]), mx[j]);
        uint64_t y3 = fp64_sub(f, fp64_mul(f, lam, fp64_sub(f, mx[i], x3)), my[i]);
        I.tx.push_back(fp64_from_mont(f, x3)); I.ty.push_back(fp64_from_mont(f, y3));
    }
    while (I.tx.size() < n_targets) {
        uint64_t xm, ym;
        if (lift(rng() % p, xm, ym)) { I.tx.push_back(fp64_from_mont(f, xm)); I.ty.push_back(fp64_from_mont(f, ym)); }
    }
    return I;
}

static RelSet to_set(const std::vector<rf_relation> &v, uint32_t n) {
    RelSet s;
    for (uint32_t i = 0; i < n; i++) s.insert(std::make_tuple(v[i].target, v[i].i, v[i].j, (int)v[i].eps, (int)v[i].tau));
    return s;
}

static int run64(const Instance &I, int backend, int lookup, uint32_t chunk, RelSet &out, rf_result &res) {
    std::vector<rf_relation> buf(I.tx.size() * 4 + 1024);
    rf_options opt{backend, lookup, chunk};
    int rc = rf64_decompose(I.p, I.fbx.data(), I.fby.data(), (uint32_t)I.fbx.size(), I.tx.data(), I.ty.data(),
                            (uint32_t)I.tx.size(), buf.data(), (uint32_t)buf.size(), &opt, &res);
    out = to_set(buf, res.n_written);
    return rc;
}

static const char *lname(int l) { return l == RF_LOOKUP_SORTED ? "sorted" : "hash"; }

static int selftest(int bits, uint32_t n_fb, uint32_t n_targets) {
    std::printf("selftest: %d-bit p, |F| = %u, %u targets (%u positive controls)\n", bits, n_fb, n_targets, n_targets / 4);
    Instance I = make_instance(bits, n_fb, n_targets, 1);
    RelSet ref; rf_result rr;
    int rc = run64(I, RF_BACKEND_EMULATE, RF_LOOKUP_SORTED, 32, ref, rr);
    if (rc != RF_OK) { std::printf("emulation failed: %s\n", rf_strerror(rc)); return 1; }
    std::printf("  emulation: %zu relations in %.3f s\n", ref.size(), rr.seconds);
    if (ref.size() < n_targets / 4) { std::printf("  FAIL: fewer relations than positive controls\n"); return 1; }
    int bad = 0;
    for (int lookup : {RF_LOOKUP_SORTED, RF_LOOKUP_HASH}) {
        for (uint32_t chunk : {16u, 32u, 64u}) {
            RelSet got; rf_result r;
            rc = run64(I, RF_BACKEND_CUDA, lookup, chunk, got, r);
            bool ok = rc == RF_OK && got == ref;
            std::printf("  cuda %-6s chunk=%-2u : %s (%zu relations, %.3f ms)%s\n", lname(lookup), chunk, ok ? "OK  " : "FAIL",
                        got.size(), r.seconds * 1e3, rc != RF_OK ? rf_strerror(rc) : "");
            if (!ok) bad++;
        }
    }
    std::printf(bad ? "SELFTEST FAILED (%d configurations)\n" : "SELFTEST PASSED\n", bad);
    return bad ? 1 : 0;
}

static int throughput(int bits, uint32_t n_fb, uint32_t n_targets) {
    std::printf("throughput: %d-bit p, |F| = %u, %u targets\n", bits, n_fb, n_targets);
    Instance I = make_instance(bits, n_fb, n_targets, 2);
    for (int lookup : {RF_LOOKUP_SORTED, RF_LOOKUP_HASH}) {
        for (uint32_t chunk : {16u, 32u, 64u}) {
            RelSet got; rf_result r;
            int rc = run64(I, RF_BACKEND_CUDA, lookup, chunk, got, r);
            if (rc != RF_OK) { std::printf("  %s chunk=%u: %s\n", lname(lookup), chunk, rf_strerror(rc)); continue; }
            double pairs = (double)n_fb * n_targets;
            std::printf("  cuda %-6s chunk=%-2u : %.3f ms  %.3g targets/s  %.3g (target,i) pairs/s\n", lname(lookup), chunk,
                        r.seconds * 1e3, n_targets / r.seconds, pairs / r.seconds);
        }
    }
    return 0;
}

static int p256(uint32_t n_fb, uint32_t n_targets) {
    std::printf("p256: |F| = %u, %u targets\n", n_fb, n_targets);
    Fp256Ctx f = fp256_ctx_p256();
    U256x64 a = f.p; a.v[0] -= 3;
    U256x64 b = {{0x3BCE3C3E27D2604Bull, 0x651D06B0CC53B0F6ull, 0xB3EBBD55769886BCull, 0x5AC635D8AA3A93E7ull}};
    U256x64 am = fp256_to_mont(f, a), bm = fp256_to_mont(f, b);
    U256x64 q;
    for (int i = 0; i < 4; i++) q.v[i] = (f.p.v[i] >> 2) | (i < 3 ? f.p.v[i + 1] << 62 : 0);
    for (int i = 0; i < 4; i++) { if (++q.v[i] != 0) break; }
    auto lift = [&](const U256x64 &xc, U256x64 &xm, U256x64 &ym) {
        xm = fp256_to_mont(f, xc);
        U256x64 rhs = fp256_add(f, fp256_add(f, fp256_mul(f, fp256_sqr(f, xm), xm), fp256_mul(f, am, xm)), bm);
        U256x64 acc = f.one;
        for (int i = 255; i >= 0; i--) { acc = fp256_sqr(f, acc); if ((q.v[i >> 6] >> (i & 63)) & 1) acc = fp256_mul(f, acc, rhs); }
        if (!u256_eq(fp256_sqr(f, acc), rhs)) return false;
        ym = acc; return true;
    };
    std::vector<U256x64> mx, my;
    std::vector<uint64_t> fbx, fby, tx, ty;
    for (uint64_t x = 1; mx.size() < n_fb; x++) {
        U256x64 xc = {{x, 0, 0, 0}}, xm, ym;
        if (!lift(xc, xm, ym)) continue;
        mx.push_back(xm); my.push_back(ym);
        U256x64 yc = fp256_from_mont(f, ym);
        for (int l = 0; l < 4; l++) { fbx.push_back(xc.v[l]); fby.push_back(yc.v[l]); }
    }
    std::mt19937_64 rng(3);
    for (uint32_t k = 0; k < n_targets; k++) {
        uint32_t i = rng() % n_fb, j = rng() % n_fb;
        if (i == j) { k--; continue; }
        U256x64 yj = (rng() & 1) ? fp256_neg(f, my[j]) : my[j];
        U256x64 lam = fp256_mul(f, fp256_sub(f, yj, my[i]), fp256_inv(f, fp256_sub(f, mx[j], mx[i])));
        U256x64 x3 = fp256_sub(f, fp256_sub(f, fp256_sqr(f, lam), mx[i]), mx[j]);
        U256x64 y3 = fp256_sub(f, fp256_mul(f, lam, fp256_sub(f, mx[i], x3)), my[i]);
        U256x64 xc = fp256_from_mont(f, x3), yc = fp256_from_mont(f, y3);
        for (int l = 0; l < 4; l++) { tx.push_back(xc.v[l]); ty.push_back(yc.v[l]); }
    }
    std::vector<rf_relation> buf(n_targets * 4 + 1024);
    rf_result rr;
    rf_options eo{RF_BACKEND_EMULATE, RF_LOOKUP_SORTED, 16};
    int rc = rf256_decompose(fbx.data(), fby.data(), n_fb, tx.data(), ty.data(), n_targets, buf.data(), (uint32_t)buf.size(), &eo, &rr);
    if (rc != RF_OK) { std::printf("emulation failed: %s\n", rf_strerror(rc)); return 1; }
    RelSet ref = to_set(buf, rr.n_written);
    std::printf("  emulation: %zu relations in %.3f s\n", ref.size(), rr.seconds);
    int bad = 0;
    for (int lookup : {RF_LOOKUP_SORTED, RF_LOOKUP_HASH}) {
        for (uint32_t chunk : {16u, 32u}) {
            rf_options opt{RF_BACKEND_CUDA, lookup, chunk};
            rf_result r;
            rc = rf256_decompose(fbx.data(), fby.data(), n_fb, tx.data(), ty.data(), n_targets, buf.data(), (uint32_t)buf.size(), &opt, &r);
            RelSet got = to_set(buf, rc == RF_OK ? r.n_written : 0);
            bool ok = rc == RF_OK && got == ref;
            std::printf("  cuda %-6s chunk=%-2u : %s (%zu relations, %.3f ms, %.3g pairs/s)%s\n", lname(lookup), chunk, ok ? "OK  " : "FAIL",
                        got.size(), r.seconds * 1e3, (double)n_fb * n_targets / r.seconds, rc != RF_OK ? rf_strerror(rc) : "");
            if (!ok) bad++;
        }
    }
    std::printf(bad ? "P256 SELFTEST FAILED (%d)\n" : "P256 SELFTEST PASSED\n", bad);
    return bad ? 1 : 0;
}

int main(int argc, char **argv) {
    std::printf("backend: %s\n", rf_backend_info());
    if (!rf_cuda_available()) { std::printf("no CUDA device; nothing to do\n"); return 2; }
    std::string cmd = argc > 1 ? argv[1] : "selftest";
    auto arg = [&](int i, long def) { return argc > i ? std::atol(argv[i]) : def; };
    if (cmd == "selftest") return selftest((int)arg(2, 61), (uint32_t)arg(3, 1024), (uint32_t)arg(4, 4096));
    if (cmd == "throughput") return throughput((int)arg(2, 61), (uint32_t)arg(3, 65536), (uint32_t)arg(4, 65536));
    if (cmd == "p256") return p256((uint32_t)arg(2, 256), (uint32_t)arg(3, 1024));
    std::printf("usage: %s selftest|throughput|p256 [args]\n", argv[0]);
    return 2;
}
