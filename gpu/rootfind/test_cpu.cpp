/* test_cpu.cpp -- host verification of the root-finding oracle headers.
 *
 * Compiles fp64.cuh / fp256x64.cuh / rootfind.cuh with a plain C++
 * compiler (no CUDA) and checks:
 *   1. Fp64 Montgomery arithmetic against unsigned __int128 reference
 *      arithmetic (mul, inv, batch_inv, to/from_mont round trip);
 *   2. P-256 Montgomery arithmetic against a schoolbook 512-bit
 *      reference reduction, and the published Montgomery constants;
 *   3. the oracle on a small synthetic a = -3 curve: positive controls
 *      R = F_i +- F_j must be found, every reported relation must re-add
 *      to R in the group, and the sorted / hash table variants and the
 *      three CHUNK sizes must produce identical relation sets;
 *   4. the C ABI (rf64_decompose) fail-closed paths: bad p, coordinate
 *      >= p, output overflow.
 * Exits non-zero on the first failure.
 *
 * The *bit-for-bit* cross-check against the Rust `prime_fast` oracle is
 * in gpu/rootfind/rs (gpu_rootfind_validate); this file is what `make test`
 * runs without cargo.
 */
#include <cstdio>
#include <cstdlib>
#include <random>
#include <set>
#include <tuple>

#include "rootfind_host.hpp"

static int failures = 0;
#define CHECK(cond, ...)                                                   \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::printf("FAIL %s:%d: ", __FILE__, __LINE__);               \
            std::printf(__VA_ARGS__);                                      \
            std::printf("\n");                                             \
            failures++;                                                    \
        }                                                                  \
    } while (0)

typedef unsigned __int128 u128;

/* --- reference helpers ---------------------------------------------- */

static uint64_t ref_mulmod(uint64_t a, uint64_t b, uint64_t p) { return (uint64_t)(((u128)a * b) % p); }
static uint64_t ref_powmod(uint64_t a, uint64_t e, uint64_t p) {
    uint64_t r = 1 % p;
    while (e) { if (e & 1) r = ref_mulmod(r, a, p); a = ref_mulmod(a, a, p); e >>= 1; }
    return r;
}
static bool is_prime_u64(uint64_t n) {
    if (n < 2) return false;
    for (uint64_t q : {2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull})
        if (n % q == 0) return n == q;
    uint64_t d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (uint64_t a : {2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull}) {
        uint64_t x = ref_powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s; r++) { x = ref_mulmod(x, x, n); if (x == n - 1) { comp = false; break; } }
        if (comp) return false;
    }
    return true;
}

/* --- 1. Fp64 -------------------------------------------------------- */

static void test_fp64(uint64_t p, std::mt19937_64 &rng) {
    Fp64Ctx f = fp64_ctx_new(p);
    CHECK(f.one == (uint64_t)((((u128)1) << 64) % p), "R mod p");
    for (int it = 0; it < 20000; it++) {
        uint64_t a = rng() % p, b = rng() % p;
        uint64_t am = fp64_to_mont(f, a), bm = fp64_to_mont(f, b);
        CHECK(fp64_from_mont(f, am) == a, "to/from mont round trip a=%llu", (unsigned long long)a);
        CHECK(fp64_from_mont(f, fp64_mul(f, am, bm)) == ref_mulmod(a, b, p), "mul %llu*%llu", (unsigned long long)a, (unsigned long long)b);
        CHECK(fp64_from_mont(f, fp64_add(f, am, bm)) == (uint64_t)(((u128)a + b) % p), "add");
        CHECK(fp64_from_mont(f, fp64_sub(f, am, bm)) == (a >= b ? a - b : a + p - b), "sub");
        if (a) {
            uint64_t inv = fp64_inv(f, am);
            CHECK(fp64_mul(f, inv, am) == f.one, "inv*a == 1 for a=%llu", (unsigned long long)a);
        }
    }
    /* edge values */
    uint64_t pm1 = fp64_to_mont(f, p - 1);
    CHECK(fp64_from_mont(f, fp64_mul(f, pm1, pm1)) == 1, "(p-1)^2 == 1");
    CHECK(fp64_inv(f, 0) == 0, "inv(0) == 0");
    /* batch inversion */
    for (int n : {1, 2, 7, 16, 32, 64}) {
        std::vector<uint64_t> x(n), s(n), orig(n);
        for (int i = 0; i < n; i++) { do { orig[i] = rng() % p; } while (orig[i] == 0); x[i] = fp64_to_mont(f, orig[i]); }
        std::vector<uint64_t> xm = x;
        fp64_batch_inv(f, x.data(), n, s.data());
        for (int i = 0; i < n; i++) CHECK(fp64_mul(f, x[i], xm[i]) == f.one, "batch_inv n=%d i=%d", n, i);
    }
}

/* --- 2. P-256 ------------------------------------------------------- */

/* 512-bit schoolbook product and reduction mod p by repeated subtraction
 * of shifted p -- slow but obviously correct. */
struct U512 { uint64_t v[8]; };
static U512 mul_wide(const U256x64 &a, const U256x64 &b) {
    U512 r{}; for (int i = 0; i < 8; i++) r.v[i] = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; j++) {
            u128 acc = (u128)a.v[j] * b.v[i] + r.v[i + j] + carry;
            r.v[i + j] = (uint64_t)acc; carry = (uint64_t)(acc >> 64);
        }
        r.v[i + 4] = carry;
    }
    return r;
}
static bool ge512(const U512 &a, const U512 &b) {
    for (int i = 7; i >= 0; i--) { if (a.v[i] != b.v[i]) return a.v[i] > b.v[i]; }
    return true;
}
static void sub512(U512 &a, const U512 &b) {
    uint64_t bw = 0;
    for (int i = 0; i < 8; i++) {
        uint64_t d = a.v[i] - b.v[i]; uint64_t b1 = a.v[i] < b.v[i];
        uint64_t t = d - bw; uint64_t b2 = d < bw;
        a.v[i] = t; bw = b1 | b2;
    }
}
static U512 shl512(const U512 &a, int bits) {
    U512 r{}; for (int i = 0; i < 8; i++) r.v[i] = 0;
    int w = bits / 64, s = bits % 64;
    for (int i = 7; i >= 0; i--) {
        if (i - w < 0) break;
        uint64_t lo = a.v[i - w] << s;
        uint64_t hi = (s && i - w - 1 >= 0) ? a.v[i - w - 1] >> (64 - s) : 0;
        r.v[i] = lo | hi;
    }
    return r;
}
static U256x64 ref_mod(U512 x, const U256x64 &p) {
    U512 P{}; for (int i = 0; i < 4; i++) { P.v[i] = p.v[i]; P.v[i + 4] = 0; }
    for (int sh = 256; sh >= 0; sh--) {
        U512 t = shl512(P, sh);
        if (ge512(x, t)) sub512(x, t);
    }
    U256x64 r; for (int i = 0; i < 4; i++) r.v[i] = x.v[i];
    return r;
}
static U256x64 ref_mulmod256(const U256x64 &a, const U256x64 &b, const U256x64 &p) { return ref_mod(mul_wide(a, b), p); }

static U256x64 rand256_below(std::mt19937_64 &rng, const U256x64 &p) {
    U256x64 r;
    for (;;) {
        for (int i = 0; i < 4; i++) r.v[i] = rng();
        r.v[3] &= p.v[3];   /* quick reject most; then exact compare */
        bool lt = false;
        for (int i = 3; i >= 0; i--) { if (r.v[i] != p.v[i]) { lt = r.v[i] < p.v[i]; break; } }
        if (lt) return r;
    }
}

static void test_fp256(std::mt19937_64 &rng) {
    Fp256Ctx f = fp256_ctx_p256();
    /* constants: R mod p and R^2 mod p from the canonical definition */
    U512 R{}; for (int i = 0; i < 8; i++) R.v[i] = 0; R.v[4] = 1;
    U256x64 r1 = ref_mod(R, f.p);
    CHECK(u256_eq(r1, f.one), "R mod p constant");
    U256x64 r2 = ref_mulmod256(r1, r1, f.p);
    CHECK(u256_eq(r2, f.r2), "R^2 mod p constant");
    CHECK(f.p.v[0] * f.pinv == (uint64_t)-1, "pinv: p0 * pinv == -1 mod 2^64");

    U256x64 one = {{1, 0, 0, 0}};
    for (int it = 0; it < 3000; it++) {
        U256x64 a = rand256_below(rng, f.p), b = rand256_below(rng, f.p);
        U256x64 am = fp256_to_mont(f, a), bm = fp256_to_mont(f, b);
        CHECK(u256_eq(fp256_from_mont(f, am), a), "p256 to/from mont");
        U256x64 prod = fp256_from_mont(f, fp256_mul(f, am, bm));
        CHECK(u256_eq(prod, ref_mulmod256(a, b, f.p)), "p256 mul it=%d", it);
        /* add / sub against reference via 512 */
        U512 s{}; for (int i = 0; i < 8; i++) s.v[i] = 0;
        uint64_t c = 0;
        for (int i = 0; i < 4; i++) { u128 t = (u128)a.v[i] + b.v[i] + c; s.v[i] = (uint64_t)t; c = (uint64_t)(t >> 64); }
        s.v[4] = c;
        CHECK(u256_eq(fp256_from_mont(f, fp256_add(f, am, bm)), ref_mod(s, f.p)), "p256 add");
        if (!u256_is_zero(a)) {
            U256x64 inv = fp256_inv(f, am);
            CHECK(u256_eq(fp256_mul(f, inv, am), f.one), "p256 inv it=%d", it);
        }
    }
    /* p-1 squared is 1 */
    U256x64 pm1 = f.p; pm1.v[0] -= 1;
    U256x64 pm1m = fp256_to_mont(f, pm1);
    CHECK(u256_eq(fp256_from_mont(f, fp256_mul(f, pm1m, pm1m)), one), "(p-1)^2 == 1 (p256)");
    for (int n : {1, 5, 16}) {
        std::vector<U256x64> x(n), s(n), orig(n);
        for (int i = 0; i < n; i++) { x[i] = fp256_to_mont(f, rand256_below(rng, f.p)); orig[i] = x[i]; }
        fp256_batch_inv(f, x.data(), n, s.data());
        for (int i = 0; i < n; i++) CHECK(u256_eq(fp256_mul(f, x[i], orig[i]), f.one), "p256 batch_inv n=%d", n);
    }
}

/* --- 3. oracle on a synthetic curve -------------------------------- */

struct Pt { uint64_t x, y; bool inf; };

struct Curve {
    Fp64Ctx f; uint64_t a, b;   /* Montgomery */
    bool on_curve(Pt P) const {
        if (P.inf) return true;
        uint64_t lhs = fp64_sqr(f, P.y);
        uint64_t rhs = fp64_add(f, fp64_add(f, fp64_mul(f, fp64_sqr(f, P.x), P.x), fp64_mul(f, a, P.x)), b);
        return lhs == rhs;
    }
    Pt neg(Pt P) const { if (!P.inf) P.y = fp64_neg(f, P.y); return P; }
    Pt add(Pt P, Pt Q) const {
        if (P.inf) return Q;
        if (Q.inf) return P;
        uint64_t lam;
        if (P.x == Q.x) {
            if (P.y != Q.y || P.y == 0) return Pt{0, 0, true};
            uint64_t x2 = fp64_sqr(f, P.x);
            uint64_t num = fp64_add(f, fp64_add(f, fp64_add(f, x2, x2), x2), a);
            lam = fp64_mul(f, num, fp64_inv(f, fp64_add(f, P.y, P.y)));
        } else {
            lam = fp64_mul(f, fp64_sub(f, Q.y, P.y), fp64_inv(f, fp64_sub(f, Q.x, P.x)));
        }
        uint64_t x3 = fp64_sub(f, fp64_sub(f, fp64_sqr(f, lam), P.x), Q.x);
        uint64_t y3 = fp64_sub(f, fp64_mul(f, lam, fp64_sub(f, P.x, x3)), P.y);
        return Pt{x3, y3, false};
    }
    /* sqrt for p = 3 mod 4 */
    bool lift_x(uint64_t xc, Pt &out) const {
        uint64_t x = fp64_to_mont(f, xc);
        uint64_t rhs = fp64_add(f, fp64_add(f, fp64_mul(f, fp64_sqr(f, x), x), fp64_mul(f, a, x)), b);
        uint64_t rc = fp64_from_mont(f, rhs);
        uint64_t yc = ref_powmod(rc, (f.p + 1) / 4, f.p);
        if (ref_mulmod(yc, yc, f.p) != rc) return false;
        uint64_t y = fp64_to_mont(f, yc);
        if (y > (f.p - 1) / 2) y = fp64_neg(f, y);
        out = Pt{x, y, false};
        return true;
    }
};

typedef std::set<std::tuple<uint32_t, uint32_t, uint32_t, int, int>> RelSet;

static RelSet run_api(uint64_t p, const std::vector<uint64_t> &fbx, const std::vector<uint64_t> &fby,
                      const std::vector<uint64_t> &tx, const std::vector<uint64_t> &ty, int lookup, uint32_t chunk,
                      rf_result &res, int &rc) {
    std::vector<rf_relation> out(4 * tx.size() + 16);
    rf_options opt{RF_BACKEND_EMULATE, lookup, chunk};
    rc = rf64_decompose(p, fbx.data(), fby.data(), (uint32_t)fbx.size(), tx.data(), ty.data(), (uint32_t)tx.size(),
                        out.data(), (uint32_t)out.size(), &opt, &res);
    RelSet s;
    for (uint32_t i = 0; i < res.n_written; i++) s.insert(std::make_tuple(out[i].target, out[i].i, out[i].j, (int)out[i].eps, (int)out[i].tau));
    return s;
}

static void test_oracle(uint64_t p, std::mt19937_64 &rng) {
    Curve c; c.f = fp64_ctx_new(p); c.a = fp64_to_mont(c.f, p - 3); c.b = 0;
    /* find b giving a non-singular curve with a point at x=1.. (order
     * irrelevant: this test only exercises the group law) */
    for (uint64_t b = 1;; b++) {
        c.b = fp64_to_mont(c.f, b);
        uint64_t disc = fp64_add(c.f, fp64_mul(c.f, fp64_to_mont(c.f, 4), fp64_mul(c.f, fp64_sqr(c.f, c.a), c.a)),
                                 fp64_mul(c.f, fp64_to_mont(c.f, 27), fp64_sqr(c.f, c.b)));
        if (disc != 0) break;
    }
    const uint32_t N = 300;
    std::vector<Pt> fb;
    std::vector<uint64_t> fbx, fby;
    for (uint64_t x = 1; fb.size() < N; x++) {
        Pt P;
        if (c.lift_x(x, P)) { fb.push_back(P); fbx.push_back(fp64_from_mont(c.f, P.x)); fby.push_back(fp64_from_mont(c.f, P.y)); }
    }
    /* targets: positive controls, then random points, then R = F_i */
    std::vector<Pt> targets;
    struct Ctl { uint32_t i, j; int ei, ej; };
    std::vector<Ctl> ctls;
    for (int k = 0; k < 60; k++) {
        uint32_t i = rng() % N, j = rng() % N;
        int ei = (rng() & 1) ? 1 : -1, ej = (rng() & 1) ? 1 : -1;
        Pt Pi = ei < 0 ? c.neg(fb[i]) : fb[i];
        Pt Pj = ej < 0 ? c.neg(fb[j]) : fb[j];
        Pt R = c.add(Pi, Pj);
        if (R.inf) { k--; continue; }
        targets.push_back(R); ctls.push_back(Ctl{i, j, ei, ej});
    }
    const size_t n_ctl = targets.size();
    for (int k = 0; k < 200; k++) {
        Pt P;
        if (c.lift_x(rng() % p, P)) targets.push_back(P); else k--;
    }
    targets.push_back(fb[7]);
    targets.push_back(c.neg(fb[11]));
    std::vector<uint64_t> tx, ty;
    for (auto &R : targets) { tx.push_back(fp64_from_mont(c.f, R.x)); ty.push_back(fp64_from_mont(c.f, R.y)); }

    rf_result res; int rc;
    RelSet base = run_api(p, fbx, fby, tx, ty, RF_LOOKUP_SORTED, 32, res, rc);
    CHECK(rc == RF_OK, "rf64_decompose rc=%d (%s)", rc, rf_strerror(rc));
    CHECK(res.backend_used == RF_BACKEND_EMULATE, "backend");

    /* every relation re-adds to its target */
    for (auto &r : base) {
        uint32_t k, i, j; int eps, tau;
        std::tie(k, i, j, eps, tau) = r;
        Pt Pi = eps < 0 ? c.neg(fb[i]) : fb[i];
        Pt sum = Pi;
        if (tau != 0) { Pt Pj = tau < 0 ? c.neg(fb[j]) : fb[j]; sum = c.add(Pi, Pj); }
        CHECK(!sum.inf && sum.x == targets[k].x && sum.y == targets[k].y, "relation (k=%u,i=%u,j=%u,%d,%d) does not re-add", k, i, j, eps, tau);
        CHECK(j >= i, "j >= i convention");
    }
    /* every positive control is found: R = ei F_i + ej F_j must appear as
     * (min, max) with matching signs (or as a degenerate when i == j and
     * the signs agree). */
    for (size_t k = 0; k < n_ctl; k++) {
        const Ctl &ct = ctls[k];
        uint32_t i = ct.i, j = ct.j; int ei = ct.ei, ej = ct.ej;
        if (i > j) { std::swap(i, j); std::swap(ei, ej); }
        bool found;
        if (i == j && ei != ej) { continue; }   /* R = O, excluded above */
        if (i == j) {
            /* R = 2 e F_i: reported as (i, e, i, e) */
            found = base.count(std::make_tuple((uint32_t)k, i, i, ei, ei)) > 0;
        } else {
            found = base.count(std::make_tuple((uint32_t)k, i, j, ei, ej)) > 0;
        }
        CHECK(found, "positive control k=%zu (i=%u,%d j=%u,%d) missing", k, i, ei, j, ej);
    }
    /* degenerate targets */
    CHECK(base.count(std::make_tuple((uint32_t)(targets.size() - 2), 7u, 7u, 1, 0)) > 0, "R = F_7 degenerate");
    CHECK(base.count(std::make_tuple((uint32_t)(targets.size() - 1), 11u, 11u, -1, 0)) > 0, "R = -F_11 degenerate");

    /* all table/chunk variants agree */
    for (int lookup : {RF_LOOKUP_SORTED, RF_LOOKUP_HASH}) {
        for (uint32_t chunk : {16u, 32u, 64u}) {
            RelSet s = run_api(p, fbx, fby, tx, ty, lookup, chunk, res, rc);
            CHECK(rc == RF_OK, "variant lookup=%d chunk=%u rc=%d", lookup, chunk, rc);
            CHECK(s == base, "variant lookup=%d chunk=%u differs (%zu vs %zu)", lookup, chunk, s.size(), base.size());
        }
    }
    std::printf("  oracle: %zu targets (%zu controls), %zu relations, all variants agree\n", targets.size(), n_ctl, base.size());

    /* 4. fail-closed API paths */
    rf_options opt{RF_BACKEND_EMULATE, RF_LOOKUP_AUTO, 0};
    std::vector<rf_relation> tiny(1);
    rc = rf64_decompose(p, fbx.data(), fby.data(), N, tx.data(), ty.data(), (uint32_t)tx.size(), tiny.data(), 1, &opt, &res);
    CHECK(rc == RF_ERR_OVERFLOW && res.n_found == base.size() && res.n_written == 1, "overflow detection rc=%d found=%u", rc, res.n_found);
    rc = rf64_decompose(p + 1, fbx.data(), fby.data(), N, tx.data(), ty.data(), 1, tiny.data(), 1, &opt, &res);
    CHECK(rc == RF_ERR_BAD_ARG, "even p rejected");
    std::vector<uint64_t> bad = tx; bad[0] = p;
    rc = rf64_decompose(p, fbx.data(), fby.data(), N, bad.data(), ty.data(), 1, tiny.data(), 1, &opt, &res);
    CHECK(rc == RF_ERR_BAD_ARG, "coordinate >= p rejected");
    opt.chunk = 24;
    rc = rf64_decompose(p, fbx.data(), fby.data(), N, tx.data(), ty.data(), 1, tiny.data(), 1, &opt, &res);
    CHECK(rc == RF_ERR_BAD_CHUNK, "bad chunk rejected");
    opt.chunk = 0; opt.backend = RF_BACKEND_CUDA;
    rc = rf64_decompose(p, fbx.data(), fby.data(), N, tx.data(), ty.data(), 1, tiny.data(), 1, &opt, &res);
    CHECK(rc == RF_ERR_NO_CUDA || rc == RF_OK || rc == RF_ERR_OVERFLOW, "cuda request without cuda -> RF_ERR_NO_CUDA (rc=%d)", rc);
}

/* --- P-256 oracle smoke: positive controls on the real curve --------- */

static void test_oracle_p256(std::mt19937_64 &rng) {
    Fp256Ctx f = fp256_ctx_p256();
    /* a = -3, b = P-256's b */
    U256x64 a = f.p; a.v[0] -= 3;
    U256x64 b = {{0x3BCE3C3E27D2604Bull, 0x651D06B0CC53B0F6ull, 0xB3EBBD55769886BCull, 0x5AC635D8AA3A93E7ull}};
    U256x64 am = fp256_to_mont(f, a), bm = fp256_to_mont(f, b);
    auto rhs = [&](const U256x64 &x) {
        return fp256_add(f, fp256_add(f, fp256_mul(f, fp256_sqr(f, x), x), fp256_mul(f, am, x)), bm);
    };
    /* sqrt via (p+1)/4: P-256 has p = 3 mod 4 */
    auto sqrt_m = [&](const U256x64 &v, U256x64 &out) {
        /* (p+1)/4 = (p >> 2) + 1 since p = 3 mod 4 */
        U256x64 q;
        for (int i = 0; i < 4; i++) q.v[i] = (f.p.v[i] >> 2) | (i < 3 ? f.p.v[i + 1] << 62 : 0);
        for (int i = 0; i < 4; i++) { if (++q.v[i] != 0) break; }   /* +1 with carry */
        U256x64 acc = f.one;
        for (int i = 255; i >= 0; i--) {
            acc = fp256_sqr(f, acc);
            if ((q.v[i >> 6] >> (i & 63)) & 1) acc = fp256_mul(f, acc, v);
        }
        if (!u256_eq(fp256_sqr(f, acc), v)) return false;
        out = acc; return true;
    };
    struct P { U256x64 x, y; };
    std::vector<P> fb;
    std::vector<uint64_t> fbx, fby;
    for (uint64_t xc = 1; fb.size() < 64; xc++) {
        U256x64 x = {{xc, 0, 0, 0}};
        U256x64 xm = fp256_to_mont(f, x), y;
        if (!sqrt_m(rhs(xm), y)) continue;
        fb.push_back(P{xm, y});
        U256x64 yc = fp256_from_mont(f, y);
        for (int l = 0; l < 4; l++) { fbx.push_back(x.v[l]); fby.push_back(yc.v[l]); }
    }
    auto add = [&](P p1, P p2, bool neg2) {
        if (neg2) p2.y = fp256_neg(f, p2.y);
        U256x64 lam = fp256_mul(f, fp256_sub(f, p2.y, p1.y), fp256_inv(f, fp256_sub(f, p2.x, p1.x)));
        U256x64 x3 = fp256_sub(f, fp256_sub(f, fp256_sqr(f, lam), p1.x), p2.x);
        U256x64 y3 = fp256_sub(f, fp256_mul(f, lam, fp256_sub(f, p1.x, x3)), p1.y);
        return P{x3, y3};
    };
    std::vector<uint64_t> tx, ty;
    struct Ctl { uint32_t i, j; int tau; };
    std::vector<Ctl> ctls;
    for (int k = 0; k < 20; k++) {
        uint32_t i = rng() % 64, j = rng() % 64;
        if (i == j) { k--; continue; }
        if (i > j) std::swap(i, j);
        bool neg = rng() & 1;
        P R = add(fb[i], fb[j], neg);
        U256x64 xc = fp256_from_mont(f, R.x), yc = fp256_from_mont(f, R.y);
        for (int l = 0; l < 4; l++) { tx.push_back(xc.v[l]); ty.push_back(yc.v[l]); }
        ctls.push_back(Ctl{i, j, neg ? -1 : 1});
    }
    std::vector<rf_relation> out(256);
    rf_result res;
    for (int lookup : {RF_LOOKUP_SORTED, RF_LOOKUP_HASH}) {
        rf_options opt{RF_BACKEND_EMULATE, lookup, 16};
        int rc = rf256_decompose(fbx.data(), fby.data(), 64, tx.data(), ty.data(), (uint32_t)ctls.size(), out.data(), 256, &opt, &res);
        CHECK(rc == RF_OK, "rf256_decompose rc=%d", rc);
        for (size_t k = 0; k < ctls.size(); k++) {
            bool found = false;
            for (uint32_t r = 0; r < res.n_written; r++)
                if (out[r].target == k && out[r].i == ctls[k].i && out[r].j == ctls[k].j && out[r].eps == 1 && out[r].tau == ctls[k].tau) found = true;
            CHECK(found, "p256 control k=%zu missing (lookup=%d)", k, lookup);
        }
    }
    std::printf("  p256 oracle: %zu controls found (sorted + hash)\n", ctls.size());
}

int main() {
    std::mt19937_64 rng(0x5EED);
    /* primes: a 61-bit p = 3 mod 4, a 40-bit one, and a tiny one */
    uint64_t p61 = (1ull << 61) - 1;   /* Mersenne prime, = 3 mod 4 */
    uint64_t p40 = (1ull << 40) - 1; while (!(p40 % 4 == 3 && is_prime_u64(p40))) p40--;
    uint64_t p20 = (1ull << 20) - 1; while (!(p20 % 4 == 3 && is_prime_u64(p20))) p20--;
    CHECK(is_prime_u64(p61), "p61 prime");

    std::printf("=== fp64 ===\n");
    for (uint64_t p : {p61, p40, p20, (uint64_t)1000003}) test_fp64(p, rng);
    std::printf("=== fp256 (P-256) ===\n");
    test_fp256(rng);
    std::printf("=== oracle (synthetic a=-3 curves) ===\n");
    for (uint64_t p : {p61, p40}) test_oracle(p, rng);
    std::printf("=== oracle (P-256 positive controls) ===\n");
    test_oracle_p256(rng);
    std::printf("backend: %s\n", rf_backend_info());
    if (failures) { std::printf("%d FAILURES\n", failures); return 1; }
    std::printf("ALL PASS\n");
    return 0;
}
