/* kangaroo_host.hpp -- host side of the kangaroo solver: jump table
 * construction, the distinguished-point table, and turning a tame/wild
 * collision into a private key.
 *
 * Also holds the Bitcoin puzzle registry.  A puzzle entry is nothing more
 * than an interval and a public key; the solver does not know or care that
 * the target is a Bitcoin key.
 */
#ifndef GPU_BTC_KANGAROO_HOST_HPP
#define GPU_BTC_KANGAROO_HOST_HPP

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include "kangaroo.cuh"

/* ---------------------------------------------------------------- *
 * small helpers on 256-bit canonical integers
 * ---------------------------------------------------------------- */
struct u256 {
    uint32_t v[8];
};

inline u256 u256_zero() { u256 r; memset(r.v, 0, sizeof r.v); return r; }

inline u256 u256_pow2(int e) {
    u256 r = u256_zero();
    if (e >= 0 && e < 256) r.v[e >> 5] = 1u << (e & 31);
    return r;
}

inline int u256_cmp(const u256 &a, const u256 &b) {
    for (int i = 7; i >= 0; i--) {
        if (a.v[i] != b.v[i]) return a.v[i] < b.v[i] ? -1 : 1;
    }
    return 0;
}

inline u256 u256_add(const u256 &a, const u256 &b) {
    u256 r;
    mp_add(r.v, a.v, b.v);
    return r;
}

inline u256 u256_sub(const u256 &a, const u256 &b) {
    u256 r;
    mp_sub(r.v, a.v, b.v);
    return r;
}

/* Parse up to 64 hex digits into a 256-bit integer.  Returns false on any
 * character that is not hex. */
inline bool u256_from_hex(const char *s, u256 &out) {
    out = u256_zero();
    size_t n = strlen(s);
    if (n == 0 || n > 64) return false;
    for (size_t i = 0; i < n; i++) {
        int d;
        char c = s[n - 1 - i];
        if (c >= '0' && c <= '9') d = c - '0';
        else if (c >= 'a' && c <= 'f') d = c - 'a' + 10;
        else if (c >= 'A' && c <= 'F') d = c - 'A' + 10;
        else return false;
        out.v[i >> 3] |= (uint32_t)d << (4 * (i & 7));
    }
    return true;
}

inline std::string u256_hex(const u256 &a) {
    char buf[80];
    int i = 7;
    while (i > 0 && a.v[i] == 0) i--;
    int p = snprintf(buf, sizeof buf, "%x", a.v[i]);
    for (int j = i - 1; j >= 0; j--)
        p += snprintf(buf + p, sizeof buf - p, "%08x", a.v[j]);
    return std::string(buf);
}

/* ---------------------------------------------------------------- *
 * SEC1 public keys
 * ---------------------------------------------------------------- */

/* Decompress a 33-byte compressed point "02<x>" / "03<x>", or read an
 * uncompressed "04<x><y>".  secp256k1 has p = 3 mod 4, so the square root
 * is a single exponentiation. */
inline bool pubkey_from_hex(const char *hex, affine_pt &out) {
    size_t n = strlen(hex);
    u256 x;
    if (n == 130 && hex[0] == '0' && hex[1] == '4') {
        char xb[65], yb[65];
        memcpy(xb, hex + 2, 64); xb[64] = 0;
        memcpy(yb, hex + 66, 64); yb[64] = 0;
        u256 y;
        if (!u256_from_hex(xb, x) || !u256_from_hex(yb, y)) return false;
        out.x = Fp::from_limbs(x.v);
        out.y = Fp::from_limbs(y.v);
        out.inf = 0;
        return Curve::affine_on_curve(out);
    }
    if (n != 66 || hex[0] != '0' || (hex[1] != '2' && hex[1] != '3')) return false;
    char xb[65];
    memcpy(xb, hex + 2, 64);
    xb[64] = 0;
    if (!u256_from_hex(xb, x)) return false;
    fp256 fx = Fp::from_limbs(x.v);
    /* y^2 = x^3 + 7 */
    fp256 rhs = Fp::add(Fp::mul(Fp::sqr(fx), fx), Curve::coeff_b());
    /* p = 3 mod 4, so the square root is y = rhs^((p+1)/4).  The exponent
     * below is (p+1)/4 = 0x3fff...bfffff0c; if it were wrong the y^2 == rhs
     * check just below would reject every key, so it is self-testing. */
    static const uint32_t exp_p1_4[8] = {
        0xbfffff0cu, 0xffffffffu, 0xffffffffu, 0xffffffffu,
        0xffffffffu, 0xffffffffu, 0xffffffffu, 0x3fffffffu};
    fp256 y = Fp::one();
    for (int bit = 255; bit >= 0; bit--) {
        y = Fp::sqr(y);
        if ((exp_p1_4[bit >> 5] >> (bit & 31)) & 1u) y = Fp::mul(y, rhs);
    }
    if (!Fp::eq(Fp::sqr(y), rhs)) return false;   /* x is not on the curve */
    /* pick the parity the prefix asks for */
    fp256 yc = Fp::to_canonical(y);
    uint32_t odd = yc.v[0] & 1u;
    uint32_t want = (uint32_t)(hex[1] - '2');
    if (odd != want) y = Fp::neg(y);
    out.x = fx;
    out.y = y;
    out.inf = 0;
    return Curve::affine_on_curve(out);
}

inline std::string pubkey_to_hex(const affine_pt &P) {
    if (P.inf) return "00";
    fp256 xc = Fp::to_canonical(P.x), yc = Fp::to_canonical(P.y);
    u256 x;
    memcpy(x.v, xc.v, sizeof x.v);
    char buf[80];
    int p = snprintf(buf, sizeof buf, "%02x", 2 + (yc.v[0] & 1u));
    for (int j = 7; j >= 0; j--) p += snprintf(buf + p, sizeof buf - p, "%08x", x.v[j]);
    return std::string(buf);
}

/* ---------------------------------------------------------------- *
 * puzzle registry
 * ---------------------------------------------------------------- */
struct Puzzle {
    int n = 0;                  /* puzzle number: key is in [2^(n-1), 2^n) */
    std::string pubkey;         /* compressed SEC1, empty if never exposed */
    std::string known_key;      /* hex, for solved puzzles used as regression tests */
    bool has_pub() const { return !pubkey.empty(); }
};

/* Reads `puzzles.txt`: one entry per line,
 *     <n> <pubkey-hex|-> [<known-key-hex>]
 * with '#' starting a comment.  Entries are validated on load: the public
 * key must decompress to a curve point, and a known key must lie in the
 * puzzle's interval and generate the stated public key. */
struct PuzzleRegistry {
    std::vector<Puzzle> entries;
    std::vector<std::string> problems;

    bool load(const char *path) {
        FILE *f = fopen(path, "r");
        if (!f) return false;
        char line[512];
        while (fgets(line, sizeof line, f)) {
            char *h = strchr(line, '#');
            if (h) *h = 0;
            Puzzle p;
            char pk[256] = {0}, kk[128] = {0};
            int got = sscanf(line, "%d %255s %127s", &p.n, pk, kk);
            if (got < 2 || p.n <= 0 || p.n > 256) continue;
            if (strcmp(pk, "-") != 0) p.pubkey = pk;
            if (got >= 3) p.known_key = kk;
            validate(p);
            entries.push_back(p);
        }
        fclose(f);
        return true;
    }

    /* Anything inconsistent is recorded rather than trusted. */
    void validate(const Puzzle &p) {
        char tag[64];
        snprintf(tag, sizeof tag, "puzzle %d: ", p.n);
        affine_pt Q;
        bool have_q = false;
        if (p.has_pub()) {
            if (!pubkey_from_hex(p.pubkey.c_str(), Q))
                problems.push_back(std::string(tag) + "public key is not a curve point");
            else have_q = true;
        }
        if (!p.known_key.empty()) {
            u256 k;
            if (!u256_from_hex(p.known_key.c_str(), k)) {
                problems.push_back(std::string(tag) + "key is not hex");
                return;
            }
            u256 lo = u256_pow2(p.n - 1), hi = u256_pow2(p.n);
            if (u256_cmp(k, lo) < 0 || u256_cmp(k, hi) >= 0)
                problems.push_back(std::string(tag) + "key is outside [2^(n-1), 2^n)");
            if (have_q) {
                affine_pt chk = Curve::to_affine(
                    Curve::scalar_mul(Curve::generator(), k.v, 0));
                if (!Curve::affine_eq(chk, Q))
                    problems.push_back(std::string(tag) + "key does not generate the public key");
            }
        }
    }

    const Puzzle *find(int n) const {
        for (const auto &p : entries) if (p.n == n) return &p;
        return nullptr;
    }
};

/* ---------------------------------------------------------------- *
 * the solver
 * ---------------------------------------------------------------- */
struct KangarooHost {
    kg_params prm;
    affine_pt Q;                  /* the target public key */
    u256 a;                       /* interval start */
    affine_pt Qshift;             /* Q - a*G */
    std::vector<kg_jump> jumps;

    /* DP table: x -> the first report seen at that point */
    std::unordered_map<std::string, kg_dp> seen;
    unsigned long long same_herd_collisions = 0;

    /* Jump scalars are pseudorandom with mean sqrt(W)/2, the two-herd
     * optimum.  Powers of two would also work but give the reachable
     * position sets an arithmetic structure worth avoiding. */
    void build_jumps() {
        uint32_t nj = 1u << prm.njump_bits;
        jumps.resize(nj);
        /* mean = 2^(w_bits/2 - 1); draw uniformly from [1, 2*mean] */
        int mean_bits = (int)prm.w_bits / 2;    /* 2*mean = 2^mean_bits */
        if (mean_bits < 1) mean_bits = 1;
        affine_pt G = Curve::generator();
        uint64_t s = 0x5DEECE66Dull ^ ((uint64_t)prm.seed << 16) ^ prm.njump_bits;
        for (uint32_t j = 0; j < nj; j++) {
            u256 sc = u256_zero();
            for (int i = 0; i < 4; i++) {
                uint64_t z = kg_splitmix64(s);
                sc.v[2 * i] = (uint32_t)z;
                sc.v[2 * i + 1] = (uint32_t)(z >> 32);
            }
            for (int l = 0; l < 8; l++) {
                int lo = 32 * l;
                if (lo >= mean_bits) sc.v[l] = 0;
                else if (lo + 32 > mean_bits) sc.v[l] &= (1u << (mean_bits - lo)) - 1u;
            }
            /* never zero: a zero jump is a fixed point */
            bool z = true;
            for (int l = 0; l < 8; l++) if (sc.v[l]) { z = false; break; }
            if (z) sc.v[0] = 1;
            memcpy(jumps[j].s, sc.v, sizeof sc.v);
            jumps[j].P = Curve::to_affine(Curve::scalar_mul(G, sc.v, 0));
        }
    }

    /* Q' = Q - centre*G, where centre = a + W/2.
     *
     * Shifting by the MIDDLE of the interval rather than its start is what
     * keeps the two herds on top of each other.  Both herds draw start
     * offsets uniformly from [0, W), so tame positions are centred on W/2;
     * measuring from the interval's centre puts the wild positions there
     * too, whatever the secret happens to be.  Shifting by `a` instead
     * leaves the wild herd centred on k' + W/2, so a key near the top of
     * the interval barely overlaps the tame herd at all. */
    u256 centre;
    void shift_target() {
        affine_pt G = Curve::generator();
        centre = u256_add(a, u256_pow2((int)prm.w_bits - 1));
        affine_pt aG = Curve::to_affine(Curve::scalar_mul(G, centre.v, 0));
        Qshift = Curve::to_affine(Curve::add(Curve::to_jac(Q),
                                             Curve::to_jac(Curve::affine_neg(aG))));
    }

    void setup(const affine_pt &target, const u256 &interval_start, uint32_t w_bits) {
        Q = target;
        a = interval_start;
        prm.w_bits = w_bits;
        build_jumps();
        shift_target();
    }

    /* Mean number of group operations a solve is expected to take. */
    double expected_steps() const { return 2.0 * ldexp(1.0, (int)prm.w_bits / 2); }

    static std::string key_of(const kg_dp &d) {
        return std::string((const char *)d.x, sizeof(d.x));
    }

    /* Feed one distinguished point.  Returns true and fills `key` when a
     * tame/wild pair pins the logarithm down. */
    bool add_dp(const kg_dp &d, u256 &key) {
        std::string k = key_of(d);
        auto it = seen.find(k);
        if (it == seen.end()) {
            seen.emplace(k, d);
            return false;
        }
        const kg_dp &e = it->second;
        if (e.herd == d.herd) {
            /* Two kangaroos of the same herd merged; nothing to learn.  Both
             * were re-seeded by the walk when they reported, so this costs
             * only the steps already spent. */
            same_herd_collisions++;
            return false;
        }
        const kg_dp &tame = (e.herd == KG_HERD_TAME) ? e : d;
        const kg_dp &wild = (e.herd == KG_HERD_TAME) ? d : e;

        /* k' = d_tame - d_wild (mod n),  k = centre + k' */
        fp256 dt, dw;
        memcpy(dt.v, tame.dist, sizeof dt.v);
        memcpy(dw.v, wild.dist, sizeof dw.v);
        fp256 kp = Fn::sub(Fn::from_canonical(dt), Fn::from_canonical(dw));
        fp256 av;
        memcpy(av.v, centre.v, sizeof av.v);
        fp256 kk = Fn::add(kp, Fn::from_canonical(av));
        fp256 kc = Fn::to_canonical(kk);
        u256 cand;
        memcpy(cand.v, kc.v, sizeof cand.v);

        affine_pt chk = Curve::to_affine(
            Curve::scalar_mul(Curve::generator(), cand.v, 0));
        if (!Curve::affine_eq(chk, Q)) {
            /* Can happen if the two kangaroos met but the reported x came
             * from opposite y (the walk does not use a negation map, so
             * this is not expected -- report it rather than hide it). */
            fprintf(stderr, "kangaroo: candidate key failed verification\n");
            return false;
        }
        key = cand;
        return true;
    }
};

#endif /* GPU_BTC_KANGAROO_HOST_HPP */
