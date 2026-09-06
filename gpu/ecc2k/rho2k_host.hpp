/* rho2k_host.hpp -- host side of the Koblitz rho: walk replay, collision
 * detection and the linear solve that turns a class collision into a
 * discrete logarithm.
 *
 * The device reports only {canonical x, walk, restart, steps}.  Everything
 * that needs arithmetic modulo r happens here, and only for the two walks
 * that actually collided.
 */
#ifndef GPU_ECC2K_RHO_HOST_HPP
#define GPU_ECC2K_RHO_HOST_HPP

#include <string>
#include <unordered_map>
#include <vector>
#include <cstdio>

#include "rho2k.cuh"
#include "scalar_ring.hpp"

struct Rho2kHost {
    rho2k_params prm;
    pt2k P, Q;
    const uint32_t *cb;
    std::unordered_map<std::string, rho2k_dp> seen;
    unsigned long long useless_collisions = 0;

    /* (1 + s^j) mod r for each j the walk can produce, in Montgomery form. */
    std::vector<sc_t> step_factor;

    void build() {
        step_factor.assign(prm.jmin + prm.nj, Sc::zero());
        sc_t s = Sc::s();
        for (uint32_t j = prm.jmin; j < prm.jmin + prm.nj; j++)
            step_factor[j] = Sc::add(Sc::one(), Sc::pow_u32(s, j));
    }

    /* Replay walk (walk, restart) for `steps` steps, tracking the
     * coefficients of the current point as a*P + b*Q.  Uses the same
     * phase_a / phase_b the kernel uses. */
    bool replay(uint32_t walk, uint32_t restart, uint32_t steps,
                sc_t &a, sc_t &b, pt2k &end) const {
        uint32_t sa[SC_WORDS], sb[SC_WORDS];
        r2k_walk_seed(walk, restart, sa, sb);
        rho2k_state st;
        st.P = Koblitz::mul2(P, sa, Q, sb);
        st.escape = 0;
        if (st.P.inf) return false;
        a = Sc::from_limbs(sa);
        b = Sc::from_limbs(sb);
        for (uint32_t i = 0; i < steps; i++) {
            f2e den;
            uint32_t j;
            if (r2k_phase_a(st, prm, cb, den, j) == R2K_MODE_INF) return false;
            r2k_phase_b(st, prm, j, F2::inv(den));
            /* one step multiplies both coefficients by (1 + s^j) */
            a = Sc::mul(a, step_factor[j]);
            b = Sc::mul(b, step_factor[j]);
        }
        end = st.P;
        return true;
    }

    static std::string key_of(const rho2k_dp &d) {
        return std::string((const char *)d.x, sizeof(d.x));
    }

    /* Find (e, sign) with A == sign * tau^e(B), searching the 2m-element
     * class.  Returns e, sets *negated; -1 if the two are not related. */
    static int relate(const pt2k &A, const pt2k &B, int *negated) {
        pt2k cur = B;
        for (int e = 0; e < F2M_M; e++) {
            if (Koblitz::eq(A, cur)) { *negated = 0; return e; }
            if (Koblitz::eq(A, Koblitz::neg(cur))) { *negated = 1; return e; }
            cur = Koblitz::frob(cur, 1);
        }
        return -1;
    }

    /* Feed one distinguished point; on a useful collision fill k (canonical
     * limbs, k*P == Q) and return true. */
    bool add_dp(const rho2k_dp &d, uint32_t k_out[SC_WORDS]) {
        std::string key = key_of(d);
        auto it = seen.find(key);
        if (it == seen.end()) {
            seen.emplace(key, d);
            return false;
        }
        const rho2k_dp &e = it->second;
        if (e.walk == d.walk && e.restart == d.restart && e.steps == d.steps)
            return false;                       /* the same report twice */

        sc_t a1, b1, a2, b2;
        pt2k P1, P2;
        if (!replay(e.walk, e.restart, e.steps, a1, b1, P1)) return false;
        if (!replay(d.walk, d.restart, d.steps, a2, b2, P2)) return false;

        int neg = 0;
        int rot = relate(P1, P2, &neg);
        if (rot < 0) {
            fprintf(stderr, "rho2k: reported collision is not a class collision\n");
            return false;
        }

        /* P1 = sigma * tau^rot(P2), and tau = s on the subgroup, so
         *     a1 + b1 k = sigma s^rot (a2 + b2 k)   (mod r)
         * =>  k = (sigma s^rot a2 - a1) / (b1 - sigma s^rot b2). */
        sc_t f = Sc::pow_u32(Sc::s(), (uint32_t)rot);
        sc_t fa = Sc::mul(f, a2), fb = Sc::mul(f, b2);
        if (neg) { fa = Sc::neg(fa); fb = Sc::neg(fb); }
        sc_t num = Sc::sub(fa, a1);
        sc_t den = Sc::sub(b1, fb);
        if (Sc::is_zero(den)) { useless_collisions++; return false; }
        sc_t k = Sc::mul(num, Sc::inv(den));
        sc_t kc = Sc::to_canonical(k);
        for (int i = 0; i < SC_WORDS; i++) k_out[i] = kc.v[i];

        pt2k chk = Koblitz::mul(P, k_out);
        if (!Koblitz::eq(chk, Q)) {
            fprintf(stderr, "rho2k: candidate k failed verification\n");
            return false;
        }
        return true;
    }
};

#endif /* GPU_ECC2K_RHO_HOST_HPP */
