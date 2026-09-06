/* rho_host.hpp -- host side of the Pollard-rho pipeline: table construction,
 * distinguished-point bookkeeping, walk replay and the final DLP solve.
 * Shared by the CUDA bench driver and the CPU test harness. */
#ifndef GPU_ECC_RHO_HOST_HPP
#define GPU_ECC_RHO_HOST_HPP

#include <string>
#include <unordered_map>
#include <vector>
#include <cstring>
#include <cstdio>

#include "rho.cuh"

struct RhoHost {
    rho_params prm;
    affine_pt P, Q;
    std::vector<affine_pt> table;   /* M[j] */
    std::vector<fp256> tc, td;      /* c_j, d_j in Fn internal form */
    std::unordered_map<std::string, rho_dp> seen;
    unsigned long long useless_collisions = 0;

    static fp256 scalar_to_fn(const uint32_t k[8]) {
        fp256 t;
        for (int l = 0; l < 8; l++) t.v[l] = k[l];
        return Fn::from_canonical(t);
    }

    void build_table() {
        uint32_t R = 1u << prm.r_bits;
        table.resize(R); tc.resize(R); td.resize(R);
        for (uint32_t j = 0; j < R; j++) {
            uint32_t c[8], d[8];
            rho_table_seed(prm.table_seed, j, c, d);
            table[j] = Curve::to_affine(Curve::double_scalar_mul(P, c, Q, d));
            tc[j] = scalar_to_fn(c);
            td[j] = scalar_to_fn(d);
        }
    }

    /* Replay walk (walk, restart) for `steps` steps, tracking the
     * coefficients of the current point as a*P + b*Q.  Uses the same
     * phase_a / phase_b pair as the kernel, so the trajectories -- cycle
     * escapes included -- cannot diverge.  Returns false if the walk hits
     * infinity earlier than expected (cannot happen for a reported DP, but
     * guards against corrupted input). */
    bool replay(uint32_t walk, uint32_t restart, uint32_t steps,
                fp256 &a, fp256 &b, affine_pt &end) const {
        uint32_t sa[8], sb[8];
        rho_walk_seed(walk, restart, sa, sb);
        rho_state st;
        st.P = Curve::to_affine(Curve::double_scalar_mul(P, sa, Q, sb));
        if (st.P.inf) return false;
        a = scalar_to_fn(sa);
        b = scalar_to_fn(sb);
        if (rho_canonical(st.P, prm)) { a = Fn::neg(a); b = Fn::neg(b); }
        rho_set_hprev(st, st.P);
        st.escape = 0;

        uint32_t done = 0;
        /* Cycle escapes consume an iteration without advancing the step
         * counter, so bound the loop generously rather than by `steps`. */
        for (uint64_t iter = 0; done < steps; iter++) {
            if (iter > (uint64_t)steps * 4 + 1024) return false;
            fp256 den;
            uint32_t j;
            int m = rho_phase_a(st, table.data(), prm, den, j);
            if (m == RHO_MODE_INF) return false;
            int neg, esc_prev;
            int advanced = rho_phase_b(st, table.data(), prm, m, j, Fp::inv(den),
                                       &neg, &esc_prev);
            fp256 na, nb;
            if (m == RHO_MODE_ESCAPE) {
                na = Fn::dbl(a); nb = Fn::dbl(b);
            } else {
                na = Fn::add(a, tc[j]); nb = Fn::add(b, td[j]);
            }
            if (neg) { na = Fn::neg(na); nb = Fn::neg(nb); }
            if (advanced) {
                a = na; b = nb;
                done++;
            } else if (esc_prev) {
                /* the walk rewound to the newly computed point */
                a = na; b = nb;
            }
        }
        end = st.P;
        return true;
    }

    static std::string key_of(const rho_dp &d) {
        return std::string((const char *)d.x, sizeof(d.x));
    }

    /* Feed one distinguished point.  Returns true and fills k_out (canonical
     * limbs, k*P = Q) when a useful collision is found. */
    bool add_dp(const rho_dp &d, uint32_t k_out[8]) {
        std::string key = key_of(d);
        auto it = seen.find(key);
        if (it == seen.end()) {
            seen.emplace(key, d);
            return false;
        }
        const rho_dp &e = it->second;
        if (e.walk == d.walk && e.restart == d.restart) return false;   /* duplicate report */
        fp256 a1, b1, a2, b2;
        affine_pt e1, e2;
        if (!replay(e.walk, e.restart, e.steps, a1, b1, e1)) return false;
        if (!replay(d.walk, d.restart, d.steps, a2, b2, e2)) return false;
        if (!Fp::eq(e1.x, e2.x)) {
            fprintf(stderr, "rho: replay mismatch -- kernel and host walks disagree\n");
            return false;
        }
        /* e1 = s * e2 with s = +-1:  a1 + b1 k = s (a2 + b2 k)
         *   =>  k = (s a2 - a1) / (b1 - s b2) */
        int same_sign = Fp::eq(e1.y, e2.y);
        fp256 sa2 = same_sign ? a2 : Fn::neg(a2);
        fp256 sb2 = same_sign ? b2 : Fn::neg(b2);
        fp256 num = Fn::sub(sa2, a1);
        fp256 den = Fn::sub(b1, sb2);
        if (Fn::is_zero(den)) { useless_collisions++; return false; }
        fp256 k = Fn::mul(num, Fn::inv(den));
        fp256 kc = Fn::to_canonical(k);
        for (int l = 0; l < 8; l++) k_out[l] = kc.v[l];
        affine_pt chk = Curve::to_affine(Curve::scalar_mul(P, k_out));
        if (!Curve::affine_eq(chk, Q)) {
            fprintf(stderr, "rho: candidate k failed verification\n");
            return false;
        }
        return true;
    }
};

#endif /* GPU_ECC_RHO_HOST_HPP */
