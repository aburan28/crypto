#!/usr/bin/env python3
"""
Thread 23, follow-up probe.

glv_hnp_phase2_eliminated.py showed that eliminating d changes sv/pv (0.388 -> 0.750
on the lam*=0.07 curve) yet leaves the K1 wall bit-for-bit identical.  This probe asks
the two questions that finding raises:

  P1. Are the outcomes identical PER SEED, or only in aggregate rate?
  P2. If the planted vector is still not lambda_1 (sv/pv < 1) after d is gone,
      what IS the shortest vector of the eliminated lattice?

Hypothesis for P2: the 2-D "lambda block"  < (n*S_K1, 0), (-lam*S_K1, S_K2) >  sits
inside BOTH lattices (rows j-1 and m+j of the eliminated basis span it), so its exact
Gauss-reduced shortest vector mu is a lattice vector in both.  If ||mu|| < ||v_planted||
then the planted vector cannot be lambda_1 no matter what is done to the d-column.

P3. Is the wall information-theoretic?  Compare against the counting bound
    m_thresh = ceil(log n / log(1/eff)),  eff = K1*K2/n.
"""

import math
import random

from glv_hnp_phase2_eliminated import (
    find_generator, lam_star, norm, scales, gen_signatures,
    build_glv_lattice, planted_vector_orig, recover_d_orig,
    build_elim_lattice, planted_vector_elim, recover_d_elim,
    reduce_matrix,
)

# ---------------------------------------------------------------------------
# exact 2-D Gauss reduction (verbatim from glv_hnp_phase2_lambda_threshold.py)
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
    def nrm2(w): return w[0] * w[0] + w[1] * w[1]
    def dot(w, z): return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v): u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0: break
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u): break
        u, v = v, u
    return u

def lambda_block_mu(n, lam, S_K1, S_K2):
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1]), w


SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

curves = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None
    curves.append((label, (p, b, n, lam, G), k1, m))


print("=" * 78)
print("Thread 23 probe — why elimination changed sv/pv but not recovery")
print("=" * 78)

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("P1: per-seed outcome agreement between orig and elim")
print("-" * 78)
print(f"{'curve':<20} {'K1':>3} {'seed':>7} {'orig':>6} {'elim':>6} {'agree':>6}")
n_agree = n_tot = 0
for label, curve, k1def, m in curves:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for k1 in [4, 6, 8]:
        for s in SEEDS:
            d = random.Random(s * 7919 + 13).randint(1, n - 1)
            sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, s)
            if len(sigs) < m: continue
            S_KAN = scales(n, k1, k2b)[3]

            Mo = build_glv_lattice(sigs, n, lam, k1, k2b)
            ro = reduce_matrix(Mo)
            oo = recover_d_orig(ro, m, n, S_KAN, d) is not None

            Me = build_elim_lattice(sigs, n, lam, k1, k2b)
            re_ = reduce_matrix(Me)
            oe = recover_d_elim(re_, m, n, sigs, lam, k1, k2b, d) is not None

            n_tot += 1
            n_agree += 1 if oo == oe else 0
            if k1 == 6:
                print(f"{label:<20} {k1:>3} {s:>7} {str(oo):>6} {str(oe):>6} "
                      f"{str(oo==oe):>6}")
print(f"\nper-seed agreement over all (curve,K1,seed): {n_agree}/{n_tot}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("P2: what is lambda_1?  ||mu(lambda-block)|| vs ||v_planted|| vs ||sv||")
print("-" * 78)
print("mu = exact shortest vector of the 2-D block <(n*S_K1,0), (-lam*S_K1,S_K2)>,")
print("which is a sublattice of BOTH the original and the eliminated lattice.")
print()
print(f"{'curve':<20} {'K1':>3} {'||mu||':>12} {'||pv||':>12} {'mu/pv':>7} "
      f"{'sv_elim':>12} {'sv=mu?':>7}")
for label, curve, k1def, m in curves:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for k1 in [k1def]:
        S_K1, _sd, S_K2, S_KAN = scales(n, k1, k2b)
        mu, wvec = lambda_block_mu(n, lam, S_K1, S_K2)
        d = random.Random(42 * 7919 + 13).randint(1, n - 1)
        sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, 42)
        pv = norm(planted_vector_elim(sigs, n, k1, k2b))
        Me = build_elim_lattice(sigs, n, lam, k1, k2b)
        re_ = reduce_matrix(Me)
        svrow = min((r for r in re_ if any(r)), key=norm)
        sv = norm(svrow)
        print(f"{label:<20} {k1:>3} {mu:>12.1f} {pv:>12.1f} {mu/pv:>7.3f} "
              f"{sv:>12.1f} {str(abs(sv-mu)<1e-6):>7}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("P2b: energy profile of the eliminated lattice's shortest vector")
print("-" * 78)
print("cols: 0..m-2 congruence | m-1 k1_0 | m..2m-1 k2_j | 2m Kannan")
print(f"{'curve':<20} {'K1':>3} {'cong':>7} {'k1_0':>7} {'k2':>7} {'kan':>7} {'nnz':>5}")
for label, curve, k1def, m in curves:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, _sd, S_K2, S_KAN = scales(n, k1, k2b)
    d = random.Random(42 * 7919 + 13).randint(1, n - 1)
    sigs = gen_signatures(G, d, m, n, lam, p, k1def, k2b, 42)
    Me = build_elim_lattice(sigs, n, lam, k1def, k2b)
    re_ = reduce_matrix(Me)
    r = min((x for x in re_ if any(x)), key=norm)
    tot = sum(x * x for x in r)
    cong = sum(r[j] ** 2 for j in range(m - 1))
    k10 = r[m - 1] ** 2
    k2e = sum(r[m + j] ** 2 for j in range(m))
    kan = r[2 * m] ** 2
    nnz = sum(1 for x in r if x)
    print(f"{label:<20} {k1def:>3} {cong/tot:>7.3f} {k10/tot:>7.3f} {k2e/tot:>7.3f} "
          f"{kan/tot:>7.3f} {nnz:>5}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("P3: is the K1 wall information-theoretic?  counting bound vs observed")
print("-" * 78)
print("eff = K1*K2/n ; m_thresh = ceil(log n / log(1/eff)) signatures needed to")
print("pin d uniquely.  Observed = success at the m actually used.")
print()
print(f"{'curve':<20} {'K1':>3} {'K2':>4} {'eff':>7} {'m_thresh':>9} {'m_used':>7} "
      f"{'enough info?':>13} {'observed':>9}")
for label, curve, k1def, m in curves[1:]:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for k1 in [2, 4, 6, 8, 12, 16]:
        eff = k1 * k2b / n
        if eff >= 1.0:
            mt = float('inf')
        else:
            mt = math.ceil(math.log(n) / math.log(1.0 / eff))
        good = 0
        for s in SEEDS:
            d = random.Random(s * 7919 + 13).randint(1, n - 1)
            sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, s)
            if len(sigs) < m: continue
            Me = build_elim_lattice(sigs, n, lam, k1, k2b)
            re_ = reduce_matrix(Me)
            good += 1 if recover_d_elim(re_, m, n, sigs, lam, k1, k2b, d) is not None else 0
        print(f"{label:<20} {k1:>3} {k2b:>4} {eff:>7.3f} {mt:>9} {m:>7} "
              f"{str(m >= mt):>13} {str(good)+'/'+str(len(SEEDS)):>9}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
