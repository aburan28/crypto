#!/usr/bin/env python3
"""
Thread 23, part 3 — the lambda-block criterion mu/pv as a PREDICTOR of the K1 wall.

From glv_hnp_phase2_elim_probe.py:
  * eliminating d is recovery-equivalent to the original lattice (45/45 per seed),
  * the counting bound m >= m_thresh holds at EVERY observed failure, so the wall is
    not information-theoretic,
  * the shortest vector of the (d-free) lattice is the exact Gauss-reduced shortest
    vector mu of the 2-D "lambda block"  < (n*S_K1, 0), (-lam*S_K1, S_K2) >,
    which is a sublattice of both constructions and is untouched by eliminating d.

That suggests a concrete criterion.  Heuristically

    ||mu||  ~  sqrt(det)     = sqrt(n * S_K1 * S_K2) ~ n^{3/2} / sqrt(K1*K2)
    ||pv||  ~  n * sqrt(2m/3 + 1)

so   mu/pv ~ sqrt( n / (K1 * K2 * (2m/3 + 1)) ),  and the planted vector can only be
lambda_1 when

    K1 * K2 * (2m/3 + 1)  <  n.                                            (CRIT)

Note (CRIT) is strictly stronger than the usual counting condition K1*K2 < n: it pays an
extra factor (2m/3 + 1), i.e. adding signatures HURTS this criterion even though it helps
the counting bound.  That would explain why T4b saw no rescue from larger m.

This script tests the criterion two ways:
  W1  in-sample : mu/pv vs success across the K1 grid on the three historical curves
  W2  m-sweep   : does mu/pv track the (absence of) m-rescue at fixed K1?
  W3  out-of-sample: fresh j=0 GLV curves, does the mu/pv = 1 crossing locate the wall?

mu is computed EXACTLY (integer Gauss reduction); pv is the exact planted-vector norm.
Nothing here is fitted -- the threshold 1.0 is predicted a priori, not tuned.
"""

import math
import random

import sympy

from glv_hnp_phase2_eliminated import (
    find_generator, glv_roots, lam_star, norm, scales, gen_signatures,
    build_elim_lattice, planted_vector_elim, recover_d_elim, reduce_matrix,
    tonelli_shanks, ec_mul,
)

# exact 2-D Gauss reduction (inlined so this script does not re-execute the probe)
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


def nu_hat(n, lam, S_K1, S_K2):
    """nu_hat = lambda_1(L2)/sqrt(det L2), the det-normalised form used by the
    2026-07-29 run (commit e845207).  det L2 = n*S_K1*S_K2 is independent of lam,
    of m, and -- unlike ||v_planted|| -- introduces no spurious m-dependence."""
    mu, _ = lambda_block_mu(n, lam, S_K1, S_K2)
    return mu / math.sqrt(n * S_K1 * S_K2)


SEEDS = [42, 1234, 9999, 555, 31337]


def mu_over_pv(curve, m, k1, seeds=SEEDS):
    """Exact mean mu/||v_planted|| over seeds, plus the LLL success count."""
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, _sd, S_K2, _sk = scales(n, k1, k2b)
    mu, _w = lambda_block_mu(n, lam, S_K1, S_K2)
    nh = nu_hat(n, lam, S_K1, S_K2)
    ratios, good, tot = [], 0, 0
    for s in seeds:
        d = random.Random(s * 7919 + 13).randint(1, n - 1)
        sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, s)
        if len(sigs) < m: continue
        tot += 1
        pv = norm(planted_vector_elim(sigs, n, k1, k2b))
        ratios.append(mu / pv)
        M = build_elim_lattice(sigs, n, lam, k1, k2b)
        red = reduce_matrix(M)
        good += 1 if recover_d_elim(red, m, n, sigs, lam, k1, k2b, d) is not None else 0
    return (sum(ratios) / len(ratios) if ratios else float('nan')), good, tot, nh


def crit_predicted_k1(n, k2, m):
    """Largest K1 with K1*K2*(2m/3+1) < n."""
    return n / (k2 * (2.0 * m / 3.0 + 1.0))


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
print("Thread 23 part 3 — does mu/pv predict the Phase-2 K1 wall?")
print("=" * 78)

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("W1: in-sample — mu/pv vs success across the K1 grid")
print("-" * 78)
print("Prediction: success collapses where mu/pv crosses 1.0 (mu shorter than planted).")
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
for label, curve, _k1, m in curves:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    print(f"\n{label}  n={n}  lam*={lam_star(lam,n):.3f}  m={m}  K2={k2b}")
    print(f"  (CRIT) predicts wall at K1 ~ {crit_predicted_k1(n, k2b, m):.1f}")
    print(f"  {'K1':>4} {'mu/pv':>8} {'mu>pv?':>7} {'nu_hat':>8} {'success':>8}")
    for k1 in K1_GRID:
        r, g, tot, nh = mu_over_pv(curve, m, k1)
        print(f"  {k1:>4} {r:>8.3f} {str(r > 1.0):>7} {nh:>8.3f} "
              f"{str(g)+'/'+str(tot):>8}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("W2: m-sweep at fixed K1 — does mu/pv explain the absence of an m-rescue?")
print("-" * 78)
print("Counting bound says more signatures help.  (CRIT) says they HURT (factor 2m/3+1).")
label, curve, _k1, _m = curves[2]
p, b, n, lam, G = curve
k2b = math.isqrt(n) + 1
for k1 in [4, 8]:
    print(f"\n{label}  K1={k1}  K2={k2b}")
    print(f"  {'m':>4} {'m_thresh':>9} {'mu/pv':>8} {'nu_hat':>8} {'success':>8}")
    eff = k1 * k2b / n
    mt = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else float('inf')
    for mm in [6, 8, 12, 16, 24, 32]:
        r, g, tot, nh = mu_over_pv(curve, mm, k1)
        print(f"  {mm:>4} {mt:>9} {r:>8.3f} {nh:>8.3f} {str(g)+'/'+str(tot):>8}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("W3: out-of-sample — fresh j=0 GLV curves, threshold fixed a priori at 1.0")
print("-" * 78)


def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                bb = num // 2
                if bb >= 0 and a * a - a * bb + bb * bb == p:
                    return (a, bb)
    return None


def j0_traces(a, b):
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]


def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None


def find_fresh(lo, hi, want):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 8 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    rt = glv_roots(nc)
                    if rt is None: continue
                    cur = build_curve(p, nc)
                    if cur is None: continue
                    out.append((f"{p}/{nc}", (p, cur[1], nc, rt[0], cur[4])))
                    break
        p = int(sympy.nextprime(p))
    return out


fresh = find_fresh(2 ** 13, 2 ** 14, 6)
print(f"found {len(fresh)} fresh 13-bit j=0 GLV curves\n")
print(f"{'curve':<14} {'lam*':>6} {'m':>3} {'K1':>4} {'mu/pv':>8} {'nu_hat':>7} "
      f"{'pred':>6} {'success':>8} {'agree':>6}")
agree = tot_cells = 0
for lab, cur in fresh:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    for m in (8, 14):
        for k1 in (2, 4, 8, 16):
            r, g, tt, nh = mu_over_pv(cur, m, k1)
            if tt == 0: continue
            pred = r > 1.0                 # predicted: recoverable
            obs = g >= (tt + 1) // 2       # observed: majority of seeds recover
            tot_cells += 1
            agree += 1 if pred == obs else 0
            print(f"{lab:<14} {lam_star(lam,n):>6.3f} {m:>3} {k1:>4} {r:>8.3f} "
                  f"{nh:>7.3f} {str(pred):>6} {str(g)+'/'+str(tt):>8} "
                  f"{str(pred==obs):>6}")
print(f"\nout-of-sample agreement (threshold fixed at mu/pv = 1.0, not tuned): "
      f"{agree}/{tot_cells} = {agree/tot_cells:.3f}" if tot_cells else "no cells")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
