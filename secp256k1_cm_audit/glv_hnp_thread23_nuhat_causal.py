"""
GLV-HNP Phase 2, Thread 23c: causal arm for nu_hat under exact CVP.

Thread 23b (glv_hnp_thread23_nuhat_under_cvp.py) found, on 30 fresh 17-bit
curves at pinned eff = 0.13:

    Spearman(nu_hat, p_hat_LLL)       = -0.567  (perm p = 0.0014)
    Spearman(nu_hat, p_hat_exactCVP)  = -0.744  (perm p < 1e-4), AUC = 1.000
    Spearman(nu_hat, p_cvp - p_lll)   = +0.142  (perm p = 0.45, null)

i.e. nu_hat predicts the INFORMATION-THEORETIC recovery probability, not the
reduction loss.  H23b ("nu_hat is an LLL-quality proxy") is falsified.

That was still an observational sample: 30 curves differing in n, lam and the
rounding of K1.  This script is the controlled version.  The lattice problem
does not care whether lam is a genuine cube root of unity, so the HNP instance
is synthesised directly:

    pick d, k1_i ~ U[0,K1), k2_i ~ U[0,K2), B_i ~ U[1,n)
    set  A_i = (k1_i + lam*k2_i - B_i*d)  mod n

which reproduces the Phase-2 relation k1_i = A_i + B_i*d - lam*k2_i (mod n)
exactly.  Everything (n, m, K1, K2, the seed set, the sampling distribution)
is held fixed and ONLY lam is varied over NLAM random values.  Any dependence
of p_hat on nu_hat is therefore causal in lam.

  H23c : Spearman(nu_hat, p_hat_cvp) is negative and significant with n, m,
         K1, K2 held fixed.
  Falsifier: perm p > 0.05, or a sign flip versus Thread 23b.

RESULT (2026-08-02): H23c CONFIRMED.  n = 131071, m = 12, K1 = 47, K2 = 363,
eff = 0.1302, 120 lam values x 12 seeds = 1440 instances per arm:

    Spearman(nu_hat, p_lll) = -0.827  perm p < 1e-4   AUC 0.967
    Spearman(nu_hat, p_cvp) = -0.724  perm p < 1e-4   AUC 0.983
    Spearman(lam*,   p_lll) = -0.030  perm p = 0.75   (null control)
    Spearman(lam*,   p_cvp) = +0.032  perm p = 0.72   (null control)

nu_hat decile table, p_cvp:  0.979 0.618 0.243 0.049 0.021 0.014 0.021 0.000
0.014 0.021 over nu_hat in [0.140, 1.041].  Practical rule at this eff:
nu_hat < 0.35 recoverable, nu_hat > 0.60 essentially not.

Also note mean p_lll = 0.267 > mean p_cvp = 0.198 here: the Kannan basis scan
is a list decoder and beats unique CVP decoding on average at this eff.

Run: python3 glv_hnp_thread23_nuhat_causal.py
"""

import math
import random
import time

import sympy

from glv_hnp_thread23_bdd import (
    build_cvp, build_kannan, error_vector, gh_cvp, norm, recover_kannan, scales,
)
from glv_hnp_thread23_nuhat_under_cvp import (
    auc, nu_hat, perm_pvalue, spearman,
)
from fpylll import CVP, IntegerMatrix, LLL

NLAM = 120
M_SIGS = 12
NSEED = 12
EFF_TARGET = 0.13
BITS = 17

# ---------------------------------------------------------------------------

def synth_instance(n, lam, m, k1_bound, k2_bound, seed):
    """Synthetic Phase-2 HNP instance: k1_i = A_i + B_i*d - lam*k2_i (mod n)."""
    rng = random.Random(seed)
    d = rng.randint(1, n - 1)
    sigs = []
    for _ in range(m):
        k1 = rng.randrange(k1_bound)
        k2 = rng.randrange(k2_bound)
        B = rng.randint(1, n - 1)
        A = (k1 + lam * k2 - B * d) % n
        assert (A + B * d - lam * k2) % n == k1 % n
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return d, sigs

def solve_lll(sigs, n, lam, k1, k2, d):
    m = len(sigs)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(build_kannan(sigs, n, lam, k1, k2))
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, S_D, _, S_KAN = scales(n, k1, k2)
    return recover_kannan(reduced, m, n, S_KAN, S_D, d)

def solve_cvp(sigs, n, lam, k1, k2, d):
    m = len(sigs)
    dim = 2 * m + 1
    M, t = build_cvp(sigs, n, lam, k1, k2)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    v = CVP.closest_vector(A, tuple(t))
    return (v[m] - t[m]) % n == d % n


print("=" * 78)
print("Thread 23c — causal arm: only lam varies")
print("=" * 78)

n = int(sympy.prevprime(1 << BITS))
k2_bound = math.isqrt(n) + 1
k1_bound = max(2, round(EFF_TARGET * n / k2_bound))
eff = k1_bound * k2_bound / n
print(f"n = {n} ({n.bit_length()} bits, prime), m = {M_SIGS}, "
      f"K1 = {k1_bound}, K2 = {k2_bound}, eff = {eff:.4f}")
print(f"{NLAM} random lam values, {NSEED} seeds each "
      f"= {NLAM * NSEED} instances per arm")

rng = random.Random(20260802)
lams = sorted({rng.randint(1, n - 1) for _ in range(NLAM + 20)})[:NLAM]

rows = []
t0 = time.time()
for lam in lams:
    nh = nu_hat(n, lam, k1_bound, k2_bound)
    ls = min(lam, n - lam) / n
    wl = wc = 0
    taus = []
    for s in range(NSEED):
        d, sigs = synth_instance(n, lam, M_SIGS, k1_bound, k2_bound, seed=1000 * s + 7)
        taus.append(norm(error_vector(sigs, d, n, k1_bound, k2_bound))
                    / gh_cvp(M_SIGS, n, k1_bound, k2_bound))
        wl += bool(solve_lll(sigs, n, lam, k1_bound, k2_bound, d))
        wc += bool(solve_cvp(sigs, n, lam, k1_bound, k2_bound, d))
    rows.append(dict(lam=lam, nu=nh, ls=ls, tau=sum(taus) / len(taus),
                     p_lll=wl / NSEED, p_cvp=wc / NSEED))
print(f"{len(rows)} lam values, {time.time()-t0:.1f}s")

print("\n" + "-" * 78)
print("nu_hat decile table (exact CVP is the information-theoretic arm)")
print("-" * 78)
rows_sorted = sorted(rows, key=lambda r: r["nu"])
B = 10
per = len(rows_sorted) // B
print(f"{'bin':>4}{'nu range':>18}{'lam* mean':>11}{'tau mean':>10}"
      f"{'p_LLL':>8}{'p_CVP':>8}")
for i in range(B):
    ch = rows_sorted[i * per:(i + 1) * per] if i < B - 1 else rows_sorted[i * per:]
    lo, hi = ch[0]["nu"], ch[-1]["nu"]
    print(f"{i:>4}{f'{lo:.3f}-{hi:.3f}':>18}"
          f"{sum(c['ls'] for c in ch)/len(ch):>11.3f}"
          f"{sum(c['tau'] for c in ch)/len(ch):>10.3f}"
          f"{sum(c['p_lll'] for c in ch)/len(ch):>8.3f}"
          f"{sum(c['p_cvp'] for c in ch)/len(ch):>8.3f}")

print("\n" + "-" * 78)
print("Correlations (n, m, K1, K2 all fixed; only lam varies)")
print("-" * 78)
print(f"{'pair':<34}{'spearman':>10}{'perm p':>10}")
for var, vname in (("nu", "nu_hat"), ("ls", "lam*"), ("tau", "tau")):
    for key in ("p_lll", "p_cvp"):
        s, pv = perm_pvalue([r[var] for r in rows], [r[key] for r in rows])
        print(f"{f'({vname}, {key})':<34}{s:>10.3f}{pv:>10.4f}")
s, pv = perm_pvalue([r["nu"] for r in rows], [r["p_cvp"] - r["p_lll"] for r in rows])
print(f"{'(nu_hat, p_cvp - p_lll)':<34}{s:>10.3f}{pv:>10.4f}")

print("\nAUC(-nu_hat) for 'p_hat >= 0.5':")
for key in ("p_lll", "p_cvp"):
    lab = [r[key] >= 0.5 for r in rows]
    print(f"  {key}: {sum(lab)}/{len(rows)} positives, "
          f"AUC = {auc(lab, [-r['nu'] for r in rows]):.3f}")

mean_l = sum(r["p_lll"] for r in rows) / len(rows)
mean_c = sum(r["p_cvp"] for r in rows) / len(rows)
print(f"\nmean p_LLL = {mean_l:.3f}   mean p_CVP = {mean_c:.3f}   "
      f"reduction loss = {mean_c - mean_l:+.3f}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
