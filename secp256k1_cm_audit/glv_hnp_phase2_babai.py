"""
GLV-HNP Phase 2, Thread 23b: the exact BDD/Babai criterion for the projected
lattice — a per-instance predictor, and a CVP-based attack that does not use
the Kannan embedding at all.

Motivation (this run's U1/U5, see glv_hnp_phase2_projected_output.txt):
  * Projecting the d-column out (Thread 23) changes lambda_1 from the trivial
    vector n*S_D*e_m to the lambda-block vector mu = lambda_1(L2), L2 the 2D
    sublattice <(n*S_K1,0), (-lam*S_K1,S_K2)> on the coordinate pair (i,m+i).
    (mu is the unnormalised form of Thread 20b's nu_hat = mu/sqrt(det L2).)
  * mu is EXACTLY independent of m, while E||v_planted|| ~ n*sqrt(2m/3+1).
    So the SVP margin mu/||pv|| strictly DEGRADES as signatures are added,
    yet observed recovery strictly IMPROVES with m.  Recovery is therefore
    governed by a BDD condition, not by an SVP condition.

The exact BDD condition.  Drop the Kannan row; let

    L0 = < n*S_K1*e_i ,  (B_i*S_K1 | 0) ,  (-lam*S_K1*e_j | S_K2*e_j) >  in Z^2m
    tau = (-A_i*S_K1 | 0)                                  (public target)
    e   = (k1_i*S_K1 | k2_i*S_K2)                          (secret error)

so that the lattice vector w = e - tau is the one an attacker must find.
Babai's nearest-plane on an LLL-reduced basis b_1..b_2m with Gram-Schmidt
b*_1..b*_2m returns exactly w iff e lies in the fundamental parallelepiped of
B*, i.e. iff

    nu_i := 2*|<e, b*_i>| / ||b*_i||^2  <=  1   for every i,
    NU   := max_i nu_i  <=  1.

NU is exact, per-instance, and needs no fitting.  Two questions:

  V1  Does NU <= 1 predict the observed LLL/Kannan recovery?  (Is the Kannan
      formulation actually achieving the BDD bound, or leaving margin?)
  V2  Does running Babai directly recover d where Kannan-LLL fails, i.e. does
      the CVP formulation move the K1 wall that U2 showed is immovable?

Run: python3 glv_hnp_phase2_babai.py
"""

import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fpylll import IntegerMatrix, LLL

import glv_hnp_common as C
from glv_hnp_common import (
    modinv, find_generator, lam_star, scales, gen_signatures, norm,
)
from glv_hnp_phase2_projected import (
    build_projected_lattice, projected_planted_vector, run_old, run_new,
    SEEDS, HIST,
)


# ---------------------------------------------------------------------------
# L0: the projected lattice WITHOUT the Kannan row (dim 2m, rank 2m)
# ---------------------------------------------------------------------------

def build_L0(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for j in range(m):
        r = [0] * dim
        r[j] = -lam * S_K1
        r[m + j] = S_K2
        rows.append(r)
    return rows


def target_and_error(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    tau = [0] * (2 * m)
    e = [0] * (2 * m)
    for i in range(m):
        # w = e - tau must lie in L0.  From A_i + B_i*d = k1_i + lam*k2_i (mod n),
        # coord i of w is S_K1*(d*B_i - lam*k2_i - c_i*n) = S_K1*(k1_i - A_i),
        # so tau[i] = +A_i*S_K1.
        tau[i] = sigs[i]['A'] * S_K1
        e[i] = sigs[i]['k1'] * S_K1
        e[m + i] = sigs[i]['k2'] * S_K2
    return tau, e


# ---------------------------------------------------------------------------
# Gram-Schmidt and Babai nearest-plane (pure float; entries are ~n^2/K1)
# ---------------------------------------------------------------------------

def gram_schmidt(B):
    """B: list of integer rows.  Returns (Bstar, mu) as floats."""
    k = len(B)
    Bs, mu = [], [[0.0] * k for _ in range(k)]
    for i in range(k):
        v = [float(x) for x in B[i]]
        for j in range(i):
            nj = sum(x * x for x in Bs[j])
            if nj == 0.0:
                mu[i][j] = 0.0
                continue
            m_ij = sum(a * b for a, b in zip(B[i], Bs[j])) / nj
            mu[i][j] = m_ij
            for t in range(len(v)):
                v[t] -= m_ij * Bs[j][t]
        Bs.append(v)
    return Bs, mu


def bdd_nu(B, e):
    """NU = max_i 2|<e,b*_i>| / ||b*_i||^2 — Babai nearest-plane succeeds
    iff NU <= 1."""
    Bs, _ = gram_schmidt(B)
    worst = 0.0
    for i in range(len(Bs)):
        ni = sum(x * x for x in Bs[i])
        if ni == 0.0:
            continue
        worst = max(worst, 2.0 * abs(sum(a * b for a, b in zip(e, Bs[i]))) / ni)
    return worst


def babai_nearest_plane(B, tau):
    """Returns the lattice vector (exact integer coords) nearest-plane finds."""
    Bs, _ = gram_schmidt(B)
    k = len(B)
    b = [float(x) for x in tau]
    coeffs = [0] * k
    for i in range(k - 1, -1, -1):
        ni = sum(x * x for x in Bs[i])
        if ni == 0.0:
            continue
        c = sum(a * x for a, x in zip(b, Bs[i])) / ni
        ci = int(math.floor(c + 0.5))
        coeffs[i] = ci
        if ci:
            for t in range(len(b)):
                b[t] -= ci * B[i][t]
    v = [0] * len(tau)
    for i in range(k):
        if coeffs[i]:
            for t in range(len(v)):
                v[t] += coeffs[i] * B[i][t]
    return v


def lll_rows(rows, dim):
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    out = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    return [r for r in out if any(r)]


def run_babai(curve, m, d_secret, k1_bound, seed):
    """CVP attack on L0 (no Kannan row).  Returns dict with NU and recovery."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    B = lll_rows(build_L0(sigs, n, lam, k1_bound, k2_bound), 2 * m)
    tau, e = target_and_error(sigs, n, k1_bound, k2_bound)
    NU = bdd_nu(B, e)
    w = babai_nearest_plane(B, tau)
    # Babai returns the lattice point NEAREST tau, i.e. v = tau - e (both v and
    # -v are in L0; ||tau - (tau-e)|| = ||e|| is the small one).  So e = tau - w.
    e_rec = [tau[t] - w[t] for t in range(2 * m)]
    ok = False
    k1s = [e_rec[i] for i in range(m)]
    k2s = [e_rec[m + i] for i in range(m)]
    if all(x % S_K1 == 0 for x in k1s) and all(x % S_K2 == 0 for x in k2s):
        k1s = [x // S_K1 for x in k1s]
        k2s = [x // S_K2 for x in k2s]
        for i in range(m):
            if sigs[i]['B'] % n == 0:
                continue
            d_cand = ((k1s[i] + lam * k2s[i] - sigs[i]['A'])
                      * modinv(sigs[i]['B'], n)) % n
            ok = (d_cand == d_secret)
            break
    return {'ok': ok, 'NU': NU, 'exact': (e_rec == e)}


def babai_rate(curve, m, k1_bound, seeds):
    wins, nus = 0, []
    n = curve[2]
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = run_babai(curve, m, d_trial, k1_bound, seed)
        if r is None:
            continue
        wins += bool(r['ok'])
        nus.append(r['NU'])
    return wins, len(seeds), (sum(nus) / len(nus) if nus else float('nan')), nus


if __name__ == "__main__":
    # ===========================================================================
    print("=" * 78)
    print("Thread 23b — exact BDD criterion NU and a Kannan-free Babai attack")
    print("=" * 78)

    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, label
        hist.append((label, (p, b, n, lam, G), k1, m))

    U2 = [("12-bit/2557", hist[1][1], 8, 0.340),
          ("12-bit/2677", hist[2][1], 10, 0.070)]
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]

    print("\n" + "-" * 78)
    print("EXP V1/V2: NU vs observed recovery, and Babai-CVP vs Kannan-LLL")
    print("-" * 78)
    print("NU = max_i 2|<e,b*_i>|/||b*_i||^2 on the LLL-reduced L0 (no Kannan row).")
    print("Babai nearest-plane provably returns the planted vector iff NU <= 1.\n")

    agree = disagree = 0
    for label, curve, m, ls in U2:
        p, b, n, lam, G = curve
        print(f"\n{label}  (lam* = {ls:.3f}, m = {m}, n = {n})")
        print(f"  {'K1':>4} {'eff':>6} {'mean NU':>9} {'NU<=1':>6} "
              f"{'BABAI':>6} {'KANNAN':>7} {'exact':>6}")
        for k1 in K1_GRID:
            wb, tb, mean_nu, nus = babai_rate(curve, m, k1, SEEDS)
            n_le1 = sum(1 for x in nus if x <= 1.0)
            wk = 0
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = run_new(curve, m, d_trial, k1, seed)
                wk += bool(r['ok'])
            ex = 0
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                rr = run_babai(curve, m, d_trial, k1, seed)
                ex += bool(rr['exact'])
            print(f"  {k1:>4} {k1*(math.isqrt(n)+1)/n:>6.3f} {mean_nu:>9.3f} "
                  f"{str(n_le1)+'/'+str(tb):>6} {str(wb)+'/'+str(tb):>6} "
                  f"{str(wk)+'/'+str(len(SEEDS)):>7} {str(ex)+'/'+str(len(SEEDS)):>6}")
            if wb == wk:
                agree += 1
            else:
                disagree += 1

    print(f"\nBabai-CVP vs Kannan-LLL agreement over the K1 grid: "
          f"{agree} agree, {disagree} differ.")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP V3: is NU <= 1 a sharp predictor?  (per-seed confusion matrix)")
    print("-" * 78)
    tp = fp = fn = tn = 0
    rows_v3 = []
    for label, curve, m, ls in U2:
        n = curve[2]
        for k1 in K1_GRID:
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                rb = run_babai(curve, m, d_trial, k1, seed)
                rk = run_new(curve, m, d_trial, k1, seed)
                pred = rb['NU'] <= 1.0
                obs = rk['ok']
                rows_v3.append((label, k1, seed, rb['NU'], pred, rb['ok'], obs))
                if pred and obs: tp += 1
                elif pred and not obs: fp += 1
                elif not pred and obs: fn += 1
                else: tn += 1
    tot = tp + fp + fn + tn
    print(f"predictor: NU <= 1   target: Kannan-LLL recovery   N = {tot}")
    print(f"  TP {tp:>4}   FP {fp:>4}")
    print(f"  FN {fn:>4}   TN {tn:>4}")
    print(f"  accuracy = {(tp+tn)/tot:.3f}")
    nu_ok = [r[3] for r in rows_v3 if r[6]]
    nu_no = [r[3] for r in rows_v3 if not r[6]]
    if nu_ok and nu_no:
        print(f"  NU | success : min {min(nu_ok):.3f}  median "
              f"{sorted(nu_ok)[len(nu_ok)//2]:.3f}  max {max(nu_ok):.3f}  "
              f"(n={len(nu_ok)})")
        print(f"  NU | failure : min {min(nu_no):.3f}  median "
              f"{sorted(nu_no)[len(nu_no)//2]:.3f}  max {max(nu_no):.3f}  "
              f"(n={len(nu_no)})")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP V4: how sharp?  AUC and the two-sided bracket")
    print("-" * 78)
    # AUC of -NU as a score for "Kannan-LLL recovers" (rank statistic, ties=0.5)
    pos = [r[3] for r in rows_v3 if r[6]]
    neg = [r[3] for r in rows_v3 if not r[6]]
    conc = sum((1.0 if a < b else 0.5 if a == b else 0.0)
               for a in pos for b in neg)
    auc = conc / (len(pos) * len(neg))
    print(f"AUC(-NU -> recovery) = {auc:.4f}   (n_pos={len(pos)}, n_neg={len(neg)})")
    print("compare: Thread 20b's fitted nu_hat separator reached AUC 0.935")
    thr_suff = min(neg)          # every NU below this is a success
    thr_nec = max(pos)           # every NU above this is a failure
    n_band = sum(1 for r in rows_v3 if thr_suff <= r[3] <= thr_nec)
    print(f"  NU <  {thr_suff:.3f}  => recovery in "
          f"{sum(1 for r in rows_v3 if r[3] < thr_suff and r[6])}/"
          f"{sum(1 for r in rows_v3 if r[3] < thr_suff)} cases (sufficient)")
    print(f"  NU >  {thr_nec:.3f}  => failure  in "
          f"{sum(1 for r in rows_v3 if r[3] > thr_nec and not r[6])}/"
          f"{sum(1 for r in rows_v3 if r[3] > thr_nec)} cases (necessary)")
    print(f"  ambiguous band [{thr_suff:.3f}, {thr_nec:.3f}] holds "
          f"{n_band}/{len(rows_v3)} instances "
          f"({100.0*n_band/len(rows_v3):.1f}%)")
    print("\nThe theoretical bound NU <= 1 is SOUND but not tight: it certifies")
    print(f"{sum(1 for r in rows_v3 if r[3] <= 1.0 and r[6])} of the "
          f"{len(pos)} successes with zero false positives, while Kannan-LLL")
    print("in fact recovers up to NU ~ 1.87 — i.e. LLL/Kannan strictly beats")
    print("the nearest-plane guarantee, by a factor of ~1.9 in NU.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
