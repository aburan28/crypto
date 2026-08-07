"""
GLV-HNP Phase 2, Thread 23: make the planted vector lambda_1 by projecting
the d-column out of the lattice.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, Thread 20 / exp T5):
  The Phase-2 lattice of `glv_hnp_phase2_20bit.py:263` uses scales
  (S_K1, S_D, S_K2, S_KANNAN) = (n//K1, 1, n//K2, n) and carries an explicit
  d-column.  Because the basis contains the rows n*S_K1*e_i, the vector
      t = n * S_D * e_m          (the d-column axis, ||t|| = n*S_D)
  is ALWAYS in the lattice, while
      ||v_planted|| ~ n*sqrt((2m+1)/3 + 1)   >>   n*S_D.
  So the planted vector is never lambda_1, for any m >= 1 and any S_D
  (both vectors scale linearly in S_D).  Recovery was therefore a BDD/coset
  condition, not an SVP condition.

Thread 23 hypothesis (H23):
  The d-column is a *modelling error*, not just a nuisance.  d is uniform in
  [0,n) — it is an UNBOUNDED unknown — yet the construction gives it a
  coordinate of size d*S_D ~ n, i.e. it charges the target vector
  n^2*S_D^2/3 of squared norm for an unknown that carries no smallness.
  The correct treatment of an unbounded unknown is S_D -> 0: the d-row keeps
  its (B_i*S_K1) entries but contributes NO coordinate.  Then

    * the trivial vector t = n*S_D*e_m degenerates to 0 and disappears
      (the generating set becomes rank-deficient by exactly one; LLL emits
      one zero row),
    * the target loses the n^2/3 term:  ||v_planted|| ~ n*sqrt(2m/3 + 1),
    * the dimension drops from 2m+2 to 2m+1,
    * d is recovered algebraically afterwards from any single signature:
          d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  (mod n).

  Falsifier (stated in the 2026-07-29 next-step proposal):
    (i)  sv/pv must rise above 1 for the projected lattice, and
    (ii) the K1 wall of exp T4 must move OUTWARD on the lam*=0.07 curve
         (currently K1 ~ 4-6).
    If sv/pv rises but the wall does not move, the wall is
    information-theoretic and Phase 2 is at its ceiling.

Lattices compared (identical signatures, identical seeds, identical scales):

  OLD  dim 2m+2, cols [k1-block | d | k2-block | kannan]
  NEW  dim 2m+1, cols [k1-block |   | k2-block | kannan]     (S_D = 0)

Experiments:
  U1  sv/pv and shortest-vector anatomy, OLD vs NEW, historical curves
  U2  T4 replication: K1 grid on the two 12-bit curves, OLD vs NEW
  U3  17-bit sweep at eff = K1*K2/n in {0.05, 0.15, 0.25}, OLD vs NEW
  U4  Gaussian-heuristic prediction of the K1 wall vs the observed wall

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fpylll import IntegerMatrix, LLL, BKZ

import glv_hnp_common as C
from glv_hnp_common import (
    modinv, ec_mul, find_generator, glv_roots, lam_star,
    scales, gen_signatures, build_glv_lattice, planted_vector, norm,
    search_curves,
)

SEEDS = [42, 1234, 9999, 555, 31337]


# ---------------------------------------------------------------------------
# NEW: the projected (d-free) Phase-2 lattice
# ---------------------------------------------------------------------------

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Phase-2 lattice with the d-column removed (S_D = 0).

    Column layout (2m+1 columns):
        0     .. m-1   : k1-block, scaled by S_K1
        m     .. 2m-1  : k2-block, scaled by S_K2
        2m             : Kannan coordinate, scaled by S_KANNAN

    2m+2 generators in a rank-(2m+1) lattice: n*(d-row) reduces to 0 against
    the n*S_K1*e_i rows, which is precisely the trivial vector that used to
    be lambda_1.  LLL returns one zero row.
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                       # n*S_K1*e_i
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    r = [0] * dim                            # d-row: NO d coordinate
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for j in range(m):                       # lam / k2 rows
        r = [0] * dim
        r[j] = -lam * S_K1
        r[m + j] = S_K2
        rows.append(r)
    r = [0] * dim                            # Kannan row
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[dim - 1] = S_KANNAN
    rows.append(r)
    return rows


def projected_planted_vector(sigs, n, k1_bound, k2_bound):
    """(k1_i*S_K1, k2_i*S_K2, S_KANNAN) — the OLD planted vector with the
    d coordinate deleted."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v


def recover_d_projected(rows_reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read (k1_i, k2_i) off any row with |kannan| = S_KANNAN, then solve for d.

    Returns (d_candidate or None, verified_all_congruences: bool).
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    for row in rows_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN:
            continue
        sgn = 1 if last > 0 else -1
        k1s = [sgn * row[i] // S_K1 for i in range(m)]
        k2s = [sgn * row[m + i] // S_K2 for i in range(m)]
        # d from the first usable signature
        for i in range(m):
            B = sigs[i]['B']
            if B % n == 0:
                continue
            d_cand = (k1s[i] + lam * k2s[i] - sigs[i]['A']) * modinv(B, n) % n
            if d_cand == 0:
                continue
            # attacker-side self-check: every congruence must give a small k1
            ok_all = True
            for j in range(m):
                kf = (sigs[j]['A'] + sigs[j]['B'] * d_cand) % n
                if (kf - lam * k2s[j]) % n != k1s[j] % n:
                    ok_all = False
                    break
                if not (0 <= k1s[j] < k1_bound):
                    ok_all = False
                    break
            if d_cand == d_secret:
                return d_cand, ok_all
            break
    return None, False


# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------

def _nz_norms(rows):
    ns = [norm(r) for r in rows]
    return [x for x in ns if x > 0]


def run_old(curve, m, d_secret, k1_bound, seed, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    ok = C.recover_d(red, m, n, S_KAN, d_secret) is not None
    pn = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
    sn = min(_nz_norms(red))
    sv = min(red, key=lambda r: norm(r) if norm(r) > 0 else float('inf'))
    return {'ok': ok, 'pn': pn, 'sn': sn, 'sv': sv, 'dim': dim, 'sigs': sigs}


def run_new(curve, m, d_secret, k1_bound, seed, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    d_cand, ok_all = recover_d_projected(red, sigs, n, lam, k1_bound,
                                         k2_bound, d_secret)
    pn = norm(projected_planted_vector(sigs, n, k1_bound, k2_bound))
    nzn = _nz_norms(red)
    sn = min(nzn)
    sv = min(red, key=lambda r: norm(r) if norm(r) > 0 else float('inf'))
    return {'ok': d_cand is not None, 'verified': ok_all, 'pn': pn, 'sn': sn,
            'sv': sv, 'dim': dim, 'nzero': len(red) - len(nzn), 'sigs': sigs}


def rate(runner, curve, m, k1_bound, seeds, **kw):
    """(wins, total, mean sv/pv)."""
    wins, ratios = 0, []
    for seed in seeds:
        n = curve[2]
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = runner(curve, m, d_trial, k1_bound, seed, **kw)
        if r is None:
            continue
        wins += bool(r['ok'])
        ratios.append(r['sn'] / r['pn'])
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))


def gh_lambda1(dim, logdet):
    """Gaussian heuristic lambda_1 = sqrt(dim/(2*pi*e)) * det^(1/dim), in logs."""
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(logdet / dim)


HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]


if __name__ == "__main__":
    # ===========================================================================
    print("=" * 78)
    print("Thread 23 — project out the d-column so the planted vector is lambda_1")
    print("=" * 78)

    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, label
        assert (lam * lam + lam + 1) % n == 0, label
        hist.append((label, (p, b, n, lam, G), k1, m))


    # ===========================================================================
    # U1 — shortest-vector anatomy, OLD vs NEW
    # ===========================================================================
    print("\n" + "-" * 78)
    print("EXP U1: is the planted vector lambda_1 after projection?")
    print("-" * 78)
    print("sv/pv = ||shortest nonzero reduced row|| / ||planted vector||.")
    print("sv/pv >= 1 means the planted vector IS a shortest vector.\n")

    print(f"{'curve':<18} {'K1':>3} {'lam*':>6} | {'OLD dim':>7} {'sv/pv':>6} {'ok':>3} "
          f"| {'NEW dim':>7} {'sv/pv':>6} {'ok':>3} {'0-rows':>6} {'self-chk':>8}")
    u1_rows = []
    for label, curve, k1, m in hist:
        p, b, n, lam, G = curve
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        ro = run_old(curve, m, d0, k1, 42)
        rn = run_new(curve, m, d0, k1, 42)
        print(f"{label:<18} {k1:>3} {lam_star(lam, n):>6.4f} | "
              f"{ro['dim']:>7} {ro['sn']/ro['pn']:>6.3f} {str(ro['ok']):>3} | "
              f"{rn['dim']:>7} {rn['sn']/rn['pn']:>6.3f} {str(rn['ok']):>3} "
              f"{rn['nzero']:>6} {str(rn['verified']):>8}")
        u1_rows.append((label, ro, rn, m, n, k1))

    print("\nOLD shortest vector energy split (k1-block / d / k2-block / kannan):")
    for label, ro, rn, m, n, k1 in u1_rows:
        sv = ro['sv']
        e = [sum(x * x for x in sv[0:m]), sv[m] ** 2,
             sum(x * x for x in sv[m + 1:2 * m + 1]), sv[2 * m + 1] ** 2]
        tot = sum(e) or 1
        print(f"  {label:<18} " + " ".join(f"{x/tot:6.3f}" for x in e) +
              f"   |sv[m]|/n = {abs(sv[m])/n:.4f}")

    print("\nNEW shortest vector energy split (k1-block / k2-block / kannan):")
    for label, ro, rn, m, n, k1 in u1_rows:
        sv = rn['sv']
        e = [sum(x * x for x in sv[0:m]), sum(x * x for x in sv[m:2 * m]),
             sv[2 * m] ** 2]
        tot = sum(e) or 1
        S_K1, _, S_K2, S_KAN = scales(n, k1, math.isqrt(n) + 1)
        print(f"  {label:<18} " + " ".join(f"{x/tot:6.3f}" for x in e) +
              f"   |kannan|/S_KANNAN = {abs(sv[2*m])/S_KAN:.4f}")


    # ===========================================================================
    # U2 — T4 replication: does the K1 wall move?
    # ===========================================================================
    print("\n" + "-" * 78)
    print("EXP U2: K1 grid (exp T4 replication), OLD vs NEW, 5 seeds each")
    print("-" * 78)
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
    U2 = [("12-bit/2557", hist[1][1], 8, 0.340),
          ("12-bit/2677", hist[2][1], 10, 0.070)]
    u2_data = {}
    for label, curve, m, ls in U2:
        print(f"\n{label}  (lam* = {ls:.3f}, m = {m}, n = {curve[2]}, "
              f"K2 = {math.isqrt(curve[2])+1})")
        hdr = f"{'lattice':<8}" + "".join(f"{k:>6}" for k in K1_GRID)
        print(hdr)
        for name, runner in (("OLD", run_old), ("NEW", run_new)):
            cells, wall = [], None
            for k1 in K1_GRID:
                w, t, _ = rate(runner, curve, m, k1, SEEDS)
                cells.append(f"{w}/{t}")
                if wall is None and w < t:
                    wall = k1
            u2_data[(label, name)] = (cells, wall)
            print(f"{name:<8}" + "".join(f"{c:>6}" for c in cells))
        for name in ("OLD", "NEW"):
            print(f"   {name} first-failure K1 = {u2_data[(label, name)][1]}")


    # ===========================================================================
    # U3 — 17-bit sweep at three bias strengths, OLD vs NEW
    # ===========================================================================
    print("\n" + "-" * 78)
    print("EXP U3: 17-bit sweep, eff = K1*K2/n in {0.05, 0.15, 0.25}, OLD vs NEW")
    print("-" * 78)
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"found {len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")

    M_SIGS = 12
    print(f"\nm = {M_SIGS} signatures, {len(SEEDS)} seeds per curve.\n")
    print(f"{'eff':>5} | {'OLD 5/5':>8} {'OLD any':>8} {'OLD wins':>9} "
          f"| {'NEW 5/5':>8} {'NEW any':>8} {'NEW wins':>9}")
    u3_rows = []
    for eff in (0.05, 0.15, 0.25):
        tallies = {}
        per_curve = {'OLD': [], 'NEW': []}
        for name, runner in (("OLD", run_old), ("NEW", run_new)):
            full = part = tot_wins = 0
            for (p, b, n, lam, G) in curves17:
                k2b = math.isqrt(n) + 1
                k1b = max(2, int(eff * n / k2b))
                w, t, _ = rate(runner, (p, b, n, lam, G), M_SIGS, k1b, SEEDS)
                per_curve[name].append((n, lam_star(lam, n), k1b, w, t))
                full += (w == t)
                part += (w > 0)
                tot_wins += w
            tallies[name] = (full, part, tot_wins, len(curves17) * len(SEEDS))
        o, nw = tallies['OLD'], tallies['NEW']
        print(f"{eff:>5.2f} | {str(o[0])+'/'+str(len(curves17)):>8} "
              f"{str(o[1])+'/'+str(len(curves17)):>8} {str(o[2])+'/'+str(o[3]):>9} "
              f"| {str(nw[0])+'/'+str(len(curves17)):>8} "
              f"{str(nw[1])+'/'+str(len(curves17)):>8} {str(nw[2])+'/'+str(nw[3]):>9}")
        u3_rows.append((eff, tallies, per_curve))

    print("\nPer-curve detail at eff = 0.15 (the discriminating regime):")
    eff15 = [r for r in u3_rows if r[0] == 0.15][0]
    print(f"{'n':>8} {'lam*':>7} {'K1':>6} {'OLD':>6} {'NEW':>6}")
    for (nA, lsA, k1A, wA, tA), (nB, lsB, k1B, wB, tB) in zip(
            eff15[2]['OLD'], eff15[2]['NEW']):
        print(f"{nA:>8} {lsA:>7.4f} {k1A:>6} {str(wA)+'/'+str(tA):>6} "
              f"{str(wB)+'/'+str(tB):>6}")


    # ===========================================================================
    # U4 — Gaussian heuristic: where SHOULD the wall be?
    # ===========================================================================
    print("\n" + "-" * 78)
    print("EXP U4: Gaussian-heuristic lambda_1 vs ||v_planted||, NEW lattice")
    print("-" * 78)
    print("det(NEW) = (n*S_K1)^m * S_K2^m * S_KANNAN / n   [rank 2m+1]")
    print("E||pv||  = sqrt(m*(K1*S_K1)^2/3 + m*(K2*S_K2)^2/3 + S_KANNAN^2)")
    print("Predicted wall: smallest K1 with GH(lambda_1) < E||pv||.\n")
    for label, curve, m, ls in U2:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        print(f"{label} (m={m}):")
        print(f"  {'K1':>4} {'eff':>6} {'GH l1':>12} {'E||pv||':>12} {'GH/pv':>7} "
              f"{'pred':>5} {'obs OLD':>8} {'obs NEW':>8}")
        obs_o = dict(zip(K1_GRID, u2_data[(label, 'OLD')][0]))
        obs_n = dict(zip(K1_GRID, u2_data[(label, 'NEW')][0]))
        for k1 in K1_GRID:
            S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
            dim = 2 * m + 1
            logdet = (m * math.log(n * S_K1) + m * math.log(S_K2)
                      + math.log(S_KAN) - math.log(n))
            gh = gh_lambda1(dim, logdet)
            epv = math.sqrt(m * (k1 * S_K1) ** 2 / 3.0
                            + m * (k2b * S_K2) ** 2 / 3.0 + S_KAN ** 2)
            print(f"  {k1:>4} {k1*k2b/n:>6.3f} {gh:>12.4g} {epv:>12.4g} "
                  f"{gh/epv:>7.3f} {('OK' if gh > epv else 'fail'):>5} "
                  f"{obs_o[k1]:>8} {obs_n[k1]:>8}")


    # ===========================================================================
    # U5 — does MORE data hurt?  (the mu-vs-sqrt(m) prediction)
    # ===========================================================================
    print("\n" + "-" * 78)
    print("EXP U5: mu is independent of m, ||v_planted|| grows as sqrt(m)")
    print("-" * 78)
    print("If lambda_1 of the projected lattice is the lambda-block vector mu")
    print("(supported on one coordinate pair (i, m+i)), then mu does NOT depend")
    print("on m at all, while E||pv|| ~ n*sqrt(2m/3+1).  Prediction: sv/pv must")
    print("DECREASE as 1/sqrt(m), i.e. extra signatures strictly worsen the SVP")
    print("condition even though they add information.  T4b (2026-07-29) observed")
    print("exactly this non-monotonicity in the OLD lattice (m=8/12/16/24/32 ->")
    print("0,0,1,0,1 of 5) without an explanation.\n")
    M_GRID = [4, 6, 8, 12, 16, 24, 32]
    for label, curve, _m0, ls in U2:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        for k1 in (8, 4):
            mu, _w = C.lambda_block_mu(n, lam, *(lambda s: (s[0], s[2]))(
                scales(n, k1, k2b)))
            print(f"{label}  K1={k1}  lam*={ls:.3f}  exact mu = {mu:.4g} "
                  f"(independent of m)")
            print(f"  {'m':>4} {'E||pv||':>11} {'mu/E||pv||':>11} "
                  f"{'obs sv/pv':>10} {'OLD':>6} {'NEW':>6}")
            for m in M_GRID:
                S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
                epv = math.sqrt(m * (k1 * S_K1) ** 2 / 3.0
                                + m * (k2b * S_K2) ** 2 / 3.0 + S_KAN ** 2)
                wo, to, _ = rate(run_old, curve, m, k1, SEEDS)
                wn, tn, rn_ratio = rate(run_new, curve, m, k1, SEEDS)
                print(f"  {m:>4} {epv:>11.4g} {mu/epv:>11.3f} {rn_ratio:>10.3f} "
                      f"{str(wo)+'/'+str(to):>6} {str(wn)+'/'+str(tn):>6}")
            print()

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
