"""
GLV-HNP Phase 2, Thread 24: derive nu_hat from the Gram-Schmidt profile of L0.

Pre-registered by the 2026-08-07 (Thread 23) log entry:

  H24: the argmax of  nu_i = 2|<e, b*_i>| / ||b*_i||^2  concentrates at the
       tail indices i ~ 2m, and  ||b*_{2m}|| ~ det(L2)/mu = lambda_2(L2)
       up to a constant.
  Falsifier: if the argmax is spread across i rather than concentrated in the
       tail, the GS-profile explanation fails and nu_hat stays empirical.

Why it matters.  Thread 20c (2026-07-29) found that
    nu_hat = lambda_1(L2)/sqrt(det L2),   L2 = <(n*S_K1,0), (-lam*S_K1,S_K2)>
separates viable from non-viable Phase-2 instances with AUC 0.935, with a
*negative* sign (skew L2 => easier), and could only motivate the sign by hand
via lambda_1*lambda_2 ~ det.  Thread 23 replaced it with the exact, unfitted
BDD certificate NU = max_i nu_i (AUC 0.978, zero false positives at NU <= 1).
If H24 holds, NU and nu_hat are the same quantity and the sign is a theorem,
not a fit.  Chaining the standard estimates,

    ||e||       ~ n*sqrt(2m/3)                     (k1*S_K1 ~ k2*S_K2 ~ n)
    lambda_2(L2) = det(L2)/lambda_1(L2) = n^3/(K1*K2*mu)     (2D, up to 2/sqrt 3)
    <e, b*_i>   ~ ||e||*||b*_i||/sqrt(dim)         (generic direction)

gives the closed form tested in W3:

    NU  ~  C * mu * K1*K2 / n^2  =  C' * nu_hat * sqrt(eff),   eff = K1*K2/n.

That is: nu_hat is the *size-free* part of NU and sqrt(eff) is the part that
moves with the bias strength.  Both signs then follow from lambda_1*lambda_2.

Numerics.  All Gram-Schmidt here is EXACT (Fraction over the integer Gram
matrix), not the float GS of glv_hnp_phase2_babai.py.  Entries of L0 are
~n^2/K1, so squared norms reach ~2^68 at 17 bits and f64 GS is not obviously
safe on dim 24; W0 quantifies the float error against the exact values.

Run: python3 glv_hnp_phase2_gsprofile.py
"""

import math
import os
import random
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import glv_hnp_common as C
from glv_hnp_common import (
    find_generator, lam_star, scales, gen_signatures, search_curves,
)
from glv_hnp_phase2_projected import SEEDS, HIST, run_new
from glv_hnp_phase2_babai import (
    build_L0, target_and_error, lll_rows, gram_schmidt, bdd_nu,
)


# ---------------------------------------------------------------------------
# The 2D lambda-block sublattice: BOTH successive minima, exactly.
# ---------------------------------------------------------------------------

def gauss_reduce_2d_full(u, v):
    """Lagrange/Gauss reduction.  Returns (u, v) with |u| = lambda_1 and
    |v| = lambda_2 (exact for rank 2)."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]

    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]

    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        q = ((2 * num + den) // (2 * den) if num >= 0
             else -((-2 * num + den) // (2 * den)))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u, v


def l2_minima(n, lam, S_K1, S_K2):
    """(lambda_1, lambda_2, det, nu_hat) of L2 = <(n*S_K1,0),(-lam*S_K1,S_K2)>."""
    u, v = gauss_reduce_2d_full((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    l1 = math.sqrt(u[0] ** 2 + u[1] ** 2)
    l2 = math.sqrt(v[0] ** 2 + v[1] ** 2)
    det = float(n) * S_K1 * S_K2
    return l1, l2, det, l1 / math.sqrt(det)


# ---------------------------------------------------------------------------
# Exact Gram-Schmidt over the integer Gram matrix
# ---------------------------------------------------------------------------

def exact_gs(B, e):
    """B: integer rows (LLL-reduced).  e: integer vector.

    Returns (Bst, ge_star, nus) where
        Bst[i]     = ||b*_i||^2                  (Fraction, exact)
        ge_star[i] = <e, b*_i>                   (Fraction, exact)
        nus[i]     = 2|<e,b*_i>| / ||b*_i||^2    (float)

    Uses the Cholesky recursion on the Gram matrix, so the vector dimension
    never enters the O(k^3) inner loop.
    """
    k = len(B)
    G = [[sum(a * b for a, b in zip(B[i], B[j])) for j in range(k)]
         for i in range(k)]
    ge = [sum(a * b for a, b in zip(e, B[i])) for i in range(k)]

    mu = [[Fraction(0)] * k for _ in range(k)]
    Bst = [Fraction(0)] * k
    for i in range(k):
        for j in range(i):
            s = Fraction(G[i][j])
            for t in range(j):
                s -= mu[i][t] * mu[j][t] * Bst[t]
            mu[i][j] = s / Bst[j] if Bst[j] != 0 else Fraction(0)
        s = Fraction(G[i][i])
        for t in range(i):
            s -= mu[i][t] * mu[i][t] * Bst[t]
        Bst[i] = s

    ge_star = [Fraction(0)] * k
    for i in range(k):
        s = Fraction(ge[i])
        for j in range(i):
            s -= mu[i][j] * ge_star[j]
        ge_star[i] = s

    nus = [(2.0 * abs(float(ge_star[i])) / float(Bst[i])) if Bst[i] != 0 else 0.0
           for i in range(k)]
    return Bst, ge_star, nus


def instance(curve, m, d_secret, k1_bound, seed, exact=True):
    """Build one L0 instance and return its full GS/nu anatomy."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    B = lll_rows(build_L0(sigs, n, lam, k1_bound, k2_bound), 2 * m)
    tau, e = target_and_error(sigs, n, k1_bound, k2_bound)
    if exact:
        Bst, _gs, nus = exact_gs(B, e)
        prof = [math.sqrt(float(x)) if x > 0 else 0.0 for x in Bst]
    else:
        Bs, _ = gram_schmidt(B)
        prof = [math.sqrt(sum(x * x for x in r)) for r in Bs]
        nus = []
        for i, r in enumerate(Bs):
            ni = sum(x * x for x in r)
            nus.append(2.0 * abs(sum(a * b for a, b in zip(e, r))) / ni
                       if ni else 0.0)
    NU = max(nus)
    arg = max(range(len(nus)), key=lambda i: nus[i])
    l1, l2, det2, nuhat = l2_minima(n, lam, S_K1, S_K2)
    return {
        'k': len(B), 'prof': prof, 'nus': nus, 'NU': NU, 'argmax': arg,
        'enorm': math.sqrt(sum(x * x for x in e)),
        'mu': l1, 'l2': l2, 'det2': det2, 'nuhat': nuhat,
        'S_K1': S_K1, 'S_K2': S_K2, 'K2': k2_bound,
    }


def auc(pos, neg):
    """AUC of the score -x for separating pos (recovery) from neg."""
    if not pos or not neg:
        return float('nan')
    conc = sum((1.0 if a < b else 0.5 if a == b else 0.0)
               for a in pos for b in neg)
    return conc / (len(pos) * len(neg))


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for t in range(i, j + 1):
                r[order[t]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float('nan')


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 24 — derive nu_hat from the GS profile of L0 (H24)")
    print("=" * 78)

    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, label
        hist.append((label, (p, b, n, lam, G), k1, m))

    U2 = [("12-bit/2557", hist[1][1], 8, 0.340),
          ("12-bit/2677", hist[2][1], 10, 0.070)]
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]

    # -----------------------------------------------------------------------
    # Collect the grid once; every experiment below reads this table.
    # -----------------------------------------------------------------------
    t0 = time.time()
    rows = []
    for label, curve, m, ls in U2:
        n = curve[2]
        for k1 in K1_GRID:
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance(curve, m, d_trial, k1, seed, exact=True)
                if r is None:
                    continue
                rk = run_new(curve, m, d_trial, k1, seed)
                r.update({'label': label, 'm': m, 'n': n, 'K1': k1,
                          'seed': seed, 'ok': bool(rk['ok']),
                          'eff': k1 * (math.isqrt(n) + 1) / n})
                rows.append(r)
    print(f"\ncollected {len(rows)} instances (exact GS) in {time.time()-t0:.1f}s")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W0: is the float GS of Thread 23b safe?  (exact vs f64)")
    print("-" * 78)
    worst_rel = 0.0
    worst_row = None
    for label, curve, m, ls in U2:
        n = curve[2]
        for k1 in K1_GRID:
            d_trial = random.Random(42 + 7777).randint(1, n - 1)
            re_ = instance(curve, m, d_trial, k1, 42, exact=True)
            rf = instance(curve, m, d_trial, k1, 42, exact=False)
            rel = abs(rf['NU'] - re_['NU']) / re_['NU']
            if rel > worst_rel:
                worst_rel, worst_row = rel, (label, k1, re_['NU'], rf['NU'])
    print(f"max |NU_float - NU_exact| / NU_exact over the 22-cell grid "
          f"(seed 42) = {worst_rel:.3e}")
    if worst_row:
        print(f"  worst cell: {worst_row[0]} K1={worst_row[1]}  "
              f"exact {worst_row[2]:.6f}  float {worst_row[3]:.6f}")
    print("(dim <= 20 here; W4 repeats the check at dim 24 / 17 bits.)")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W1: where is the argmax of nu_i?  (H24, clause 1)")
    print("-" * 78)
    print("index reported as the 1-based position from the TAIL: 1 = b*_{2m}.\n")
    from collections import Counter
    for label, curve, m, ls in U2:
        sub = [r for r in rows if r['label'] == label]
        tail = Counter(r['k'] - r['argmax'] for r in sub)
        k = sub[0]['k']
        print(f"{label}  (m={m}, dim={k}, N={len(sub)})")
        print("  tail-position of argmax : " +
              "  ".join(f"{pos}:{cnt}" for pos, cnt in sorted(tail.items())))
        frac = sum(r['argmax'] / (r['k'] - 1) for r in sub) / len(sub)
        in_last_quarter = sum(1 for r in sub if r['argmax'] >= 0.75 * (r['k'] - 1))
        print(f"  mean argmax/(dim-1) = {frac:.3f}   "
              f"argmax in last quarter: {in_last_quarter}/{len(sub)}")
    allsub = rows
    print(f"\nOVERALL: argmax == last index in "
          f"{sum(1 for r in allsub if r['argmax'] == r['k']-1)}/{len(allsub)}; "
          f"in last two indices in "
          f"{sum(1 for r in allsub if r['argmax'] >= r['k']-2)}/{len(allsub)}")
    print("(uniform-null expectation for 'last index' = "
          f"{100.0/rows[0]['k']:.1f}%)")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W1b: the GS profile itself, log2 ||b*_i|| (seed 42)")
    print("-" * 78)
    for label, curve, m, ls in U2:
        n = curve[2]
        for k1 in (4, 8, 16, 32):
            d_trial = random.Random(42 + 7777).randint(1, n - 1)
            r = instance(curve, m, d_trial, k1, 42, exact=True)
            if r is None:
                continue
            prof = " ".join(f"{math.log2(x):5.1f}" if x > 0 else "  -inf"
                            for x in r['prof'])
            print(f"{label} K1={k1:<3} log2||b*_i||: {prof}")
            print(f"{'':>{len(label)+8}} log2 lambda_1(L2)={math.log2(r['mu']):.1f}"
                  f"  log2 lambda_2(L2)={math.log2(r['l2']):.1f}"
                  f"  argmax i={r['argmax']} (tail pos {r['k']-r['argmax']})"
                  f"  NU={r['NU']:.3f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W2: is ||b*_last|| ~ lambda_2(L2)?  (H24, clause 2)")
    print("-" * 78)
    print(f"{'curve':<14} {'K1':>4} {'lam_1(L2)':>11} {'lam_2(L2)':>11} "
          f"{'||b*_last||':>12} {'ratio':>7} {'l1*l2/det':>10}")
    ratios = []
    for label, curve, m, ls in U2:
        for k1 in K1_GRID:
            sub = [r for r in rows if r['label'] == label and r['K1'] == k1]
            if not sub:
                continue
            bl = sum(r['prof'][-1] for r in sub) / len(sub)
            r0 = sub[0]
            rat = bl / r0['l2']
            ratios.append(rat)
            print(f"{label:<14} {k1:>4} {r0['mu']:>11.4g} {r0['l2']:>11.4g} "
                  f"{bl:>12.4g} {rat:>7.3f} "
                  f"{r0['mu']*r0['l2']/r0['det2']:>10.4f}")
    print(f"\nratio ||b*_last|| / lambda_2(L2): mean {sum(ratios)/len(ratios):.3f}  "
          f"min {min(ratios):.3f}  max {max(ratios):.3f}  "
          f"spread {max(ratios)/min(ratios):.2f}x")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W3: the closed form  NU ~ C * nu_hat * sqrt(eff)")
    print("-" * 78)
    print("From ||e|| ~ n*sqrt(2m/3), lambda_2 = det(L2)/mu, <e,b*> ~ "
          "||e|| ||b*||/sqrt(dim):")
    print("   NU ~ (2/sqrt 3) * mu * K1*K2 / n^2  =  (2/sqrt 3) * nu_hat "
          "* sqrt(K1*K2/n)\n")
    pred_raw, obs = [], []
    for r in rows:
        eff = r['K1'] * r['K2'] / r['n']
        r['pred'] = r['nuhat'] * math.sqrt(eff)
        pred_raw.append(r['pred'])
        obs.append(r['NU'])
    Cfit = (sum(o / p for o, p in zip(obs, pred_raw) if p > 0)
            / sum(1 for p in pred_raw if p > 0))
    resid = [o / (Cfit * p) for o, p in zip(obs, pred_raw) if p > 0]
    print(f"fitted C = {Cfit:.4f}   (2/sqrt 3 = {2/math.sqrt(3):.4f})")
    print(f"NU / (C * nu_hat * sqrt(eff)):  mean {sum(resid)/len(resid):.3f}  "
          f"min {min(resid):.3f}  max {max(resid):.3f}  "
          f"spread {max(resid)/min(resid):.2f}x")
    print(f"Spearman(pred, NU) = {spearman(pred_raw, obs):.4f}   N={len(obs)}")

    print("\nper-cell means:")
    print(f"{'curve':<14} {'K1':>4} {'eff':>6} {'nu_hat':>8} {'pred':>9} "
          f"{'C*pred':>9} {'NU obs':>8} {'obs/pred':>9} {'rec':>5}")
    for label, curve, m, ls in U2:
        for k1 in K1_GRID:
            sub = [r for r in rows if r['label'] == label and r['K1'] == k1]
            if not sub:
                continue
            nu_m = sum(r['NU'] for r in sub) / len(sub)
            pr = sub[0]['pred']
            wins = sum(1 for r in sub if r['ok'])
            print(f"{label:<14} {k1:>4} {sub[0]['eff']:>6.3f} "
                  f"{sub[0]['nuhat']:>8.4f} {pr:>9.4f} {Cfit*pr:>9.4f} "
                  f"{nu_m:>8.4f} {nu_m/pr:>9.4f} "
                  f"{str(wins)+'/'+str(len(sub)):>5}")

    pos = [r['pred'] for r in rows if r['ok']]
    neg = [r['pred'] for r in rows if not r['ok']]
    posNU = [r['NU'] for r in rows if r['ok']]
    negNU = [r['NU'] for r in rows if not r['ok']]
    print(f"\nAUC(-pred -> recovery)  = {auc(pos, neg):.4f}   "
          f"[closed form, NO per-instance lattice work]")
    print(f"AUC(-NU   -> recovery)  = {auc(posNU, negNU):.4f}   "
          f"[exact BDD certificate, Thread 23b]")
    print("Thread 20b's fitted nu_hat separator reached AUC 0.935.")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W4: is the NU bracket n-independent?  (17-bit re-measure)")
    print("-" * 78)
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"found {len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")
    M17 = 12
    rows17 = []
    t0 = time.time()
    for eff in (0.05, 0.15, 0.25):
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), M17, d_trial, k1b, seed,
                             exact=True)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), M17, d_trial, k1b, seed)
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff})
                rows17.append(r)
    print(f"{len(rows17)} 17-bit instances (dim {rows17[0]['k']}, exact GS) "
          f"in {time.time()-t0:.1f}s")

    p17 = [r['NU'] for r in rows17 if r['ok']]
    n17 = [r['NU'] for r in rows17 if not r['ok']]
    print(f"\nNU | success : min {min(p17):.3f}  median "
          f"{sorted(p17)[len(p17)//2]:.3f}  max {max(p17):.3f}  (n={len(p17)})")
    print(f"NU | failure : min {min(n17):.3f}  median "
          f"{sorted(n17)[len(n17)//2]:.3f}  max {max(n17):.3f}  (n={len(n17)})")
    tp = sum(1 for r in rows17 if r['NU'] <= 1.0 and r['ok'])
    fp = sum(1 for r in rows17 if r['NU'] <= 1.0 and not r['ok'])
    print(f"NU <= 1 certificate:  TP {tp}  FP {fp}   "
          f"(FP must be 0 if nearest-plane theory holds)")
    print(f"bracket at 17 bits : sufficient NU < {min(n17):.3f} , "
          f"necessary NU > {max(p17):.3f}")
    print("compare 12 bits (Thread 23b): sufficient NU < 1.188 , "
          "necessary NU > 1.869")
    print(f"AUC(-NU -> recovery) at 17 bits = {auc(p17, n17):.4f}")
    pr17p = [r['nuhat'] * math.sqrt(r['eff']) for r in rows17 if r['ok']]
    pr17n = [r['nuhat'] * math.sqrt(r['eff']) for r in rows17 if not r['ok']]
    print(f"AUC(-nu_hat*sqrt(eff) -> recovery) at 17 bits = "
          f"{auc(pr17p, pr17n):.4f}")
    allp = [r['nuhat'] * math.sqrt(r['eff']) for r in rows17]
    allo = [r['NU'] for r in rows17]
    Cf17 = sum(o / p for o, p in zip(allo, allp) if p > 0) / len(allp)
    print(f"fitted C at 17 bits = {Cf17:.4f}   (12 bits: {Cfit:.4f})")

    print(f"\nargmax at last index in "
          f"{sum(1 for r in rows17 if r['argmax'] == r['k']-1)}/{len(rows17)}; "
          f"last two in "
          f"{sum(1 for r in rows17 if r['argmax'] >= r['k']-2)}/{len(rows17)}")
    r17 = [r['prof'][-1] / r['l2'] for r in rows17]
    print(f"||b*_last||/lambda_2(L2) at 17 bits: mean "
          f"{sum(r17)/len(r17):.3f}  min {min(r17):.3f}  max {max(r17):.3f}")

    # float-vs-exact at dim 24
    worst = 0.0
    for r0 in rows17[:20]:
        pass
    chk = 0
    worst24 = 0.0
    for eff in (0.15,):
        for (p, b, n, lam, G) in curves17[:6]:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            d_trial = random.Random(42 + 7777).randint(1, n - 1)
            re_ = instance((p, b, n, lam, G), M17, d_trial, k1b, 42, exact=True)
            rf = instance((p, b, n, lam, G), M17, d_trial, k1b, 42, exact=False)
            if re_ and rf:
                chk += 1
                worst24 = max(worst24, abs(rf['NU'] - re_['NU']) / re_['NU'])
    print(f"\nW0 repeat at dim {M17*2}: max relative NU error (float vs exact) "
          f"over {chk} cells = {worst24:.3e}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
