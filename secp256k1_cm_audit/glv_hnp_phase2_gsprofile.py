"""
GLV-HNP Phase 2, Thread 24: derive nu_hat from the Gram-Schmidt profile of L0.

Pre-registered in the 2026-08-07 (Thread 23) log entry.

Background.  Thread 23b established that Phase-2 recovery is a BDD condition,
certified per-instance by

    nu_i := 2*|<e, b*_i>| / ||b*_i||^2 ,      NU := max_i nu_i ,

on the LLL-reduced basis of the Kannan-free projected lattice

    L0 = < n*S_K1*e_i , (B_i*S_K1 | 0) , (-lam*S_K1*e_j | S_K2*e_j) >  in Z^2m
    e  = (k1_i*S_K1 | k2_i*S_K2)              (secret error)

NU <= 1 was a SOUND certificate (0 FP in 110 instances, AUC 0.978) but loose by
~1.9x against what Kannan-LLL actually achieves.  NU is a `max_i`, so its value
is decided by ONE Gram-Schmidt index.  Thread 23 also showed the *unnormalised*
lambda-block shortest vector

    mu = lambda_1(L2),   L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >  in Z^2

is exactly the object Thread 20b fitted as nu_hat = mu/sqrt(det L2), and that
its correlation with recovery is NEGATIVE (bigger mu => more failures) — an
empirical fit with no derivation.

Hypothesis under test (H24, pre-registered verbatim):

    H24: the argmax of nu_i concentrates at the tail indices i ~ 2m, and
         ||b*_{2m}|| ~ det(L2)/mu = lambda_2(L2) up to a constant.

If H24 holds, the negative sign of nu_hat is DERIVED rather than fitted:
small mu => skew L2 => large lambda_2(L2) => large tail GS norms => small
tail nu_i => small NU => BDD succeeds.

Falsifier (pre-registered): if argmax_i is spread across i rather than
concentrated in the tail, the GS-profile explanation fails and nu_hat stays
empirical.

Experiments
    W1  full GS profile log2||b*_i|| and nu_i profile over the U2 grid;
        histogram of argmax_i by relative position.
    W2  ||b*_{2m}|| vs lambda_2(L2) = det(L2)/mu — the quantitative half of H24.
    W3  is NU decided by the tail?  NU_tail (last quarter) vs NU, and how much
        predictive power survives if only the tail is inspected.
    W4  does the GS profile explain the SIGN of nu_hat?  mu vs tail GS norms
        vs recovery, at matched K1.

Run: python3 glv_hnp_phase2_gsprofile.py
"""

import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import find_generator, scales, gen_signatures
from glv_hnp_phase2_projected import run_new, SEEDS, HIST
from glv_hnp_phase2_babai import (
    build_L0, target_and_error, gram_schmidt, lll_rows,
)


# ---------------------------------------------------------------------------
# exact 2D reduction returning BOTH successive minima
# ---------------------------------------------------------------------------

def lagrange_2d(u, v):
    """Lagrange-reduced basis of a rank-2 integer lattice.  For a reduced
    basis |u| = lambda_1 and |v| = lambda_2 (Lagrange/Gauss is optimal in 2D)."""
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]
    if n2(u) > n2(v):
        u, v = v, u
    while True:
        num, den = dot(v, u), n2(u)
        if den == 0:
            break
        q = ((2 * num + den) // (2 * den) if num >= 0
             else -((-2 * num + den) // (2 * den)))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if n2(v) >= n2(u):
            break
        u, v = v, u
    return u, v


def l2_minima(n, lam, S_K1, S_K2):
    """(lambda_1, lambda_2, det) of the lambda-block sublattice L2."""
    u, v = lagrange_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    l1 = math.sqrt(u[0] * u[0] + u[1] * u[1])
    l2 = math.sqrt(v[0] * v[0] + v[1] * v[1])
    return l1, l2, float(n) * S_K1 * S_K2


# ---------------------------------------------------------------------------
# per-instance GS + nu profile
# ---------------------------------------------------------------------------

def profile(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    B = lll_rows(build_L0(sigs, n, lam, k1_bound, k2_bound), 2 * m)
    _tau, e = target_and_error(sigs, n, k1_bound, k2_bound)
    Bs, _ = gram_schmidt(B)
    gs, nus = [], []
    for i in range(len(Bs)):
        ni = sum(x * x for x in Bs[i])
        gs.append(math.sqrt(ni))
        nus.append(0.0 if ni == 0.0
                   else 2.0 * abs(sum(a * c for a, c in zip(e, Bs[i]))) / ni)
    NU = max(nus)
    arg = max(range(len(nus)), key=lambda i: nus[i])
    return {'gs': gs, 'nus': nus, 'NU': NU, 'arg': arg, 'dim': len(Bs)}


def bucket(arg, dim, nb=4):
    """Which quarter of the GS index range the argmax sits in (0 = head)."""
    return min(nb - 1, (arg * nb) // dim)


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 24 — is NU decided by the tail of the GS profile, and is the")
    print("            tail governed by lambda_2(L2) = det(L2)/mu?")
    print("=" * 78)

    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, label
        hist.append((label, (p, b, n, lam, G), k1, m))

    U2 = [("12-bit/2557", hist[1][1], 8, 0.340),
          ("12-bit/2677", hist[2][1], 10, 0.070)]
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
    U2_by_label = [(lab, cur) for lab, cur, _m, _ls in U2]

    ALL = []   # (label, K1, seed, prof, recovered, mu, lam2, detL2, m)

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W1: GS profile and argmax_i of nu_i over the U2 grid")
    print("-" * 78)
    print("dim(L0) = 2m.  'argmax' = index i (1-based) attaining NU.")
    print("'q' = quarter of the index range holding the argmax (1=head, 4=tail).\n")

    for label, curve, m, ls in U2:
        p, b, n, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound=1, k2_bound=k2_bound)
        print(f"\n{label}  (lam* = {ls:.3f}, m = {m}, dim = {2*m})")
        print(f"  {'K1':>4} {'mu':>10} {'lam2(L2)':>10} {'argmax(1b)':>11} "
              f"{'mean NU':>9} {'q-hist(1..4)':>14} {'recov':>6}")
        for k1 in K1_GRID:
            S_K1, _S_D, S_K2, _S_KAN = scales(n, k1, k2_bound)
            mu, lam2, detL2 = l2_minima(n, lam, S_K1, S_K2)
            args, nus, qh, rec = [], [], [0] * 4, 0
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                pr = profile(curve, m, d_trial, k1, seed)
                rk = run_new(curve, m, d_trial, k1, seed)
                if pr is None:
                    continue
                args.append(pr['arg'] + 1)
                nus.append(pr['NU'])
                qh[bucket(pr['arg'], pr['dim'])] += 1
                rec += bool(rk['ok'])
                ALL.append((label, k1, seed, pr, bool(rk['ok']),
                            mu, lam2, detL2, m))
            print(f"  {k1:>4} {mu:>10.0f} {lam2:>10.0f} "
                  f"{str(args):>11} {sum(nus)/len(nus):>9.3f} "
                  f"{str(qh):>14} {str(rec)+'/'+str(len(SEEDS)):>6}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W1b: aggregate argmax position — the PRE-REGISTERED FALSIFIER")
    print("-" * 78)
    qh = [0] * 4
    tail_frac_num = 0
    for (_l, _k, _s, pr, _r, _mu, _l2, _d, _m) in ALL:
        qh[bucket(pr['arg'], pr['dim'])] += 1
        if pr['arg'] >= pr['dim'] - 2:      # last three indices
            tail_frac_num += 1
    tot = len(ALL)
    print(f"N = {tot} instances")
    for q in range(4):
        print(f"  quarter {q+1} ({'head' if q == 0 else 'tail' if q == 3 else '  '}): "
              f"{qh[q]:>4}  ({100.0*qh[q]/tot:5.1f}%)")
    print(f"  argmax in the LAST THREE indices: {tail_frac_num}/{tot} "
          f"({100.0*tail_frac_num/tot:.1f}%)")
    print("H24 predicts quarter 4 dominant.  If the mass is spread, H24 fails.")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W2: is ||b*_{2m}|| ~ det(L2)/mu = lambda_2(L2)?")
    print("-" * 78)
    print("Ratios ||b*_last|| / lambda_2(L2) and ||b*_last|| / (det L2 / mu).")
    print("H24 wants these ~constant across K1 and across curves.\n")
    print(f"  {'curve':<12} {'K1':>4} {'||b*_last||':>12} {'lam2(L2)':>10} "
          f"{'ratio':>8} {'detL2/mu':>10} {'ratio':>8}")
    by_key = {}
    for (l, k, s, pr, r, mu, lam2, detL2, m) in ALL:
        by_key.setdefault((l, k), []).append((pr, mu, lam2, detL2))
    ratios_l2, ratios_dm = [], []
    for (l, k), v in sorted(by_key.items(), key=lambda t: (t[0][0], t[0][1])):
        last = sum(pr['gs'][-1] for pr, _, _, _ in v) / len(v)
        _pr, mu, lam2, detL2 = v[0]
        r1 = last / lam2
        r2 = last / (detL2 / mu)
        ratios_l2.append(r1)
        ratios_dm.append(r2)
        print(f"  {l:<12} {k:>4} {last:>12.1f} {lam2:>10.1f} {r1:>8.3f} "
              f"{detL2/mu:>10.1f} {r2:>8.3f}")
    def spread(xs):
        mn, mx = min(xs), max(xs)
        return mn, mx, mx / mn if mn > 0 else float('inf')
    a, b_, c = spread(ratios_l2)
    print(f"\n  ||b*_last||/lambda_2 : min {a:.3f} max {b_:.3f} "
          f"spread {c:.2f}x  over {len(ratios_l2)} cells")
    a, b_, c = spread(ratios_dm)
    print(f"  ||b*_last||/(detL2/mu): min {a:.3f} max {b_:.3f} spread {c:.2f}x")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W3: does a TAIL-ONLY certificate keep NU's predictive power?")
    print("-" * 78)
    for frac_name, lo in (("last quarter", 0.75), ("last half", 0.5),
                          ("last index only", None)):
        pos, neg, ident = [], [], 0
        for (_l, _k, _s, pr, rec, _mu, _l2, _d, _m) in ALL:
            d = pr['dim']
            if lo is None:
                sub = pr['nus'][-1:]
            else:
                sub = pr['nus'][int(lo * d):]
            nt = max(sub)
            if abs(nt - pr['NU']) < 1e-9:
                ident += 1
            (pos if rec else neg).append(nt)
        conc = sum((1.0 if x < y else 0.5 if x == y else 0.0)
                   for x in pos for y in neg)
        auc = conc / (len(pos) * len(neg))
        print(f"  {frac_name:<16}  NU_tail == NU in {ident:>3}/{tot} "
              f"({100.0*ident/tot:5.1f}%)   AUC(-NU_tail) = {auc:.4f}")
    pos = [pr['NU'] for (_l, _k, _s, pr, rec, *_x) in ALL if rec]
    neg = [pr['NU'] for (_l, _k, _s, pr, rec, *_x) in ALL if not rec]
    conc = sum((1.0 if x < y else 0.5 if x == y else 0.0)
               for x in pos for y in neg)
    print(f"  {'full NU (ref)':<16}  {'':<24}   "
          f"AUC(-NU)      = {conc/(len(pos)*len(neg)):.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W4: does the GS profile explain the SIGN of nu_hat?")
    print("-" * 78)
    print("nu_hat = mu/sqrt(det L2) correlated NEGATIVELY with recovery (Thread")
    print("20b).  H24's mechanism: mu up => lambda_2 down => tail GS norms down")
    print("=> nu_i up => NU up => BDD fails.  Check each link at matched K1.\n")
    print(f"  {'K1':>4} | {'2557: mu':>9} {'lam2':>8} {'b*_last':>9} {'NU':>7} "
          f"{'rec':>4} | {'2677: mu':>9} {'lam2':>8} {'b*_last':>9} {'NU':>7} {'rec':>4}")
    for k in K1_GRID:
        cells = []
        for l in ("12-bit/2557", "12-bit/2677"):
            v = by_key.get((l, k))
            rec = sum(1 for (ll, kk, _s, _p, r, *_x) in ALL
                      if ll == l and kk == k and r)
            last = sum(pr['gs'][-1] for pr, _, _, _ in v) / len(v)
            nu = sum(pr['NU'] for pr, _, _, _ in v) / len(v)
            _pr, mu, lam2, _d = v[0]
            cells.append((mu, lam2, last, nu, rec, len(v)))
        (m1, l1, b1, n1, r1, t1), (m2, l2_, b2, n2, r2, t2) = cells
        print(f"  {k:>4} | {m1:>9.0f} {l1:>8.0f} {b1:>9.1f} {n1:>7.3f} "
              f"{str(r1)+'/'+str(t1):>4} | {m2:>9.0f} {l2_:>8.0f} {b2:>9.1f} "
              f"{n2:>7.3f} {str(r2)+'/'+str(t2):>4}")

    # rank correlations over all cells (Spearman via ranks, ties averaged)
    def spearman(xs, ys):
        def rank(a):
            order = sorted(range(len(a)), key=lambda i: a[i])
            r = [0.0] * len(a)
            i = 0
            while i < len(order):
                j = i
                while j + 1 < len(order) and a[order[j + 1]] == a[order[i]]:
                    j += 1
                avg = (i + j) / 2.0 + 1.0
                for t in range(i, j + 1):
                    r[order[t]] = avg
                i = j + 1
            return r
        rx, ry = rank(xs), rank(ys)
        mx = sum(rx) / len(rx); my = sum(ry) / len(ry)
        num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
        den = math.sqrt(sum((a - mx) ** 2 for a in rx)
                        * sum((b - my) ** 2 for b in ry))
        return num / den if den else float('nan')

    mus = [x[5] for x in ALL]
    lams = [x[6] for x in ALL]
    lasts = [x[3]['gs'][-1] for x in ALL]
    nusv = [x[3]['NU'] for x in ALL]
    recs = [1.0 if x[4] else 0.0 for x in ALL]
    print(f"\n  Spearman over all {tot} instances:")
    print(f"    rho(mu, lambda_2)     = {spearman(mus, lams):+.3f}   "
          f"(exact duality: mu*lambda_2 = det L2 up to <=2/sqrt3)")
    print(f"    rho(lambda_2, b*_last)= {spearman(lams, lasts):+.3f}")
    print(f"    rho(b*_last, NU)      = {spearman(lasts, nusv):+.3f}")
    print(f"    rho(NU, recovery)     = {spearman(nusv, recs):+.3f}")
    print(f"    rho(mu, recovery)     = {spearman(mus, recs):+.3f}   "
          f"(Thread 20b's nu_hat sign)")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W4b: the nu_hat sign, WITHIN matched K1 (removes the K1 confound)")
    print("-" * 78)
    print("Pooling over K1 confounds every correlation: K1 drives det(L2) and")
    print("hence mu, lambda_2 and recovery simultaneously.  Thread 20b's nu_hat")
    print("comparison was CROSS-CURVE AT MATCHED K1, so redo it that way: for")
    print("each K1, sign of (mu_2677 - mu_2557) vs (rec_2677 - rec_2557).\n")
    print(f"  {'K1':>4} {'d(mu)':>9} {'d(lam2)':>9} {'d(b*last)':>10} "
          f"{'d(NU)':>8} {'d(rec)':>7}  {'nu_hat sign ok?':>16}")
    agree_sign = ties = 0
    for k in K1_GRID:
        va = by_key[("12-bit/2557", k)]
        vb = by_key[("12-bit/2677", k)]
        ra = sum(1 for (ll, kk, _s, _p, r, *_x) in ALL
                 if ll == "12-bit/2557" and kk == k and r)
        rb = sum(1 for (ll, kk, _s, _p, r, *_x) in ALL
                 if ll == "12-bit/2677" and kk == k and r)
        dmu = vb[0][1] - va[0][1]
        dl2 = vb[0][2] - va[0][2]
        dbl = (sum(pr['gs'][-1] for pr, *_ in vb) / len(vb)
               - sum(pr['gs'][-1] for pr, *_ in va) / len(va))
        dnu = (sum(pr['NU'] for pr, *_ in vb) / len(vb)
               - sum(pr['NU'] for pr, *_ in va) / len(va))
        drec = rb - ra
        if drec == 0:
            ok = "tie"
            ties += 1
        else:
            # nu_hat's claimed sign: larger mu => worse recovery
            ok = "YES" if (dmu > 0) == (drec < 0) else "no"
            agree_sign += (ok == "YES")
        print(f"  {k:>4} {dmu:>9.0f} {dl2:>9.0f} {dbl:>10.1f} {dnu:>8.3f} "
              f"{drec:>7} {ok:>16}")
    print(f"\n  sign of nu_hat confirmed in {agree_sign}/{11-ties} "
          f"non-tied K1 cells ({ties} ties).")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W5: NU_half as a separator — sharpness, and is it still SOUND?")
    print("-" * 78)
    print("W3 found AUC(-NU_half) = 0.996 > AUC(-NU) = 0.978.  NU_half drops")
    print("the head indices, so it can no longer be a nearest-plane certificate;")
    print("check directly whether NU_half <= 1 still implies Babai success.\n")

    from glv_hnp_phase2_babai import run_babai

    rowsW5 = []
    for (l, k, s, pr, rec, mu, lam2, detL2, m) in ALL:
        d = pr['dim']
        nh = max(pr['nus'][d // 2:])
        rowsW5.append((l, k, s, pr['NU'], nh, rec))

    for name, idx in (("NU (full)", 3), ("NU_half", 4)):
        tp = fp = fn = tn = 0
        for r in rowsW5:
            pred, obs = r[idx] <= 1.0, r[5]
            if pred and obs: tp += 1
            elif pred and not obs: fp += 1
            elif not pred and obs: fn += 1
            else: tn += 1
        pos = [r[idx] for r in rowsW5 if r[5]]
        neg = [r[idx] for r in rowsW5 if not r[5]]
        thr_suff, thr_nec = min(neg), max(pos)
        band = sum(1 for r in rowsW5 if thr_suff <= r[idx] <= thr_nec)
        print(f"  {name}  vs Kannan-LLL recovery   N = {len(rowsW5)}")
        print(f"    TP {tp:>3}  FP {fp:>3}  FN {fn:>3}  TN {tn:>3}   "
              f"accuracy {(tp+tn)/len(rowsW5):.3f}")
        print(f"    sufficient below {thr_suff:.3f} | necessary above "
              f"{thr_nec:.3f} | ambiguous band {band}/{len(rowsW5)} "
              f"({100.0*band/len(rowsW5):.1f}%)")

    # soundness check: NU_half <= 1 must NOT be claimed as a Babai certificate
    sound_fp = sound_n = 0
    for (l, k, s, pr, rec, mu, lam2, detL2, m) in ALL:
        d = pr['dim']
        nh = max(pr['nus'][d // 2:])
        if nh <= 1.0:
            curve = dict(U2_by_label)[l]
            d_trial = random.Random(s + 7777).randint(1, curve[2] - 1)
            rb = run_babai(curve, m, d_trial, k, s)
            sound_n += 1
            sound_fp += (not rb['ok'])
    print(f"\n  soundness: of {sound_n} instances with NU_half <= 1, "
          f"{sound_fp} FAIL Babai nearest-plane.")
    print("  (NU <= 1 had 0 such failures by construction — it is the theorem.)")
    print("  => NU_half is a sharper SCORE but not a certificate; the pair")
    print("     (NU for soundness, NU_half for ranking) is the useful package.")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W6: where does the 0.47 constant come from?  det(L0) vs det(L2)^m")
    print("-" * 78)
    print("L0 contains L2^m (one copy of L2 on each coordinate pair (i, m+i))")
    print("plus the B-row, so det(L0) = det(L2)^m / [L0 : L2^m], and the whole")
    print("GS profile is pulled down by that index.  Measure the index.\n")
    print(f"  {'curve':<12} {'K1':>4} {'log2 detL0':>11} {'m*log2 detL2':>13} "
          f"{'index':>12} {'index<=n?':>10} {'I^(-1/2m)':>10} {'b*last/lam2':>12}")
    for (l, k), v in sorted(by_key.items(), key=lambda t: (t[0][0], t[0][1])):
        pr0, mu, lam2, detL2 = v[0]
        mm = len(pr0['gs']) // 2
        ld0 = sum(math.log2(x) for x in pr0['gs'])
        ld2 = mm * math.log2(detL2)
        idx = 2.0 ** (ld2 - ld0)
        n_l = 2659 if l == "12-bit/2557" else 2647
        last = sum(p['gs'][-1] for p, *_ in v) / len(v)
        print(f"  {l:<12} {k:>4} {ld0:>11.2f} {ld2:>13.2f} {idx:>12.1f} "
              f"{str(idx <= n_l * 1.001):>10} {idx**(-1.0/(2*mm)):>10.3f} "
              f"{last/lam2:>12.3f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W7: OUT-OF-SAMPLE — do NU / NU_half transfer to 17-bit curves?")
    print("-" * 78)
    print("Everything above is fitted on two 12-bit curves.  The 12-bit brackets")
    print("were NU in [1.188, 1.869] and NU_half in [1.068, 1.252].  Re-measure")
    print("on the 20 independent 17-bit curves of Thread 23 exp U3 (m = 12,")
    print("eff = K1*K2/n in {0.05, 0.15, 0.25}, 5 seeds) — 300 fresh instances.")
    print("If the brackets are n-independent, NU gives a size-free viability")
    print("test for Phase 2.\n")

    import time as _time
    from glv_hnp_common import search_curves

    t0 = _time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"found {len(curves17)} 17-bit j=0 GLV curves in {_time.time()-t0:.1f}s")

    M17 = 12
    rows17 = []
    for eff in (0.05, 0.15, 0.25):
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                pr = profile((p, b, n, lam, G), M17, d_trial, k1b, seed)
                if pr is None:
                    continue
                rk = run_new((p, b, n, lam, G), M17, d_trial, k1b, seed)
                nh = max(pr['nus'][pr['dim'] // 2:])
                rows17.append((eff, n, k1b, seed, pr['NU'], nh,
                               bool(rk['ok']), pr['arg'], pr['dim']))
    print(f"{len(rows17)} instances in {_time.time()-t0:.1f}s")

    # reproduction guard: Thread 23 exp U3 measured NEW-lattice wins
    # 99/100, 21/100, 9/100 at eff = 0.05, 0.15, 0.25.
    REF_U3 = {0.05: 99, 0.15: 21, 0.25: 9}
    print("  reproduction check vs Thread 23 exp U3 (NEW-lattice wins):")
    for eff in (0.05, 0.15, 0.25):
        w = sum(1 for r in rows17 if r[0] == eff and r[6])
        t = sum(1 for r in rows17 if r[0] == eff)
        flag = "OK" if w == REF_U3[eff] else f"DRIFT (ref {REF_U3[eff]})"
        print(f"    eff {eff:.2f}: {w}/{t}   {flag}")
    print()

    for name, idx, band12 in (("NU (full)", 4, (1.188, 1.869)),
                              ("NU_half", 5, (1.068, 1.252))):
        tp = fp = fn = tn = 0
        for r in rows17:
            pred, obs = r[idx] <= 1.0, r[6]
            if pred and obs: tp += 1
            elif pred and not obs: fp += 1
            elif not pred and obs: fn += 1
            else: tn += 1
        pos = [r[idx] for r in rows17 if r[6]]
        neg = [r[idx] for r in rows17 if not r[6]]
        conc = sum((1.0 if x < y else 0.5 if x == y else 0.0)
                   for x in pos for y in neg)
        auc = conc / (len(pos) * len(neg)) if pos and neg else float('nan')
        ts, tn_ = (min(neg) if neg else float('nan')), (max(pos) if pos else float('nan'))
        print(f"  {name}   N = {len(rows17)}   AUC = {auc:.4f}")
        print(f"    TP {tp:>4}  FP {fp:>4}  FN {fn:>4}  TN {tn:>4}   "
              f"accuracy {(tp+tn)/len(rows17):.3f}")
        print(f"    17-bit bracket [{ts:.3f}, {tn_:.3f}]   "
              f"12-bit bracket [{band12[0]:.3f}, {band12[1]:.3f}]")
        # does the 12-bit bracket still hold out-of-sample?
        v_lo = sum(1 for r in rows17 if r[idx] < band12[0] and not r[6])
        v_hi = sum(1 for r in rows17 if r[idx] > band12[1] and r[6])
        n_lo = sum(1 for r in rows17 if r[idx] < band12[0])
        n_hi = sum(1 for r in rows17 if r[idx] > band12[1])
        print(f"    transfer: {v_lo}/{n_lo} violations below the 12-bit lower "
              f"threshold, {v_hi}/{n_hi} above the upper\n")

    qh17 = [0] * 4
    for r in rows17:
        qh17[bucket(r[7], r[8])] += 1
    print(f"  argmax quarter histogram at 17 bits: {qh17}  "
          f"(tail = {100.0*qh17[3]/len(rows17):.1f}%)")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
