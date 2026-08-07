"""
GLV-HNP Phase 2, Thread 25: is mu a genuine SECOND coordinate, or is its
apparent power mediated by NU after all?

Pre-registered by the 2026-08-07 run #2 (Thread 24) log entry:

  H25: within the ambiguous band  1.040 <= NU <= 2.199  (where the
       nearest-plane guarantee gives no answer), AUC(-mu -> Kannan-LLL
       recovery) stays >= 0.8.
  If yes  : mu is a genuine second coordinate and (NU, mu) is the
            2-parameter viability test Phase 2 has been missing; deliverable
            is a logistic fit on (log NU, log mu) with its decision boundary.
  If no   : mu's apparent power is entirely mediated by NU, the W5 result is
            a stratification artifact, and the closed form should be retired.

Background.  Thread 24 (W5/W6) found that nu_hat = mu/sqrt(det L2) predicts
recovery BETTER than the exact BDD certificate NU when the bias strength
eff = K1*K2/n is held fixed (pooled AUC 0.935 vs 0.800), and that NU and
nu_hat*sqrt(eff) are mutually UNCORRELATED inside every eff stratum
(Spearman -0.28 .. +0.16).  Two uncorrelated quantities both predicting
recovery is the signature of two mechanisms.  This script tests that directly
by conditioning on NU instead of on eff.

Secondary (also pre-registered, from W1b):  the GS profile of L0 is two
blocks -- the first m norms sit at lambda_1(L2) and only the second block
moves with K1 -- and the step between blocks VANISHES right as the K1 wall is
crossed.  Define

    step = log2 ||b*_{m+1}||  -  log2 ||b*_1||

and test whether step -> 0 predicts the wall better than either NU or mu.

The instance table is the same 5 strata x 20 curves x 5 seeds at 17 bits used
by glv_hnp_phase2_gsprofile_strat.py, regenerated identically (every draw is
seeded, so the rows are reproducible), plus the `step` column.  --dump-json
writes the table so future runs can re-analyse without recomputing.

Gram-Schmidt is float, justified by W0 of the parent script (max relative NU
error vs exact Fractions ~6.6e-16 at dim 24 / 17 bits).

Run: python3 glv_hnp_phase2_nuband.py [--dump-json FILE]
"""

import json
import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

# Pre-registered band, verbatim from Thread 24 W4 (17 bits, 300 instances):
#   sufficient  NU < 1.040   (zero false positives observed at NU <= 1)
#   necessary   NU > 2.199   (no recovery ever observed above)
NU_LO = 1.040
NU_HI = 2.199

M17 = 12
EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)


# ---------------------------------------------------------------------------
# Logistic regression (IRLS with a small ridge, no scipy/sklearn)
# ---------------------------------------------------------------------------

def _solve(A, b):
    """Gaussian elimination with partial pivoting.  A is n x n, b is n."""
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c]) < 1e-14:
            return None
        M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r == c:
                continue
            f = M[r][c] / M[c][c]
            for k in range(c, n + 1):
                M[r][k] -= f * M[c][k]
    return [M[i][n] / M[i][i] for i in range(n)]


def logistic_fit(X, y, ridge=1e-6, iters=60):
    """X: list of feature rows (no intercept column).  y: 0/1.
    Returns (beta, deviance) with beta[0] the intercept."""
    n, p = len(X), len(X[0]) + 1
    Z = [[1.0] + list(row) for row in X]
    beta = [0.0] * p
    for _ in range(iters):
        eta = [sum(b * z for b, z in zip(beta, row)) for row in Z]
        pr = [1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, e)))) for e in eta]
        W = [max(pi * (1 - pi), 1e-9) for pi in pr]
        H = [[sum(Z[i][a] * W[i] * Z[i][b] for i in range(n))
              + (ridge if a == b else 0.0) for b in range(p)] for a in range(p)]
        g = [sum(Z[i][a] * (y[i] - pr[i]) for i in range(n))
             - ridge * beta[a] for a in range(p)]
        step = _solve(H, g)
        if step is None:
            break
        beta = [b + s for b, s in zip(beta, step)]
        if max(abs(s) for s in step) < 1e-10:
            break
    eta = [sum(b * z for b, z in zip(beta, row)) for row in Z]
    pr = [1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, e)))) for e in eta]
    dev = -2.0 * sum(math.log(max(1e-12, pi if yi else 1 - pi))
                     for pi, yi in zip(pr, y))
    return beta, dev


def logistic_score(beta, row):
    return beta[0] + sum(b * x for b, x in zip(beta[1:], row))


def standardize(X):
    p = len(X[0])
    mus = [sum(r[j] for r in X) / len(X) for j in range(p)]
    sds = [math.sqrt(max(1e-18, sum((r[j] - mus[j]) ** 2 for r in X) / len(X)))
           for j in range(p)]
    return [[(r[j] - mus[j]) / sds[j] for j in range(p)] for r in X], mus, sds


def auc_from_scores(scores, y):
    """AUC for a score where LARGER -> more likely y = 1."""
    pos = [s for s, t in zip(scores, y) if t]
    neg = [s for s, t in zip(scores, y) if not t]
    if not pos or not neg:
        return float('nan')
    conc = sum((1.0 if a > b else 0.5 if a == b else 0.0)
               for a in pos for b in neg)
    return conc / (len(pos) * len(neg))


# ---------------------------------------------------------------------------
# Table
# ---------------------------------------------------------------------------

def build_table():
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")

    t0 = time.time()
    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), M17, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), M17, d_trial, k1b, seed)
                prof = r['prof']
                # W1b two-block step: head sits at lambda_1(L2), second block
                # moves with K1; the step vanishes as the wall is crossed.
                step = (math.log2(prof[M17]) - math.log2(prof[0])
                        if prof[0] > 0 and prof[M17] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                r.pop('prof', None)
                r.pop('nus', None)
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")
    return rows


def band_report(rows, lo, hi, label):
    sub = [r for r in rows if lo <= r['NU'] <= hi]
    pos = [r for r in sub if r['ok']]
    neg = [r for r in sub if not r['ok']]
    print(f"\n{label}: {lo:.3f} <= NU <= {hi:.3f}   "
          f"N={len(sub)}  recovered {len(pos)}/{len(sub)}")
    if not pos or not neg:
        print("  (degenerate — one class empty, no AUC)")
        return None
    scores = [('mu', 'mu'), ('nu_hat', 'nuhat'), ('NU', 'NU'),
              ('lam*', 'lamstar'), ('eff', 'eff'), ('step', 'step')]
    print(f"  {'score':>8} {'AUC(-x)':>9}")
    out = {}
    for name, key in scores:
        a = auc([r[key] for r in pos], [r[key] for r in neg])
        out[name] = a
        print(f"  {name:>8} {a:>9.4f}")
    return out


if __name__ == "__main__":
    dump = None
    if "--dump-json" in sys.argv:
        dump = sys.argv[sys.argv.index("--dump-json") + 1]

    print("=" * 78)
    print("Thread 25 — is mu a second coordinate? condition on NU, not on eff")
    print("=" * 78)
    print(f"pre-registered band: {NU_LO} <= NU <= {NU_HI}")
    print(f"pre-registered H25 : AUC(-mu -> recovery) >= 0.8 inside the band")

    rows = build_table()

    if dump:
        with open(dump, "w") as fh:
            json.dump(rows, fh)
        print(f"table dumped to {dump} ({os.path.getsize(dump)} bytes)")

    npos = sum(1 for r in rows if r['ok'])
    print(f"\noverall recovery {npos}/{len(rows)}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1: how the pre-registered NU bracket partitions this table")
    print("-" * 78)
    lowb = [r for r in rows if r['NU'] < NU_LO]
    midb = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    hib = [r for r in rows if r['NU'] > NU_HI]
    for nm, g in (("NU < 1.040  (sufficient)", lowb),
                  ("1.040 <= NU <= 2.199 (ambiguous)", midb),
                  ("NU > 2.199  (necessary violated)", hib)):
        w = sum(1 for r in g if r['ok'])
        rate = f"{w}/{len(g)}" if g else "0/0"
        pct = f"{100.0*w/len(g):.1f}%" if g else "--"
        print(f"  {nm:<34} N={len(g):>4}  recovered {rate:>8}  ({pct})")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X2: THE TEST — does mu still separate INSIDE the NU band?")
    print("-" * 78)
    print("AUC > 0.5 means SMALLER score -> more likely to recover.")
    main = band_report(rows, NU_LO, NU_HI, "pre-registered band")

    if main is not None:
        verdict = "CONFIRMED" if main['mu'] >= 0.8 else "FALSIFIED"
        print(f"\n  H25 (AUC(-mu) >= 0.8 inside the band): {verdict} "
              f"(AUC = {main['mu']:.4f})")

    print("\nrobustness — same test on shifted/narrowed bands:")
    for lo, hi, lab in ((1.0, 2.0, "band [1.0, 2.0]"),
                        (1.2, 1.8, "narrow [1.2, 1.8]"),
                        (0.8, 1.5, "low    [0.8, 1.5]"),
                        (1.5, 3.0, "high   [1.5, 3.0]")):
        band_report(rows, lo, hi, lab)

    # data-driven band from THIS table
    okNU = [r['NU'] for r in rows if r['ok']]
    bad = [r['NU'] for r in rows if not r['ok']]
    if okNU and bad:
        print(f"\nthis table's own bracket: sufficient NU < {min(bad):.4f} , "
              f"necessary NU > {max(okNU):.4f}")
        band_report(rows, min(bad), max(okNU), "self-consistent band")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X3: two-parameter logistic fit on (log NU, log mu)")
    print("-" * 78)
    y = [1 if r['ok'] else 0 for r in rows]
    f_nu = [[math.log(r['NU'])] for r in rows]
    f_mu = [[math.log(r['mu'])] for r in rows]
    f_both = [[math.log(r['NU']), math.log(r['mu'])] for r in rows]

    for name, F in (("log NU only", f_nu), ("log mu only", f_mu),
                    ("log NU + log mu", f_both)):
        Z, mus, sds = standardize(F)
        beta, dev = logistic_fit(Z, y)
        sc = [logistic_score(beta, z) for z in Z]
        a = auc_from_scores(sc, y)
        acc = sum(1 for s, t in zip(sc, y) if (s > 0) == bool(t)) / len(y)
        coefs = "  ".join(f"{b:+.3f}" for b in beta[1:])
        print(f"  {name:<16} deviance {dev:>8.2f}  AUC {a:.4f}  "
              f"acc {acc:.3f}  std-coef {coefs}")

    Z1, _, _ = standardize(f_nu)
    _, dev1 = logistic_fit(Z1, y)
    Z2, m2, s2 = standardize(f_both)
    beta2, dev2 = logistic_fit(Z2, y)
    lr = dev1 - dev2
    print(f"\n  likelihood-ratio (adding log mu to log NU): {lr:.2f} on 1 df")
    print(f"  chi2 95% = 3.84, 99.9% = 10.83 -> "
          f"{'SIGNIFICANT' if lr > 10.83 else 'significant' if lr > 3.84 else 'not significant'}")

    # decision boundary in raw units
    b0, bn, bm = beta2
    print("\n  decision boundary (p = 1/2), raw units:")
    print(f"    {bn/s2[0]:+.4f}*log(NU) {bm/s2[1]:+.4f}*log(mu) "
          f"{b0 - bn*m2[0]/s2[0] - bm*m2[1]/s2[1]:+.4f} = 0")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X4: leave-one-curve-out CV (group by n) — is the pair real,")
    print("        or is it memorising the 20 curves?")
    print("-" * 78)
    groups = sorted({r['n'] for r in rows})
    for name, feat in (("log NU", lambda r: [math.log(r['NU'])]),
                       ("log mu", lambda r: [math.log(r['mu'])]),
                       ("log NU + log mu",
                        lambda r: [math.log(r['NU']), math.log(r['mu'])]),
                       ("log NU + log nu_hat",
                        lambda r: [math.log(r['NU']), math.log(r['nuhat'])])):
        preds, ys = [], []
        for gname in groups:
            tr = [r for r in rows if r['n'] != gname]
            te = [r for r in rows if r['n'] == gname]
            Xtr = [feat(r) for r in tr]
            ytr = [1 if r['ok'] else 0 for r in tr]
            _, mus, sds = standardize(Xtr)
            Ztr = [[(v - mus[j]) / sds[j] for j, v in enumerate(x)]
                   for x in Xtr]
            beta, _ = logistic_fit(Ztr, ytr)
            for r in te:
                z = [(v - mus[j]) / sds[j] for j, v in enumerate(feat(r))]
                preds.append(logistic_score(beta, z))
                ys.append(1 if r['ok'] else 0)
        a = auc_from_scores(preds, ys)
        acc = sum(1 for s, t in zip(preds, ys) if (s > 0) == bool(t)) / len(ys)
        print(f"  {name:<22} CV-AUC {a:.4f}   CV-acc {acc:.3f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X5: the W1b two-block step  step = log2||b*_{m+1}|| - log2||b*_1||")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>8} | {'mean step':>10} {'step|rec':>10} "
          f"{'step|fail':>10} | {'AUC step':>9} {'AUC -step':>10}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff
               and not math.isnan(r['step'])]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        ms = sum(r['step'] for r in sub) / len(sub)
        mp = sum(r['step'] for r in pos) / len(pos) if pos else float('nan')
        mn = sum(r['step'] for r in neg) / len(neg) if neg else float('nan')
        if pos and neg:
            a_small = auc([r['step'] for r in pos], [r['step'] for r in neg])
            a_big = 1.0 - a_small
        else:
            a_small = a_big = float('nan')
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(pos))+'/'+str(len(sub)):>8} | {ms:>10.4f} "
              f"{mp:>10.4f} {mn:>10.4f} | {a_small:>9.4f} {a_big:>10.4f}")

    allr = [r for r in rows if not math.isnan(r['step'])]
    p_ = [r for r in allr if r['ok']]
    n_ = [r for r in allr if not r['ok']]
    print(f"\npooled AUC(-step) = {auc([r['step'] for r in p_], [r['step'] for r in n_]):.4f}"
          f"   AUC(+step) = {1.0-auc([r['step'] for r in p_], [r['step'] for r in n_]):.4f}")
    print(f"|step| test: AUC(-|step|) = "
          f"{auc([abs(r['step']) for r in p_], [abs(r['step']) for r in n_]):.4f}"
          f"  (W1b predicts step -> 0 at the wall)")
    print(f"Spearman(step, NU)     = {spearman([r['step'] for r in allr], [r['NU'] for r in allr]):+.4f}")
    print(f"Spearman(step, mu)     = {spearman([r['step'] for r in allr], [r['mu'] for r in allr]):+.4f}")
    print(f"Spearman(step, nu_hat) = {spearman([r['step'] for r in allr], [r['nuhat'] for r in allr]):+.4f}")
    print(f"Spearman(step, eff)    = {spearman([r['step'] for r in allr], [r['eff'] for r in allr]):+.4f}")

    if main is not None:
        print("\nstep inside the pre-registered NU band:")
        band_report(rows, NU_LO, NU_HI, "  (repeat)")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X6: is mu's power mediated by NU?  Spearman(mu, NU) overall and")
    print("        within each NU band.")
    print("-" * 78)
    print(f"overall  Spearman(mu, NU) = "
          f"{spearman([r['mu'] for r in rows], [r['NU'] for r in rows]):+.4f}")
    for nm, g in (("NU < 1.040", lowb), ("band", midb), ("NU > 2.199", hib)):
        if len(g) > 3:
            print(f"  {nm:<12} N={len(g):>4}  Spearman(mu, NU) = "
                  f"{spearman([r['mu'] for r in g], [r['NU'] for r in g]):+.4f}"
                  f"   Spearman(nu_hat, NU) = "
                  f"{spearman([r['nuhat'] for r in g], [r['NU'] for r in g]):+.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
