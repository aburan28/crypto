"""
GLV-HNP Phase 2, Thread 25: is mu a SECOND coordinate, or is its apparent
power mediated by NU after all?

Pre-registered by the 2026-08-07 (Thread 24) log entry:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where the nearest-plane
       certificate gives no answer), AUC(-mu -> Kannan-LLL recovery) stays
       >= 0.8.
  Falsifier: if AUC inside the band collapses to ~0.5, mu's apparent power is
       entirely mediated by NU, the W5 result is a stratification artifact,
       and the closed form should be retired.

Secondary, from W1b: the GS profile of L0 is m exact copies of lambda_1(L2)
followed by a second block, and the step between the blocks *vanishes* right
as the K1 wall is crossed.  Define

    step = log2 ||b*_{m+1}||  -  log2 ||b*_1||          (0-indexed: prof[m]/prof[0])

and test whether step predicts the wall better than either NU or mu.

Tertiary: if H25 holds, fit a logistic decision boundary on the pair and
report cross-validated AUC against each coordinate alone.

Data: the same 5 x 20 x 5 = 500-instance 17-bit stratified design as
glv_hnp_phase2_gsprofile_strat.py (identical curves, bounds and seeds), plus
`step` and a --dump-json flag so the table survives the run.

Run: python3 glv_hnp_phase2_conditional.py [--dump-json PATH]
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

M17 = 12
EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)

# The Thread 24 W4 bracket at 17 bits: NU < 1.040 was sufficient (0 FP over
# 300 instances), NU > 2.199 was necessary.  In between, nearest-plane says
# nothing.
BAND_LO, BAND_HI = 1.040, 2.199


# ---------------------------------------------------------------------------

def build_rows():
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")

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
                step = (math.log2(prof[M17]) - math.log2(prof[0])
                        if prof[0] > 0 and prof[M17] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step,
                          'prof0': prof[0], 'profm': prof[M17]})
                # the full profile / nu vector is large and not needed again
                r.pop('nus', None)
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")
    return rows


def split(rows):
    return ([r for r in rows if r['ok']], [r for r in rows if not r['ok']])


def auc_of(rows, key, sign=1.0):
    """AUC for 'smaller sign*key -> recovery'.  nan if a class is empty."""
    pos, neg = split(rows)
    if not pos or not neg:
        return float('nan')
    return auc([sign * r[key] for r in pos], [sign * r[key] for r in neg])


def fmt(x, w=8, p=4):
    return f"{'--':>{w}}" if x != x else f"{x:>{w}.{p}f}"


# ---------------------------------------------------------------------------
# plain logistic regression (standardised features, gradient ascent on LL)
# ---------------------------------------------------------------------------

def logistic_fit(X, y, iters=20000, lr=0.2):
    d = len(X[0])
    n = len(X)
    w = [0.0] * d
    b = 0.0
    for _ in range(iters):
        gw = [0.0] * d
        gb = 0.0
        for xi, yi in zip(X, y):
            z = b + sum(wj * xj for wj, xj in zip(w, xi))
            pz = 1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, z))))
            e = yi - pz
            for j in range(d):
                gw[j] += e * xi[j]
            gb += e
        for j in range(d):
            w[j] += lr * gw[j] / n
        b += lr * gb / n
    return w, b


def logistic_score(w, b, xi):
    return b + sum(wj * xj for wj, xj in zip(w, xi))


def standardise(X):
    d = len(X[0])
    mu = [sum(r[j] for r in X) / len(X) for j in range(d)]
    sd = [math.sqrt(sum((r[j] - mu[j]) ** 2 for r in X) / len(X)) or 1.0
          for j in range(d)]
    return [[(r[j] - mu[j]) / sd[j] for j in range(d)] for r in X], mu, sd


def cv_auc(X, y, folds=5, seed=11):
    """Stratified-by-order k-fold CV AUC of the logistic score."""
    idx = list(range(len(X)))
    random.Random(seed).shuffle(idx)
    scores = [0.0] * len(X)
    for f in range(folds):
        te = [i for k, i in enumerate(idx) if k % folds == f]
        tr = [i for k, i in enumerate(idx) if k % folds != f]
        if len({y[i] for i in tr}) < 2:
            return float('nan')
        Xtr, mu, sd = standardise([X[i] for i in tr])
        w, b = logistic_fit(Xtr, [y[i] for i in tr])
        for i in te:
            xi = [(X[i][j] - mu[j]) / sd[j] for j in range(len(mu))]
            scores[i] = logistic_score(w, b, xi)
    pos = [-scores[i] for i in range(len(X)) if y[i] == 1]
    neg = [-scores[i] for i in range(len(X)) if y[i] == 0]
    if not pos or not neg:
        return float('nan')
    return auc(pos, neg)


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    dump = None
    if "--dump-json" in sys.argv:
        dump = sys.argv[sys.argv.index("--dump-json") + 1]

    print("=" * 78)
    print("Thread 25 — does mu separate INSIDE an NU band?")
    print("=" * 78)
    print()

    rows = build_rows()

    if dump:
        with open(dump, "w") as fh:
            json.dump([{k: v for k, v in r.items() if k != 'prof'}
                       for r in rows], fh)
        print(f"dumped {len(rows)} rows -> {dump}")

    pos, neg = split(rows)
    print(f"\noverall recovery {len(pos)}/{len(rows)}")

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W8: the ambiguous band.  NU in [%.3f, %.3f] is exactly where"
          % (BAND_LO, BAND_HI))
    print("        nearest-plane gives no answer.  Does anything else speak?")
    print("-" * 78)

    band = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    bpos, bneg = split(band)
    print(f"\nband occupancy {len(band)}/{len(rows)}   recovery "
          f"{len(bpos)}/{len(band)}")
    print(f"\n{'score':>12} {'AUC in band':>12} {'AUC all':>9}   note")
    for key, sign, note in (
        ('mu',      1.0, 'H25 target: >= 0.80 confirms a 2nd coordinate'),
        ('nuhat',   1.0, 'size-free mu'),
        ('step',   -1.0, 'W1b: larger step -> recovery'),
        ('NU',      1.0, 'residual power of NU inside its own band'),
        ('eff',     1.0, 'CONTROL: band may just re-select bias strength'),
        ('lamstar', 1.0, 'CONTROL: falsified by Thread 20'),
    ):
        print(f"{key:>12} {fmt(auc_of(band, key, sign), 12)} "
              f"{fmt(auc_of(rows, key, sign), 9)}   {note}")

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W8b: same, but WITHIN each eff stratum, so the band cannot be")
    print("         re-selecting bias strength.")
    print("-" * 78)
    print(f"\n{'eff':>5} {'N band':>7} {'rec':>8} | {'AUC mu':>8} "
          f"{'AUC step':>9} {'AUC NU':>8} {'AUC lam*':>9}")
    for eff in EFFS:
        sub = [r for r in rows
               if r['effq'] == eff and BAND_LO <= r['NU'] <= BAND_HI]
        p, ng = split(sub)
        print(f"{eff:>5.2f} {len(sub):>7} {str(len(p))+'/'+str(len(sub)):>8} | "
              f"{fmt(auc_of(sub, 'mu'), 8)} {fmt(auc_of(sub, 'step', -1.0), 9)} "
              f"{fmt(auc_of(sub, 'NU'), 8)} {fmt(auc_of(sub, 'lamstar'), 9)}")

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W9: `step` as a predictor in its own right, per stratum.")
    print("        step = log2||b*_{m+1}|| - log2||b*_1||,  m = %d" % M17)
    print("-" * 78)
    print(f"\n{'eff':>5} {'N':>5} {'rec':>8} | {'AUC step':>9} {'AUC mu':>8} "
          f"{'AUC NU':>8} | {'step|ok':>9} {'step|no':>9}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        p, ng = split(sub)
        mp = sum(r['step'] for r in p) / len(p) if p else float('nan')
        mn = sum(r['step'] for r in ng) / len(ng) if ng else float('nan')
        print(f"{eff:>5.2f} {len(sub):>5} {str(len(p))+'/'+str(len(sub)):>8} | "
              f"{fmt(auc_of(sub, 'step', -1.0), 9)} {fmt(auc_of(sub, 'mu'), 8)} "
              f"{fmt(auc_of(sub, 'NU'), 8)} | {fmt(mp, 9, 3)} {fmt(mn, 9, 3)}")

    print("\nrank correlations (pooled, N=%d):" % len(rows))
    for a, bq in (('step', 'NU'), ('step', 'mu'), ('step', 'eff'),
                  ('mu', 'NU'), ('nuhat', 'NU')):
        print(f"  Spearman({a:>5}, {bq:<3}) = "
              f"{spearman([r[a] for r in rows], [r[bq] for r in rows]):+.4f}")

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W10: 2-parameter viability test.  Cross-validated (5-fold) AUC")
    print("         of a logistic fit, vs each coordinate alone.")
    print("-" * 78)

    def logf(r, k):
        v = r[k]
        return math.log(v) if v > 0 else -60.0

    y = [1 if r['ok'] else 0 for r in rows]
    combos = (
        (('NU',),            'NU alone'),
        (('mu',),            'mu alone'),
        (('nuhat',),         'nu_hat alone'),
        (('NU', 'mu'),       'the pair proposed by Thread 24'),
        (('NU', 'nuhat'),    'size-free variant of the pair'),
        (('NU', 'nuhat', 'step'), 'plus the GS step'),
    )
    print(f"\n{'features':>26} {'CV AUC':>8}   note")
    for keys, note in combos:
        X = [[logf(r, k) if k != 'step' else r['step'] for k in keys]
             for r in rows]
        print(f"{'(' + ', '.join(keys) + ')':>26} "
              f"{fmt(cv_auc(X, y), 8)}   {note}")

    print("\nin-sample decision boundary on (log NU, log nu_hat), "
          "standardised features:")
    X = [[logf(r, 'NU'), logf(r, 'nuhat')] for r in rows]
    Xs, mus, sds = standardise(X)
    w, b0 = logistic_fit(Xs, y)
    print(f"  z = {b0:+.4f} {w[0]:+.4f}*u1 {w[1]:+.4f}*u2   "
          f"(recover when z > 0)")
    print(f"  u1 = (log NU     - {mus[0]:+.4f}) / {sds[0]:.4f}")
    print(f"  u2 = (log nu_hat - {mus[1]:+.4f}) / {sds[1]:.4f}")
    acc = sum(1 for xi, yi in zip(Xs, y)
              if (logistic_score(w, b0, xi) > 0) == (yi == 1)) / len(y)
    print(f"  in-sample accuracy {acc:.4f}  "
          f"(base rate {max(sum(y), len(y)-sum(y))/len(y):.4f})")

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W14: the two coefficients above are nearly equal once the")
    print("         standardisation is undone.  If the ratio is 1 the")
    print("         boundary is not a fitted plane but the single product")
    print("         PI = NU * nu_hat.")
    print("-" * 78)
    c1, c2 = w[0] / sds[0], w[1] / sds[1]
    inter = b0 - w[0] * mus[0] / sds[0] - w[1] * mus[1] / sds[1]
    print(f"\n  z = {inter:+.4f} {c1:+.4f}*log NU {c2:+.4f}*log nu_hat")
    print(f"  exponent ratio c1/c2 = {c1/c2:.4f}     "
          f"(1.0 => a pure product)")
    cut = math.exp(inter / (-(c1 + c2) / 2.0))
    print(f"  equal-exponent boundary:  PI = NU*nu_hat < {cut:.4f}")

    PI = [r['NU'] * r['nuhat'] for r in rows]
    ppos = [p for p, r in zip(PI, rows) if r['ok']]
    pneg = [p for p, r in zip(PI, rows) if not r['ok']]
    print(f"\n  AUC(-PI) pooled, N={len(rows)}: {auc(ppos, pneg):.4f}"
          f"   [NU alone {auc_of(rows,'NU'):.4f}, "
          f"nu_hat alone {auc_of(rows,'nuhat'):.4f}]")
    print("  PI is threshold-free and unfitted, so this AUC is directly")
    print("  comparable to the 5-fold CV numbers of W10.")

    print(f"\n  fixed cut sweep (pooled):")
    print(f"  {'cut':>6} {'TP':>4} {'FP':>4} {'FN':>4} {'TN':>4} "
          f"{'acc':>7} {'Youden J':>9}")
    for t in (0.80, 0.90, cut, 1.00, 1.10, 1.20):
        tp = sum(1 for p, r in zip(PI, rows) if p < t and r['ok'])
        fp = sum(1 for p, r in zip(PI, rows) if p < t and not r['ok'])
        fn = sum(1 for p, r in zip(PI, rows) if p >= t and r['ok'])
        tn = sum(1 for p, r in zip(PI, rows) if p >= t and not r['ok'])
        print(f"  {t:>6.2f} {tp:>4} {fp:>4} {fn:>4} {tn:>4} "
              f"{(tp+tn)/len(rows):>7.4f} "
              f"{tp/(tp+fn) - fp/(fp+tn):>9.3f}")

    print(f"\n  per-stratum, at the fixed cut PI < {cut:.2f}:")
    print(f"  {'eff':>5} {'N':>5} {'AUC(-PI)':>9} | {'TP':>4} {'FP':>4} "
          f"{'FN':>4} {'TN':>4} {'acc':>7}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        sp = [r['NU'] * r['nuhat'] for r in sub if r['ok']]
        sn = [r['NU'] * r['nuhat'] for r in sub if not r['ok']]
        a = auc(sp, sn) if sp and sn else float('nan')
        tp = sum(1 for r in sub if r['NU'] * r['nuhat'] < cut and r['ok'])
        fp = sum(1 for r in sub if r['NU'] * r['nuhat'] < cut and not r['ok'])
        fn = sum(1 for r in sub if r['NU'] * r['nuhat'] >= cut and r['ok'])
        tn = sum(1 for r in sub
                 if r['NU'] * r['nuhat'] >= cut and not r['ok'])
        print(f"  {eff:>5.2f} {len(sub):>5} {fmt(a, 9)} | {tp:>4} {fp:>4} "
              f"{fn:>4} {tn:>4} {(tp+tn)/len(sub):>7.3f}")

    print("\n  CAVEAT: the exponent ratio and the cut are both read off THESE")
    print("  500 rows.  The product FORM is the structural claim; the cut")
    print("  value is in-sample and must be re-measured at another bit")
    print("  length before it is quoted as a constant.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
