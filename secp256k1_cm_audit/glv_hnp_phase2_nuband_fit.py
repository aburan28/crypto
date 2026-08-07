"""
GLV-HNP Phase 2, Thread 25b: the 2-parameter viability test.

Thread 25 (glv_hnp_phase2_nuband.py) conditioned recovery on NU instead of on
eff and found that the pre-registered H25 splits:

    inside the band 1.040 <= NU <= 2.199   AUC(-mu)     = 0.693   (H25 FAILS)
                                           AUC(-nu_hat) = 0.840   (clears 0.8)
                                           AUC(+step)   = 0.820   (clears 0.8)

so the second coordinate is real but it is the SIZE-NORMALISED nu_hat =
mu/sqrt(det L2), not the raw lambda_1(L2) = mu that Thread 24's W5 nominated.
W5 could not see the difference because it held eff fixed, which pins det L2
within a stratum; conditioning on NU lets eff vary and separates the two.

This script re-analyses the dumped table (no recomputation) to produce the
deliverable the pre-registration asked for: a logistic fit on (log NU,
log nu_hat) with its decision boundary in raw units, its grouped-CV score,
and an operating point with a confusion table.

Run: python3 glv_hnp_phase2_nuband_fit.py [table.json]
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_nuband import (
    NU_LO, NU_HI, logistic_fit, logistic_score, standardize, auc_from_scores,
)
from glv_hnp_phase2_gsprofile import auc, spearman

FEATS = {
    'log NU': lambda r: [math.log(r['NU'])],
    'log nu_hat': lambda r: [math.log(r['nuhat'])],
    'log mu': lambda r: [math.log(r['mu'])],
    'step': lambda r: [r['step']],
    'log NU + log mu': lambda r: [math.log(r['NU']), math.log(r['mu'])],
    'log NU + log nu_hat': lambda r: [math.log(r['NU']), math.log(r['nuhat'])],
    'log NU + step': lambda r: [math.log(r['NU']), r['step']],
    'log NU + log nu_hat + step':
        lambda r: [math.log(r['NU']), math.log(r['nuhat']), r['step']],
}


def grouped_cv(rows, feat):
    """Leave-one-curve-out CV, grouping by n."""
    preds, ys = [], []
    for gname in sorted({r['n'] for r in rows}):
        tr = [r for r in rows if r['n'] != gname]
        te = [r for r in rows if r['n'] == gname]
        Xtr = [feat(r) for r in tr]
        ytr = [1 if r['ok'] else 0 for r in tr]
        _, mus, sds = standardize(Xtr)
        Ztr = [[(v - mus[j]) / sds[j] for j, v in enumerate(x)] for x in Xtr]
        beta, _ = logistic_fit(Ztr, ytr)
        for r in te:
            z = [(v - mus[j]) / sds[j] for j, v in enumerate(feat(r))]
            preds.append(logistic_score(beta, z))
            ys.append(1 if r['ok'] else 0)
    return preds, ys


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "glv_hnp_phase2_nuband_table.json"
    with open(path) as fh:
        rows = json.load(fh)
    rows = [r for r in rows if not math.isnan(r.get('step', float('nan')))]
    y = [1 if r['ok'] else 0 for r in rows]

    print("=" * 78)
    print("Thread 25b — the 2-parameter viability test (NU, nu_hat)")
    print("=" * 78)
    print(f"table {path}: {len(rows)} instances, {sum(y)} recovered")

    print("\n" + "-" * 78)
    print("EXP Y1: in-sample fit and leave-one-curve-out CV, all candidates")
    print("-" * 78)
    print(f"  {'features':<28} {'dev':>8} {'AUC':>7} {'acc':>6} | "
          f"{'CV-AUC':>7} {'CV-acc':>7}")
    for name, feat in FEATS.items():
        X = [feat(r) for r in rows]
        Z, mus, sds = standardize(X)
        beta, dev = logistic_fit(Z, y)
        sc = [logistic_score(beta, z) for z in Z]
        a = auc_from_scores(sc, y)
        acc = sum(1 for s, t in zip(sc, y) if (s > 0) == bool(t)) / len(y)
        cp, cy = grouped_cv(rows, feat)
        ca = auc_from_scores(cp, cy)
        cacc = sum(1 for s, t in zip(cp, cy) if (s > 0) == bool(t)) / len(cy)
        print(f"  {name:<28} {dev:>8.2f} {a:>7.4f} {acc:>6.3f} | "
              f"{ca:>7.4f} {cacc:>7.3f}")

    print("\n" + "-" * 78)
    print("EXP Y2: likelihood-ratio tests against the log-NU-only baseline")
    print("-" * 78)
    Z0, _, _ = standardize([FEATS['log NU'](r) for r in rows])
    _, dev0 = logistic_fit(Z0, y)
    print(f"  baseline (log NU) deviance = {dev0:.2f}")
    for name in ('log NU + log mu', 'log NU + log nu_hat', 'log NU + step',
                 'log NU + log nu_hat + step'):
        X = [FEATS[name](r) for r in rows]
        Z, _, _ = standardize(X)
        _, dev = logistic_fit(Z, y)
        df = len(X[0]) - 1
        lr = dev0 - dev
        print(f"  + {name:<28} LR {lr:>7.2f} on {df} df   "
              f"({'SIGNIFICANT' if lr > 10.83 else 'significant' if lr > 3.84 else 'ns'})")

    print("\n" + "-" * 78)
    print("EXP Y3: the boundary, in raw units")
    print("-" * 78)
    X = [FEATS['log NU + log nu_hat'](r) for r in rows]
    Z, mus, sds = standardize(X)
    beta, dev = logistic_fit(Z, y)
    b0, bn, bh = beta
    cn, ch = bn / sds[0], bh / sds[1]
    c0 = b0 - bn * mus[0] / sds[0] - bh * mus[1] / sds[1]
    print(f"  logit p(recover) = {cn:+.4f}*log(NU) {ch:+.4f}*log(nu_hat) {c0:+.4f}")
    print(f"  boundary p=1/2   : log(NU) {ch/cn:+.4f}*log(nu_hat) {c0/cn:+.4f} = 0")
    print(f"  i.e.  NU * nu_hat^{ch/cn:.3f} = {math.exp(-c0/cn):.5g}")
    print(f"  equivalently: recover iff  NU * nu_hat^{ch/cn:.3f} < "
          f"{math.exp(-c0/cn):.5g}")

    T = [math.log(r['NU']) + (ch / cn) * math.log(r['nuhat']) for r in rows]
    thr = -c0 / cn
    tp = sum(1 for t, yy in zip(T, y) if t < thr and yy)
    fp = sum(1 for t, yy in zip(T, y) if t < thr and not yy)
    fn = sum(1 for t, yy in zip(T, y) if t >= thr and yy)
    tn = sum(1 for t, yy in zip(T, y) if t >= thr and not yy)
    print(f"\n  confusion at p=1/2:  TP {tp}  FP {fp}  FN {fn}  TN {tn}"
          f"   acc {(tp+tn)/len(y):.3f}")
    print(f"  AUC of the combined statistic T = log NU + {ch/cn:.3f}*log nu_hat"
          f"  : {auc([t for t, yy in zip(T, y) if yy], [t for t, yy in zip(T, y) if not yy]):.4f}")

    print("\n  conservative operating points on T (sweep):")
    print(f"  {'thr':>8} {'TP':>5} {'FP':>5} {'FN':>5} {'TN':>5} {'prec':>7} {'rec':>7}")
    lo, hi = min(T), max(T)
    for i in range(11):
        t0 = lo + (hi - lo) * i / 10.0
        tp = sum(1 for t, yy in zip(T, y) if t < t0 and yy)
        fp = sum(1 for t, yy in zip(T, y) if t < t0 and not yy)
        fn = sum(1 for t, yy in zip(T, y) if t >= t0 and yy)
        tn = sum(1 for t, yy in zip(T, y) if t >= t0 and not yy)
        prec = tp / (tp + fp) if tp + fp else float('nan')
        rec = tp / (tp + fn) if tp + fn else float('nan')
        print(f"  {t0:>8.4f} {tp:>5} {fp:>5} {fn:>5} {tn:>5} "
              f"{prec:>7.3f} {rec:>7.3f}")

    print("\n" + "-" * 78)
    print("EXP Y4: orthogonality of the two coordinates (Spearman)")
    print("-" * 78)
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    for nm, g in (("all", rows), ("inside NU band", band)):
        print(f"  {nm:<16} N={len(g):>4}  "
              f"rho(NU, nu_hat) = {spearman([r['NU'] for r in g], [r['nuhat'] for r in g]):+.4f}   "
              f"rho(NU, mu) = {spearman([r['NU'] for r in g], [r['mu'] for r in g]):+.4f}   "
              f"rho(NU, step) = {spearman([r['NU'] for r in g], [r['step'] for r in g]):+.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
