"""
GLV-HNP Phase 2, Thread 25b: the lambda_2 law.

Thread 25's conditioning experiment (glv_hnp_phase2_nu_conditioned.py) ended
somewhere its own pre-registered hypothesis did not point.  H25 asked whether
mu = lambda_1(L2) is a second coordinate of the LLL wall alongside NU.  It is
a second coordinate, but it is the WRONG one: mu's pooled rank-AUC over the
whole 1800-instance table is 0.357, i.e. INVERTED, while inside every fixed-eff
stratum it runs 0.70-0.92.  That is a Simpson reversal, and it means mu is
uninterpretable unless it is normalised across curves.

nu_hat = mu / sqrt(det L2) is Thread 20's normalisation.  Fitting the exponent
  score(alpha) = mu / det(L2)^alpha
freely gives alpha ~ 0.88-0.89 on train seeds, and alpha = 1 does better still
on held-out seeds.  Since Thread 24 W2 established det(L2) = mu * lambda_2 to
within 4.3%, alpha = 1 says the statistic is simply

        1 / lambda_2(L2)        i.e.  LARGER lambda_2 -> recovery

This script re-derives that from the dumped table and runs the controls that
decide whether it is real:

  Y1  alpha sweep, held-out
  Y2  lambda_2 as a standalone score, pooled and per-stratum
  Y3  CONTROL — is lambda_2 just a proxy for the curve order n, or for det L2
      alone?  If AUC(n) or AUC(det2) matches AUC(lambda_2), the law is vacuous.
  Y4  the two-coordinate model (NU, lambda_2) and what anything else adds
  Y5  variance decomposition: which coordinate lives at which level

Reads the JSON dumped by glv_hnp_phase2_nu_conditioned.py --dump-json.
Run: python3 glv_hnp_phase2_lambda2_law.py [table.json]
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman
from glv_hnp_phase2_nu_conditioned import logistic_fit, auc_p, EFFS

DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "glv_hnp_phase2_nu_conditioned_table.json")


def feat(r, k):
    if k == 'step':
        return r['step']
    if k == 'eff':
        return r['eff']
    return math.log(r[k])


def split(rows):
    seeds = sorted({r['seed'] for r in rows})
    tr = set(seeds[0::2])
    return ([r for r in rows if r['seed'] in tr],
            [r for r in rows if r['seed'] not in tr])


def model(rows, A, B, keys):
    y = [1.0 if r['ok'] else 0.0 for r in rows]
    beta, ll = logistic_fit([[feat(r, k) for k in keys] for r in rows], y)
    bt, _ = logistic_fit([[feat(r, k) for k in keys] for r in A],
                         [1.0 if r['ok'] else 0.0 for r in A])
    et = [bt[0] + sum(x * v for x, v in zip(bt[1:], [feat(r, k) for k in keys]))
          for r in B]
    ha = auc([-e for e, r in zip(et, B) if r['ok']],
             [-e for e, r in zip(et, B) if not r['ok']])
    return beta, ll, ha


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    rows = json.load(open(path))
    A, B = split(rows)
    print("=" * 78)
    print("Thread 25b — the lambda_2 law")
    print("=" * 78)
    print(f"{len(rows)} instances from {os.path.basename(path)}  "
          f"(train {len(A)} / held-out {len(B)} by seed parity)")

    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Y1: exponent sweep for score(alpha) = mu / det(L2)^alpha")
    print("    alpha fitted on TRAIN seeds only; AUC on HELD-OUT seeds.")
    print("-" * 78)
    yA = [1.0 if r['ok'] else 0.0 for r in A]
    b, _ = logistic_fit(
        [[math.log(r['NU']), math.log(r['mu']), math.log(r['det2'])] for r in A],
        yA)
    a_fit = -b[3] / b[2]
    print(f"  free-fit alpha on train = {a_fit:.4f}   (nu_hat assumes 0.5)")
    print(f"\n{'alpha':>8} {'held-out AUC':>13} {'pooled AUC':>12}   note")
    for a in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, a_fit, 0.9, 1.0, 1.1, 1.25):
        hp = [r['mu'] / r['det2'] ** a for r in B if r['ok']]
        hn = [r['mu'] / r['det2'] ** a for r in B if not r['ok']]
        pp = [r['mu'] / r['det2'] ** a for r in rows if r['ok']]
        pn = [r['mu'] / r['det2'] ** a for r in rows if not r['ok']]
        note = ("<- nu_hat" if a == 0.5 else
                "<- free fit" if a == a_fit else
                "<- 1/lambda_2" if a == 1.0 else "")
        print(f"{a:>8.4f} {auc(hp, hn):>13.4f} {auc(pp, pn):>12.4f}   {note}")

    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Y2: lambda_2(L2) as a standalone score (larger -> recover)")
    print("-" * 78)
    pos = [-r['l2'] for r in rows if r['ok']]
    neg = [-r['l2'] for r in rows if not r['ok']]
    a, p = auc_p(pos, neg)
    print(f"  pooled    AUC(+lambda_2) = {a:.4f}   p = {p:.2e}  (N={len(rows)})")
    hp = [-r['l2'] for r in B if r['ok']]
    hn = [-r['l2'] for r in B if not r['ok']]
    print(f"  held-out  AUC(+lambda_2) = {auc(hp, hn):.4f}          (N={len(B)})")
    print(f"\n{'eff':>6} {'N':>6} {'rec':>10} {'AUC +l2':>9} {'AUC -NU':>9} "
          f"{'AUC -mu':>9} {'AUC -nuhat':>11}")
    for e in EFFS:
        s = [r for r in rows if r['effq'] == e]
        ps = [r for r in s if r['ok']]
        ns = [r for r in s if not r['ok']]
        if not ps or not ns:
            print(f"{e:>6.2f} {len(s):>6} "
                  f"{str(len(ps)) + '/' + str(len(s)):>10}   (degenerate)")
            continue
        print(f"{e:>6.2f} {len(s):>6} "
              f"{str(len(ps)) + '/' + str(len(s)):>10} "
              f"{auc([-r['l2'] for r in ps], [-r['l2'] for r in ns]):>9.4f} "
              f"{auc([r['NU'] for r in ps], [r['NU'] for r in ns]):>9.4f} "
              f"{auc([r['mu'] for r in ps], [r['mu'] for r in ns]):>9.4f} "
              f"{auc([r['nuhat'] for r in ps], [r['nuhat'] for r in ns]):>11.4f}")

    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Y3: CONTROL — is lambda_2 just n, or just det(L2)?")
    print("    If a trivial size proxy matches it, the law is vacuous.")
    print("-" * 78)
    print(f"{'score':>16} {'pooled AUC':>11} {'held-out':>10}   "
          f"Spearman vs log lambda_2")
    ll2 = [math.log(r['l2']) for r in rows]
    for k, s, lab in ((('l2'), -1, 'lambda_2'), ('det2', -1, 'det(L2)'),
                      ('n', -1, 'n'), ('mu', 1, 'mu'), ('K1', 1, 'K1'),
                      ('nuhat', 1, 'nu_hat')):
        pp = [s * r[k] for r in rows if r['ok']]
        pn = [s * r[k] for r in rows if not r['ok']]
        hp = [s * r[k] for r in B if r['ok']]
        hn = [s * r[k] for r in B if not r['ok']]
        rho = spearman([math.log(r[k]) for r in rows], ll2)
        print(f"{lab:>16} {auc(pp, pn):>11.4f} {auc(hp, hn):>10.4f}   {rho:>+.4f}")
    print("\n  n is the curve order; all 20 curves are 17-bit, so n spans only")
    print("  [2^16, 2^17].  A high AUC(n) would mean the effect is size, not")
    print("  geometry.")

    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Y4: two-coordinate model and marginal value of everything else")
    print("-" * 78)
    print(f"{'model':>26} {'loglik':>10} {'held-out AUC':>13}   coefficients")
    for keys in (('NU',), ('l2',), ('nuhat',),
                 ('NU', 'nuhat'), ('NU', 'l2'),
                 ('NU', 'l2', 'mu'), ('NU', 'l2', 'step'),
                 ('NU', 'l2', 'mu', 'step', 'eff')):
        beta, ll, ha = model(rows, A, B, keys)
        print(f"{'+'.join(keys):>26} {ll:>10.2f} {ha:>13.4f}   "
              f"{' '.join(f'{x:+.3f}' for x in beta)}")
    beta, _, _ = model(rows, A, B, ('NU', 'l2'))
    if abs(beta[2]) > 1e-9:
        print(f"\n  (NU, lambda_2) 50% contour:  "
              f"lambda_2 = {math.exp(-beta[0]/beta[2]):.4g} * "
              f"NU^{-beta[1]/beta[2]:+.4f}")

    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Y5: at which level does each coordinate live?")
    print("    Cells are (curve, eff); mu, lambda_2, det2 are EXACTLY constant")
    print("    within a cell, NU and step vary with the signature seed.")
    print("-" * 78)
    cells = {}
    for r in rows:
        cells.setdefault((r['n'], r['effq']), []).append(r)
    print(f"{'quantity':>12} {'between-cell':>13} {'within-cell':>12} {'ICC':>8}")
    for k in ('l2', 'mu', 'det2', 'NU', 'step'):
        v = [feat(r, k) for r in rows]
        gm = sum(v) / len(v)
        sb = sum(len(c) * (sum(feat(x, k) for x in c) / len(c) - gm) ** 2
                 for c in cells.values())
        sw = sum((feat(x, k) - sum(feat(q, k) for q in c) / len(c)) ** 2
                 for c in cells.values() for x in c)
        print(f"{k:>12} {sb:>13.3f} {sw:>12.3f} "
              f"{sb / (sb + sw) if sb + sw else float('nan'):>8.4f}")
    live = [c for c in cells.values() if 0 < sum(1 for r in c if r['ok']) < len(c)]
    tp = 0
    cn = 0.0
    for c in live:
        ps = [r for r in c if r['ok']]
        ns = [r for r in c if not r['ok']]
        w = len(ps) * len(ns)
        tp += w
        cn += auc([r['NU'] for r in ps], [r['NU'] for r in ns]) * w
    print(f"\n  {len(live)} mixed cells, {tp} pos/neg pairs")
    print(f"  within-cell pooled AUC(-NU) = {cn/tp:.4f}  "
          f"<- NU is the seed-level coordinate")
    print(f"  lambda_2 has no within-cell variation at all "
          f"<- it is the curve-level coordinate")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
