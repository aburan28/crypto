"""
GLV-HNP Phase 2, Thread 25: is mu a genuine SECOND coordinate once NU is
held fixed, or is its power in W5 entirely mediated by NU?

Pre-registered by the 2026-08-07 (Thread 24) log entry:

  H25: within the ambiguous NU band [1.04, 2.20] (the 17-bit V3/V4 bracket,
       glv_hnp_phase2_babai_output.txt -- sufficient NU<1.040, necessary
       NU>2.199 -- where nearest-plane gives no answer), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.

Falsifier: if AUC(-mu) inside the band collapses toward 0.5, W5's apparent
mu-power was entirely mediated by NU (curves with small mu also tend to have
small NU within a stratum) and the closed form should be retired as a
stratification artifact.

Secondary (W1b follow-up): step = log2(||b*_{m+1}||) - log2(||b*_1||), the
two-block GS-profile step Thread 24/W1b found vanishing right at the K1
wall. Test whether step predicts recovery better than NU or mu, globally
and inside the ambiguous band.

Data: reuses the 500-instance 17-bit table from
glv_hnp_phase2_gsprofile_strat.py (5 eff strata x 20 curves x 5 seeds, dim
24, float GS -- W0 established float GS is safe to ~1e-15 relative error at
this dimension). No new curve search or lattice reduction.

Run: python3 glv_hnp_phase2_thread25.py [--from-json]
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman
from glv_hnp_phase2_gsprofile_strat import collect_rows

JSON_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "glv_hnp_phase2_gsprofile_strat_rows.json")

NU_LO, NU_HI = 1.040, 2.199  # ambiguous band, from V3/V4 (17-bit, Thread 23)


def auc_ci_lo(pos, neg, n_boot=2000, seed=0):
    """Crude bootstrap 5th-percentile lower bound on AUC (percentile method)."""
    import random
    rng = random.Random(seed)
    vals = []
    for _ in range(n_boot):
        p = [rng.choice(pos) for _ in pos]
        q = [rng.choice(neg) for _ in neg]
        vals.append(auc(p, q))
    vals.sort()
    return vals[int(0.05 * len(vals))]


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 -- does mu survive conditioning on NU?  (H25)")
    print("=" * 78)

    if "--from-json" in sys.argv and os.path.exists(JSON_PATH):
        rows = json.load(open(JSON_PATH))
        print(f"\nloaded {len(rows)} rows from {JSON_PATH}")
    else:
        rows, curves17 = collect_rows()
        print(f"\ncollected {len(rows)} instances fresh "
              f"(dim {rows[0]['k']})")

    m = rows[0]['k'] // 2
    for r in rows:
        r['step'] = math.log2(r['prof'][m]) - math.log2(r['prof'][0])

    print("\n" + "-" * 78)
    print(f"EXP X1: AUC(-mu -> recovery) inside the ambiguous NU band "
          f"[{NU_LO}, {NU_HI}]")
    print("-" * 78)
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    pos = [r['mu'] for r in band if r['ok']]
    neg = [r['mu'] for r in band if not r['ok']]
    a_mu = auc(pos, neg)
    lo_mu = auc_ci_lo(pos, neg)
    print(f"band size N={len(band)} (of {len(rows)} total), "
          f"recovery {len(pos)}/{len(band)}")
    print(f"AUC(-mu -> recovery) in-band   = {a_mu:.4f}  "
          f"(bootstrap 5th pct >= {lo_mu:.4f})")

    pos_nh = [r['nuhat'] for r in band if r['ok']]
    neg_nh = [r['nuhat'] for r in band if not r['ok']]
    print(f"AUC(-nu_hat -> recovery) in-band = {auc(pos_nh, neg_nh):.4f}")
    pos_st = [r['step'] for r in band if r['ok']]
    neg_st = [r['step'] for r in band if not r['ok']]
    print(f"AUC(-step -> recovery) in-band   = {auc(pos_st, neg_st):.4f}")
    pos_NU = [r['NU'] for r in band if r['ok']]
    neg_NU = [r['NU'] for r in band if not r['ok']]
    print(f"AUC(-NU -> recovery) in-band     = {auc(pos_NU, neg_NU):.4f}  "
          f"(sanity check: NU is ~constant-range by construction of the "
          f"band, so this should be near 0.5)")

    print(f"\nH25 verdict: AUC(-mu) in-band = {a_mu:.4f} "
          f"{'>= 0.8 -> H25 HOLDS' if a_mu >= 0.8 else '< 0.8 -> H25 FALSIFIED'}")

    print("\n" + "-" * 78)
    print("EXP X2: sub-band scan -- does mu's power hold at finer NU cuts?")
    print("-" * 78)
    cuts = [(1.04, 1.4), (1.4, 1.8), (1.8, 2.199)]
    print(f"{'NU range':>14} {'N':>5} {'rec':>7} {'AUC mu':>8} "
          f"{'AUC nu_hat':>10} {'AUC step':>9}")
    for lo, hi in cuts:
        sub = [r for r in rows if lo <= r['NU'] < hi]
        p = [r['mu'] for r in sub if r['ok']]
        q = [r['mu'] for r in sub if not r['ok']]
        pn = [r['nuhat'] for r in sub if r['ok']]
        qn = [r['nuhat'] for r in sub if not r['ok']]
        ps = [r['step'] for r in sub if r['ok']]
        qs = [r['step'] for r in sub if not r['ok']]
        if not p or not q:
            print(f"[{lo:>5.2f},{hi:>5.2f}) {len(sub):>5} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {'(degenerate)':>8}")
            continue
        print(f"[{lo:>5.2f},{hi:>5.2f}) {len(sub):>5} "
              f"{str(len(p))+'/'+str(len(sub)):>7} {auc(p, q):>8.4f} "
              f"{auc(pn, qn):>10.4f} {auc(ps, qs):>9.4f}")

    print("\n" + "-" * 78)
    print("EXP X3 (secondary): global comparison of NU, mu, and step as "
          "predictors")
    print("-" * 78)
    allpos = [r for r in rows if r['ok']]
    allneg = [r for r in rows if not r['ok']]
    for key, label in (('NU', 'NU'), ('mu', 'mu'), ('step', 'step'),
                        ('nuhat', 'nu_hat')):
        p = [r[key] for r in allpos]
        q = [r[key] for r in allneg]
        print(f"AUC(-{label:<7} -> recovery) global = {auc(p, q):.4f}")
    print(f"\nSpearman(step, NU)     = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"Spearman(step, mu)     = "
          f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")
    print(f"Spearman(step, log mu) = "
          f"{spearman([r['step'] for r in rows], [math.log(r['mu']) for r in rows]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
