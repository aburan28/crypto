"""
GLV-HNP Phase 2, Thread 25: is mu = lambda_1(L2) a SECOND coordinate, or is
its apparent power in W5 entirely mediated by NU?

Pre-registered by the 2026-08-07 (Thread 24, run #2) log entry.  W5/W6 found
that NU (exact BDD certificate, sound nearest-plane guarantee) and
nu_hat = mu/sqrt(det L2) (closed form, no lattice work) are mutually
uncorrelated at fixed bias strength eff, yet both predict Kannan-LLL recovery
cross-curve.  The open question is whether that is two independent
mechanisms or a single one that NU only partially captures.

  H25: within the ambiguous NU band 1.04 <= NU <= 2.20 (17-bit bracket from
       W4, where the nearest-plane sufficient/necessary certificate gives no
       answer either way), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.
  Falsifier: if AUC(-mu) inside the band collapses towards 0.5, mu's power
       is mediated by NU and the W5 result is a stratification artifact of
       pooling across eff.

Secondary: W1b showed the GS profile of L0 is m exact copies of
lambda_1(L2), then a second block that moves with K1, and the step between
them shrinks as the K1 wall is crossed.  Define
    step = log2(||b*_m||) - log2(||b*_0||)     (0-indexed; ||b*_0|| = mu)
and test whether step predicts recovery, overall and inside the NU band.

Run: python3 glv_hnp_phase2_thread25.py [--dump-json FILE]
"""

import argparse
import json
import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

NU_BAND = (1.040, 2.199)   # 17-bit sufficient/necessary bracket, W4 (2026-08-07 #2)

FIELDS = ('n', 'K1', 'eff', 'effq', 'seed', 'NU', 'mu', 'nuhat', 'step', 'ok')


def collect(dump_path=None):
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")
    M17 = 12
    EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)

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
                step = math.log2(r['prof'][M17]) - math.log2(r['prof'][0])
                rows.append({
                    'n': n, 'K1': k1b, 'eff': k1b * k2b / n, 'effq': eff,
                    'seed': seed, 'NU': r['NU'], 'mu': r['mu'],
                    'nuhat': r['nuhat'], 'step': step, 'ok': bool(rk['ok']),
                })
    print(f"{len(rows)} instances (float GS, dim {2*M17}) in "
          f"{time.time()-t0:.1f}s")

    if dump_path:
        with open(dump_path, 'w') as f:
            json.dump(rows, f)
        print(f"dumped {len(rows)} rows to {dump_path}")
    return rows


def analyze(rows):
    print("\n" + "=" * 78)
    print(f"EXP T25a: H25 — does mu separate INSIDE the ambiguous NU band "
          f"{NU_BAND}?")
    print("=" * 78)

    band = [r for r in rows if NU_BAND[0] <= r['NU'] <= NU_BAND[1]]
    below = [r for r in rows if r['NU'] < NU_BAND[0]]
    above = [r for r in rows if r['NU'] > NU_BAND[1]]
    print(f"rows total {len(rows)}  |  NU<band {len(below)}  "
          f"in-band {len(band)}  NU>band {len(above)}")
    if below:
        print(f"  NU<band: recovery rate "
              f"{sum(r['ok'] for r in below)}/{len(below)} "
              f"(expect ~all recover; nearest-plane sufficient region)")
    if above:
        print(f"  NU>band: recovery rate "
              f"{sum(r['ok'] for r in above)}/{len(above)} "
              f"(expect ~none; nearest-plane necessary region)")

    for name, sub in (('in-band', band), ('NU<band', below), ('NU>band', above),
                       ('pooled', rows)):
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{name:>9}  N={len(sub):>4}  rec {len(pos)}/{len(sub)}  "
                  f"(degenerate, AUC undefined)")
            continue
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"{name:>9}  N={len(sub):>4}  rec {len(pos)}/{len(sub):<4}  "
              f"AUC(-mu) {a_mu:.4f}  AUC(-nu_hat) {a_nh:.4f}  "
              f"AUC(-NU) {a_nu:.4f}  AUC(-step) {a_st:.4f}")

    print("\nVerdict: H25 holds iff in-band AUC(-mu) >= 0.80.")

    print("\n" + "-" * 78)
    print("EXP T25a-control: within-band AUC(-mu) per eff stratum "
          "(is it just re-reading eff again?)")
    print("-" * 78)
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in band if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"  eff={eff:.2f}  in-band N={len(sub)}  "
                  f"rec {len(pos)}/{len(sub)}  (degenerate)")
            continue
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        print(f"  eff={eff:.2f}  in-band N={len(sub)}  "
              f"rec {len(pos)}/{len(sub)}  AUC(-mu) {a_mu:.4f}")

    print("\n" + "=" * 78)
    print("EXP T25b: does the GS-profile step predict recovery?")
    print(" step = log2||b*_m|| - log2||b*_0||  (0 at the exact wall, by W1b)")
    print("=" * 78)
    pos = [r for r in rows if r['ok']]
    neg = [r for r in rows if not r['ok']]
    print(f"pooled (N={len(rows)}): AUC(-step -> recovery) = "
          f"{auc([r['step'] for r in pos], [r['step'] for r in neg]):.4f}")
    print(f"pooled: Spearman(step, NU) = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"pooled: Spearman(step, mu) = "
          f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")
    print(f"in-band (N={len(band)}): AUC(-step -> recovery) = "
          f"{auc([r['step'] for r in band if r['ok']], [r['step'] for r in band if not r['ok']]):.4f}"
          if band and any(r['ok'] for r in band) and any(not r['ok'] for r in band)
          else "in-band: degenerate")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None)
    ap.add_argument("--load-json", default=None,
                     help="skip generation, re-analyze a prior --dump-json")
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — conditioning on NU: is mu a second coordinate?")
    print("=" * 78)

    if args.load_json:
        with open(args.load_json) as f:
            rows = json.load(f)
        print(f"loaded {len(rows)} rows from {args.load_json}")
    else:
        rows = collect(args.dump_json)

    analyze(rows)
