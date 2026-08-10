"""
GLV-HNP Phase 2, Thread 25: is mu a second coordinate, or is its power in W5
entirely mediated by NU?

Pre-registered by the 2026-08-07 (Thread 24) log entry.  W5/W6 established
recovery = f(NU, X) with X ~ mu-driven and X uncorrelated with NU (Spearman
in [-0.28, +0.16] across eff-strata).  That is necessary but not sufficient
for mu to be a genuine second parameter: it is also consistent with mu doing
all its work through NU inside a stratum and the two only decorrelating
because eff moves them differently.

  H25: within the ambiguous NU band [1.04, 2.20] (Thread 23b/W4's ---
       nearest-plane gives no verdict there), AUC(-mu -> recovery) >= 0.80.
  Falsifier: AUC inside the band collapses towards 0.5 (mu's apparent power
       in W5 was a stratification artifact of eff, not real conditional
       information on top of NU).

Secondary (W1b follow-up): the GS profile of L0 is m exact copies of
lambda_1(L2) followed by a second block that starts moving with K1 right at
the point the wall crosses (Thread 24 log, W1b).  Define

    step = log2(||b*_{m+1}||) - log2(||b*_1||)

(0-indexed: prof[m] vs prof[0]) and test whether step predicts recovery
better than NU or mu alone.

Data: same generator as glv_hnp_phase2_gsprofile_strat.py (float GS, per
W0/W4 of the parent script max relative NU error ~1e-15 at dim 24), same
20 curves x 5 strata x 5 seeds = 500-instance table, so this is a
re-analysis, not a new sweep -- deterministic given (M17, EFFS, SEEDS).

Run: python3 glv_hnp_phase2_thread25.py [--dump-json out.json]
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

# From the Thread 24 log (W4, 17-bit bracket):
#   sufficient NU < 1.040 , necessary NU > 2.199
NU_BAND = (1.040, 2.199)


def collect(dump_path=None):
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    M17 = 12
    EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)

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
                r['step'] = (math.log2(r['prof'][M17]) - math.log2(r['prof'][0])
                             if r['prof'][0] > 0 and r['prof'][M17] > 0
                             else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n)})
                rows.append(r)

    if dump_path:
        slim = [{k: v for k, v in r.items() if k not in ('prof', 'nus')}
                for r in rows]
        with open(dump_path, 'w') as f:
            json.dump(slim, f)
        print(f"[dumped {len(slim)} rows -> {dump_path}]")
    return rows


def auc_ci_note(pos, neg):
    """AUC plus a crude count so a thin band doesn't read as a strong AUC."""
    return auc(pos, neg), len(pos), len(neg)


if __name__ == "__main__":
    dump_path = None
    if "--dump-json" in sys.argv:
        i = sys.argv.index("--dump-json")
        dump_path = sys.argv[i + 1] if i + 1 < len(sys.argv) else \
            "glv_hnp_phase2_thread25_rows.json"

    print("=" * 78)
    print("Thread 25 — is mu a 2nd coordinate, or is it mediated by NU?  (H25)")
    print("=" * 78)

    t0 = time.time()
    rows = collect(dump_path)
    print(f"\n{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print(f"EXP H25: AUC(-mu -> recovery) inside the ambiguous NU band "
          f"[{NU_BAND[0]}, {NU_BAND[1]}]")
    print("-" * 78)
    band = [r for r in rows if NU_BAND[0] <= r['NU'] <= NU_BAND[1]]
    print(f"band population: {len(band)}/{len(rows)} instances "
          f"({len(band)/len(rows)*100:.1f}%)")
    pos_b = [r for r in band if r['ok']]
    neg_b = [r for r in band if not r['ok']]
    if pos_b and neg_b:
        a_mu, np_, nn_ = auc_ci_note([r['mu'] for r in pos_b],
                                      [r['mu'] for r in neg_b])
        a_nh, _, _ = auc_ci_note([r['nuhat'] for r in pos_b],
                                  [r['nuhat'] for r in neg_b])
        a_nu, _, _ = auc_ci_note([r['NU'] for r in pos_b],
                                  [r['NU'] for r in neg_b])
        print(f"  rec {len(pos_b)}/{len(band)}   "
              f"AUC(-mu) = {a_mu:.4f} (n+={np_}, n-={nn_})   "
              f"AUC(-nu_hat) = {a_nh:.4f}   AUC(-NU, sanity) = {a_nu:.4f}")
        print(f"  H25 threshold 0.80: {'PASS' if a_mu >= 0.80 else 'FAIL'}")
    else:
        print("  degenerate band (all-success or all-failure); H25 untestable "
              "as stated")

    # Robustness: same test per eff-stratum *within* the band, and per
    # narrower NU sub-bands, so a single lucky split doesn't carry H25.
    print("\nWithin-band AUC(-mu), split by eff stratum (band population only):")
    for eff in (0.05, 0.10, 0.15, 0.20, 0.25):
        sub = [r for r in band if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if pos and neg:
            a, _, _ = auc_ci_note([r['mu'] for r in pos], [r['mu'] for r in neg])
            print(f"  eff={eff:.2f}  N={len(sub):>3}  rec {len(pos)}/{len(sub)}  "
                  f"AUC(-mu)={a:.4f}")
        else:
            print(f"  eff={eff:.2f}  N={len(sub):>3}  degenerate")

    print("\nNU sub-bands (finer than the single H25 band):")
    edges = [0.5, 1.04, 1.4, 1.8, 2.199, 3.5]
    for lo, hi in zip(edges, edges[1:]):
        sub = [r for r in rows if lo <= r['NU'] < hi]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if pos and neg:
            a, _, _ = auc_ci_note([r['mu'] for r in pos], [r['mu'] for r in neg])
            print(f"  NU in [{lo:.2f},{hi:.2f})  N={len(sub):>3}  "
                  f"rec {len(pos)}/{len(sub)}  AUC(-mu)={a:.4f}")
        else:
            print(f"  NU in [{lo:.2f},{hi:.2f})  N={len(sub):>3}  degenerate")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SECONDARY: does 'step' (block-transition height) beat NU/mu?")
    print("-" * 78)
    print("step = log2(||b*_{m+1}||) - log2(||b*_1||)   (0 at the K1 wall,")
    print("       per W1b the head sits flat at log2 lambda_1(L2) until the")
    print("       wall is crossed, so step should fall as K1 grows)\n")
    pos_all = [r for r in rows if r['ok']]
    neg_all = [r for r in rows if not r['ok']]
    a_step, _, _ = auc_ci_note([r['step'] for r in pos_all],
                                [r['step'] for r in neg_all])
    a_mu_all, _, _ = auc_ci_note([r['mu'] for r in pos_all],
                                  [r['mu'] for r in neg_all])
    a_nu_all, _, _ = auc_ci_note([r['NU'] for r in pos_all],
                                  [r['NU'] for r in neg_all])
    print(f"pooled (N={len(rows)}):")
    print(f"  AUC(-step) = {a_step:.4f}   AUC(-mu) = {a_mu_all:.4f}   "
          f"AUC(-NU) = {a_nu_all:.4f}")
    print(f"  Spearman(step, NU) = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"  Spearman(step, mu) = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\nper-eff-stratum AUC(-step -> recovery), eff held fixed "
          "(honest test, same discipline as W5):")
    for eff in (0.05, 0.10, 0.15, 0.20, 0.25):
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if pos and neg:
            a, _, _ = auc_ci_note([r['step'] for r in pos], [r['step'] for r in neg])
            print(f"  eff={eff:.2f}  N={len(sub):>3}  rec {len(pos)}/{len(sub)}  "
                  f"AUC(-step)={a:.4f}")
        else:
            print(f"  eff={eff:.2f}  N={len(sub):>3}  degenerate")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
