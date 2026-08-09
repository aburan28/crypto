"""
GLV-HNP Phase 2, Thread 25: does mu carry information NU does not?

Pre-registered by the 2026-08-07 (Thread 24) log entry, after W5/W6 found
that NU (exact BDD certificate) and nu_hat = lambda_1(L2)/sqrt(det L2) are
uncorrelated at fixed eff (Spearman ~ [-0.28, +0.16] within every stratum)
yet both individually separate recovery (AUC 0.75-0.93 for nu_hat, but NU
degrades from 0.978 at 12 bits to 0.860 at 17 bits and is even ANTI-
predictive at eff=0.15 alone).

  H25: within the ambiguous NU band [1.04, 2.20] (17 bits, ex-Thread-24 W4 —
       where nearest-plane gives no clean answer either way), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.

Falsifier: if AUC(-mu) inside the band collapses to ~0.5, mu's apparent
power in W5 was entirely mediated by NU (a stratification artifact of
conditioning on eff instead of NU), and the closed form should be retired
as a *predictor* (it can remain a sound bound via nu_hat*sqrt(eff) ~ NU/C
for extreme cases only).

Secondary (also proposed in the Thread 24 log entry): does the GS-profile
"step" from the first block to the second,

    step = log2(||b*_{m+1}||) - log2(||b*_1||)     (0-indexed: prof[m] vs prof[0])

predict the wall better than NU or mu?  W1b (Thread 24) showed the profile
head is m exact copies of lambda_1(L2) and that this step vanishes exactly
as the K1 wall is crossed, so step -> 0 is a natural third candidate.

Reuses the exact 500-instance 17-bit generation of glv_hnp_phase2_gsprofile_strat.py
(same search_curves call, same SEEDS, same EFFS) so results are directly
comparable to the W5/W6 table already in the log.

Run: python3 glv_hnp_phase2_nu_band.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

# Bracket measured by Thread 24 / W4 at 17 bits, dim 24:
#   sufficient NU < 1.040 (all successes below this),
#   necessary  NU > 2.199 (all failures above this).
NU_LO, NU_HI = 1.040, 2.199

if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — does mu separate INSIDE the NU-ambiguous band?")
    print("=" * 78)

    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")
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
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n)})
                # step diagnostic: log2||b*_{m+1}|| - log2||b*_1||, 0-indexed
                m = M17
                b1, bm1 = r['prof'][0], r['prof'][m]
                r['step'] = (math.log2(bm1) - math.log2(b1)
                             if b1 > 0 and bm1 > 0 else float('nan'))
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print(f"EXP H25: AUC(-mu -> recovery) inside the ambiguous NU band "
          f"[{NU_LO}, {NU_HI}]")
    print("-" * 78)
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    below = [r for r in rows if r['NU'] < NU_LO]
    above = [r for r in rows if r['NU'] > NU_HI]
    print(f"below band (NU < {NU_LO}): {len(below)}  "
          f"({sum(1 for r in below if r['ok'])} recovered)")
    print(f"in   band  [{NU_LO},{NU_HI}]: {len(band)}  "
          f"({sum(1 for r in band if r['ok'])} recovered)")
    print(f"above band (NU > {NU_HI}): {len(above)}  "
          f"({sum(1 for r in above if r['ok'])} recovered)")

    pos = [r['mu'] for r in band if r['ok']]
    neg = [r['mu'] for r in band if not r['ok']]
    if pos and neg:
        a_mu_band = auc(pos, neg)
        print(f"\nAUC(-mu -> recovery) within band, N={len(band)}: "
              f"{a_mu_band:.4f}   (H25 threshold: >= 0.80)")
        print(f"H25 {'HOLDS' if a_mu_band >= 0.80 else 'FALSIFIED'}")
    else:
        a_mu_band = float('nan')
        print(f"\nband is degenerate (pos={len(pos)}, neg={len(neg)}) "
              "-- cannot compute AUC")

    # control: NU itself inside the band should be near-useless by
    # construction (that's what "ambiguous" means) -- sanity check
    if pos and neg:
        posNU = [r['NU'] for r in band if r['ok']]
        negNU = [r['NU'] for r in band if not r['ok']]
        print(f"AUC(-NU -> recovery) within band (sanity, should be ~0.5): "
              f"{auc(posNU, negNU):.4f}")
        posNH = [r['nuhat'] for r in band if r['ok']]
        negNH = [r['nuhat'] for r in band if not r['ok']]
        print(f"AUC(-nu_hat -> recovery) within band: "
              f"{auc(posNH, negNH):.4f}")
        posLS = [r['lamstar'] for r in band if r['ok']]
        negLS = [r['lamstar'] for r in band if not r['ok']]
        print(f"AUC(-lam*   -> recovery) within band (control): "
              f"{auc(posLS, negLS):.4f}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP H25b: same test, per-eff-stratum (does mu need eff pooled in?)")
    print("-" * 78)
    print(f"{'eff':>5} {'N band':>7} {'rec':>7} {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p = [r['mu'] for r in sub if r['ok']]
        ng = [r['mu'] for r in sub if not r['ok']]
        if p and ng:
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {auc(p, ng):>8.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {'(degenerate)':>8}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP S1 (secondary): does step = log2||b*_{m+1}|| - log2||b*_1||")
    print("        predict recovery better than NU or mu, pooled?")
    print("-" * 78)
    valid = [r for r in rows if r['step'] == r['step']]  # drop nan
    pos_s = [r['step'] for r in valid if r['ok']]
    neg_s = [r['step'] for r in valid if not r['ok']]
    print(f"N valid = {len(valid)}/{len(rows)}")
    print(f"AUC(-step -> recovery), pooled = {auc(pos_s, neg_s):.4f}")
    posNU_all = [r['NU'] for r in rows if r['ok']]
    negNU_all = [r['NU'] for r in rows if not r['ok']]
    posMU_all = [r['mu'] for r in rows if r['ok']]
    negMU_all = [r['mu'] for r in rows if not r['ok']]
    print(f"AUC(-NU   -> recovery), pooled = "
          f"{auc(posNU_all, negNU_all):.4f}   (reference)")
    print(f"AUC(-mu   -> recovery), pooled = "
          f"{auc(posMU_all, negMU_all):.4f}   (reference)")
    print(f"step | success: mean {sum(pos_s)/len(pos_s):.3f}  "
          f"median {sorted(pos_s)[len(pos_s)//2]:.3f}")
    print(f"step | failure: mean {sum(neg_s)/len(neg_s):.3f}  "
          f"median {sorted(neg_s)[len(neg_s)//2]:.3f}")
    print(f"Spearman(step, NU)  = {spearman([r['step'] for r in valid], [r['NU'] for r in valid]):.4f}")
    print(f"Spearman(step, mu)  = {spearman([r['step'] for r in valid], [r['mu'] for r in valid]):.4f}")

    print("\nEXP S1b: AUC(-step) within the NU-ambiguous band")
    band_s = [r for r in band if r['step'] == r['step']]
    pos_sb = [r['step'] for r in band_s if r['ok']]
    neg_sb = [r['step'] for r in band_s if not r['ok']]
    if pos_sb and neg_sb:
        print(f"N={len(band_s)}  AUC(-step -> recovery) = "
              f"{auc(pos_sb, neg_sb):.4f}")
    else:
        print("degenerate")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
