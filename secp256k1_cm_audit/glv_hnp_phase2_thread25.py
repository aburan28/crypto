"""
GLV-HNP Phase 2, Thread 25: does mu separate recovery INSIDE the NU
ambiguous band, and does the GS-profile "step" do better than either?

Pre-registered by the 2026-08-07 (autolab run #2) log entry (Thread 24,
EXP W5/W6): W5 found AUC(-mu -> recovery) = 0.75-0.93 in every eff-fixed
17-bit stratum, W6 found NU and nu_hat*sqrt(eff) are UNCORRELATED at fixed
eff (Spearman in [-0.28, +0.16], sign-flipping across strata). Two mutually
uncorrelated quantities that both predict recovery govern different things:
NU is the exact BDD nearest-plane certificate (sound, AUC 0.978, zero false
positives at NU <= 1, Thread 23b) but loose by ~1.9x in the necessary
direction (empirical wall NU ~ 1.87-2.20 vs the NU<=1 guarantee). The open
question is whether that slack band has a second, NU-independent cause.

H25: within the ambiguous band [1.040, 2.199] (the 17-bit two-sided bracket
     measured in Thread 24 EXP W4 -- NU < 1.040 always recovers, NU > 2.199
     always fails), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

Falsifier: if AUC(-mu) inside the band is ~0.5 (chance), mu's apparent power
in W5 is entirely mediated by NU/eff correlation (a stratification artifact),
and nu_hat should be retired as a *causal* coordinate (it can remain useful
as a lattice-free proxy outside the band).

Secondary (W1b-motivated): step = log2(||b*_{m+1}||) - log2(||b*_1||), the
GS-profile jump from the "m copies of lambda_1(L2)" head (W1b) to the second
block. Test AUC(-step -> recovery), band-conditioned, against NU and mu.

Numerics: float GS (justified by Thread 24 EXP W0/W4: max relative NU error
vs exact Fractions ~1e-15 at dim 20 and dim 24), matching
glv_hnp_phase2_gsprofile_strat.py so this table is a superset/rerun of that
one -- same construction, same 500-instance grid (20 curves x 5 seeds x 5
eff strata, dim 24), with the ambiguous-band analysis added and a
--dump-json flag so the table survives the run for future threads.

Run: python3 glv_hnp_phase2_thread25.py [--dump-json OUT.json]
"""

import argparse
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

# Two-sided bracket measured at 17 bits, Thread 24 EXP W4 (log line ~6410-6417).
NU_SUFFICIENT = 1.040   # NU < this  -> 33/33 recovered in that run
NU_NECESSARY = 2.199    # NU > this  -> 50/50 failed in that run


def band_report(rows, lo, hi):
    band = [r for r in rows if lo <= r['NU'] <= hi]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band NU in [{lo:.3f}, {hi:.3f}]: N={len(band)}  "
          f"recovered {len(pos)}/{len(band)}")
    if not pos or not neg:
        print("  degenerate band (one class empty) -- cannot compute AUC")
        return band, None
    a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
    a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
    a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
    a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
    a_ef = auc([r['eff'] for r in pos], [r['eff'] for r in neg])
    print(f"  AUC(-mu)     = {a_mu:.4f}")
    print(f"  AUC(-nu_hat) = {a_nh:.4f}")
    print(f"  AUC(-NU)     = {a_nu:.4f}   (should be near 0.5: NU is ~constant in-band)")
    print(f"  AUC(-step)   = {a_st:.4f}")
    print(f"  AUC(-eff)    = {a_ef:.4f}   (control: is this just re-reading eff?)")
    return band, dict(mu=a_mu, nuhat=a_nh, NU=a_nu, step=a_st, eff=a_ef, N=len(band))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None)
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 -- does mu separate recovery INSIDE the NU ambiguous band?")
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
                m = r['k'] // 2
                step = math.log2(r['prof'][m]) - math.log2(r['prof'][0]) \
                    if r['prof'][0] > 0 and r['prof'][m] > 0 else 0.0
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        slim = [{k: v for k, v in r.items() if k not in ('prof', 'nus')}
                for r in rows]
        with open(args.dump_json, "w") as f:
            json.dump(slim, f)
        print(f"\ndumped {len(slim)} rows to {args.dump_json}")

    print("\n" + "-" * 78)
    print("SANITY: pooled AUCs, reproducing Thread 24 W6 pooled numbers")
    print("-" * 78)
    pos_all = [r for r in rows if r['ok']]
    neg_all = [r for r in rows if not r['ok']]
    print(f"pooled AUC(-NU)     = "
          f"{auc([r['NU'] for r in pos_all], [r['NU'] for r in neg_all]):.4f}")
    print(f"pooled AUC(-mu)     = "
          f"{auc([r['mu'] for r in pos_all], [r['mu'] for r in neg_all]):.4f}")
    print(f"pooled AUC(-nu_hat) = "
          f"{auc([r['nuhat'] for r in pos_all], [r['nuhat'] for r in neg_all]):.4f}")
    print(f"pooled AUC(-step)   = "
          f"{auc([r['step'] for r in pos_all], [r['step'] for r in neg_all]):.4f}")

    print("\n" + "-" * 78)
    print(f"H25 TEST: the NU ambiguous band [{NU_SUFFICIENT}, {NU_NECESSARY}]")
    print("-" * 78)
    band, scores = band_report(rows, NU_SUFFICIENT, NU_NECESSARY)

    print("\n" + "-" * 78)
    print("H25, per-eff-stratum inside the band (is it just re-reading eff?)")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC mu':>8} {'AUC step':>9}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | (degenerate)")
            continue
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_mu:>8.4f} {a_st:>9.4f}")

    print("\n" + "-" * 78)
    print("SECONDARY: does step predict the wall better than NU or mu, globally?")
    print("-" * 78)
    print(f"AUC(-step -> recovery), pooled            = "
          f"{auc([r['step'] for r in pos_all], [r['step'] for r in neg_all]):.4f}")
    print(f"Spearman(step, NU)                        = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"Spearman(step, mu)                        = "
          f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    if scores is not None:
        if scores['mu'] >= 0.8:
            print(f"H25 SURVIVES: AUC(-mu) = {scores['mu']:.4f} >= 0.8 inside the "
                  f"ambiguous band (N={scores['N']}). mu is a genuine second "
                  f"coordinate independent of NU.")
        elif scores['mu'] <= 0.6:
            print(f"H25 FALSIFIED: AUC(-mu) = {scores['mu']:.4f}, near chance, "
                  f"inside the ambiguous band (N={scores['N']}). mu's power in "
                  f"W5 is mediated by eff/NU; retire it as a causal coordinate.")
        else:
            print(f"H25 INCONCLUSIVE: AUC(-mu) = {scores['mu']:.4f} inside the "
                  f"band (N={scores['N']}) -- above chance but below the 0.8 bar.")
        if scores['step'] > scores['mu'] and scores['step'] > scores['NU']:
            print(f"step ({scores['step']:.4f}) beats both mu ({scores['mu']:.4f}) "
                  f"and NU ({scores['NU']:.4f}) inside the band.")
    else:
        print("band degenerate at this grid -- cannot evaluate H25 as stated.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
