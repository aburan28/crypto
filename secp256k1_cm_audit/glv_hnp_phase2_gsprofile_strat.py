"""
GLV-HNP Phase 2, Thread 24b: is the closed-form separator actually doing
cross-curve work, or is it just re-reading the bias strength?

W3 (glv_hnp_phase2_gsprofile.py) found AUC(-nu_hat*sqrt(eff) -> recovery)
= 0.992 on the 22-cell U2 grid, beating the exact BDD certificate NU (0.978).
That grid varies K1 over 11 values on only 2 curves, and recovery is monotone
in K1 within a curve, so a predictor monotone in K1 scores high almost for
free.  The honest test holds eff FIXED and asks whether nu_hat still ranks
curves correctly.

W5  per-eff-stratum AUC at 17 bits (20 curves x 5 seeds per stratum), for
    three scores: nu_hat (no lattice work), NU (exact BDD certificate),
    and lam* (the quantity Thread 20 falsified).
W6  is C = NU/(nu_hat*sqrt(eff)) stable across strata, i.e. does the closed
    form carry an absolute scale or only a ranking?

Gram-Schmidt is float here, justified by W0/W4 of the parent script
(max relative NU error vs exact Fractions ~1e-15 at dim 20 and dim 24).

Run: python3 glv_hnp_phase2_gsprofile_strat.py
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

DUMP_JSON = "--dump-json" in sys.argv

if __name__ == "__main__":
    print("=" * 78)
    print("Thread 24b — cross-curve test of the closed-form separator (eff fixed)")
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
                m_ = r['k'] // 2
                step = (math.log2(r['prof'][m_]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][m_] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if DUMP_JSON:
        outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "glv_hnp_phase2_gsprofile_strat_rows.json")
        slim = [{k: v for k, v in r.items() if k not in ('prof', 'nus')}
                for r in rows]
        for r, s in zip(slim, rows):
            r['prof_head'] = s['prof'][0]
            r['prof_mid'] = s['prof'][s['k'] // 2]
            r['prof_tail'] = s['prof'][-1]
        with open(outpath, "w") as f:
            json.dump(slim, f)
        print(f"dumped {len(slim)} rows to {outpath}")

    print("\n" + "-" * 78)
    print("EXP W5: AUC within each eff stratum — eff is CONSTANT, so the only")
    print("        signal left is the cross-curve geometry.")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC nu_hat':>11} {'AUC NU':>8} "
          f"{'AUC lam*':>9} | {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | "
                  f"{'(degenerate)':>11}")
            continue
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_ls = auc([r['lamstar'] for r in pos], [r['lamstar'] for r in neg])
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_nh:>11.4f} "
              f"{a_nu:>8.4f} {a_ls:>9.4f} | {a_mu:>8.4f}")

    print("\nAUC > 0.5 means SMALLER score -> more likely to recover.")
    print("Thread 20 falsified lam* as a predictor; it is the control column.")

    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    print(f"\npooled over all strata (N={len(rows)}):")
    print(f"  AUC(-nu_hat*sqrt(eff)) = "
          f"{auc([r['nuhat']*math.sqrt(r['eff']) for r in pooled_pos], [r['nuhat']*math.sqrt(r['eff']) for r in pooled_neg]):.4f}")
    print(f"  AUC(-NU)               = "
          f"{auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg]):.4f}")
    print(f"  AUC(-nu_hat alone)     = "
          f"{auc([r['nuhat'] for r in pooled_pos], [r['nuhat'] for r in pooled_neg]):.4f}")
    print(f"  AUC(-eff alone)        = "
          f"{auc([r['eff'] for r in pooled_pos], [r['eff'] for r in pooled_neg]):.4f}")

    print("\n" + "-" * 78)
    print("EXP W6: is C = NU / (nu_hat*sqrt(eff)) an absolute constant?")
    print("-" * 78)
    print(f"{'eff':>5} {'mean C':>9} {'min':>8} {'max':>8} {'spread':>8} "
          f"{'Spearman(pred,NU)':>19}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        cs = [r['NU'] / (r['nuhat'] * math.sqrt(r['eff'])) for r in sub
              if r['nuhat'] > 0]
        pr = [r['nuhat'] * math.sqrt(r['eff']) for r in sub]
        nu = [r['NU'] for r in sub]
        print(f"{eff:>5.2f} {sum(cs)/len(cs):>9.3f} {min(cs):>8.3f} "
              f"{max(cs):>8.3f} {max(cs)/min(cs):>8.2f}x "
              f"{spearman(pr, nu):>19.4f}")

    print("\n" + "-" * 78)
    print("EXP W7: per-curve detail at the discriminating stratum")
    print("-" * 78)
    disc = None
    best = -1
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        w = sum(1 for r in sub if r['ok'])
        bal = min(w, len(sub) - w)
        if bal > best:
            best, disc = bal, eff
    sub = [r for r in rows if r['effq'] == disc]
    print(f"stratum eff = {disc:.2f}  ({best} of the minority class)\n")
    print(f"{'n':>8} {'lam*':>7} {'nu_hat':>8} {'mean NU':>9} {'rec':>6}")
    bycurve = {}
    for r in sub:
        bycurve.setdefault(r['n'], []).append(r)
    for n in sorted(bycurve, key=lambda x: bycurve[x][0]['nuhat']):
        g = bycurve[n]
        print(f"{n:>8} {g[0]['lamstar']:>7.4f} {g[0]['nuhat']:>8.4f} "
              f"{sum(x['NU'] for x in g)/len(g):>9.4f} "
              f"{str(sum(1 for x in g if x['ok']))+'/'+str(len(g)):>6}")

    print("\n" + "-" * 78)
    print("EXP W8 (Thread 25 / H25): does mu still separate INSIDE the")
    print("        ambiguous NU band, where nearest-plane gives no answer?")
    print("-" * 78)
    NU_LO, NU_HI = 1.040, 2.199
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    posb = [r for r in band if r['ok']]
    negb = [r for r in band if not r['ok']]
    print(f"band [{NU_LO}, {NU_HI}]: N={len(band)}  rec={len(posb)}/{len(band)}")
    if posb and negb:
        a_mu_band = auc([r['mu'] for r in posb], [r['mu'] for r in negb])
        a_nh_band = auc([r['nuhat'] for r in posb], [r['nuhat'] for r in negb])
        print(f"  AUC(-mu -> recovery)     inside band = {a_mu_band:.4f}   "
              f"(H25 predicts >= 0.8)")
        print(f"  AUC(-nu_hat -> recovery) inside band = {a_nh_band:.4f}")
        print(f"  H25 verdict: {'HOLDS' if a_mu_band >= 0.8 else 'FALSIFIED'}")
    else:
        print("  degenerate: one class empty inside the band")

    print("\n" + "-" * 78)
    print("EXP W9: does the head->mid GS-profile step predict the wall?")
    print("  step = log2||b*_{m+1}|| - log2||b*_1||  (W1b: vanishes at the wall)")
    print("-" * 78)
    posS = [r['step'] for r in rows if r['ok'] and not math.isnan(r['step'])]
    negS = [r['step'] for r in rows if not r['ok'] and not math.isnan(r['step'])]
    print(f"pooled AUC(-step -> recovery) = {auc(posS, negS):.4f}   N={len(posS)+len(negS)}")
    print(f"pooled AUC(-NU   -> recovery) = "
          f"{auc([r['NU'] for r in rows if r['ok']], [r['NU'] for r in rows if not r['ok']]):.4f}")
    print(f"pooled AUC(-mu   -> recovery) = "
          f"{auc([r['mu'] for r in rows if r['ok']], [r['mu'] for r in rows if not r['ok']]):.4f}")
    print(f"\n{'eff':>5} {'N':>5} {'rec':>7} | {'AUC step':>9} {'AUC NU':>8} {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            continue
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_step:>9.4f} "
              f"{a_nu:>8.4f} {a_mu:>8.4f}")
    print(f"\nSpearman(step, NU) pooled = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"Spearman(step, log(mu)) pooled = "
          f"{spearman([r['step'] for r in rows], [math.log(r['mu']) for r in rows]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
