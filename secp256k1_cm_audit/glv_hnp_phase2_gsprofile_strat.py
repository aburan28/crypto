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

Thread 25 (2026-08-07 run #2 proposal) adds two more experiments on the SAME
rows table (no new data collection):

W8  H25 — does mu = lambda_1(L2) separate recovery *inside* the NU ambiguous
    band [1.04, 2.20] (the band where the nearest-plane certificate itself
    gives no answer, per Thread 23b/24 W4)?  If AUC(-mu) >= 0.8 there, mu is
    a genuine second coordinate, not something NU already mediates.
W9  the W1b GS-profile "step" observation, made quantitative:
    step = log2(||b*_{m+1}||) - log2(||b*_1||), the jump from the head block
    (m copies of lambda_1(L2)) to the second block.  Tests whether step is a
    better wall predictor than NU or mu.

Gram-Schmidt is float here, justified by W0/W4 of the parent script
(max relative NU error vs exact Fractions ~1e-15 at dim 20 and dim 24).

Run: python3 glv_hnp_phase2_gsprofile_strat.py [--dump-json PATH]
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

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None,
                     help="write the collected rows table to this path")
    args = ap.parse_args()

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
                step = (math.log2(r['prof'][M17]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][M17] > 0
                        else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        dump_fields = ('n', 'K1', 'ok', 'eff', 'effq', 'lamstar', 'step',
                        'NU', 'nuhat', 'mu', 'l2', 'det2', 'k', 'argmax')
        with open(args.dump_json, 'w') as f:
            json.dump([{k: r[k] for k in dump_fields} for r in rows], f)
        print(f"dumped {len(rows)} rows to {args.dump_json}")

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

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W8 (Thread 25, H25): does mu separate INSIDE the NU ambiguous")
    print("band [1.04, 2.20] (17-bit bracket, Thread 24 W4)?")
    print("-" * 78)
    NU_LO, NU_HI = 1.040, 2.199
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    pos_b = [r for r in band if r['ok']]
    neg_b = [r for r in band if not r['ok']]
    print(f"band population: {len(band)}/{len(rows)}  "
          f"({len(pos_b)} recover, {len(neg_b)} fail)")
    if pos_b and neg_b:
        a_mu_band = auc([r['mu'] for r in pos_b], [r['mu'] for r in neg_b])
        a_nh_band = auc([r['nuhat'] for r in pos_b], [r['nuhat'] for r in neg_b])
        a_nu_band = auc([r['NU'] for r in pos_b], [r['NU'] for r in neg_b])
        print(f"AUC(-mu     -> recovery) inside band = {a_mu_band:.4f}")
        print(f"AUC(-nu_hat -> recovery) inside band = {a_nh_band:.4f}")
        print(f"AUC(-NU     -> recovery) inside band = {a_nu_band:.4f}  "
              f"(sanity check: NU is near-constant here by construction)")
        print(f"\nH25 {'HOLDS' if a_mu_band >= 0.8 else 'FALSIFIED'} "
              f"(threshold 0.8, observed {a_mu_band:.4f})")
    else:
        print("band is empty on one side — cannot evaluate H25 on this table "
              "(need a wider eff sweep or more seeds).")

    # per-eff breakdown of the band, since eff still varies inside it
    print("\nband breakdown by eff stratum:")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if p and ng:
            print(f"  eff={eff:.2f}  N={len(sub):>3}  "
                  f"AUC(-mu)={auc([r['mu'] for r in p], [r['mu'] for r in ng]):.4f}")
        else:
            print(f"  eff={eff:.2f}  N={len(sub):>3}  (degenerate, "
                  f"{len(p)} pos / {len(ng)} neg)")

    # eff-restricted pooled AUC: the naive pooled AUC above re-introduces the
    # same cross-scale confound W3/W6 diagnosed for the closed form (mu's
    # absolute scale differs across eff strata), and eff=0.05 is a near-
    # degenerate class (99/100 recover overall, per W5) that adds noise
    # without information. Report the pooled AUC over the four
    # non-degenerate strata only, which is the honest reading of "does mu
    # separate inside the band, once eff is controlled".
    band_ctrl = [r for r in band if r['effq'] != 0.05]
    pos_c = [r for r in band_ctrl if r['ok']]
    neg_c = [r for r in band_ctrl if not r['ok']]
    if pos_c and neg_c:
        a_mu_ctrl = auc([r['mu'] for r in pos_c], [r['mu'] for r in neg_c])
        print(f"\nAUC(-mu -> recovery), band restricted to eff in "
              f"{{0.10,0.15,0.20,0.25}} (N={len(band_ctrl)}): {a_mu_ctrl:.4f}")
        print(f"H25 (eff-controlled) "
              f"{'HOLDS' if a_mu_ctrl >= 0.8 else 'FALSIFIED'} "
              f"(threshold 0.8, observed {a_mu_ctrl:.4f})")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP W9: step = log2||b*_{m+1}|| - log2||b*_1||  (W1b, quantified)")
    print("-" * 78)
    print("W1b: the head of the GS profile is m copies of lambda_1(L2); the")
    print("step to the second block vanishes near the recovery wall. Test")
    print("whether step beats NU/mu as a wall predictor, pooled and per-eff.")
    have_step = [r for r in rows if not math.isnan(r['step'])]
    print(f"\n{len(have_step)}/{len(rows)} rows have a finite step value")
    pos_s = [r['step'] for r in have_step if r['ok']]
    neg_s = [r['step'] for r in have_step if not r['ok']]
    if pos_s and neg_s:
        print(f"step | success : mean {sum(pos_s)/len(pos_s):.3f}  "
              f"min {min(pos_s):.3f}  max {max(pos_s):.3f}")
        print(f"step | failure : mean {sum(neg_s)/len(neg_s):.3f}  "
              f"min {min(neg_s):.3f}  max {max(neg_s):.3f}")
        print(f"AUC(step -> recovery), SMALLER step -> more likely to fail "
              f"(H24/W1b direction), so score as -(-step)=step for 'larger "
              f"step is safer':")
        a_step_pos = auc([-x for x in pos_s], [-x for x in neg_s])
        print(f"  AUC(-step -> recovery) = {a_step_pos:.4f}  "
              f"(pooled, compare NU pooled AUC above)")
    print("\nper-eff AUC(-step):")
    for eff in EFFS:
        sub = [r for r in have_step if r['effq'] == eff]
        p = [r['step'] for r in sub if r['ok']]
        ng = [r['step'] for r in sub if not r['ok']]
        if p and ng:
            print(f"  eff={eff:.2f}  N={len(sub):>3}  "
                  f"AUC(-step)={auc([-x for x in p], [-x for x in ng]):.4f}  "
                  f"Spearman(step,NU)={spearman([r['step'] for r in sub], [r['NU'] for r in sub]):.4f}")
        else:
            print(f"  eff={eff:.2f}  N={len(sub):>3}  (degenerate)")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
