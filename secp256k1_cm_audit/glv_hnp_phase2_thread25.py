"""
GLV-HNP Phase 2, Thread 25: find the second mechanism by conditioning on NU.

W5/W6 (glv_hnp_phase2_gsprofile_strat.py, 2026-08-07 run #2) established that
recovery = f(NU, X) with X ~ mu-driven (mu = lambda_1(L2)) and X essentially
independent of NU (Spearman(pred, NU) in [-0.28, +0.16] across strata). W4
found NU is a sound zero-FP certificate at NU <= 1.04 but the ambiguous band
widens to [1.040, 2.199] at 17 bits.

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where the exact
       nearest-plane certificate gives no answer either way), mu still
       separates viable from non-viable instances: AUC(-mu -> recovery) >= 0.8.

Falsifier: if AUC(-mu -> recovery) inside the band collapses to ~0.5, mu's
apparent power in W5 was entirely mediated by NU (a stratification artifact:
within an eff-stratum, low-mu curves also tend to have low NU), and the
closed form should be retired as a genuine second coordinate.

Secondary (W1b follow-up): step = log2(||b*_{m+1}||) - log2(||b*_1||), the
jump in the GS profile from the first (m-fold-repeated lambda_1(L2)) block to
the second block. W1b observed this step visually vanishes as K1 crosses the
wall. Test whether step predicts recovery, standalone and inside the NU band.

Run: python3 glv_hnp_phase2_thread25.py [--dump-json OUT.json]
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

NU_BAND_LO = 1.040
NU_BAND_HI = 2.199


def collect_rows():
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
                m = M17
                step = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][m] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                # 'prof' is a list, not JSON-trivial to keep compact; drop it
                # from the dumped record but keep everything else.
                r.pop('prof', None)
                r.pop('nus', None)
                rows.append(r)
    return rows, len(curves17)


if __name__ == "__main__":
    dump_path = None
    if "--dump-json" in sys.argv:
        i = sys.argv.index("--dump-json")
        dump_path = sys.argv[i + 1] if i + 1 < len(sys.argv) else "thread25_rows.json"

    print("=" * 78)
    print("Thread 25 — conditioning on NU: is mu a genuine second coordinate?")
    print("=" * 78)

    t0 = time.time()
    rows, ncurves = collect_rows()
    print(f"\n{ncurves} 17-bit j=0 GLV curves, {len(rows)} instances "
          f"(float GS, dim {rows[0]['k']}) in {time.time()-t0:.1f}s")

    if dump_path:
        with open(dump_path, "w") as f:
            json.dump(rows, f)
        print(f"dumped {len(rows)} rows to {dump_path}")

    print("\n" + "-" * 78)
    print(f"H25: AUC(-mu -> recovery) inside the ambiguous NU band "
          f"[{NU_BAND_LO}, {NU_BAND_HI}]")
    print("-" * 78)

    band = [r for r in rows if NU_BAND_LO <= r['NU'] <= NU_BAND_HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"N in band = {len(band)} / {len(rows)} total "
          f"({len(pos)} recover, {len(neg)} fail)")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_ls = auc([r['lamstar'] for r in pos], [r['lamstar'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"  AUC(-mu     -> recovery) = {a_mu:.4f}   (H25 threshold 0.80)")
        print(f"  AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"  AUC(-lam*   -> recovery) = {a_ls:.4f}")
        print(f"  AUC(-step   -> recovery) = {a_step:.4f}")
        print(f"  verdict: H25 {'HOLDS' if a_mu >= 0.80 else 'FALSIFIED'} "
              f"(AUC={a_mu:.4f} {'>=' if a_mu >= 0.80 else '<'} 0.80)")
    else:
        print("  degenerate band (all-pos or all-neg) -- cannot compute AUC")

    print("\nper-eff-stratum breakdown, band-only (AUC step reported as")
    print("+step, i.e. 1 - auc(-step,...), since step and recovery are")
    print("POSITIVELY correlated -- opposite convention from mu/NU/nu_hat):")
    print(f"{'eff':>5} {'N_band':>7} {'rec':>7} {'AUC mu':>8} "
          f"{'AUC nu_hat':>10} {'AUC +step':>10}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if p and ng:
            am = auc([r['mu'] for r in p], [r['mu'] for r in ng])
            anh = auc([r['nuhat'] for r in p], [r['nuhat'] for r in ng])
            astp = 1.0 - auc([r['step'] for r in p], [r['step'] for r in ng])
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {am:>8.4f} "
                  f"{anh:>10.4f} {astp:>10.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {'(degen)':>8}")

    print("\n" + "-" * 78)
    print("Secondary: does the GS-profile step predict recovery on its own?")
    print("-" * 78)
    allpos = [r for r in rows if r['ok']]
    allneg = [r for r in rows if not r['ok']]
    a_step_all = auc([r['step'] for r in allpos], [r['step'] for r in allneg])
    a_nu_all = auc([r['NU'] for r in allpos], [r['NU'] for r in allneg])
    print(f"AUC(-step -> recovery), pooled N={len(rows)}: {a_step_all:.4f}  "
          f"(i.e. AUC(+step)={1-a_step_all:.4f} -- step is POSITIVELY "
          f"correlated with recovery, opposite sign from mu/NU/nu_hat)")
    print(f"AUC(-NU   -> recovery), pooled N={len(rows)}: {a_nu_all:.4f}  (reference)")
    print(f"Spearman(step, NU) = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"Spearman(step, mu) = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\nstep by eff stratum (mean, recover vs fail):")
    print(f"{'eff':>5} {'mean step|rec':>14} {'mean step|fail':>15}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in rows if r['effq'] == eff]
        p = [r['step'] for r in sub if r['ok']]
        ng = [r['step'] for r in sub if not r['ok']]
        if p and ng:
            print(f"{eff:>5.2f} {sum(p)/len(p):>14.3f} {sum(ng)/len(ng):>15.3f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
