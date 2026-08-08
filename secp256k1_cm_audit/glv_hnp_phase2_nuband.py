"""
GLV-HNP Phase 2, Thread 25: does mu carry a second, NU-independent signal?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry.  W5/W6 found that
NU (the exact BDD certificate) and nu_hat*sqrt(eff) (the closed-form
predictor) are uncorrelated at fixed eff (Spearman ~ -0.2..+0.2), yet both
separate recovery well in isolation.  If they are genuinely two coordinates,
mu should still separate recovery INSIDE the band where NU alone gives no
answer.

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (Thread 23b's two-sided
       empirical bracket, widened to the 17-bit re-measure), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.
  Falsifier: if AUC inside the band drops to ~0.5, mu's apparent power in W5
       is entirely mediated by NU and the W5 result is a stratification
       artifact.

Secondary (W1b follow-up): step = log2(||b*_{m+1}||) - log2(||b*_1||), i.e.
the jump from the flat lambda_1(L2)-repeated head to the second (K1-moving)
block.  Test whether step predicts recovery, and whether it is just a
monotone reparametrisation of NU/mu or carries independent information.

Run: python3 glv_hnp_phase2_nuband.py
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

DUMP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "glv_hnp_phase2_nuband_table.json")


def build_table():
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
                m = r['k'] // 2
                head = r['prof'][0]
                second = r['prof'][m]
                step = (math.log2(second) - math.log2(head)
                        if head > 0 and second > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'm': m, 'step': step})
                rows.append(r)
    return rows


def to_jsonable(rows):
    keys = ('k', 'NU', 'argmax', 'enorm', 'mu', 'l2', 'det2', 'nuhat',
            'S_K1', 'S_K2', 'K2', 'n', 'K1', 'ok', 'eff', 'effq',
            'lamstar', 'm', 'step')
    return [{key: r[key] for key in keys} for r in rows]


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — is mu a second coordinate, independent of NU? (H25)")
    print("=" * 78)

    if "--from-json" in sys.argv and os.path.exists(DUMP_PATH):
        with open(DUMP_PATH) as f:
            rows = json.load(f)
        print(f"\nloaded {len(rows)} instances from {DUMP_PATH}")
    else:
        t0 = time.time()
        rows = build_table()
        print(f"\n{len(rows)} instances (float GS, dim {rows[0]['k']}) "
              f"in {time.time()-t0:.1f}s")
        if "--dump-json" in sys.argv:
            with open(DUMP_PATH, "w") as f:
                json.dump(to_jsonable(rows), f)
            print(f"dumped table to {DUMP_PATH}")

    print("\n" + "-" * 78)
    print("EXP X1: H25 — AUC(-mu -> recovery) inside the NU-ambiguous band")
    print("-" * 78)
    LO, HI = 1.04, 2.20
    band = [r for r in rows if LO <= r['NU'] <= HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band NU in [{LO}, {HI}]: N={len(band)}  "
          f"rec={len(pos)}/{len(band)}")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nuhat = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        a_nu_inband = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        print(f"AUC(-mu     -> recovery | band) = {a_mu:.4f}")
        print(f"AUC(-nu_hat -> recovery | band) = {a_nuhat:.4f}")
        print(f"AUC(-step   -> recovery | band) = {a_step:.4f}")
        print(f"AUC(-NU     -> recovery | band) = {a_nu_inband:.4f}  "
              "(sanity: NU restricted to its own ambiguous band)")
        verdict = "H25 HOLDS" if a_mu >= 0.8 else "H25 FALSIFIED"
        print(f"\n{verdict}: AUC(-mu) inside band = {a_mu:.4f}"
              f" ({'>=' if a_mu >= 0.8 else '<'} 0.8 threshold)")
    else:
        print("degenerate band (all-pos or all-neg) -- cannot compute AUC")

    print("\nper-band-decile detail (sorted by mu):")
    band_sorted = sorted(band, key=lambda r: r['mu'])
    ndec = max(1, len(band_sorted) // 5)
    for i in range(0, len(band_sorted), ndec):
        chunk = band_sorted[i:i + ndec]
        rec = sum(1 for r in chunk if r['ok'])
        print(f"  mu in [{chunk[0]['mu']:.3g}, {chunk[-1]['mu']:.3g}]  "
              f"rec {rec}/{len(chunk)}")

    print("\n" + "-" * 78)
    print("EXP X2: pooled AUC of mu and step vs NU, for context")
    print("-" * 78)
    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    print(f"N={len(rows)}, rec={len(pooled_pos)}/{len(rows)}")
    print(f"AUC(-mu)   pooled (mixes eff, expect inflated) = "
          f"{auc([r['mu'] for r in pooled_pos], [r['mu'] for r in pooled_neg]):.4f}")
    print(f"AUC(-step) pooled = "
          f"{auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg]):.4f}")
    print(f"AUC(-NU)   pooled = "
          f"{auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg]):.4f}")

    print("\n" + "-" * 78)
    print("EXP X3: per-eff-stratum AUC(-step), same design as W5")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC step':>9} {'AUC mu':>8} "
          f"{'Spearman(step,NU)':>19} {'Spearman(step,mu)':>19}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | (degenerate)")
            continue
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        sp_nu = spearman([r['step'] for r in sub], [r['NU'] for r in sub])
        sp_mu = spearman([r['step'] for r in sub], [r['mu'] for r in sub])
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_step:>9.4f} "
              f"{a_mu:>8.4f} {sp_nu:>19.4f} {sp_mu:>19.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
