"""
GLV-HNP Phase 2, Thread 25: does mu (equivalently nu_hat) carry independent
signal INSIDE the NU-ambiguous band, and does the GS-profile "step" (block-1
to block-2 jump) predict the wall?

Pre-registered by the 2026-08-07 (Thread 24) log entry:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where the exact BDD
       certificate NU gives no answer), AUC(-mu -> Kannan-LLL recovery)
       stays >= 0.8.

  Secondary (from W1b): step = log2(||b*_{m+1}||) - log2(||b*_1||); does
       step -> 0 predict the wall better than NU or mu?

Reuses the exact generation code from glv_hnp_phase2_gsprofile_strat.py
(same curves17/EFFS/SEEDS) so results are directly comparable to W5/W6/W7.

Run: python3 glv_hnp_thread25_nu_band.py
"""

import math
import random
import time

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

NU_LO, NU_HI = 1.040, 2.199   # W4's ambiguous bracket at 17 bits


def report(name, sub):
    pos = [r for r in sub if r['ok']]
    neg = [r for r in sub if not r['ok']]
    if not pos or not neg:
        print(f"  {name:<12} N={len(sub):>4} rec={len(pos)}/{len(sub)}  (degenerate)")
        return
    a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
    a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
    a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
    print(f"  {name:<12} N={len(sub):>4} rec={len(pos):>3}/{len(sub):<3} "
          f"AUC(-step)={a_step:>7.4f} AUC(-mu)={a_mu:>7.4f} AUC(-NU)={a_nu:>7.4f}")


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — H25 (mu inside the NU-ambiguous band) + step falsifier")
    print("=" * 78)

    t0 = time.time()
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
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n)})
                m = M17
                r['step'] = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                             if r['prof'][0] > 0 and r['prof'][m] > 0 else None)
                rows.append(r)
    rows = [r for r in rows if r['step'] is not None]
    print(f"{len(rows)} instances (dim {rows[0]['k']}, m={M17}) "
          f"in {time.time()-t0:.1f}s\n")

    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]

    print(f"H25 test, band NU in [{NU_LO},{NU_HI}]:")
    print("-" * 78)
    print("  pooled across eff (naive, confounded like W3):")
    report("band", band)
    print("\n  per-eff inside band (honest, eff held fixed):")
    for effq in EFFS:
        report(f"eff={effq:.2f}", [r for r in band if r['effq'] == effq])

    print("\nStep falsifier, unconditional per-eff stratum (mirrors W5, adds step):")
    print("-" * 78)
    for effq in EFFS:
        report(f"eff={effq:.2f}", [r for r in rows if r['effq'] == effq])

    print(f"\nSpearman(step, NU) pooled = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"Spearman(step, mu) pooled = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
