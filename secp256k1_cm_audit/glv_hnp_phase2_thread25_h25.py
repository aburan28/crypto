"""
GLV-HNP Phase 2, Thread 25: does mu carry information NU does not?

W5/W6 (glv_hnp_phase2_gsprofile_strat.py, 2026-08-07) established that at
fixed eff, nu_hat (equivalently mu = lambda_1(L2)) separates recovery with
AUC 0.75-0.93 while NU (the exact BDD certificate) is uncorrelated with it
(Spearman ~0, sign varies by stratum) and separates worse at larger sizes
(AUC 0.978 -> 0.860, 12->17 bits).  Interpretation offered: NU governs
Babai nearest-plane, mu governs a second, distinct mechanism (the LLL wall
Kannan-LLL crosses that nearest-plane does not).

H25 (pre-registered 2026-08-07): within the ambiguous NU band where the
sound nearest-plane certificate gives no answer, mu still separates
recovery.

  H25: for 1.04 <= NU <= 2.20 (17-bit ambiguous band from W4), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.

Falsifier: if AUC(-mu) inside the band drops near 0.5, mu's W5 power was
mediated entirely by NU (a stratification artifact) and the closed-form
statistic should be retired as a second coordinate.

Secondary (W1b follow-up): step = log2(prof[m]) - log2(prof[0]) is the gap
between the head-block plateau (m exact copies of lambda_1(L2), per W1b)
and the first entry of the second GS block.  Test whether step -> 0 predicts
the wall at least as well as NU or mu.

This script reuses the exact instance-generation loop of
glv_hnp_phase2_gsprofile_strat.py (17-bit curves, EFFS grid, dim 24, float
GS) so the row table is numerically the same population W5/W6 analysed.
Adds --dump-json so the table survives the run for further reanalysis
without regenerating instances (curve search + LLL is the expensive part).

Run: python3 glv_hnp_phase2_thread25_h25.py [--dump-json out.json]
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

# W4 (2026-08-07): 17-bit ambiguous band bracket.
NU_BAND_LO, NU_BAND_HI = 1.040, 2.199


def collect_rows(m17=12, effs=(0.05, 0.10, 0.15, 0.20, 0.25)):
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    rows = []
    for eff in effs:
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), m17, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), m17, d_trial, k1b, seed)
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n)})
                # Secondary statistic: head-plateau -> second-block step.
                prof = r['prof']
                r['step'] = (math.log2(prof[m17]) - math.log2(prof[0])
                             if prof[0] > 0 and prof[m17] > 0 else float('nan'))
                rows.append(r)
    return rows, len(curves17)


def report(rows, ncurves):
    print("=" * 78)
    print("Thread 25 — H25: does mu separate inside the NU-ambiguous band?")
    print("=" * 78)
    print(f"\n{ncurves} 17-bit j=0 GLV curves, {len(rows)} instances "
          f"(float GS, dim {rows[0]['k']})")

    band = [r for r in rows if NU_BAND_LO <= r['NU'] <= NU_BAND_HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print("\n" + "-" * 78)
    print(f"H25: ambiguous band {NU_BAND_LO} <= NU <= {NU_BAND_HI}  "
          f"(N={len(band)}, {len(pos)} recover / {len(neg)} fail)")
    print("-" * 78)
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"  AUC(-mu)     = {a_mu:.4f}   (H25 threshold: >= 0.80)")
        print(f"  AUC(-nu_hat) = {a_nh:.4f}")
        print(f"  AUC(-NU)     = {a_nu:.4f}   (expect ~0.5: no info left "
              f"inside its own ambiguous band, by construction)")
        print(f"  AUC(-step)   = {a_st:.4f}")
        verdict = "CONFIRMED" if a_mu >= 0.80 else "FALSIFIED"
        print(f"\n  H25 verdict: {verdict} (AUC(-mu) = {a_mu:.4f} "
              f"{'>=' if a_mu >= 0.80 else '<'} 0.80)")
    else:
        print("  degenerate: one class empty in the band")

    print("\n" + "-" * 78)
    print("Per-eff-stratum breakdown of the band (does mu's power inside "
          "the band come from one eff value or hold across all?)")
    print("-" * 78)
    print(f"{'eff':>5} {'N_band':>7} {'rec':>7} {'AUC -mu':>9} "
          f"{'AUC -step':>10}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if p and ng:
            a_mu = auc([r['mu'] for r in p], [r['mu'] for r in ng])
            a_st = auc([r['step'] for r in p], [r['step'] for r in ng])
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {a_mu:>9.4f} "
                  f"{a_st:>10.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {'(degen)':>9}")

    print("\n" + "-" * 78)
    print("Secondary: step statistic vs NU/mu over the FULL population "
          "(all strata, not just the band)")
    print("-" * 78)
    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    a_step = auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg])
    a_mu_full = auc([r['mu'] for r in pooled_pos], [r['mu'] for r in pooled_neg])
    a_nu_full = auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg])
    print(f"  AUC(-step) pooled = {a_step:.4f}")
    print(f"  AUC(-mu)   pooled = {a_mu_full:.4f}")
    print(f"  AUC(-NU)   pooled = {a_nu_full:.4f}")
    sp_step_nu = spearman([r['step'] for r in rows], [r['NU'] for r in rows])
    sp_step_mu = spearman([r['step'] for r in rows], [r['mu'] for r in rows])
    print(f"  Spearman(step, NU) = {sp_step_nu:.4f}")
    print(f"  Spearman(step, mu) = {sp_step_mu:.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    dump_path = None
    if "--dump-json" in sys.argv:
        i = sys.argv.index("--dump-json")
        dump_path = sys.argv[i + 1] if i + 1 < len(sys.argv) else "thread25_rows.json"

    t0 = time.time()
    rows, ncurves = collect_rows()
    print(f"collected {len(rows)} instances in {time.time()-t0:.1f}s",
          file=sys.stderr)

    if dump_path:
        slim = [{k: v for k, v in r.items() if k not in ('prof', 'nus')}
                for r in rows]
        with open(dump_path, "w") as f:
            json.dump(slim, f)
        print(f"dumped {len(slim)} rows to {dump_path}", file=sys.stderr)

    report(rows, ncurves)
