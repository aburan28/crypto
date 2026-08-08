"""
GLV-HNP Phase 2, Thread 25: is mu a genuine SECOND coordinate alongside NU,
or was its apparent power in W5 entirely mediated by NU?

Background (2026-08-07 #2 log entry, Thread 24):
  W4  NU (exact BDD certificate) is a sound but size-degrading separator at
      17 bits: NU < 1.040 -> 33/33 recover, NU > 2.199 -> 50/50 fail, and the
      *ambiguous band* [1.040, 2.199] holds the rest.
  W5  nu_hat = mu/sqrt(det L2) is a strong cross-curve separator (AUC
      0.75-0.93) even with eff = K1*K2/n held fixed.
  W6  NU and nu_hat*sqrt(eff) are mutually UNCORRELATED at fixed eff
      (Spearman ~ 0, sign varies by stratum) -- so if both predict recovery
      they must be reading different structure, not the same one twice.

H25 (this run): within the ambiguous NU band [1.040, 2.199], AUC(-mu ->
    Kannan-LLL recovery) stays >= 0.8.
Falsifier: AUC(-mu) collapses to ~0.5 inside the band -> mu's W5 power was a
    stratification artifact (mu correlates with eff/NU well enough outside
    the band to look predictive pooled, but carries nothing once NU is
    controlled) -> retire the closed form.

Secondary (from Thread 24's W1b): the GS profile head is m exact copies of
    lambda_1(L2); quantify step = log2(||b*_{m}||) - log2(||b*_1||) (0-indexed:
    prof[m] vs prof[0], the first vector of the second block vs the first
    vector of the head) and test whether step -> 0 predicts recovery better
    than NU or mu alone.

Reuses the exact instance()/curve-search/generation machinery of
glv_hnp_phase2_gsprofile_strat.py so the table is a byte-for-byte re-roll of
the same experiment (same curves17, EFFS, SEEDS) -- only the analysis is new.
Float GS (justified by W0/W4 of the parent script: exact-vs-float NU relative
error ~1e-15 at this dimension).

Run: python3 glv_hnp_phase2_gsprofile_stratnu.py [--dump-json OUT.json]
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

NU_BAND = (1.040, 2.199)  # ambiguous bracket, W4 @ 17 bits


def build_rows(m=12):
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), m, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), m, d_trial, k1b, seed)
                step = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][m] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    return rows, curves17


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 -- does mu separate INSIDE the ambiguous NU band?")
    print("=" * 78)

    t0 = time.time()
    rows, curves17 = build_rows(m=12)
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves, {len(rows)} instances "
          f"(dim {rows[0]['k']}) in {time.time()-t0:.1f}s")

    dump_path = None
    if "--dump-json" in sys.argv:
        i = sys.argv.index("--dump-json")
        dump_path = sys.argv[i + 1] if i + 1 < len(sys.argv) else \
            "glv_hnp_phase2_gsprofile_stratnu_rows.json"
        with open(dump_path, "w") as f:
            json.dump([{k: v for k, v in r.items()
                        if k not in ('prof', 'nus')} for r in rows], f)
        print(f"dumped {len(rows)} rows -> {dump_path}")

    print("\n" + "-" * 78)
    print("H25: AUC(-mu -> recovery) INSIDE the ambiguous NU band "
          f"[{NU_BAND[0]}, {NU_BAND[1]}]")
    print("-" * 78)
    band = [r for r in rows if NU_BAND[0] <= r['NU'] <= NU_BAND[1]]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band N = {len(band)}  (recover {len(pos)} / fail {len(neg)})")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"  AUC(-mu)     = {a_mu:.4f}   (H25 threshold: >= 0.80)")
        print(f"  AUC(-nu_hat) = {a_nh:.4f}")
        print(f"  AUC(-NU)     = {a_nu:.4f}   (near 0.5 expected: NU is ~constant in-band)")
        print(f"  AUC(-step)   = {a_step:.4f}")
        verdict = "CONFIRMED" if a_mu >= 0.80 else "FALSIFIED"
        print(f"\n  H25 verdict (naive pool over the whole band): {verdict} "
              f"(AUC(-mu) = {a_mu:.4f} vs threshold 0.80)")
        print("  NOTE: W5/W6 already showed raw mu is confounded with eff "
              "pooled; see the eff-conditioned re-test below before trusting this.")
    else:
        print("  degenerate band (all one class) -- H25 untestable on this draw")

    print("\n" + "-" * 78)
    print("Per-eff-stratum breakdown of the band -- W5/W6 already showed mu is")
    print("confounded with eff pooled, so pool-then-AUC over the whole band is")
    print("not a fair test; the honest test holds eff fixed WITHIN the band too.")
    print("-" * 78)
    print(f"{'eff':>5} {'band N':>7} {'rec':>7} {'AUC(-mu)':>9} {'AUC(-nuhat)':>12}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ne = [r for r in sub if not r['ok']]
        a_mu = auc([r['mu'] for r in p], [r['mu'] for r in ne]) if p and ne else float('nan')
        a_nh = auc([r['nuhat'] for r in p], [r['nuhat'] for r in ne]) if p and ne else float('nan')
        flag = "  (degenerate: <2 in minority class)" if p and ne and min(len(p), len(ne)) < 5 else ""
        print(f"{eff:>5.2f} {len(sub):>7} {str(len(p))+'/'+str(len(sub)):>7} "
              f"{a_mu:>9.4f} {a_nh:>12.4f}{flag}")

    non_degen = [eff for eff in sorted(set(r['effq'] for r in rows))
                 if min(sum(1 for r in band if r['effq'] == eff and r['ok']),
                        sum(1 for r in band if r['effq'] == eff and not r['ok'])) >= 5]
    band_ex = [r for r in band if r['effq'] in non_degen]
    p, ne = [r for r in band_ex if r['ok']], [r for r in band_ex if not r['ok']]
    if p and ne:
        a_mu = auc([r['mu'] for r in p], [r['mu'] for r in ne])
        a_nh = auc([r['nuhat'] for r in p], [r['nuhat'] for r in ne])
        print(f"\nband, non-degenerate strata only (eff in {non_degen}): "
              f"N={len(band_ex)} rec={len(p)}/{len(band_ex)}")
        print(f"  AUC(-mu)     = {a_mu:.4f}   (H25 threshold: >= 0.80)")
        print(f"  AUC(-nu_hat) = {a_nh:.4f}")
        verdict2 = "CONFIRMED" if a_mu >= 0.80 else "FALSIFIED"
        print(f"  H25 verdict (eff-conditioned): {verdict2}")

    print("\n" + "-" * 78)
    print("SECONDARY: step = log2||b*_m|| - log2||b*_1|| (0-indexed: second "
          "block first entry minus head first entry)")
    print("-" * 78)
    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    print(f"pooled (N={len(rows)}):")
    print(f"  AUC(-step)   = {auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg]):.4f}")
    print(f"  AUC(-NU)     = {auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg]):.4f}")
    print(f"  AUC(-mu)     = {auc([r['mu'] for r in pooled_pos], [r['mu'] for r in pooled_neg]):.4f}")
    print(f"  Spearman(step, NU) = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"  Spearman(step, mu) = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")
    print(f"  step range: min {min(r['step'] for r in rows):.3f}  "
          f"max {max(r['step'] for r in rows):.3f}  "
          f"mean {sum(r['step'] for r in rows)/len(rows):.3f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
