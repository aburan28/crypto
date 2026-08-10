"""
GLV-HNP Phase 2, Thread 25: does mu carry information NU does not, or is its
apparent power in W5/W6 (glv_hnp_phase2_gsprofile_strat.py) a stratification
artifact?

W5 found: within an eff-stratum, AUC(-nu_hat*sqrt(eff) -> recovery) = 0.75-0.93
and AUC(-mu -> recovery) tracks it closely, while AUC(-NU -> recovery) is
weak (0.35-0.73) and W6 showed NU and nu_hat*sqrt(eff) are RANK-UNCORRELATED
(Spearman in [-0.28, +0.16] per stratum). Since NU is a sound BDD certificate
(zero false positives, W4) and mu/nu_hat is a different, stronger cross-curve
predictor, the open question is whether they are the same information seen
through two different-strength lenses, or two genuinely different
coordinates.

H25 (primary): within NU's own ambiguous band (1.04 <= NU <= 2.20, where the
nearest-plane certificate alone cannot decide recovery either way -- the
band is read off Thread 24's W4 bracket), mu (equivalently nu_hat) still
separates recovery with AUC >= 0.8.
  - If TRUE: mu is a genuine second coordinate; NU is not simply a noisier
    version of it, and (NU, mu) is a real 2-parameter viability test.
  - If FALSE (AUC ~ 0.5 inside the band): mu's power in W5 is mediated
    entirely through NU (i.e. mu and NU are both monotone in some latent
    variable, and conditioning on NU should have removed the mu signal but
    it didn't survive re-testing), and the closed form should be retired
    as "redundant with NU" rather than "a second mechanism".

Secondary (W1b follow-up): the GS profile head is m exact copies of
lambda_1(L2) and the step to the second block vanishes right at the K1 wall
(Thread 24, exp W1b). Define step = log2(prof[m]) - log2(prof[0]) (prof has
dim 2m for the un-projected L0 lattice used by instance()/gsprofile.py, so
index m is the first coordinate of the second block) and test whether
AUC(-step -> recovery) beats NU or mu as a standalone predictor.

Data: identical generation to glv_hnp_phase2_gsprofile_strat.py (same
search_curves call, same SEEDS, same EFFS, same instance()/run_new() calls)
so this is a re-analysis of the same 500-instance population, not a new
experiment -- curve search and gen_signatures are seeded deterministically,
so re-running reproduces the exact same rows.

Run: python3 glv_hnp_phase2_thread25.py
"""

import math
import random
import time

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — does mu survive conditioning on NU?")
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
                    if r['prof'][0] > 0 and r['prof'][m] > 0 else float('nan')
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    assert len(rows) == 500, f"expected 500 instances (reproduction check), got {len(rows)}"

    print("\n" + "-" * 78)
    print("H25: AUC(-mu -> recovery) inside NU's ambiguous band [1.04, 2.20]")
    print("     (band = Thread 24 W4's 17-bit sufficient/necessary bracket)")
    print("-" * 78)
    band = [r for r in rows if 1.040 <= r['NU'] <= 2.199]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"N in band = {len(band)} / {len(rows)}  "
          f"(pos={len(pos)}, neg={len(neg)})")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"  AUC(-mu)      = {a_mu:.4f}")
        print(f"  AUC(-nu_hat)  = {a_nh:.4f}")
        print(f"  AUC(-NU)      = {a_nu:.4f}  (sanity: should be ~0.5, band is NU-flat by construction)")
        print(f"  AUC(-step)    = {a_step:.4f}")
        print(f"  H25 verdict: {'CONFIRMED' if a_mu >= 0.8 else 'FALSIFIED'} "
              f"(threshold 0.8, observed {a_mu:.4f})")
    else:
        print("  degenerate band (all one class) -- cannot compute AUC")

    print("\n" + "-" * 78)
    print("Same test, per eff-stratum (band membership varies in size per stratum)")
    print("-" * 78)
    print(f"{'eff':>5} {'N band':>7} {'rec':>7} | {'AUC mu':>8} {'AUC nu_hat':>10} "
          f"{'AUC step':>9}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p2 = [r for r in sub if r['ok']]
        n2 = [r for r in sub if not r['ok']]
        if not p2 or not n2:
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{(str(len(p2))+'/'+str(len(sub))):>7} | (degenerate)")
            continue
        print(f"{eff:>5.2f} {len(sub):>7} "
              f"{(str(len(p2))+'/'+str(len(sub))):>7} | "
              f"{auc([r['mu'] for r in p2], [r['mu'] for r in n2]):>8.4f} "
              f"{auc([r['nuhat'] for r in p2], [r['nuhat'] for r in n2]):>10.4f} "
              f"{auc([r['step'] for r in p2], [r['step'] for r in n2]):>9.4f}")

    print("\n" + "-" * 78)
    print("Secondary: step = log2(prof[m]) - log2(prof[0]) as a standalone")
    print("           predictor, pooled and per-stratum, vs NU and mu")
    print("-" * 78)
    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    print(f"pooled (N={len(rows)}):")
    print(f"  AUC(-step) = "
          f"{auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg]):.4f}")
    print(f"  AUC(-NU)   = "
          f"{auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg]):.4f}")
    print(f"  AUC(-mu)   = "
          f"{auc([r['mu'] for r in pooled_pos], [r['mu'] for r in pooled_neg]):.4f}")
    print(f"  Spearman(step, NU)     = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"  Spearman(step, mu)     = "
          f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print(f"\n{'eff':>5} {'N':>5} {'rec':>7} | {'AUC step':>9} {'AUC NU':>8} {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        p2 = [r for r in sub if r['ok']]
        n2 = [r for r in sub if not r['ok']]
        if not p2 or not n2:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{(str(len(p2))+'/'+str(len(sub))):>7} | (degenerate)")
            continue
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{(str(len(p2))+'/'+str(len(sub))):>7} | "
              f"{auc([r['step'] for r in p2], [r['step'] for r in n2]):>9.4f} "
              f"{auc([r['NU'] for r in p2], [r['NU'] for r in n2]):>8.4f} "
              f"{auc([r['mu'] for r in p2], [r['mu'] for r in n2]):>8.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
