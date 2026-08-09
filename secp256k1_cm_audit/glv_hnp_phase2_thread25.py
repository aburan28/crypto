"""
GLV-HNP Phase 2, Thread 25: does mu carry a SECOND coordinate on top of NU,
and does the W1b two-block GS-profile step predict the K1 wall directly?

Proposed by the 2026-08-07 #2 (Thread 24) log entry, from W5/W6:
  - NU (exact BDD certificate) is a sound size-free certificate but a
    size-degrading separator (AUC 0.978 at 12 bits -> 0.860 at 17 bits).
  - mu = lambda_1(L2) (no lattice work) tracks recovery with AUC 0.75-0.93
    in every non-degenerate eff-stratum, HOLDING eff fixed.
  - NU and nu_hat*sqrt(eff) are mutually uncorrelated at fixed eff
    (Spearman in [-0.28, +0.16] across strata) -> they measure different
    mechanisms, not one.

H25: within the ambiguous NU band [1.04, 2.20] (nearest-plane gives no
     verdict there), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.
Falsifier: AUC in-band drops toward 0.5 -> mu's apparent power was entirely
     mediated by NU (a stratification artifact) and should be retired.

Secondary (from W1b): the GS profile head is m copies of lambda_1(L2); the
step to the second block vanishes as the K1 wall is crossed. Define
     step = log2(||b*_{m}||) - log2(||b*_0||)     (0-indexed, k=2m)
and test whether step predicts recovery, and whether it is just a restated
NU/mu or new information (Spearman against both).

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


def stratified_auc(rows, key, effs):
    """Concordance pooled ONLY within each fixed-eff stratum, then summed.

    Plain pooling across strata mixes scales that are each individually
    monotone in eff (Thread 24 W6's C-drift artifact); this is the honest
    number when eff is a nuisance variable.
    """
    num, den = 0.0, 0
    for eff in effs:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r[key] for r in sub if r['ok']]
        neg = [r[key] for r in sub if not r['ok']]
        if not pos or not neg:
            continue
        conc = sum((1.0 if a < b else 0.5 if a == b else 0.0)
                   for a in pos for b in neg)
        num += conc
        den += len(pos) * len(neg)
    return num / den if den else float('nan')

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None,
                     help="write the raw instance table to this path")
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — does mu survive conditioning on NU? (H25) + GS step")
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
                step = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][m] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'm': m, 'step': step})
                del r['prof']  # keep the dump small; not needed downstream
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        with open(args.dump_json, "w") as f:
            json.dump(rows, f)
        print(f"wrote {len(rows)} rows to {args.dump_json}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("H25: AUC(-mu -> recovery) WITHIN the ambiguous NU band [1.04, 2.20]")
    print("-" * 78)
    LO, HI = 1.04, 2.20
    band = [r for r in rows if LO <= r['NU'] <= HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band size N={len(band)}  (pos={len(pos)}, neg={len(neg)}) "
          f"out of {len(rows)} total")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        print(f"AUC(-mu     -> recovery) in-band, POOLED (naive)   = {a_mu:.4f}")
        print(f"AUC(-nu_hat -> recovery) in-band, POOLED (naive)   = {a_nh:.4f}")
        print(f"AUC(-NU     -> recovery) in-band, POOLED (naive)   = {a_nu:.4f}  "
              f"(expect ~0.5: NU is constant-ish inside its own band)")
        print("(naive pooling mixes strata that are each individually "
              "monotone in eff -- Thread 24 W6's confound. STRATIFIED below "
              "is the honest number.)")
        effs_nondeg = tuple(e for e in EFFS if e != 0.05)  # 0.05 is 99/100
        sa_mu = stratified_auc(band, 'mu', effs_nondeg)
        sa_nu = stratified_auc(band, 'NU', effs_nondeg)
        print(f"\nAUC(-mu -> recovery) in-band, STRATIFIED (excl eff=0.05, "
              f"99/100 degenerate) = {sa_mu:.4f}")
        print(f"AUC(-NU -> recovery) in-band, STRATIFIED (same strata)      "
              f"          = {sa_nu:.4f}")
        print(f"\nVerdict: H25 {'HOLDS' if sa_mu >= 0.8 else 'FALSIFIED'} "
              f"(threshold 0.8, stratified statistic)")
    else:
        print("band is degenerate (all-pos or all-neg) -- cannot test H25 "
              "with a single band; falling back to per-eff-stratum bands.")

    print("\nper-eff-stratum in-band AUC(-mu):")
    print(f"{'eff':>5} {'band N':>7} {'pos':>5} {'neg':>5} {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p_ = [r for r in sub if r['ok']]
        n_ = [r for r in sub if not r['ok']]
        if p_ and n_:
            print(f"{eff:>5.2f} {len(sub):>7} {len(p_):>5} {len(n_):>5} "
                  f"{auc([r['mu'] for r in p_], [r['mu'] for r in n_]):>8.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>7} {len(p_):>5} {len(n_):>5} "
                  f"{'(degenerate)':>8}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Secondary: does the GS-profile step predict recovery / the wall?")
    print("step = log2||b*_m|| - log2||b*_0||  (block-boundary jump, W1b)")
    print("-" * 78)
    valid = [r for r in rows if not math.isnan(r['step'])]
    print(f"valid step values: {len(valid)}/{len(rows)}")
    pos = [r for r in valid if r['ok']]
    neg = [r for r in valid if not r['ok']]
    a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
    print(f"AUC(-step -> recovery), naive pooled = {a_step:.4f}  "
          f"(step -> 0 predicted to favour recovery if W1b's wall-crossing "
          f"picture is right)")
    print(f"AUC(+step -> recovery), naive pooled = {1 - a_step:.4f}")
    effs_nondeg = tuple(e for e in EFFS if e != 0.05)
    sa_step = stratified_auc(valid, 'step', effs_nondeg)
    print(f"AUC(-step -> recovery), STRATIFIED (excl eff=0.05) = {sa_step:.4f}")
    print(f"AUC(+step -> recovery), STRATIFIED (excl eff=0.05) = "
          f"{1 - sa_step:.4f}   <- sign that actually holds, see per-stratum "
          f"table below")

    print("\nper-eff-stratum AUC(step) [reporting whichever sign is >= 0.5]:")
    print(f"{'eff':>5} {'N':>5} {'rec':>7} {'AUC -step':>10} {'AUC +step':>10}")
    for eff in EFFS:
        sub = [r for r in valid if r['effq'] == eff]
        p_ = [r for r in sub if r['ok']]
        n_ = [r for r in sub if not r['ok']]
        if p_ and n_:
            a = auc([r['step'] for r in p_], [r['step'] for r in n_])
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p_))+'/'+str(len(sub)):>7} {a:>10.4f} {1-a:>10.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p_))+'/'+str(len(sub)):>7} {'(degenerate)':>10}")

    print(f"\nSpearman(step, NU) pooled  = {spearman([r['step'] for r in valid], [r['NU'] for r in valid]):.4f}")
    print(f"Spearman(step, mu) pooled  = {spearman([r['step'] for r in valid], [r['mu'] for r in valid]):.4f}")
    print("(if |Spearman| is near 1, step is a restatement of an existing "
          "quantity, not new information)")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
