"""
GLV-HNP Phase 2, Thread 25: is mu a genuine second coordinate alongside NU,
or is its apparent power (W5, 2026-08-07) entirely mediated by NU?

Pre-registered by the 2026-08-07 (Thread 24) log entry:

  H25: within the ambiguous NU band 1.04 <= NU <= 2.20 (where the nearest-
       plane certificate gives no answer either way), AUC(-mu -> Kannan-LLL
       recovery) stays >= 0.8.

  If yes: mu is a genuine second coordinate and (NU, mu) is the 2-parameter
          viability test Phase 2 has been missing.
  If no:  W5's stratified-by-eff result was an artifact of eff itself (mu
          tracks eff-driven bias strength, not a NU-orthogonal signal), and
          the closed form should be retired as a predictor.

Secondary (also pre-registered): W1b found the GS profile is m exact copies
of lambda_1(L2) followed by a K1-sensitive second block, with the step
between blocks vanishing right at the K1 wall.  Test whether

    step = log2(||b*_{m+1}||) - log2(||b*_1||)

(0-based: prof[m] is the first vector of the second block, prof[0] the
first of the head) predicts recovery better than NU or mu alone.

Data: identical collection to glv_hnp_phase2_gsprofile_strat.py (W5/W6) --
same 20 curves x 5 eff-strata x 5 seeds at 17 bits, dim 24, float GS
(justified by W0/W4: float-vs-exact NU relative error <1e-14 at this
dimension).  Re-collected here rather than reusing the prior run's stdout
because the prior run did not persist the row table; --dump-json fixes that
for future threads.

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

EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
M17 = 12

# Bracket from the 2026-08-07 W4 measurement (500 inst -> now re-derived below
# on this run's own data so the band is not hard-coded from a different seed
# set; the pre-registered H25 band is quoted for comparison).
H25_BAND_PREREG = (1.04, 2.20)


def collect(curves17):
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
                prof = r['prof']
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n),
                          'step': (math.log2(prof[M17]) - math.log2(prof[0]))
                                  if prof[0] > 0 and prof[M17] > 0 else None})
                rows.append(r)
    return rows


def dump_json(rows, path):
    keys = ('n', 'K1', 'eff', 'effq', 'ok', 'NU', 'mu', 'nuhat', 'lamstar',
            'step', 'argmax', 'k', 'l2', 'det2')
    with open(path, 'w') as f:
        json.dump([{k: r[k] for k in keys} for r in rows], f)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None)
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — does mu separate inside the ambiguous NU band? (H25)")
    print("=" * 78)

    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")

    t0 = time.time()
    rows = collect(curves17)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        dump_json(rows, args.dump_json)
        print(f"dumped {len(rows)} rows to {args.dump_json}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP H25: AUC(-mu -> recovery) restricted to the ambiguous NU band")
    print("-" * 78)

    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    lo, hi = H25_BAND_PREREG
    print(f"pre-registered band (2026-08-07, N=500): NU in [{lo}, {hi}]")

    band = [r for r in rows if lo <= r['NU'] <= hi]
    bpos = [r for r in band if r['ok']]
    bneg = [r for r in band if not r['ok']]
    print(f"this run's band population: N={len(band)}  pos={len(bpos)}  "
          f"neg={len(bneg)}")
    if bpos and bneg:
        a_mu_band = auc([r['mu'] for r in bpos], [r['mu'] for r in bneg])
        a_nh_band = auc([r['nuhat'] for r in bpos], [r['nuhat'] for r in bneg])
        a_nu_band = auc([r['NU'] for r in bpos], [r['NU'] for r in bneg])
        print(f"  AUC(-mu     -> recovery) in-band = {a_mu_band:.4f}")
        print(f"  AUC(-nu_hat -> recovery) in-band = {a_nh_band:.4f}")
        print(f"  AUC(-NU     -> recovery) in-band = {a_nu_band:.4f}  "
              f"(sanity: should be ~0.5, NU is the stratifier)")
        verdict = "HOLDS" if a_mu_band >= 0.8 else "FALSIFIED"
        print(f"\nH25 (AUC(-mu) >= 0.8 in-band): {verdict}  "
              f"(observed {a_mu_band:.4f})")
    else:
        print("  degenerate: band is one-class on this run's data, cannot "
              "score H25 as pre-registered.")

    if bpos and bneg:
        print("\nis nu_hat's in-band AUC (0.84 pooled) a real cross-eff "
              "effect or a W6-style averaging artifact? per-stratum, "
              "in-band:")
        print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC -mu':>8} "
              f"{'AUC -nu_hat':>11}")
        for eff in EFFS:
            sub = [r for r in band if r['effq'] == eff]
            pos = [r for r in sub if r['ok']]
            neg = [r for r in sub if not r['ok']]
            if not pos or not neg:
                print(f"{eff:>5.2f} {len(sub):>5} "
                      f"{str(len(pos))+'/'+str(len(sub)):>7} | (degenerate)")
                continue
            a_mu_s = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
            a_nh_s = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_mu_s:>8.4f} "
                  f"{a_nh_s:>11.4f}")

    # Repeat using the band re-derived from THIS run's own success/failure
    # extremes, in case the 2026-08-07 bracket doesn't reproduce with a
    # fresh curve draw (search_curves is randomized via sympy.nextprime scan
    # but not seeded, so the 20 curves differ run to run).
    p_nu = [r['NU'] for r in pooled_pos]
    n_nu = [r['NU'] for r in pooled_neg]
    if p_nu and n_nu:
        lo2, hi2 = min(n_nu), max(p_nu)
        print(f"\nthis run's own bracket: sufficient NU < {lo2:.3f}, "
              f"necessary NU > {hi2:.3f}")
        if hi2 > lo2:
            band2 = [r for r in rows if lo2 <= r['NU'] <= hi2]
            b2p = [r for r in band2 if r['ok']]
            b2n = [r for r in band2 if not r['ok']]
            if b2p and b2n:
                a_mu2 = auc([r['mu'] for r in b2p], [r['mu'] for r in b2n])
                print(f"own-bracket band: N={len(band2)} pos={len(b2p)} "
                      f"neg={len(b2n)}  AUC(-mu) = {a_mu2:.4f}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP H25b (secondary): step = log2||b*_{m+1}|| - log2||b*_1|| "
          "as a wall predictor")
    print("-" * 78)
    valid = [r for r in rows if r['step'] is not None]
    print(f"valid rows: {len(valid)}/{len(rows)}")
    vp = [r for r in valid if r['ok']]
    vn = [r for r in valid if not r['ok']]
    if vp and vn:
        a_step = auc([r['step'] for r in vp], [r['step'] for r in vn])
        print(f"AUC(-step -> recovery), pooled = {a_step:.4f}  "
              f"(step measures HOW FAR the second block has risen above the "
              f"lambda_1(L2) plateau; larger step -> harder instance, so "
              f"score sign matches NU/mu convention: -step)")
        print(f"AUC(+step -> recovery), pooled = {1-a_step:.4f}")
        print(f"Spearman(step, NU) = {spearman([r['step'] for r in valid], [r['NU'] for r in valid]):.4f}")
        print(f"Spearman(step, mu) = {spearman([r['step'] for r in valid], [r['mu'] for r in valid]):.4f}")

        print(f"\n{'eff':>5} {'N':>5} {'rec':>7} | {'AUC -step':>10} "
              f"{'AUC -mu':>8} {'AUC -NU':>8}")
        for eff in EFFS:
            sub = [r for r in valid if r['effq'] == eff]
            pos = [r for r in sub if r['ok']]
            neg = [r for r in sub if not r['ok']]
            if not pos or not neg:
                print(f"{eff:>5.2f} {len(sub):>5} "
                      f"{str(len(pos))+'/'+str(len(sub)):>7} | (degenerate)")
                continue
            a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
            a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
            a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_st:>10.4f} "
                  f"{a_mu:>8.4f} {a_nu:>8.4f}")
    else:
        print("degenerate: one class only, cannot score.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
