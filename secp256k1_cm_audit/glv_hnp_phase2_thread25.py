"""
GLV-HNP Phase 2, Thread 25: find the second mechanism by conditioning on NU.

Pre-registered by the 2026-08-07 (Thread 24) log entry:

  H25: within the ambiguous NU band 1.04 <= NU <= 2.20 (where nearest-plane
       gives no answer -- W4's 17-bit bracket), AUC(-mu -> Kannan-LLL
       recovery) stays >= 0.8.
  Falsifier: if AUC(-mu) drops to ~0.5 inside the band, mu's apparent power
       (W5, AUC 0.75-0.93 per eff-stratum) was entirely mediated by NU, and
       the closed-form nu_hat/mu route should be retired.

Secondary (also pre-registered): W1b showed the GS profile head is m exact
copies of lambda_1(L2) and the step to the second block vanishes right as
the K1 wall is crossed. Define
    step = log2(||b*_{m+1}||) - log2(||b*_1||)
(0-indexed: prof[m] is the first index of the second block, prof[0] the
first index of the first block) and test whether step predicts recovery
better than NU or mu.

Data: same recipe as glv_hnp_phase2_gsprofile_strat.py W5/W6 (17-bit j=0 GLV
curves, m=12, 5 eff strata x 20 curves x 5 seeds), regenerated here rather
than reusing a cached table -- per the 2026-08-07 log's own proposal, this
run also writes the table to JSON via --dump-json so future runs can reuse
it without paying the ~generation cost again.

Run: python3 glv_hnp_phase2_thread25.py [--dump-json out.json]
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


def collect(curves17, m, effs):
    rows = []
    for eff in effs:
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
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n)})
                r['step'] = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                             if r['prof'][m] > 0 and r['prof'][0] > 0
                             else float('nan'))
                # keep the JSON dump small: drop the per-index profile/nus
                r.pop('prof', None)
                r.pop('nus', None)
                rows.append(r)
    return rows


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — conditioning on NU: is mu a genuine second coordinate?")
    print("=" * 78)

    dump_path = None
    if "--dump-json" in sys.argv:
        i = sys.argv.index("--dump-json")
        dump_path = sys.argv[i + 1] if i + 1 < len(sys.argv) else \
            "glv_hnp_phase2_thread25_table.json"

    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")
    M17 = 12
    EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)

    t0 = time.time()
    rows = collect(curves17, M17, EFFS)
    print(f"{len(rows)} instances (float GS, dim {2*M17}) "
          f"in {time.time()-t0:.1f}s")

    if dump_path:
        with open(dump_path, "w") as f:
            json.dump(rows, f)
        print(f"table dumped to {dump_path} ({os.path.getsize(dump_path)} bytes)")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1 (H25): AUC(-mu -> recovery) INSIDE the ambiguous NU band")
    print("-" * 78)
    BAND_LO, BAND_HI = 1.040, 2.199  # W4's 17-bit bracket
    band = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    pos_b = [r for r in band if r['ok']]
    neg_b = [r for r in band if not r['ok']]
    print(f"band [{BAND_LO}, {BAND_HI}]: N={len(band)}  "
          f"pos={len(pos_b)}  neg={len(neg_b)}")
    if pos_b and neg_b:
        a_mu_band = auc([r['mu'] for r in pos_b], [r['mu'] for r in neg_b])
        a_nh_band = auc([r['nuhat'] for r in pos_b], [r['nuhat'] for r in neg_b])
        a_nu_band = auc([r['NU'] for r in pos_b], [r['NU'] for r in neg_b])
        print(f"  AUC(-mu     -> recovery) INSIDE band = {a_mu_band:.4f}")
        print(f"  AUC(-nu_hat -> recovery) INSIDE band = {a_nh_band:.4f}")
        print(f"  AUC(-NU     -> recovery) INSIDE band = {a_nu_band:.4f}  "
              f"[should be ~0.5: band is where NU alone can't decide]")
        print(f"\n  H25 threshold 0.8: mu {'HOLDS' if a_mu_band >= 0.8 else 'FAILS'} "
              f"({a_mu_band:.4f})")
    else:
        print("  degenerate band (all-pos or all-neg) -- cannot evaluate H25")

    print("\nfor comparison, whole-table AUCs (unconditioned):")
    pos_a = [r for r in rows if r['ok']]
    neg_a = [r for r in rows if not r['ok']]
    print(f"  AUC(-mu)     = {auc([r['mu'] for r in pos_a], [r['mu'] for r in neg_a]):.4f}")
    print(f"  AUC(-nu_hat) = {auc([r['nuhat'] for r in pos_a], [r['nuhat'] for r in neg_a]):.4f}")
    print(f"  AUC(-NU)     = {auc([r['NU'] for r in pos_a], [r['NU'] for r in neg_a]):.4f}")

    print("\nband AUC(-mu) by eff stratum (band may be thin per-stratum):")
    print(f"{'eff':>5} {'N':>5} {'pos':>4} {'neg':>4} {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ne = [r for r in sub if not r['ok']]
        if p and ne:
            print(f"{eff:>5.2f} {len(sub):>5} {len(p):>4} {len(ne):>4} "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in ne]):>8.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>5} {len(p):>4} {len(ne):>4} "
                  f"{'(degen)':>8}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X2: does step = log2||b*_{m+1}|| - log2||b*_1|| predict the wall?")
    print("-" * 78)
    valid = [r for r in rows if not math.isnan(r['step'])]
    pos_s = [r for r in valid if r['ok']]
    neg_s = [r for r in valid if not r['ok']]
    print(f"N valid = {len(valid)} (of {len(rows)})")
    if pos_s and neg_s:
        a_step = auc([r['step'] for r in pos_s], [r['step'] for r in neg_s])
        # step should be SMALL (near 0) when the wall is crossed -> failure,
        # so recovery should correlate with LARGE step. Report AUC of +step
        # (not -step) for recovery, i.e. 1 - the "-x" convention AUC.
        a_step_pos = 1.0 - a_step
        print(f"  AUC(-step -> recovery) = {a_step:.4f}")
        print(f"  AUC(+step -> recovery) = {a_step_pos:.4f}  "
              f"[step measured to vanish AT the wall, so large step = healthy]")
        print(f"  step | success: mean {sum(r['step'] for r in pos_s)/len(pos_s):.3f}  "
              f"median {sorted(r['step'] for r in pos_s)[len(pos_s)//2]:.3f}")
        print(f"  step | failure: mean {sum(r['step'] for r in neg_s)/len(neg_s):.3f}  "
              f"median {sorted(r['step'] for r in neg_s)[len(neg_s)//2]:.3f}")
        print(f"  Spearman(step, NU)     = "
              f"{spearman([r['step'] for r in valid], [r['NU'] for r in valid]):.4f}")
        print(f"  Spearman(step, mu)     = "
              f"{spearman([r['step'] for r in valid], [r['mu'] for r in valid]):.4f}")
        print(f"  Spearman(step, nu_hat) = "
              f"{spearman([r['step'] for r in valid], [r['nuhat'] for r in valid]):.4f}")
        best = max(a_step_pos, 1 - a_step)  # already handled sign; keep simple
        print(f"\n  step vs NU (0.860) vs mu ({auc([r['mu'] for r in pos_a], [r['mu'] for r in neg_a]):.4f}) "
              f"pooled AUC comparison: step={max(a_step, a_step_pos):.4f}")
    else:
        print("  degenerate: cannot evaluate step AUC")

    print("\nstep by eff stratum:")
    print(f"{'eff':>5} {'N':>5} {'AUC(+step)':>10}")
    for eff in EFFS:
        sub = [r for r in valid if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ne = [r for r in sub if not r['ok']]
        if p and ne:
            a = 1.0 - auc([r['step'] for r in p], [r['step'] for r in ne])
            print(f"{eff:>5.2f} {len(sub):>5} {a:>10.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>5} {'(degen)':>10}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X3: logistic decision boundary on (log NU, log nu_hat)")
    print("-" * 78)
    print("X1 found nu_hat (AUC 0.840 in-band) succeeds where raw mu (0.693)")
    print("needs per-stratum conditioning, so nu_hat -- not mu -- is the")
    print("second coordinate. Fit w0 + w1*logNU + w2*log(nu_hat) by gradient")
    print("descent on the FULL 500-row table (logistic regression, no deps).\n")

    xs = [(math.log(r['NU']), math.log(r['nuhat'])) for r in rows]
    ys = [1.0 if r['ok'] else 0.0 for r in rows]
    mx0 = sum(x[0] for x in xs) / len(xs)
    mx1 = sum(x[1] for x in xs) / len(xs)
    sx0 = math.sqrt(sum((x[0] - mx0) ** 2 for x in xs) / len(xs))
    sx1 = math.sqrt(sum((x[1] - mx1) ** 2 for x in xs) / len(xs))
    zs = [((x[0] - mx0) / sx0, (x[1] - mx1) / sx1) for x in xs]

    w = [0.0, 0.0, 0.0]  # bias, w_logNU, w_lognuhat (standardized units)
    lr, n = 0.3, len(zs)
    for epoch in range(3000):
        g = [0.0, 0.0, 0.0]
        for (z0, z1), y in zip(zs, ys):
            p = 1.0 / (1.0 + math.exp(-(w[0] + w[1] * z0 + w[2] * z1)))
            err = p - y
            g[0] += err
            g[1] += err * z0
            g[2] += err * z1
        w = [w[i] - lr * g[i] / n for i in range(3)]

    def score(z0, z1):
        return w[0] + w[1] * z0 + w[2] * z1

    scores = [score(z0, z1) for z0, z1 in zs]
    pos_sc = [s for s, y in zip(scores, ys) if y == 1.0]
    neg_sc = [s for s, y in zip(scores, ys) if y == 0.0]
    # auc() expects "-x separates pos from neg"; the logit is oriented so
    # LARGER score -> more likely recovery, so feed -scores.
    a_joint = auc([-s for s in pos_sc], [-s for s in neg_sc])
    print(f"standardized fit: logit = {w[0]:.3f} + {w[1]:.3f}*z(logNU) "
          f"+ {w[2]:.3f}*z(log nu_hat)")
    print(f"  z(logNU)     mean={mx0:.4f} sd={sx0:.4f}")
    print(f"  z(lognu_hat) mean={mx1:.4f} sd={sx1:.4f}")
    print(f"training AUC (joint, full table) = {a_joint:.4f}")
    print(f"  vs AUC(-NU) alone (full table)     = "
          f"{auc([r['NU'] for r in pos_a], [r['NU'] for r in neg_a]):.4f}")
    print(f"  vs AUC(-nu_hat) alone (full table) = "
          f"{auc([r['nuhat'] for r in pos_a], [r['nuhat'] for r in neg_a]):.4f}")

    band_z = [((math.log(r['NU']) - mx0) / sx0,
                (math.log(r['nuhat']) - mx1) / sx1) for r in band]
    band_sc = [score(z0, z1) for z0, z1 in band_z]
    band_pos_sc = [s for s, r in zip(band_sc, band) if r['ok']]
    band_neg_sc = [s for s, r in zip(band_sc, band) if not r['ok']]
    if band_pos_sc and band_neg_sc:
        a_joint_band = auc([-s for s in band_pos_sc], [-s for s in band_neg_sc])
        print(f"\njoint-score AUC INSIDE the NU-ambiguous band = "
              f"{a_joint_band:.4f}  (nu_hat alone in-band: 0.8403)")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
