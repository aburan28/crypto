"""
GLV-HNP Phase 2, Thread 25: does mu carry a second, NU-independent signal?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry.  W5/W6 showed
recovery = f(NU, X) with X ~ mu-driven and X uncorrelated with NU (Spearman
in [-0.28, +0.16] across strata).  This session asks the concrete question:

  H25: within the ambiguous NU band [1.04, 2.20] (17 bits, ambiguous meaning
       neither W4's sufficient-NU<1.04 nor necessary-NU>2.20 bracket fires),
       AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

If yes: mu is a genuine second coordinate, distinct from the NU nearest-plane
certificate, and (NU, mu) is a 2-parameter viability test.  Fit and report a
logistic decision boundary on (log NU, log mu).

Secondary (also pre-registered): W1b found the GS profile head is m exact
copies of lambda_1(L2), and the step to the second block *vanishes* exactly
as the K1 wall is crossed.  Define
    step = log2(||b*_{m+1}||) - log2(||b*_1||)     (0-indexed: prof[m]-prof[0])
and test whether step predicts recovery better than NU or mu alone.

Data: identical generation to glv_hnp_phase2_gsprofile_strat.py (17-bit j=0
GLV curves, dim 24, 5 eff strata x SEEDS), so results are directly comparable
to the W4/W5/W6 tables.  Adds --dump-json so the row table survives the run
(requested by the 2026-08-07 #2 log entry's cost estimate).

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

# W4's 17-bit bracket: sufficient NU < 1.040, necessary NU > 2.199.
NU_LO, NU_HI = 1.040, 2.199


def build_rows():
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
                step = (math.log2(r['prof'][r['k'] // 2]) -
                        math.log2(r['prof'][0])) if r['prof'][0] > 0 else float('nan')
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    return rows, len(curves17), M17


def logistic_fit_2d(xs, ys, labels, iters=4000, lr=0.5):
    """Plain-Python gradient-descent logistic regression on standardized
    (x, y) -> label in {0, 1}.  Returns (w0, w1, w2, mx, sx, my, sy) for
        p = sigma(w0 + w1*(x-mx)/sx + w2*(y-my)/sy)."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n) or 1.0
    sy = math.sqrt(sum((y - my) ** 2 for y in ys) / n) or 1.0
    zx = [(x - mx) / sx for x in xs]
    zy = [(y - my) / sy for y in ys]
    w0 = w1 = w2 = 0.0
    for _ in range(iters):
        g0 = g1 = g2 = 0.0
        for x, y, lab in zip(zx, zy, labels):
            z = w0 + w1 * x + w2 * y
            p = 1.0 / (1.0 + math.exp(-z)) if z > -50 else 0.0
            err = p - lab
            g0 += err
            g1 += err * x
            g2 += err * y
        w0 -= lr * g0 / n
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
    return w0, w1, w2, mx, sx, my, sy


def logistic_predict(w, x, y):
    w0, w1, w2, mx, sx, my, sy = w
    z = w0 + w1 * (x - mx) / sx + w2 * (y - my) / sy
    return 1.0 / (1.0 + math.exp(-z)) if z > -50 else 0.0


if __name__ == "__main__":
    dump_path = None
    if "--dump-json" in sys.argv:
        i = sys.argv.index("--dump-json")
        dump_path = sys.argv[i + 1] if i + 1 < len(sys.argv) else "thread25_rows.json"

    print("=" * 78)
    print("Thread 25 — does mu separate inside the NU-ambiguous band?")
    print("=" * 78)

    t0 = time.time()
    rows, ncurves, M17 = build_rows()
    print(f"\n{ncurves} 17-bit curves, {len(rows)} instances (dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if dump_path:
        with open(dump_path, "w") as f:
            json.dump([{k: v for k, v in r.items() if k not in ('prof', 'nus')}
                       for r in rows], f)
        print(f"dumped {len(rows)} rows -> {dump_path}")

    print("\n" + "-" * 78)
    print(f"H25: AUC(-mu -> recovery) inside ambiguous band NU in [{NU_LO}, {NU_HI}]")
    print("-" * 78)
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    pos_all = [r for r in rows if r['ok']]
    neg_all = [r for r in rows if not r['ok']]
    print(f"full table: N={len(rows)}  rec={len(pos_all)}/{len(rows)}  "
          f"AUC(-NU)={auc([r['NU'] for r in pos_all], [r['NU'] for r in neg_all]):.4f}  "
          f"AUC(-mu)={auc([r['mu'] for r in pos_all], [r['mu'] for r in neg_all]):.4f}")

    pos_b = [r for r in band if r['ok']]
    neg_b = [r for r in band if not r['ok']]
    print(f"\nband [{NU_LO},{NU_HI}]: N={len(band)}  rec={len(pos_b)}/{len(band)}")
    if pos_b and neg_b:
        a_mu = auc([r['mu'] for r in pos_b], [r['mu'] for r in neg_b])
        a_nu = auc([r['NU'] for r in pos_b], [r['NU'] for r in neg_b])
        a_step = auc([r['step'] for r in pos_b], [r['step'] for r in neg_b])
        print(f"  AUC(-mu -> rec)   = {a_mu:.4f}   (H25 threshold: >= 0.80)")
        print(f"  AUC(-NU -> rec)   = {a_nu:.4f}   (sanity: should be ~0.5, band is where NU is uninformative)")
        print(f"  AUC(-step -> rec) = {a_step:.4f}")
        print(f"\n  H25 verdict: {'HOLDS' if a_mu >= 0.80 else 'FALSIFIED'} "
              f"(AUC(-mu)={a_mu:.4f} {'>=' if a_mu >= 0.80 else '<'} 0.80)")
    else:
        print("  degenerate band (all-pos or all-neg); cannot compute AUC")

    print("\nper-stratum band breakdown:")
    print(f"{'eff':>5} {'N_band':>7} {'rec':>7} {'AUC(-mu)':>9} {'AUC(-NU)':>9}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in band if r['effq'] == eff]
        if not sub:
            continue
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if p and ng:
            print(f"{eff:>5.2f} {len(sub):>7} {str(len(p))+'/'+str(len(sub)):>7} "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>9.4f} "
                  f"{auc([r['NU'] for r in p], [r['NU'] for r in ng]):>9.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>7} {str(len(p))+'/'+str(len(sub)):>7} "
                  f"{'(degen)':>9} {'(degen)':>9}")

    print("\n" + "-" * 78)
    print("Logistic fit on (log NU, log mu) -> recovery, full table")
    print("-" * 78)
    xs = [math.log(r['NU']) for r in rows]
    ys = [math.log(r['mu']) for r in rows if r['mu'] > 0]
    # guard: mu should always be > 0 (lambda_1 of a nontrivial lattice)
    assert len(ys) == len(rows), "mu <= 0 encountered, investigate"
    labels = [1.0 if r['ok'] else 0.0 for r in rows]
    w = logistic_fit_2d(xs, ys, labels)
    w0, w1, w2, mx, sx, my, sy = w
    preds = [logistic_predict(w, x, y) for x, y in zip(xs, ys)]
    # auc() convention (matches every other score in this codebase): SMALLER
    # score -> recovery.  Predicted probability p is the opposite orientation
    # (larger p -> recovery), so negate it before calling auc().
    a_joint = auc([-p for p, l in zip(preds, labels) if l == 1.0],
                   [-p for p, l in zip(preds, labels) if l == 0.0])
    print(f"  standardized weights: w0={w0:.4f}  w_logNU={w1:.4f}  w_logmu={w2:.4f}")
    print(f"  (logNU mean/std = {mx:.4f}/{sx:.4f}, logmu mean/std = {my:.4f}/{sy:.4f})")
    print(f"  AUC(joint predicted prob -> recovery) = {a_joint:.4f}")
    print(f"  decision boundary (p=0.5): "
          f"{w1/sx:.4f}*logNU + {w2/sy:.4f}*logmu = "
          f"{w1*mx/sx + w2*my/sy - w0:.4f}")

    print("\n" + "-" * 78)
    print("Secondary: does step = log2(prof[m]) - log2(prof[0]) predict the wall?")
    print("-" * 78)
    valid = [r for r in rows if not math.isnan(r['step'])]
    p = [r for r in valid if r['ok']]
    ng = [r for r in valid if not r['ok']]
    a_step_full = auc([r['step'] for r in p], [r['step'] for r in ng])
    a_nu_full = auc([r['NU'] for r in p], [r['NU'] for r in ng])
    a_mu_full = auc([r['mu'] for r in p], [r['mu'] for r in ng])
    print(f"  full table AUC(-step -> rec) = {a_step_full:.4f}   "
          f"(cf. AUC(-NU)={a_nu_full:.4f}, AUC(-mu)={a_mu_full:.4f})")
    print(f"  Spearman(step, NU)  = {spearman([r['step'] for r in valid], [r['NU'] for r in valid]):.4f}")
    print(f"  Spearman(step, mu)  = {spearman([r['step'] for r in valid], [r['mu'] for r in valid]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
