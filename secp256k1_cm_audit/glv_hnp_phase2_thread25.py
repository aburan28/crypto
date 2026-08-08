"""
GLV-HNP Phase 2, Thread 25: find the second mechanism by conditioning on NU.

Pre-registered by the 2026-08-07 (Thread 24) log entry.  W5/W6 established
that Kannan-LLL recovery = f(NU, X) where X is mu-driven (mu = lambda_1(L2))
and X is uncorrelated with NU at fixed eff (Spearman in [-0.28, +0.16] across
strata).  NU is a sound nearest-plane certificate (0 FP at NU <= 1, 410
instances) but a size-degrading separator (AUC 0.978 -> 0.860, 12 -> 17
bits).  The question this script answers: inside the band where NU gives no
answer, does mu (or the GS-profile "step") still separate recoverable from
non-recoverable instances?

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (17-bit bracket from
       Thread 24 W4), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.
  Falsifier: AUC collapses to ~0.5 inside the band -> mu's apparent power in
       W5 was mediated entirely by NU (a stratification artifact), and the
       closed form should be retired as a genuine second coordinate.

Secondary (also proposed in the Thread 24 W1b note): the GS profile of L0 is
m exact copies of lambda_1(L2) followed by a second block that moves with
K1; the step from block 1 to block 2,

    step = log2(||b*_{m+1}||) - log2(||b*_1||)

vanishes right as the K1 wall is crossed (W1b, 12-bit qualitative).  Test
whether step is a better cross-curve separator than NU or mu, quantitatively,
at 17 bits with eff held fixed (same design as Thread 24 W5).

Same instance generation as glv_hnp_phase2_gsprofile_strat.py (float GS,
justified by W0/W4 there): 17-bit curves, M=12 (dim 24), 5 eff strata x 20
curves x 5 seeds.  --dump-json writes the row table so a later thread does
not have to regenerate it.

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

NU_BAND = (1.040, 2.199)  # Thread 24 W4, 17-bit ambiguous bracket


def logistic_fit(xs, ys, labels, lr=0.5, iters=4000):
    """2-feature logistic regression by plain gradient descent (no numpy/
    sklearn dependency).  xs: list of (x1, x2), ys: list of 0/1 labels.
    Returns (w1, w2, b)."""
    n = len(xs)
    mx1 = sum(x[0] for x in xs) / n
    mx2 = sum(x[1] for x in xs) / n
    sx1 = (sum((x[0] - mx1) ** 2 for x in xs) / n) ** 0.5 or 1.0
    sx2 = (sum((x[1] - mx2) ** 2 for x in xs) / n) ** 0.5 or 1.0
    zs = [((x[0] - mx1) / sx1, (x[1] - mx2) / sx2) for x in xs]

    w1 = w2 = b = 0.0
    for _ in range(iters):
        g1 = g2 = gb = 0.0
        for (z1, z2), y in zip(zs, ys):
            s = w1 * z1 + w2 * z2 + b
            p = 1.0 / (1.0 + math.exp(-s)) if s > -30 else 0.0
            err = p - y
            g1 += err * z1
            g2 += err * z2
            gb += err
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
        b -= lr * gb / n
    # un-standardize: score = w1*(x1-mx1)/sx1 + w2*(x2-mx2)/sx2 + b
    #               = (w1/sx1)*x1 + (w2/sx2)*x2 + (b - w1*mx1/sx1 - w2*mx2/sx2)
    W1 = w1 / sx1
    W2 = w2 / sx2
    B = b - w1 * mx1 / sx1 - w2 * mx2 / sx2

    preds = []
    for (x1, x2) in xs:
        s = W1 * x1 + W2 * x2 + B
        preds.append(1.0 / (1.0 + math.exp(-s)) if s > -30 else 0.0)
    correct = sum(1 for p, y in zip(preds, ys) if (p >= 0.5) == (y >= 0.5))
    return W1, W2, B, correct / n


if __name__ == "__main__":
    dump_path = None
    if "--dump-json" in sys.argv:
        i = sys.argv.index("--dump-json")
        dump_path = sys.argv[i + 1] if i + 1 < len(sys.argv) else "thread25_rows.json"

    print("=" * 78)
    print("Thread 25 — conditioning on NU: is mu a genuine second coordinate?")
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
                step = math.log2(r['prof'][M17]) - math.log2(r['prof'][0]) \
                    if r['prof'][0] > 0 and r['prof'][M17] > 0 else float('nan')
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if dump_path:
        keep = ['n', 'K1', 'eff', 'effq', 'ok', 'NU', 'mu', 'l2', 'det2',
                'nuhat', 'lamstar', 'step', 'k', 'argmax']
        with open(dump_path, 'w') as f:
            json.dump([{k: r[k] for k in keep} for r in rows], f)
        print(f"wrote {len(rows)} rows to {dump_path}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print(f"H25: inside the ambiguous NU band [{NU_BAND[0]}, {NU_BAND[1]}], "
          f"does mu still separate?")
    print("-" * 78)
    band = [r for r in rows if NU_BAND[0] <= r['NU'] <= NU_BAND[1]]
    outside = [r for r in rows if not (NU_BAND[0] <= r['NU'] <= NU_BAND[1])]
    print(f"in band: {len(band)}/{len(rows)}   outside: {len(outside)}")
    bpos = [r for r in band if r['ok']]
    bneg = [r for r in band if not r['ok']]
    print(f"band recovery: {len(bpos)}/{len(band)}")
    if bpos and bneg:
        a_mu = auc([r['mu'] for r in bpos], [r['mu'] for r in bneg])
        a_nh = auc([r['nuhat'] for r in bpos], [r['nuhat'] for r in bneg])
        a_nu = auc([r['NU'] for r in bpos], [r['NU'] for r in bneg])
        a_st = auc([r['step'] for r in bpos], [r['step'] for r in bneg])
        a_ls = auc([r['lamstar'] for r in bpos], [r['lamstar'] for r in bneg])
        print(f"  AUC(-mu     -> recovery) = {a_mu:.4f}")
        print(f"  AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"  AUC(-NU     -> recovery) = {a_nu:.4f}  (control: should be "
              f"~0.5, NU is constant-ish inside its own band)")
        print(f"  AUC(-step   -> recovery) = {a_st:.4f}")
        print(f"  AUC(-lam*   -> recovery) = {a_ls:.4f}  (control)")
        verdict = "HOLDS" if a_mu >= 0.8 or a_mu <= 0.2 else "FALSIFIED"
        print(f"\nH25 verdict: {verdict} (threshold AUC>=0.8 or <=0.2 for a "
              f"clean second coordinate; observed {a_mu:.4f})")
    else:
        print("degenerate band (all-pos or all-neg); cannot compute AUC")

    print("\nper-eff-stratum, restricted to the band:")
    print(f"{'eff':>5} {'N':>4} {'rec':>7} | {'AUC mu':>8} {'AUC step':>9}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if len(sub) < 4 or not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>4} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | (too few / degenerate)")
            continue
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"{eff:>5.2f} {len(sub):>4} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_mu:>8.4f} {a_st:>9.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Secondary: step = log2||b*_{m+1}|| - log2||b*_1|| as a wall predictor")
    print("-" * 78)
    print("Pooled over all rows (5 strata x 20 curves x 5 seeds, N="
          f"{len(rows)}):")
    pos = [r for r in rows if r['ok']]
    neg = [r for r in rows if not r['ok']]
    print(f"  AUC(-step -> recovery) = {auc([r['step'] for r in pos], [r['step'] for r in neg]):.4f}")
    print(f"  AUC(-mu   -> recovery) = {auc([r['mu'] for r in pos], [r['mu'] for r in neg]):.4f}")
    print(f"  AUC(-NU   -> recovery) = {auc([r['NU'] for r in pos], [r['NU'] for r in neg]):.4f}")
    print(f"  step | success: mean {sum(r['step'] for r in pos)/len(pos):.3f}"
          f"  step | failure: mean {sum(r['step'] for r in neg)/len(neg):.3f}")
    print(f"  Spearman(step, NU) = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"  Spearman(step, mu) = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\nper-eff-stratum (full table, not just the band):")
    print(f"{'eff':>5} {'N':>4} {'rec':>7} | {'AUC step':>9}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            continue
        print(f"{eff:>5.2f} {len(sub):>4} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | "
              f"{auc([r['step'] for r in pos], [r['step'] for r in neg]):>9.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Logistic fit on (log NU, log mu) -> recovery (plain gradient descent)")
    print("-" * 78)
    xs = [(math.log(r['NU']), math.log(r['mu'])) for r in rows if r['mu'] > 0]
    ys = [1.0 if r['ok'] else 0.0 for r in rows if r['mu'] > 0]
    W1, W2, B, acc = logistic_fit(xs, ys, None)
    print(f"score = {W1:.4f}*log(NU) + {W2:.4f}*log(mu) + {B:.4f}")
    print(f"training accuracy (0.5 threshold) = {acc:.4f}  (N={len(xs)}, "
          f"base rate {sum(ys)/len(ys):.3f})")
    scores = [W1 * x1 + W2 * x2 + B for x1, x2 in xs]
    spos = [s for s, y in zip(scores, ys) if y == 1.0]
    sneg = [s for s, y in zip(scores, ys) if y == 0.0]
    print(f"AUC(score -> recovery) = {auc(sneg, spos):.4f}  "
          f"(score is POSITIVE-oriented here, so pos/neg swapped vs auc()'s "
          f"usual '-x' convention)")
    print(f"compare: AUC(-NU alone) pooled = "
          f"{auc([r['NU'] for r in rows if r['ok']], [r['NU'] for r in rows if not r['ok']]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
