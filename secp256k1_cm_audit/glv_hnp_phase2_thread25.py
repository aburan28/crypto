"""
GLV-HNP Phase 2, Thread 25: is there a second coordinate beyond NU?

W5/W6 (2026-08-07, glv_hnp_phase2_gsprofile_strat.py) established:
  - NU (exact BDD certificate) is a sound size-free CERTIFICATE (0 FP at
    NU <= 1 over 410 instances) but a size-DEGRADING separator: the
    ambiguous band widens from [1.19, 1.87] at 12 bits to [1.04, 2.20] at
    17 bits.
  - mu = lambda_1(L2) (equivalently nu_hat) does real cross-curve work with
    eff held fixed (AUC 0.75-0.93), and is UNCORRELATED with NU at fixed eff
    (Spearman -0.28..+0.16 across strata; W6). So mu and NU are not the same
    quantity wearing different clothes.

H25 (pre-registered by the 2026-08-07 #2 log entry):
    Within the ambiguous band 1.04 <= NU <= 2.20 (where the NU <= 1
    certificate gives no answer either way), AUC(-mu -> Kannan-LLL recovery)
    stays >= 0.8.
  If yes: (NU, mu) is a genuine 2-parameter viability test -- NU screens the
    clear cases, mu resolves the ambiguous middle.
  If no: mu's apparent power in W5 was entirely mediated by NU (i.e. an eff
    stratum with more recoveries also has smaller mu on average, and W5's
    "eff fixed" control did not fully kill that channel), and the closed
    form should be retired.

Secondary (W1b): the GS profile's head is m exact copies of lambda_1(L2);
the step to the second block vanishes as the K1 wall is crossed. Define

    step_i = log2(prof[m]) - log2(prof[0])     (block-1 head vs block-2 head)

and test whether step predicts recovery at least as well as NU or mu.

Run: python3 glv_hnp_phase2_thread25.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman


def logistic_fit_2d(xs, ys, labels, iters=4000, lr=0.5):
    """Minimal 2-feature logistic regression by gradient descent, no deps.
    xs, ys: standardised features. labels: 1=recovered, 0=not."""
    mx, sx = sum(xs) / len(xs), (sum((v - sum(xs) / len(xs)) ** 2 for v in xs) / len(xs)) ** 0.5
    my, sy = sum(ys) / len(ys), (sum((v - sum(ys) / len(ys)) ** 2 for v in ys) / len(ys)) ** 0.5
    sx = sx or 1.0
    sy = sy or 1.0
    zx = [(v - mx) / sx for v in xs]
    zy = [(v - my) / sy for v in ys]
    w0 = w1 = w2 = 0.0
    n = len(labels)
    for _ in range(iters):
        g0 = g1 = g2 = 0.0
        for a, b, t in zip(zx, zy, labels):
            z = w0 + w1 * a + w2 * b
            p = 1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0
            err = p - t
            g0 += err
            g1 += err * a
            g2 += err * b
        w0 -= lr * g0 / n
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
    return w0, w1, w2, (mx, sx, my, sy)


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — is mu a second coordinate beyond NU? (H25)")
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
                step = math.log2(r['prof'][M17]) - math.log2(r['prof'][0])
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    NU_LO, NU_HI = 1.040, 2.199  # 17-bit bracket from W4 (2026-08-07 #2)

    print("\n" + "-" * 78)
    print(f"H25: within the ambiguous NU band [{NU_LO:.3f}, {NU_HI:.3f}], "
          "does mu still separate?")
    print("-" * 78)

    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band size: {len(band)}/{len(rows)}  "
          f"({len(pos)} recovered, {len(neg)} not)")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"AUC(-mu -> recovery)      within band = {a_mu:.4f}")
        print(f"AUC(-nu_hat -> recovery)  within band = {a_nh:.4f}")
        print(f"AUC(-NU -> recovery)      within band = {a_nu:.4f}  "
              "(expected ~0.5: NU is the stratifying variable)")
        print(f"AUC(-step -> recovery)    within band = {a_step:.4f}")
        print(f"\nH25 verdict: {'HOLDS' if a_mu >= 0.8 else 'FALSIFIED'} "
              f"(threshold 0.8, observed {a_mu:.4f})")
    else:
        print("degenerate band (all-recovered or all-failed): cannot test H25")

    print("\n" + "-" * 78)
    print("H25 per-eff-stratum breakdown (band membership varies by stratum)")
    print("-" * 78)
    print(f"{'eff':>5} {'band N':>7} {'rec':>7} {'AUC mu':>8} {'AUC step':>9}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        q = [r for r in sub if not r['ok']]
        if not p or not q:
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {'(degen)':>8}")
            continue
        a_mu = auc([r['mu'] for r in p], [r['mu'] for r in q])
        a_step = auc([r['step'] for r in p], [r['step'] for r in q])
        print(f"{eff:>5.2f} {len(sub):>7} "
              f"{str(len(p))+'/'+str(len(sub)):>7} {a_mu:>8.4f} {a_step:>9.4f}")

    print("\n" + "-" * 78)
    print("Secondary: step = log2(||b*_{m+1}||) - log2(||b*_1||), pooled")
    print("-" * 78)
    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    a_step_all = auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg])
    a_nu_all = auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg])
    a_mu_all = auc([r['mu'] for r in pooled_pos], [r['mu'] for r in pooled_neg])
    print(f"AUC(-step -> recovery) pooled (N={len(rows)}) = {a_step_all:.4f}")
    print(f"AUC(-NU   -> recovery) pooled                 = {a_nu_all:.4f}")
    print(f"AUC(-mu   -> recovery) pooled                 = {a_mu_all:.4f}")
    print(f"Spearman(step, NU) pooled = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"step | success : mean {sum(r['step'] for r in pooled_pos)/len(pooled_pos):.3f}  "
          f"min {min(r['step'] for r in pooled_pos):.3f}  "
          f"max {max(r['step'] for r in pooled_pos):.3f}")
    print(f"step | failure : mean {sum(r['step'] for r in pooled_neg)/len(pooled_neg):.3f}  "
          f"min {min(r['step'] for r in pooled_neg):.3f}  "
          f"max {max(r['step'] for r in pooled_neg):.3f}")

    print("\n" + "-" * 78)
    print("2-parameter logistic fit on (log NU, log mu), full 500-instance set")
    print("-" * 78)
    log_nu = [math.log(r['NU']) for r in rows]
    log_mu = [math.log(r['mu']) for r in rows]
    labels = [1.0 if r['ok'] else 0.0 for r in rows]
    w0, w1, w2, (mx, sx, my, sy) = logistic_fit_2d(log_nu, log_mu, labels)
    preds = []
    for a, b in zip(log_nu, log_mu):
        z = w0 + w1 * (a - mx) / sx + w2 * (b - my) / sy
        preds.append(1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0)
    fit_pos = [p for p, t in zip(preds, labels) if t == 1.0]
    fit_neg = [p for p, t in zip(preds, labels) if t == 0.0]
    a_fit = auc(fit_neg, fit_pos)  # higher predicted p -> recovery, so flip
    print(f"weights (standardised): w0={w0:.3f} w_logNU={w1:.3f} w_logmu={w2:.3f}")
    print(f"AUC of fitted p(recovery) = {a_fit:.4f}  "
          f"(vs AUC(-NU) alone = {a_nu_all:.4f}, AUC(-mu) alone = {a_mu_all:.4f})")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
