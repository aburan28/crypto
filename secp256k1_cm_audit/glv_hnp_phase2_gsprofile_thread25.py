"""
GLV-HNP Phase 2, Thread 25: does mu carry a second, NU-independent signal?

W5/W6 of glv_hnp_phase2_gsprofile_strat.py (Thread 24b, 2026-08-07) found that
at fixed eff, nu_hat (equivalently mu = lambda_1(L2)) separates recovery with
AUC 0.75-0.93 while the exact BDD certificate NU does NOT (AUC 0.35-0.73,
sometimes anti-predictive), and the two scores are uncorrelated within a
stratum (Spearman ~ -0.28..+0.16). Interpretation offered there: NU governs
Babai nearest-plane viability (a sound but separately-acting certificate);
mu/nu_hat tracks a second, distinct mechanism.

Pre-registered by the 2026-08-07 #2 log entry:

  H25: within the ambiguous NU band [1.04, 2.20] (17-bit bracket from W4,
       where nearest-plane theory gives no verdict either way), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.

  Falsifier: if AUC(-mu) inside the band collapses toward 0.5, mu's apparent
       power in W5 was mediated by NU after all (a stratification artifact:
       eff-matched strata still let NU vary, and mu was silently reading
       that), and the closed form should be retired as a redundant reading
       of NU.

Secondary (W1b of the parent Thread 24 script): the profile head is m exact
copies of lambda_1(L2); the step to the second GS block vanishes right as the
K1 wall is crossed. Define

    step = log2(||b*_{m}||) - log2(||b*_{m-1}||)     (0-indexed: prof[m]/prof[m-1])

(the jump from the last head entry to the first tail entry) and test whether
`step` predicts recovery at least as well as NU or mu.

Uses the same instance-generation and eff-stratified 17-bit grid as
glv_hnp_phase2_gsprofile_strat.py (search_curves(2**16, 2**17, per_bin=2,
nbins=10), M=12, EFFS = (0.05, 0.10, 0.15, 0.20, 0.25), SEEDS from
glv_hnp_phase2_projected), so this is a fresh re-collection of the same
population Thread 24b analysed, not a replay of stored numbers.

Run: python3 glv_hnp_phase2_gsprofile_thread25.py
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

NU_BAND = (1.040, 2.199)   # W4 (17-bit) ambiguous bracket, ./RESEARCH_AUTOLAB_LOG.md 2026-08-07 #2

if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — does mu separate INSIDE the NU-ambiguous band? (H25)")
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
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n),
                          'step': math.log2(r['prof'][m]) - math.log2(r['prof'][m - 1])
                          if r['prof'][m] > 0 and r['prof'][m - 1] > 0 else float('nan')})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    print("\n" + "-" * 78)
    print(f"H25: stratify by NU band [{NU_BAND[0]}, {NU_BAND[1]}] "
          "(the ambiguous zone from W4)")
    print("-" * 78)
    below = [r for r in rows if r['NU'] < NU_BAND[0]]
    inside = [r for r in rows if NU_BAND[0] <= r['NU'] <= NU_BAND[1]]
    above = [r for r in rows if r['NU'] > NU_BAND[1]]
    print(f"below band  (NU < {NU_BAND[0]}): N={len(below):4d}  "
          f"rec={sum(1 for r in below if r['ok'])}/{len(below)}  "
          "(NU already predicts recovery here)")
    print(f"ABOVE band  (NU > {NU_BAND[1]}): N={len(above):4d}  "
          f"rec={sum(1 for r in above if r['ok'])}/{len(above)}  "
          "(NU already predicts failure here)")
    print(f"INSIDE band ({NU_BAND[0]} <= NU <= {NU_BAND[1]}): N={len(inside):4d}  "
          f"rec={sum(1 for r in inside if r['ok'])}/{len(inside)}  "
          "<- nearest-plane theory gives no verdict here")

    pos = [r for r in inside if r['ok']]
    neg = [r for r in inside if not r['ok']]
    print(f"\ninside band: {len(pos)} recovered / {len(neg)} failed")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        a_ls = auc([r['lamstar'] for r in pos], [r['lamstar'] for r in neg])
        print(f"  AUC(-mu    -> recovery) = {a_mu:.4f}   <- H25 target (need >= 0.8)")
        print(f"  AUC(-nuhat -> recovery) = {a_nh:.4f}")
        print(f"  AUC(-NU    -> recovery) = {a_nu:.4f}   (near 0.5 expected: NU is uninformative here BY CONSTRUCTION)")
        print(f"  AUC(-step  -> recovery) = {a_step:.4f}")
        print(f"  AUC(-lam*  -> recovery) = {a_ls:.4f}   (control, Thread 20 falsified)")
        verdict = "CONFIRMED" if a_mu >= 0.8 else "FALSIFIED"
        print(f"\nH25 verdict: {verdict} (AUC(-mu) = {a_mu:.4f} vs threshold 0.8)")
    else:
        print("  band is empty on one side; H25 cannot be evaluated on this draw")

    print("\n" + "-" * 78)
    print("Cross-check: does mu separate INSIDE each eff stratum intersected")
    print("with the band? (rules out eff itself being smuggled back in)")
    print("-" * 78)
    print(f"{'eff':>5} {'N in band':>10} {'rec':>7} | {'AUC mu':>8} {'AUC NU':>8} "
          f"{'AUC step':>9}")
    for eff in EFFS:
        sub = [r for r in inside if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>10} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | {'(degenerate)':>8}")
            continue
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"{eff:>5.2f} {len(sub):>10} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_mu:>8.4f} {a_nu:>8.4f} "
              f"{a_step:>9.4f}")

    print("\n" + "-" * 78)
    print("Secondary: does `step` (jump from head block to tail block) predict")
    print("recovery pooled over all strata, unconditionally?")
    print("-" * 78)
    pooled_pos = [r for r in rows if r['ok'] and not math.isnan(r['step'])]
    pooled_neg = [r for r in rows if not r['ok'] and not math.isnan(r['step'])]
    print(f"AUC(-step -> recovery), pooled (N={len(pooled_pos)+len(pooled_neg)}) = "
          f"{auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg]):.4f}")
    print(f"AUC(-mu   -> recovery), pooled = "
          f"{auc([r['mu'] for r in pooled_pos], [r['mu'] for r in pooled_neg]):.4f}")
    print(f"AUC(-NU   -> recovery), pooled = "
          f"{auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg]):.4f}")
    print(f"Spearman(step, NU) pooled = "
          f"{spearman([r['step'] for r in rows if not math.isnan(r['step'])], [r['NU'] for r in rows if not math.isnan(r['step'])]):.4f}")
    print(f"Spearman(step, mu) pooled = "
          f"{spearman([r['step'] for r in rows if not math.isnan(r['step'])], [r['mu'] for r in rows if not math.isnan(r['step'])]):.4f}")

    print("\n" + "-" * 78)
    print("Joint fit: logistic regression on (log NU, log nuhat), pooled over")
    print("all 500 instances. Deliverable requested by the 2026-08-07 #2 log")
    print("entry if a second coordinate survives H25.")
    print("-" * 78)
    xs = [(math.log(r['NU']), math.log(r['nuhat'])) for r in rows if r['nuhat'] > 0]
    ys = [1.0 if r['ok'] else 0.0 for r in rows if r['nuhat'] > 0]
    mx0 = sum(x[0] for x in xs) / len(xs)
    mx1 = sum(x[1] for x in xs) / len(xs)
    sx0 = math.sqrt(sum((x[0] - mx0) ** 2 for x in xs) / len(xs))
    sx1 = math.sqrt(sum((x[1] - mx1) ** 2 for x in xs) / len(xs))
    zs = [((x[0] - mx0) / sx0, (x[1] - mx1) / sx1) for x in xs]

    w0, w1, b = 0.0, 0.0, 0.0
    lr = 0.5
    for it in range(3000):
        g0 = g1 = gb = 0.0
        for (z0, z1), y in zip(zs, ys):
            p = 1.0 / (1.0 + math.exp(-(w0 * z0 + w1 * z1 + b)))
            err = p - y
            g0 += err * z0
            g1 += err * z1
            gb += err
        n = len(zs)
        w0 -= lr * g0 / n
        w1 -= lr * g1 / n
        b -= lr * gb / n

    print(f"standardized coefficients: w(log NU)={w0:.4f}  "
          f"w(log nuhat)={w1:.4f}  bias={b:.4f}")
    print(f"(both should be negative: larger NU or larger nuhat -> less likely "
          "to recover)")
    scores = [-(w0 * z0 + w1 * z1) for z0, z1 in zs]
    pos_s = [s for s, y in zip(scores, ys) if y == 1.0]
    neg_s = [s for s, y in zip(scores, ys) if y == 0.0]
    print(f"AUC(joint logistic score -> recovery) = {auc(pos_s, neg_s):.4f}  "
          f"(N={len(zs)})")
    print(f"  vs AUC(-NU alone)    pooled = "
          f"{auc([z0 for z0, y in zip([z[0] for z in zs], ys) if y==1.0], [z0 for z0, y in zip([z[0] for z in zs], ys) if y==0.0]):.4f}")
    print(f"  vs AUC(-nuhat alone) pooled = "
          f"{auc([z1 for z1, y in zip([z[1] for z in zs], ys) if y==1.0], [z1 for z1, y in zip([z[1] for z in zs], ys) if y==0.0]):.4f}")
    decision_slope = -w0 / w1 if w1 != 0 else float('inf')
    print(f"\ndecision boundary in standardized space: "
          f"z(log NU) = {decision_slope:.4f} * z(log nuhat) + const")
    print("(a nonzero, non-infinite slope means the boundary is genuinely 2D, "
          "not axis-aligned on either variable alone)")

    print("\n" + "-" * 78)
    print("Held-out check: fit on 10 curves, evaluate joint AUC on the other 10")
    print("(curve-level split, not instance-level, so no leakage via seeds)")
    print("-" * 78)
    curve_ns = sorted({r['n'] for r in rows})
    train_ns = set(curve_ns[0::2])
    test_ns = set(curve_ns[1::2])
    tr = [r for r in rows if r['n'] in train_ns and r['nuhat'] > 0]
    te = [r for r in rows if r['n'] in test_ns and r['nuhat'] > 0]

    def fit_logistic(rs):
        xs = [(math.log(r['NU']), math.log(r['nuhat'])) for r in rs]
        ys = [1.0 if r['ok'] else 0.0 for r in rs]
        m0 = sum(x[0] for x in xs) / len(xs)
        m1 = sum(x[1] for x in xs) / len(xs)
        s0 = math.sqrt(sum((x[0] - m0) ** 2 for x in xs) / len(xs))
        s1 = math.sqrt(sum((x[1] - m1) ** 2 for x in xs) / len(xs))
        zz = [((x[0] - m0) / s0, (x[1] - m1) / s1) for x in xs]
        ww0, ww1, bb = 0.0, 0.0, 0.0
        for it in range(3000):
            g0 = g1 = gb = 0.0
            for (z0, z1), y in zip(zz, ys):
                p = 1.0 / (1.0 + math.exp(-(ww0 * z0 + ww1 * z1 + bb)))
                err = p - y
                g0 += err * z0
                g1 += err * z1
                gb += err
            n = len(zz)
            ww0 -= lr * g0 / n
            ww1 -= lr * g1 / n
            bb -= lr * gb / n
        return (ww0, ww1, bb, m0, m1, s0, s1)

    ww0, ww1, bb, m0, m1, s0, s1 = fit_logistic(tr)
    print(f"trained on {len(train_ns)} curves (N={len(tr)}): "
          f"w(log NU)={ww0:.4f}  w(log nuhat)={ww1:.4f}")
    te_scores = [-(ww0 * (math.log(r['NU']) - m0) / s0
                    + ww1 * (math.log(r['nuhat']) - m1) / s1) for r in te]
    te_pos = [s for s, r in zip(te_scores, te) if r['ok']]
    te_neg = [s for s, r in zip(te_scores, te) if not r['ok']]
    print(f"held-out on {len(test_ns)} unseen curves (N={len(te)}): "
          f"AUC(joint) = {auc(te_pos, te_neg):.4f}")
    te_nu_pos = [r['NU'] for r in te if r['ok']]
    te_nu_neg = [r['NU'] for r in te if not r['ok']]
    print(f"  vs held-out AUC(-NU alone)    = {auc(te_nu_pos, te_nu_neg):.4f}")
    te_nh_pos = [r['nuhat'] for r in te if r['ok']]
    te_nh_neg = [r['nuhat'] for r in te if not r['ok']]
    print(f"  vs held-out AUC(-nuhat alone) = {auc(te_nh_pos, te_nh_neg):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
