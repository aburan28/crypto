"""
GLV-HNP Phase 2, Thread 25: does mu carry a *second* coordinate once NU is
fixed, or is its apparent power in W5 entirely mediated by NU?

W5/W6 (glv_hnp_phase2_gsprofile_strat.py, Thread 24b) found NU and
nu_hat*sqrt(eff) are uncorrelated at fixed eff (Spearman -0.28..+0.16) yet
both predict recovery.  NU is the exact BDD nearest-plane certificate
(sound, zero false positives at NU<=1, but a size-degrading separator:
AUC 0.978 at 12 bits -> 0.860 at 17 bits).  The ambiguous band where NU
alone gives no answer is [1.04, 2.20] (17 bits, dim 24).

  H25: within that band, AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.
  Falsifier: AUC collapses to ~0.5 inside the band, i.e. mu's pooled AUC
       (0.94, W5-ish) was entirely riding on the NU-eff correlation and
       adds nothing once NU is already known.

Secondary (W1b): step = log2(||b*_{m+1}||) - log2(||b*_1||) is the GS-profile
gap between the m identical lambda_1(L2) copies (block 1) and the second
block; W1b observed step -> 0 exactly as the recovery wall is crossed.
Test whether -step alone, or (NU, step) jointly, separates better than
(NU, mu).

Reuses the exact same 500-instance collection as Thread 24b (20 17-bit
j=0 GLV curves x 5 eff strata x 5 seeds), float GS (justified by W0/W4 of
glv_hnp_phase2_gsprofile.py: relative error vs exact Fractions ~1e-15 at
dim 24).

Run: python3 glv_hnp_phase2_nu_band.py [--dump-json PATH]
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

M17 = 12
EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
NU_BAND = (1.040, 2.199)  # Thread 24b W4, 17-bit sufficient/necessary bracket


def collect():
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
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
                step = (math.log2(r['prof'][M17]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][M17] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    return rows


def logistic_fit_2d(rows, key1, key2, iters=2000, lr=0.05):
    """Minimal hand-rolled logistic regression on (log key1, log key2) ->
    ok, via gradient descent. Returns (w1, w2, b, final_loss)."""
    xs = []
    ys = []
    for r in rows:
        v1, v2 = r[key1], r[key2]
        if v1 is None or v2 is None:
            continue
        if isinstance(v1, float) and (math.isnan(v1) or v1 <= 0):
            continue
        x1 = math.log(v1)
        x2 = v2  # step is already log-scale / signed
        xs.append((x1, x2))
        ys.append(1.0 if r['ok'] else 0.0)
    # standardize
    m1 = sum(x[0] for x in xs) / len(xs)
    m2 = sum(x[1] for x in xs) / len(xs)
    s1 = math.sqrt(sum((x[0] - m1) ** 2 for x in xs) / len(xs)) or 1.0
    s2 = math.sqrt(sum((x[1] - m2) ** 2 for x in xs) / len(xs)) or 1.0
    Z = [((x[0] - m1) / s1, (x[1] - m2) / s2) for x in xs]
    w1 = w2 = b = 0.0
    n = len(Z)
    for _ in range(iters):
        g1 = g2 = gb = 0.0
        loss = 0.0
        for (z1, z2), y in zip(Z, ys):
            zscore = w1 * z1 + w2 * z2 + b
            p = 1.0 / (1.0 + math.exp(-max(-30, min(30, zscore))))
            err = p - y
            g1 += err * z1
            g2 += err * z2
            gb += err
            loss += -(y * math.log(max(p, 1e-12)) +
                       (1 - y) * math.log(max(1 - p, 1e-12)))
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
        b -= lr * gb / n
    return w1, w2, b, loss / n, (m1, s1, m2, s2)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None)
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — is mu a second coordinate, conditional on NU?")
    print("=" * 78)

    t0 = time.time()
    rows = collect()
    print(f"\n{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        with open(args.dump_json, "w") as f:
            json.dump([{k: v for k, v in r.items() if k != 'nus' and k != 'prof'}
                       for r in rows], f)
        print(f"dumped {len(rows)} rows to {args.dump_json}")

    print("\n" + "-" * 78)
    print(f"H25: within NU band [{NU_BAND[0]}, {NU_BAND[1]}], "
          f"AUC(-mu -> recovery) >= 0.8 ?")
    print("-" * 78)
    band = [r for r in rows if NU_BAND[0] <= r['NU'] <= NU_BAND[1]]
    print(f"band population: {len(band)} / {len(rows)} instances "
          f"({100*len(band)/len(rows):.1f}%)")
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"  recoveries in band: {len(pos)}/{len(band)}")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_step = auc([-r['step'] for r in pos], [-r['step'] for r in neg])
        a_nu_in_band = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        print(f"  AUC(-mu    -> recovery), in-band  = {a_mu:.4f}")
        print(f"  AUC(-nuhat -> recovery), in-band  = {a_nh:.4f}")
        print(f"  AUC(step   -> recovery), in-band  = {a_step:.4f}  "
              f"(step>0 predicted to fail per W1b)")
        print(f"  AUC(-NU    -> recovery), in-band  = {a_nu_in_band:.4f}  "
              f"(sanity: NU should be near-random here by construction)")
        verdict = "SURVIVES" if a_mu >= 0.8 else "FALSIFIED"
        print(f"\n  H25 verdict: {verdict} (threshold 0.8, observed {a_mu:.4f})")
    else:
        print("  degenerate band (all one class) -- cannot compute AUC")

    print("\n" + "-" * 78)
    print("Per-eff-stratum in-band AUC(-mu), to check H25 isn't an artifact "
          "of pooling across eff")
    print("-" * 78)
    print(f"{'eff':>5} {'N in band':>10} {'rec':>7} {'AUC mu':>8} {'AUC step':>9}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not p or not ng:
            print(f"{eff:>5.2f} {len(sub):>10} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {'(degenerate)':>8}")
            continue
        a_mu = auc([r['mu'] for r in p], [r['mu'] for r in ng])
        a_st = auc([-r['step'] for r in p], [-r['step'] for r in ng])
        print(f"{eff:>5.2f} {len(sub):>10} "
              f"{str(len(p))+'/'+str(len(sub)):>7} {a_mu:>8.4f} {a_st:>9.4f}")

    print("\n" + "-" * 78)
    print("W1b secondary: step = log2||b*_{m+1}|| - log2||b*_1|| vs recovery, "
          "ALL instances (not just in-band)")
    print("-" * 78)
    valid = [r for r in rows if not math.isnan(r['step'])]
    pos_a = [r for r in valid if r['ok']]
    neg_a = [r for r in valid if not r['ok']]
    a_step_all = auc([-r['step'] for r in pos_a], [-r['step'] for r in neg_a])
    print(f"  AUC(-step -> recovery), all {len(valid)} instances = {a_step_all:.4f}")
    print(f"  AUC(-NU   -> recovery), all instances (Thread 24b) = "
          f"{auc([r['NU'] for r in pos_a], [r['NU'] for r in neg_a]):.4f}")
    print(f"  Spearman(step, NU) = {spearman([r['step'] for r in valid], [r['NU'] for r in valid]):.4f}")

    print("\n" + "-" * 78)
    print("Joint logistic fit on (log NU, step) -> recovery")
    print("-" * 78)
    w1, w2, b, loss, (m1, s1, m2, s2) = logistic_fit_2d(valid, 'NU', 'step')
    print(f"  standardized weights: w_logNU={w1:.3f}  w_step={w2:.3f}  b={b:.3f}")
    print(f"  final mean cross-entropy = {loss:.4f}")
    print(f"  (logNU mean={m1:.3f} sd={s1:.3f} | step mean={m2:.3f} sd={s2:.3f})")
    # AUC of the joint score as a sanity check against NU alone / mu alone
    scores = []
    for r in valid:
        z1 = (math.log(r['NU']) - m1) / s1
        z2 = (r['step'] - m2) / s2
        scores.append(w1 * z1 + w2 * z2 + b)
    pos_s = [s for s, r in zip(scores, valid) if r['ok']]
    neg_s = [s for s, r in zip(scores, valid) if not r['ok']]
    print(f"  AUC(joint logistic score -> recovery) = "
          f"{auc([-s for s in pos_s], [-s for s in neg_s]):.4f}  "
          f"(oriented so higher score = more likely recovery)")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
