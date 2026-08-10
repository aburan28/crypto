"""
GLV-HNP Phase 2, Thread 25: find the second mechanism by conditioning on NU.

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry:

  W5/W6 established recovery = f(NU, X) with X ~ mu-driven (mu = lambda_1(L2))
  and X independent of NU (Spearman(pred, NU) ~ 0 or negative in every
  eff-stratum). NU is a sound BDD certificate (0 FP at NU<=1) but a
  size-degrading separator (AUC 0.978 -> 0.860, 12->17 bits); mu/nu_hat is a
  better cross-curve separator (AUC 0.75-0.93 at fixed eff) but is not a
  certificate for anything on its own.

  H25: within the ambiguous NU band [1.04, 2.20] (17-bit bracket, where the
       nearest-plane certificate gives no verdict either way), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.
  Falsifier: if mu's apparent power inside the band collapses towards 0.5,
       mu's W5 signal was entirely mediated by NU (via the shared eff trend)
       and the "two mechanisms" reading of W5/W6 is a stratification
       artifact -- retire the closed form.

Secondary (from W1b): the GS profile of L0 is m exact copies of
lambda_1(L2) followed by a block that moves with K1; the step from block 1
to block 2 vanishes right as the K1 wall is crossed. Quantify
  step_i = log2(||b*_{m+i}||) - log2(||b*_i||)   averaged over i=1..m
and test it as a third predictor, alongside NU and mu, of Kannan-LLL
recovery.

Data: identical generation to glv_hnp_phase2_gsprofile_strat.py (17-bit
curves, M=12, 5 eff strata x 20 curves x 5 seeds = 500 instances target),
but with exact=True so `prof` (needed for `step`) is on the Fraction path
used by W0/W4 of the parent script, not re-derived. Run once, dump to JSON
so re-analysis needs no new lattice work (per the Thread 24 cost note).

Run: python3 glv_hnp_phase2_thread25.py [--dump-json FILE] [--from-json FILE]
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
# 17-bit bracket from the 2026-08-07 #2 (Thread 24) log entry, W4.
NU_BAND_LO, NU_BAND_HI = 1.040, 2.199


def collect():
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")

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
                step = sum(math.log2(r['prof'][m + i]) - math.log2(r['prof'][i])
                           for i in range(m) if r['prof'][i] > 0
                           and r['prof'][m + i] > 0) / m
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                # 'prof' holds python floats but is only needed to derive
                # 'step'; drop it so the JSON dump stays small.
                r.pop('prof')
                r.pop('nus')
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k'] if rows else '?'}) "
          f"in {time.time()-t0:.1f}s")
    return rows


def load(path):
    with open(path) as f:
        return json.load(f)


def save(rows, path):
    with open(path, 'w') as f:
        json.dump(rows, f)


def logistic_fit_1d(x, y, iters=8000, lr=0.5):
    n = len(y)
    m, s = sum(x) / n, (sum((v - sum(x) / n) ** 2 for v in x) / n) ** 0.5
    z = [(v - m) / s for v in x]
    w = b = 0.0
    for _ in range(iters):
        gw = gb = 0.0
        for a, t in zip(z, y):
            p = 1.0 / (1.0 + math.exp(-(w * a + b)))
            err = p - t
            gw += err * a
            gb += err
        w -= lr * gw / n
        b -= lr * gb / n
    correct = sum((1.0 / (1.0 + math.exp(-(w * a + b))) >= 0.5) == (t == 1)
                  for a, t in zip(z, y))
    return {'w': w, 'b': b, 'train_acc': correct / n}


def logistic_fit_2d(x1, x2, y, iters=20000, lr=0.5):
    """Plain gradient-descent logistic regression, y in {0,1}, standardised
    features. Returns (w1, w2, b, decision_boundary_fn_in_original_units)."""
    n = len(y)
    m1, s1 = sum(x1) / n, (sum((v - sum(x1) / n) ** 2 for v in x1) / n) ** 0.5
    m2, s2 = sum(x2) / n, (sum((v - sum(x2) / n) ** 2 for v in x2) / n) ** 0.5
    z1 = [(v - m1) / s1 for v in x1]
    z2 = [(v - m2) / s2 for v in x2]
    w1 = w2 = b = 0.0
    for _ in range(iters):
        g1 = g2 = gb = 0.0
        for a, c, t in zip(z1, z2, y):
            p = 1.0 / (1.0 + math.exp(-(w1 * a + w2 * c + b)))
            err = p - t
            g1 += err * a
            g2 += err * c
            gb += err
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
        b -= lr * gb / n
    # accuracy at p=0.5
    correct = 0
    for a, c, t in zip(z1, z2, y):
        p = 1.0 / (1.0 + math.exp(-(w1 * a + w2 * c + b)))
        correct += (p >= 0.5) == (t == 1)
    return {'w1': w1, 'w2': w2, 'b': b, 'm1': m1, 's1': s1, 'm2': m2, 's2': s2,
            'train_acc': correct / n}


def report(rows):
    print(f"\nloaded {len(rows)} instances")
    print("\n" + "-" * 78)
    print(f"EXP H25: AUC(-mu -> recovery) inside vs outside the NU band "
          f"[{NU_BAND_LO}, {NU_BAND_HI}]")
    print("-" * 78)
    band = [r for r in rows if NU_BAND_LO <= r['NU'] <= NU_BAND_HI]
    out = [r for r in rows if not (NU_BAND_LO <= r['NU'] <= NU_BAND_HI)]
    print(f"in-band: {len(band)}/{len(rows)} instances "
          f"({sum(1 for r in band if r['ok'])}/{len(band)} recover)")
    print(f"out-of-band: {len(out)}/{len(rows)} instances "
          f"({sum(1 for r in out if r['ok'])}/{len(out)} recover, "
          f"should be near-deterministic by construction of the bracket)")

    def auc_mu(sub, label):
        pos = [r['mu'] for r in sub if r['ok']]
        neg = [r['mu'] for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"  {label}: degenerate (pos={len(pos)} neg={len(neg)})")
            return None
        a = auc(pos, neg)
        print(f"  {label}: AUC(-mu -> recovery) = {a:.4f}  "
              f"(N={len(sub)}, {len(pos)} pos / {len(neg)} neg)")
        return a

    a_in = auc_mu(band, "IN-BAND ")
    a_out = auc_mu(out, "OUT-BAND")
    a_all = auc_mu(rows, "ALL     ")

    print(f"\nH25 target: AUC(-mu) in-band >= 0.80.  "
          f"Observed: {a_in if a_in is not None else float('nan'):.4f}  "
          f"=> {'HOLDS' if (a_in is not None and a_in >= 0.80) else 'FALSIFIED'}")

    # per-eff-stratum in-band breakdown (mu and eff both vary with the
    # stratum design, so also check within-band-within-stratum)
    print("\nin-band AUC(-mu), stratified by effq (checks eff isn't leaking "
          "back in through the band selection):")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        pos = [r['mu'] for r in sub if r['ok']]
        neg = [r['mu'] for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"  eff={eff:.2f}  N={len(sub):3d}  degenerate "
                  f"(pos={len(pos)} neg={len(neg)})")
            continue
        print(f"  eff={eff:.2f}  N={len(sub):3d}  AUC(-mu) = "
              f"{auc(pos, neg):.4f}")

    print("\n" + "-" * 78)
    print("EXP secondary: step = mean_i[log2||b*_{m+i}|| - log2||b*_i||] "
          "as a third predictor")
    print("-" * 78)
    # auc(pos, neg) = P(pos < neg); step is expected LARGER for success
    # (W1b: step -> 0 right at the wall), the OPPOSITE sign convention from
    # mu/NU/nu_hat, so the "smaller score -> success" AUC here is 1-auc(pos,neg).
    pos = [r['step'] for r in rows if r['ok']]
    neg = [r['step'] for r in rows if not r['ok']]
    step_auc_pooled = 1 - auc(pos, neg)
    print(f"AUC(step -> recovery), LARGER step -> success = "
          f"{step_auc_pooled:.4f}  (pooled, all 500)")
    print(f"step | success: mean {sum(pos)/len(pos):.4f}  min {min(pos):.4f}  "
          f"max {max(pos):.4f}")
    print(f"step | failure: mean {sum(neg)/len(neg):.4f}  min {min(neg):.4f}  "
          f"max {max(neg):.4f}")
    print("in-band only (where NU gives no verdict):")
    posb = [r['step'] for r in band if r['ok']]
    negb = [r['step'] for r in band if not r['ok']]
    if posb and negb:
        step_auc_band = 1 - auc(posb, negb)
        print(f"  AUC(step -> recovery) in-band = {step_auc_band:.4f}  "
              f"(N={len(band)}, no eff-stratification needed)")
    sp = spearman([r['step'] for r in rows], [r['NU'] for r in rows])
    sp_mu = spearman([r['step'] for r in rows], [math.log(r['mu']) for r in rows])
    print(f"Spearman(step, NU)        = {sp:.4f}")
    print(f"Spearman(step, log mu)    = {sp_mu:.4f}")

    print("\n" + "-" * 78)
    print("EXP: logistic fit on (log NU, log mu) -> recovery, full pooled set")
    print("-" * 78)
    x1 = [math.log(r['NU']) for r in rows]
    x2 = [math.log(r['mu']) for r in rows]
    y = [1 if r['ok'] else 0 for r in rows]
    fit = logistic_fit_2d(x1, x2, y)
    print(f"standardised weights: w(logNU)={fit['w1']:.4f}  "
          f"w(logmu)={fit['w2']:.4f}  b={fit['b']:.4f}  "
          f"train accuracy={fit['train_acc']:.4f}")
    print(f"(feature standardisation: logNU ~ N({fit['m1']:.3f},{fit['s1']:.3f}), "
          f"logmu ~ N({fit['m2']:.3f},{fit['s2']:.3f}))")
    only_nu = logistic_fit_1d(x1, y)
    only_mu = logistic_fit_1d(x2, y)
    print(f"logNU-only  train accuracy = {only_nu['train_acc']:.4f}")
    print(f"logmu-only  train accuracy = {only_mu['train_acc']:.4f}")
    print(f"joint (logNU,logmu) train accuracy = {fit['train_acc']:.4f}")

    x3 = [r['step'] for r in rows]
    fit_step = logistic_fit_2d(x1, x3, y)
    only_step = logistic_fit_1d(x3, y)
    print(f"\nstep-only   train accuracy = {only_step['train_acc']:.4f}")
    print(f"joint (logNU,step) train accuracy = {fit_step['train_acc']:.4f}  "
          f"weights: w(logNU)={fit_step['w1']:.4f}  w(step)={fit_step['w2']:.4f}")

    print("\n" + "-" * 78)
    print("EXP: does the pooled in-band H25 test reintroduce the eff-confound?")
    print("-" * 78)
    print("in-band-within-stratum, EXCLUDING the degenerate eff=0.05 stratum "
          "(24 obs, only 3 negatives):")
    band = [r for r in rows if NU_BAND_LO <= r['NU'] <= NU_BAND_HI]
    clean = [r for r in band if r['effq'] != 0.05]
    pos = [r['mu'] for r in clean if r['ok']]
    neg = [r['mu'] for r in clean if not r['ok']]
    print(f"  N={len(clean)}  AUC(-mu -> recovery) = {auc(pos, neg):.4f}")
    print(f"  vs pooled-including-eff=0.05: see H25 block above (0.69)")
    print(f"  vs raw per-stratum range: 0.86-0.93 (eff 0.10-0.25)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None)
    ap.add_argument("--from-json", default=None)
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — conditioning on NU: is mu a genuine second coordinate?")
    print("=" * 78)

    if args.from_json:
        rows = load(args.from_json)
    else:
        rows = collect()
        if args.dump_json:
            save(rows, args.dump_json)
            print(f"dumped {len(rows)} rows to {args.dump_json}")

    report(rows)

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
