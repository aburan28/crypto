"""
GLV-HNP Phase 2, Thread 25b: is the product rule NU*nu_hat size-stable?

Thread 25a (glv_hnp_phase2_nu_conditional.py) found that NU and nu_hat are
two genuine coordinates -- each separates recovery inside bands of the other
-- and that the P=1/2 boundary of the logistic fit on (log NU, log nu_hat)
collapses to a near-perfect product:

    PR = NU * nu_hat^(1/1.052) ~ NU * nu_hat  <  0.96   =>  recover

at 17 bits: in-sample AUC 0.9397, leave-one-curve-out AUC 0.9286, accuracy
0.866, against 0.7564 for NU alone and 0.6648 for nu_hat alone.

That was fitted and tested at ONE size.  The whole reason to want a second
coordinate is that NU alone degrades with size -- Thread 24 W4 measured
AUC(-NU) falling 0.978 -> 0.860 from 12 to 17 bits, and the ambiguous bracket
widening from [1.188, 1.869] to [1.040, 2.199].  So the question that decides
whether the product rule is a result or a curve-fit is:

  H25b: AUC(-NU*nu_hat) stays >= 0.90 at 12, 14, 17 and 20 bits, and the
        fitted P=1/2 threshold stays inside a factor of 2 across those sizes.
  Falsifier: if the threshold drifts monotonically with n, or the AUC decays
        the way AUC(-NU) does, then the product rule is a 17-bit fit and the
        pair (NU, nu_hat) needs a size-dependent normalisation.

The lattice dimension is held FIXED at 2m = 24 across all four sizes, so any
drift is a size effect and not a dimension effect.

Run: python3 glv_hnp_phase2_product_rule.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc

M = 12                       # dim = 2m = 24 at every size
EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
SIZES = (12, 14, 17, 20)
EXPO = 1.0 / 1.052           # from the 17-bit logistic boundary


def collect(bits, ncurves=12):
    curves = search_curves(1 << (bits - 1), 1 << bits, per_bin=2, nbins=10)
    curves = curves[:ncurves]
    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), M, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), M, d_trial, k1b, seed)
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'bits': bits})
                r['PR'] = r['NU'] * (r['nuhat'] ** EXPO)
                rows.append(r)
    return curves, rows


def split_auc(rows, key):
    pos = [r[key] for r in rows if r['ok']]
    neg = [r[key] for r in rows if not r['ok']]
    return auc(pos, neg)


def best_threshold(rows, key):
    """Threshold maximising balanced accuracy; returns (thr, bal_acc)."""
    vals = sorted({r[key] for r in rows})
    P = sum(1 for r in rows if r['ok'])
    N = len(rows) - P
    if P == 0 or N == 0:
        return float('nan'), float('nan')
    best = (float('nan'), -1.0)
    for i, v in enumerate(vals):
        thr = v if i + 1 == len(vals) else 0.5 * (v + vals[i + 1])
        tp = sum(1 for r in rows if r[key] < thr and r['ok'])
        fp = sum(1 for r in rows if r[key] < thr and not r['ok'])
        ba = 0.5 * (tp / P + (N - fp) / N)
        if ba > best[1]:
            best = (thr, ba)
    return best


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25b — is the product rule NU*nu_hat size-stable? (H25b)")
    print("=" * 78)
    print(f"\ndim = 2m = {2*M} held FIXED at every size; "
          f"eff strata {EFFS}; seeds {SEEDS}")
    print(f"product score PR = NU * nu_hat^{EXPO:.4f}  "
          f"(exponent fitted at 17 bits)")

    data = {}
    for bits in SIZES:
        t0 = time.time()
        curves, rows = collect(bits)
        data[bits] = rows
        w = sum(1 for r in rows if r['ok'])
        print(f"  {bits:>2} bits: {len(curves):>2} curves, {len(rows):>3} "
              f"instances, {w:>3} recovered, {time.time()-t0:>5.1f}s")

    print("\n" + "-" * 78)
    print("EXP B1: AUC by size.  NU alone is the column that is known to")
    print("        decay; the question is whether the product follows it.")
    print("-" * 78)
    print(f"{'bits':>5} {'N':>5} {'rec':>9} {'AUC PR':>8} {'AUC NU*nuh':>11} "
          f"{'AUC NU':>8} {'AUC nuhat':>10} {'AUC mu':>8}")
    for bits in SIZES:
        rows = data[bits]
        w = sum(1 for r in rows if r['ok'])
        plain = [r['NU'] * r['nuhat'] for r in rows]
        a_pl = auc([v for v, r in zip(plain, rows) if r['ok']],
                   [v for v, r in zip(plain, rows) if not r['ok']])
        print(f"{bits:>5} {len(rows):>5} "
              f"{str(w) + '/' + str(len(rows)):>9} "
              f"{split_auc(rows, 'PR'):>8.4f} {a_pl:>11.4f} "
              f"{split_auc(rows, 'NU'):>8.4f} "
              f"{split_auc(rows, 'nuhat'):>10.4f} "
              f"{split_auc(rows, 'mu'):>8.4f}")

    print("\n" + "-" * 78)
    print("EXP B2: does the P=1/2 threshold move?  Best-balanced-accuracy cut")
    print("        on each score, fitted independently at each size.")
    print("-" * 78)
    print(f"{'bits':>5} | {'PR thr':>8} {'bal acc':>8} | {'NU thr':>8} "
          f"{'bal acc':>8} | {'nuhat thr':>10} {'bal acc':>8}")
    thrs = {}
    for bits in SIZES:
        rows = data[bits]
        tp, bp = best_threshold(rows, 'PR')
        tn_, bn = best_threshold(rows, 'NU')
        th, bh = best_threshold(rows, 'nuhat')
        thrs[bits] = tp
        print(f"{bits:>5} | {tp:>8.4f} {bp:>8.4f} | {tn_:>8.4f} {bn:>8.4f} "
              f"| {th:>10.4f} {bh:>8.4f}")
    ok = [t for t in thrs.values() if t == t]
    if ok:
        print(f"\nPR threshold spread across sizes: "
              f"{min(ok):.4f} .. {max(ok):.4f}  "
              f"= {max(ok)/min(ok):.2f}x   (H25b allows 2x)")
        print(f"17-bit logistic P=1/2 boundary was 0.9607")

    print("\n" + "-" * 78)
    print("EXP B3: transfer.  Apply the 17-bit rule PR < 0.9607 UNCHANGED at")
    print("        every other size -- no refitting at all.")
    print("-" * 78)
    THR = 0.9607
    print(f"{'bits':>5} {'TP':>5} {'FP':>5} {'FN':>5} {'TN':>5} {'acc':>8} "
          f"{'bal acc':>8} {'FP rate':>8}")
    for bits in SIZES:
        rows = data[bits]
        tp = sum(1 for r in rows if r['PR'] < THR and r['ok'])
        fp = sum(1 for r in rows if r['PR'] < THR and not r['ok'])
        fn = sum(1 for r in rows if r['PR'] >= THR and r['ok'])
        tn = len(rows) - tp - fp - fn
        P, N = tp + fn, fp + tn
        ba = 0.5 * (tp / P + tn / N) if P and N else float('nan')
        print(f"{bits:>5} {tp:>5} {fp:>5} {fn:>5} {tn:>5} "
              f"{(tp+tn)/len(rows):>8.4f} {ba:>8.4f} "
              f"{fp/N if N else float('nan'):>8.4f}")

    print("\n" + "-" * 78)
    print("EXP B4: the sound-certificate question.  NU <= 1 is a THEOREM for")
    print("        nearest-plane.  Is there a PR cut with zero false")
    print("        positives, and how much more does it certify than NU<=1?")
    print("-" * 78)
    print(f"{'bits':>5} {'NU<=1 TP':>9} {'NU<=1 FP':>9} | "
          f"{'max safe PR':>12} {'TP there':>9} {'gain':>7}")
    for bits in SIZES:
        rows = data[bits]
        tp1 = sum(1 for r in rows if r['NU'] <= 1.0 and r['ok'])
        fp1 = sum(1 for r in rows if r['NU'] <= 1.0 and not r['ok'])
        fails = [r['PR'] for r in rows if not r['ok']]
        safe = min(fails) if fails else float('inf')
        tp2 = sum(1 for r in rows if r['PR'] < safe and r['ok'])
        print(f"{bits:>5} {tp1:>9} {fp1:>9} | {safe:>12.4f} {tp2:>9} "
              f"{tp2/tp1 if tp1 else float('nan'):>6.2f}x")
    print("\n(max safe PR = largest cut with 0 FP in this sample; it is an")
    print(" EMPIRICAL bound, not a theorem, and will shrink with more data.)")

    print("\n" + "-" * 78)
    print("EXP B5: per-size within-eff-stratum AUC, the Thread 24 W5 control")
    print("        repeated at every size.")
    print("-" * 78)
    print(f"{'bits':>5} {'eff':>5} {'N':>4} {'rec':>7} {'AUC PR':>8} "
          f"{'AUC NU':>8} {'AUC nuhat':>10}")
    for bits in SIZES:
        for eff in EFFS:
            sub = [r for r in data[bits] if r['effq'] == eff]
            if not sub:
                continue
            w = sum(1 for r in sub if r['ok'])
            if w == 0 or w == len(sub):
                print(f"{bits:>5} {eff:>5.2f} {len(sub):>4} "
                      f"{str(w) + '/' + str(len(sub)):>7} {'(degenerate)':>8}")
                continue
            print(f"{bits:>5} {eff:>5.2f} {len(sub):>4} "
                  f"{str(w) + '/' + str(len(sub)):>7} "
                  f"{split_auc(sub, 'PR'):>8.4f} {split_auc(sub, 'NU'):>8.4f} "
                  f"{split_auc(sub, 'nuhat'):>10.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
