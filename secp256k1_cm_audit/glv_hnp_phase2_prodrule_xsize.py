"""
GLV-HNP Phase 2, Thread 25b: does the product rule NU * nu_hat <= 1 hold at
other bit-lengths, or is ~1 an artifact of the 17-bit table it was read off?

X5 of glv_hnp_phase2_nu_conditioned.py fitted a logistic on (log NU,
log nu_hat) over 500 17-bit instances and got

    logit = -0.323  -8.479*log NU  -8.055*log nu_hat

whose two slopes agree to 5%, i.e. the fitted boundary is the product rule
NU * nu_hat = exp(-0.323/8.27) = 0.962.  X7 then tested the product directly:
pooled AUC 0.9395, eff-stratified 0.9162, against 0.7996 for the exact BDD
certificate NU alone.

This script is the out-of-sample check.  It is a VALIDATION, not a discovery
run -- the rule was read off the 17-bit table before these sizes were touched.

  H25b: the product NU*nu_hat separates at AUC >= 0.85 at every bit-length,
        and the accuracy-optimal threshold stays within [0.8, 1.2] as n grows
        by 2^10.
  Falsifier: an AUC that decays with n, or a threshold that drifts
        monotonically out of the bracket, means the product is a 17-bit fit
        and only the size-free ranking survives.

Run: python3 glv_hnp_phase2_prodrule_xsize.py
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
from glv_hnp_phase2_gsprofile_strat import EFFS
from glv_hnp_phase2_nu_conditioned import stratified_auc

SIZES = (12, 14, 17, 20, 22)
M = 12


def table(bits, m=M, per_bin=2, nbins=10):
    curves = search_curves(1 << (bits - 1), 1 << bits, per_bin=per_bin,
                           nbins=nbins)
    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves:
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
                          'eff': k1b * k2b / n, 'effq': eff, 'bits': bits,
                          'lamstar': lam_star(lam, n),
                          'prod': r['NU'] * r['nuhat']})
                rows.append(r)
    return curves, rows


def best_threshold(rows, key='prod'):
    """Accuracy-optimal cut, Youden-optimal cut, and the largest zero-FP cut.

    The accuracy-optimal cut moves with class balance, and the recovery rate
    falls from 362/500 to 156/500 across SIZES, so it cannot be read as a
    drift in the rule itself.  Youden's J = TPR - FPR is balance-free and is
    the column to compare across sizes.
    """
    vals = sorted({r[key] for r in rows})
    npos = sum(1 for r in rows if r['ok'])
    nneg = len(rows) - npos
    best_t, best_acc = float('nan'), -1.0
    j_t, best_j = float('nan'), -2.0
    for v in vals:
        tp = sum(1 for r in rows if r['ok'] and r[key] <= v)
        fp = sum(1 for r in rows if not r['ok'] and r[key] <= v)
        acc = (tp + (nneg - fp)) / len(rows)
        if acc > best_acc:
            best_acc, best_t = acc, v
        j = (tp / npos if npos else 0) - (fp / nneg if nneg else 0)
        if j > best_j:
            best_j, j_t = j, v
    zfp = [r[key] for r in rows if not r['ok']]
    zfp_t = min(zfp) if zfp else float('inf')
    tp0 = sum(1 for r in rows if r['ok'] and r[key] < zfp_t)
    return best_t, best_acc, j_t, best_j, zfp_t, tp0


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25b — cross-size validation of the product rule NU*nu_hat")
    print("=" * 78)

    tables = {}
    for bits in SIZES:
        t0 = time.time()
        curves, rows = table(bits)
        tables[bits] = rows
        rec = sum(1 for r in rows if r['ok'])
        print(f"  {bits:>2} bits: {len(curves):>2} curves, {len(rows):>4} rows, "
              f"{rec:>3} recoveries, dim {rows[0]['k']}, "
              f"n in [{min(r['n'] for r in rows)}, {max(r['n'] for r in rows)}]"
              f"  ({time.time()-t0:.1f}s)")

    print("\n" + "-" * 78)
    print("EXP Y1: AUC by bit-length.  'strat' = pooled over eff cells,")
    print("        i.e. bias strength held fixed inside every comparison.")
    print("-" * 78)
    print(f"{'bits':>5} {'N':>5} {'rec':>9} | {'prod':>7} {'NU':>7} "
          f"{'nu_hat':>7} {'eff':>7} | {'prod strat':>11} {'NU strat':>9} "
          f"{'nuh strat':>10}")
    for bits in SIZES:
        rows = tables[bits]
        cells = [[r for r in rows if r['effq'] == e] for e in EFFS]

        def a(key):
            pos = [r[key] for r in rows if r['ok']]
            neg = [r[key] for r in rows if not r['ok']]
            return auc(pos, neg) if pos and neg else float('nan')
        print(f"{bits:>5} {len(rows):>5} "
              f"{str(sum(1 for r in rows if r['ok']))+'/'+str(len(rows)):>9} | "
              f"{a('prod'):>7.4f} {a('NU'):>7.4f} {a('nuhat'):>7.4f} "
              f"{a('eff'):>7.4f} | {stratified_auc(cells,'prod')[0]:>11.4f} "
              f"{stratified_auc(cells,'NU')[0]:>9.4f} "
              f"{stratified_auc(cells,'nuhat')[0]:>10.4f}")

    print("\n" + "-" * 78)
    print("EXP Y2: does the threshold move?  (17-bit fit put it at 0.962)")
    print("-" * 78)
    print(f"{'bits':>5} {'rec rate':>9} {'acc-opt t':>10} {'acc':>7} "
          f"{'J-opt t':>8} {'J':>6} {'acc@1.0':>8} {'max 0-FP t':>11} "
          f"{'TP there':>9}")
    for bits in SIZES:
        rows = tables[bits]
        t, acc, jt, j, zfp, tp = best_threshold(rows)
        npos = sum(1 for r in rows if r['ok'])
        acc1 = sum(1 for r in rows if (r['prod'] <= 1.0) == r['ok']) / len(rows)
        print(f"{bits:>5} {npos/len(rows):>9.3f} {t:>10.3f} {acc:>7.3f} "
              f"{jt:>8.3f} {j:>6.3f} {acc1:>8.3f} {zfp:>11.3f} "
              f"{str(tp)+'/'+str(npos):>9}")
    print("\n  acc-opt t moves with class balance (rec rate column); J-opt t")
    print("  does not.  Compare the J-opt column across sizes, not acc-opt.")

    print("\n" + "-" * 78)
    print("EXP Y3: control — the same two columns for the sound certificate")
    print("        NU <= 1, which theory says has no false positives ever.")
    print("-" * 78)
    print(f"{'bits':>5} {'TP@NU<=1':>9} {'FP@NU<=1':>9} {'recall':>8}")
    for bits in SIZES:
        rows = tables[bits]
        tp = sum(1 for r in rows if r['ok'] and r['NU'] <= 1.0)
        fp = sum(1 for r in rows if not r['ok'] and r['NU'] <= 1.0)
        pos = sum(1 for r in rows if r['ok'])
        print(f"{bits:>5} {tp:>9} {fp:>9} {tp/pos if pos else float('nan'):>8.3f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
