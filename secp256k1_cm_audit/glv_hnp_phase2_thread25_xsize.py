"""
GLV-HNP Phase 2, Thread 25b: is the product criterion  P = NU * nuhat < 1
SIZE-FREE?

Thread 25 X5 found, on 1000 instances at 17 bits:

    free logistic fit    NU * nuhat^0.9363 = 0.9801
    constrained fit      NU * nuhat        = 0.9579   (LR vs free 1.47, 1 df)
    unfitted cut P < 1   acc 0.8510,  AUC(-P) 0.9350
    per-eff-stratum AUC(-P) = 0.89 .. 0.94  (vs NU alone 0.42 .. 0.85)

Both factors are dimensionless, so the threshold ought to transport.  This
script re-measures the SAME unfitted cut at 12, 14, 17 and 20 bits with the
lattice dimension held at dim = 2m = 24, refits the threshold independently
per size, and reports how far it moves.

Falsifier: if the per-size refitted threshold moves by more than ~2x across
an 8-bit range, or if the fixed cut P<1 loses accuracy monotonically with n,
the product is a 17-bit coincidence and not a criterion.

Run:  python3 glv_hnp_phase2_thread25_xsize.py [--rebuild]
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc
from glv_hnp_phase2_thread25 import (
    HERE, EFFS, load_table, logistic_fit, auc_score, acc_at_half,
)

M = 12                       # dim 2m = 24 at every size
SIZES = (12, 14, 17, 20)
SEEDS_X = [42, 1234, 9999, 555, 31337]


def stats(rows, label):
    y = [1 if r['ok'] else 0 for r in rows]
    P = [r['NU'] * r['nuhat'] for r in rows]
    pos = [p for p, t in zip(P, y) if t]
    neg = [p for p, t in zip(P, y) if not t]
    if not pos or not neg:
        print(f"{label:>10} {len(rows):>5} {sum(y):>5}   (degenerate)")
        return None
    a = auc(pos, neg)
    tp = sum(1 for p, t in zip(P, y) if p < 1 and t)
    fp = sum(1 for p, t in zip(P, y) if p < 1 and not t)
    fn = sum(1 for p, t in zip(P, y) if p >= 1 and t)
    tn = sum(1 for p, t in zip(P, y) if p >= 1 and not t)
    acc = (tp + tn) / len(y)
    beta, _ = logistic_fit([[math.log(p)] for p in P], y)
    thr = math.exp(-beta[0] / beta[1]) if beta[1] else float('nan')
    aNU = auc([r['NU'] for r in rows if r['ok']],
              [r['NU'] for r in rows if not r['ok']])
    anh = auc([r['nuhat'] for r in rows if r['ok']],
              [r['nuhat'] for r in rows if not r['ok']])
    print(f"{label:>10} {len(rows):>5} {sum(y):>5} | {a:>7.4f} {aNU:>7.4f} "
          f"{anh:>8.4f} | {acc:>6.4f} {tp:>5} {fp:>4} {fn:>4} {tn:>5} | "
          f"{thr:>7.4f} | {min(neg):>7.4f} {max(pos):>7.4f}")
    return {'bits': label, 'auc': a, 'acc': acc, 'thr': thr,
            'suff': min(neg), 'nec': max(pos), 'N': len(rows)}


if __name__ == "__main__":
    print("=" * 110)
    print("Thread 25b — is  P = NU * nuhat < 1  size-free?  dim 2m = 24 at "
          "every size.")
    print("=" * 110)
    rebuild = "--rebuild" in sys.argv

    tables = {}
    for bits in SIZES:
        cache = os.path.join(HERE, f"glv_hnp_phase2_thread25_t{bits}.json")
        print(f"\n--- {bits} bits ---")
        tables[bits] = load_table(rebuild=rebuild, cache=cache, bits=bits,
                                  m=M, effs=EFFS, seeds=SEEDS_X)

    print("\n" + "-" * 110)
    print("PER-SIZE: the SAME unfitted cut P < 1, plus an independent refit.")
    print("-" * 110)
    print(f"{'bits':>10} {'N':>5} {'rec':>5} | {'AUC P':>7} {'AUC NU':>7} "
          f"{'AUC nuh':>8} | {'acc':>6} {'TP':>5} {'FP':>4} {'FN':>4} "
          f"{'TN':>5} | {'refit':>7} | {'suff<':>7} {'nec>':>7}")
    summary = []
    for bits in SIZES:
        s = stats(tables[bits], f"{bits}-bit")
        if s:
            summary.append(s)

    allrows = [r for bits in SIZES for r in tables[bits]]
    s = stats(allrows, "POOLED")

    print("\n'suff<' = largest P that never fails is below this (min P over")
    print("failures); 'nec>' = max P over successes.  Between them P is")
    print("ambiguous.  refit = per-size maximum-likelihood threshold.")

    if len(summary) >= 2:
        ts = [x['thr'] for x in summary]
        print(f"\nrefitted threshold across {len(summary)} sizes: "
              f"min {min(ts):.4f}  max {max(ts):.4f}  spread {max(ts)/min(ts):.3f}x")
        print(f"pre-registered falsifier: spread > 2x => not a criterion.  "
              f"VERDICT: {'SURVIVES' if max(ts)/min(ts) <= 2 else 'FALSIFIED'}")

    print("\n" + "-" * 110)
    print("PER-SIZE x PER-EFF AUC(-P): the criterion must also survive with")
    print("the bias strength held fixed.")
    print("-" * 110)
    print(f"{'bits':>6} | " + " ".join(f"{f'eff={e:.2f}':>14}" for e in EFFS))
    for bits in SIZES:
        cells = []
        for eff in EFFS:
            sub = [r for r in tables[bits] if r['effq'] == eff]
            pos = [r['NU'] * r['nuhat'] for r in sub if r['ok']]
            neg = [r['NU'] * r['nuhat'] for r in sub if not r['ok']]
            if not pos or not neg:
                cells.append(f"{'.':>7} ({len(sub):>3}) ")
            else:
                cells.append(f"{auc(pos, neg):>7.4f} ({len(sub):>3}) ")
        print(f"{bits:>6} | " + " ".join(f"{c:>14}" for c in cells))

    print("\n" + "-" * 110)
    print("CROSS-SIZE TRANSFER: fit on ONE size, score the others.  This is")
    print("the honest out-of-sample test of the threshold.")
    print("-" * 110)
    print(f"{'fit on':>8} {'coefficients':>26} | " +
          " ".join(f"{f'{b}b acc':>9}" for b in SIZES))
    for fb in SIZES:
        tr = tables[fb]
        y = [1 if r['ok'] else 0 for r in tr]
        X = [[math.log(r['NU'] * r['nuhat'])] for r in tr]
        beta, _ = logistic_fit(X, y)
        accs = []
        for tb in SIZES:
            te = tables[tb]
            yt = [1 if r['ok'] else 0 for r in te]
            Xt = [[math.log(r['NU'] * r['nuhat'])] for r in te]
            z = [beta[0] + beta[1] * x[0] for x in Xt]
            pr = [1.0 / (1.0 + math.exp(-max(-30, min(30, v)))) for v in z]
            accs.append(acc_at_half(pr, yt))
        cs = f"[{beta[0]:+.3f} {beta[1]:+.3f}]"
        print(f"{fb:>8} {cs:>26} | " +
              " ".join(f"{a:>9.4f}" for a in accs))

    print("\n" + "=" * 110)
    print("done")
    print("=" * 110)
