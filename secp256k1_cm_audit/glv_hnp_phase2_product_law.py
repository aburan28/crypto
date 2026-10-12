"""
GLV-HNP Phase 2, Thread 25b: is  PI = NU * nu_hat ~ 1  the viability curve?

POST-HOC.  This hypothesis was NOT pre-registered.  It fell out of EXP X4 of
`glv_hnp_phase2_nu_conditioned.py`, where a logistic fit on (log NU,
log nu_hat) over 1000 17-bit instances returned

    logit p = -0.2332 - 8.2056*log(NU) - 7.5446*log(nu_hat)
    => decision boundary   NU * nu_hat^0.9194 = 0.9720

i.e. an exponent within 8% of 1 and a constant within 3% of 1.  A fit on one
size is worth nothing on its own — two free parameters landing near (1, 1) is
exactly what an overfit looks like.  This script does the only test that can
distinguish a law from a coincidence: refit independently at three sizes and
score the FIXED, unfitted rule out of sample.

  H25b: the exponent alpha = -w2/w1 and the constant c = exp(-w0/w1) of the
        (log NU, log nu_hat) logistic are size-stable, and the parameter-free
        rule PI = NU*nu_hat < 1 separates recovery from failure at every size.
  Falsifier: if alpha or c drifts materially with bit-length, or if a 12-bit
        fit applied out of sample to 20 bits loses accuracy relative to the
        in-sample fit, PI is a per-size fit and not a law.

Note on what PI can and cannot be.  NU <= 1 is a SUFFICIENT condition for
Babai nearest-plane (Thread 23b/24: 96 TP / 0 FP over 410 instances) and that
is a theorem.  PI is a claim about KANNAN-LLL recovery, which is strictly
stronger than nearest-plane, so PI < 1 is not implied by, and does not imply,
the nearest-plane bound.  Any threshold behaviour here is empirical.

Experiments:
  Y1  independent (log NU, log nu_hat) logistic fits at 12 / 17 / 20 bits;
      compare alpha and c.
  Y2  the parameter-free rule PI < 1: AUC, accuracy, TP/FP, and the
      sufficient/necessary bracket at each size.
  Y3  out-of-sample transfer: fit at one size, score at the others.
  Y4  double conditioning: AUC(-PI) inside NU quintiles and inside eff
      strata.  PI must beat BOTH parents in BOTH conditionings, or it is
      just whichever parent dominates that slice.

Run: python3 glv_hnp_phase2_product_law.py [--dump-json FILE]
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
from glv_hnp_phase2_nu_conditioned import (
    ALL_SEEDS, EFFS, auc_boot, quantile_bins, bin_of, logistic, score_of,
)

# (bit-length, m).  m = 12 matches Thread 24's 17-bit runs; 12-bit uses m = 8
# as Thread 23b/24 did, so the dimension is 16 there and 24 at 17/20 bits.
SIZES = [(12, 8), (17, 12), (20, 12)]


def table_for(bits, m, seeds, verbose=True):
    t0 = time.time()
    curves = search_curves(1 << (bits - 1), 1 << bits, per_bin=2, nbins=10)
    if verbose:
        print(f"  {bits} bits: {len(curves)} curves in {time.time()-t0:.1f}s", end="")
    t0 = time.time()
    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in seeds:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), m, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), m, d_trial, k1b, seed)
                r.pop('prof', None)
                r.pop('nus', None)
                r.update({'bits': bits, 'm': m, 'n': n, 'K1': k1b,
                          'seed': seed, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n),
                          'PI': r['NU'] * r['nuhat']})
                rows.append(r)
    if verbose:
        print(f",  {len(rows)} instances (dim {rows[0]['k']}) "
              f"in {time.time()-t0:.1f}s")
    return rows


def feat(r):
    return [1.0, math.log(r['NU']), math.log(r['nuhat'])]


def fit_report(rows):
    y = [1.0 if r['ok'] else 0.0 for r in rows]
    X = [feat(r) for r in rows]
    w = logistic(X, y)
    s = [score_of(w, xi) for xi in X]
    a = auc([-s[i] for i in range(len(rows)) if rows[i]['ok']],
            [-s[i] for i in range(len(rows)) if not rows[i]['ok']])
    acc = sum(1 for i in range(len(rows)) if (s[i] > 0) == rows[i]['ok']) / len(rows)
    alpha = w[2] / w[1] if abs(w[1]) > 1e-12 else float('nan')
    c = math.exp(-w[0] / w[1]) if abs(w[1]) > 1e-12 else float('nan')
    return w, alpha, c, a, acc


def bracket(rows, key):
    pos = [r[key] for r in rows if r['ok']]
    neg = [r[key] for r in rows if not r['ok']]
    return (min(neg) if neg else float('nan'),
            max(pos) if pos else float('nan'))


if __name__ == "__main__":
    dump = None
    if "--dump-json" in sys.argv:
        i = sys.argv.index("--dump-json")
        dump = sys.argv[i + 1] if i + 1 < len(sys.argv) else "product_law.json"

    print("=" * 78)
    print("Thread 25b — is PI = NU * nu_hat ~ 1 the viability curve?  (H25b)")
    print("POST-HOC hypothesis from X4; this is its first out-of-sample test.")
    print("=" * 78)

    print("\nbuilding tables:")
    tables = {}
    for bits, m in SIZES:
        tables[bits] = table_for(bits, m, ALL_SEEDS)
    allrows = [r for bits, _ in SIZES for r in tables[bits]]
    if dump:
        with open(dump, "w") as fh:
            json.dump(allrows, fh)
        print(f"tables dumped to {dump} ({os.path.getsize(dump)} bytes)")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP Y1: independent logistic fits per size")
    print("        logit p = w0 + w1*log(NU) + w2*log(nu_hat)")
    print("        boundary:  NU * nu_hat^alpha = c ,  alpha = w2/w1 , "
          "c = exp(-w0/w1)")
    print("-" * 78)
    print(f"{'bits':>5} {'dim':>4} {'N':>5} {'rec':>9} {'w1':>9} {'w2':>9} "
          f"{'alpha':>8} {'c':>8} {'AUC':>8} {'acc':>7}")
    fits = {}
    for bits, m in SIZES:
        rows = tables[bits]
        w, alpha, c, a, acc = fit_report(rows)
        fits[bits] = w
        rec = sum(1 for r in rows if r['ok'])
        print(f"{bits:>5} {rows[0]['k']:>4} {len(rows):>5} "
              f"{str(rec)+'/'+str(len(rows)):>9} {w[1]:>9.4f} {w[2]:>9.4f} "
              f"{alpha:>8.4f} {c:>8.4f} {a:>8.4f} {acc:>7.4f}")
    print("\nH25b predicts alpha ~ 1 and c ~ 1 at every size.")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP Y2: the PARAMETER-FREE rule  PI = NU * nu_hat  vs threshold 1")
    print("-" * 78)
    print(f"{'bits':>5} {'N':>5} {'rec':>9} | {'AUC(-PI)':>9} {'AUC(-NU)':>9} "
          f"{'AUC(-nuh)':>10} | {'acc@PI<1':>9} {'base':>6} | "
          f"{'TP':>4} {'FP':>4} {'FN':>4} {'TN':>4}")
    for bits, m in SIZES:
        rows = tables[bits]
        pos = [r for r in rows if r['ok']]
        neg = [r for r in rows if not r['ok']]
        aPI = auc([r['PI'] for r in pos], [r['PI'] for r in neg])
        aNU = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        anh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        tp = sum(1 for r in rows if r['PI'] < 1 and r['ok'])
        fp = sum(1 for r in rows if r['PI'] < 1 and not r['ok'])
        fn = sum(1 for r in rows if r['PI'] >= 1 and r['ok'])
        tn = sum(1 for r in rows if r['PI'] >= 1 and not r['ok'])
        acc = (tp + tn) / len(rows)
        base = max(len(pos), len(neg)) / len(rows)
        print(f"{bits:>5} {len(rows):>5} "
              f"{str(len(pos))+'/'+str(len(rows)):>9} | {aPI:>9.4f} "
              f"{aNU:>9.4f} {anh:>10.4f} | {acc:>9.4f} {base:>6.4f} | "
              f"{tp:>4} {fp:>4} {fn:>4} {tn:>4}")

    print(f"\n{'bits':>5} {'PI | recovered':>34} {'PI | failed':>34}")
    for bits, m in SIZES:
        rows = tables[bits]
        pos = sorted(r['PI'] for r in rows if r['ok'])
        neg = sorted(r['PI'] for r in rows if not r['ok'])
        f = lambda v: (f"min {v[0]:.3f} med {v[len(v)//2]:.3f} max {v[-1]:.3f}"
                       if v else "  (none)")
        print(f"{bits:>5} {f(pos):>34} {f(neg):>34}")
    print(f"\n{'bits':>5} {'sufficient PI <':>17} {'necessary PI >':>16} "
          f"{'ambiguous width':>16}")
    for bits, m in SIZES:
        suf, nec = bracket(tables[bits], 'PI')
        print(f"{bits:>5} {suf:>17.4f} {nec:>16.4f} {nec/suf:>15.2f}x")
    print("compare the NU bracket (Thread 24 W4): 12 bits [1.188, 1.869] "
          "= 1.57x , 17 bits [1.040, 2.199] = 2.11x")
    print(f"\n{'bits':>5} {'sufficient NU <':>17} {'necessary NU >':>16} "
          f"{'ambiguous width':>16}")
    for bits, m in SIZES:
        suf, nec = bracket(tables[bits], 'NU')
        print(f"{bits:>5} {suf:>17.4f} {nec:>16.4f} {nec/suf:>15.2f}x")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP Y3: out-of-sample transfer — fit at one size, score at another")
    print("-" * 78)
    print("accuracy of the fitted logistic (rows = fit size, cols = test size)")
    hdr0 = "fit / test"
    print(f"{hdr0:>10} " + " ".join(f"{str(b)+'b':>9}" for b, _ in SIZES))
    for fb, _ in SIZES:
        w = fits[fb]
        cells = []
        for tb, _ in SIZES:
            rows = tables[tb]
            s = [score_of(w, feat(r)) for r in rows]
            acc = sum(1 for i in range(len(rows))
                      if (s[i] > 0) == rows[i]['ok']) / len(rows)
            cells.append(f"{acc:>9.4f}")
        print(f"{str(fb)+' bits':>10} " + " ".join(cells))
    prow = []
    for tb, _ in SIZES:
        rows = tables[tb]
        acc = sum(1 for r in rows if (r['PI'] < 1) == r['ok']) / len(rows)
        prow.append(f"{acc:>9.4f}")
    print(f"{'PI<1 (free)':>10} " + " ".join(prow))
    print("\nThe last row uses NO fitted parameter.  If it matches the fitted")
    print("rows, the two degrees of freedom were never being used.")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP Y4: double conditioning — does PI beat both parents inside")
    print("        NU quintiles AND inside eff strata?")
    print("-" * 78)
    for bits, m in SIZES:
        rows = tables[bits]
        print(f"\n{bits} bits, conditioned on NU quintile:")
        edges = quantile_bins([r['NU'] for r in rows], 5)
        print(f"{'NU bin':>18} {'N':>5} {'rec':>8} | {'AUC PI':>8} "
              f"{'AUC NU':>8} {'AUC nuhat':>10}")
        lo = float('-inf')
        for bi in range(5):
            hi = edges[bi] if bi < len(edges) else float('inf')
            sub = [r for r in rows if bin_of(r['NU'], edges) == bi]
            pos = [r for r in sub if r['ok']]
            neg = [r for r in sub if not r['ok']]
            rng = (f"[{lo:.3f},{hi:.3f})" if bi < 4 else f"[{lo:.3f},inf)")
            if bi == 0:
                rng = f"(-inf,{hi:.3f})"
            if not pos or not neg:
                print(f"{rng:>18} {len(sub):>5} "
                      f"{str(len(pos))+'/'+str(len(sub)):>8} |   (degenerate)")
            else:
                cells = [auc([r[k] for r in pos], [r[k] for r in neg])
                         for k in ('PI', 'NU', 'nuhat')]
                print(f"{rng:>18} {len(sub):>5} "
                      f"{str(len(pos))+'/'+str(len(sub)):>8} | {cells[0]:>8.4f} "
                      f"{cells[1]:>8.4f} {cells[2]:>10.4f}")
            lo = hi

        print(f"{bits} bits, conditioned on eff stratum:")
        print(f"{'eff':>18} {'N':>5} {'rec':>8} | {'AUC PI':>8} "
              f"{'AUC NU':>8} {'AUC nuhat':>10}")
        for eff in EFFS:
            sub = [r for r in rows if r['effq'] == eff]
            pos = [r for r in sub if r['ok']]
            neg = [r for r in sub if not r['ok']]
            if not pos or not neg:
                print(f"{eff:>18.2f} {len(sub):>5} "
                      f"{str(len(pos))+'/'+str(len(sub)):>8} |   (degenerate)")
                continue
            cells = [auc([r[k] for r in pos], [r[k] for r in neg])
                     for k in ('PI', 'NU', 'nuhat')]
            print(f"{eff:>18.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>8} | {cells[0]:>8.4f} "
                  f"{cells[1]:>8.4f} {cells[2]:>10.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP Y5: pooled across all three sizes")
    print("-" * 78)
    pos = [r for r in allrows if r['ok']]
    neg = [r for r in allrows if not r['ok']]
    print(f"N = {len(allrows)}  ({len(pos)} recovered)")
    for k in ('PI', 'NU', 'nuhat', 'eff', 'lamstar'):
        a, lo, hi = auc_boot([r[k] for r in pos], [r[k] for r in neg])
        print(f"  AUC(-{k:<8}) = {a:.4f}  [{lo:.4f}, {hi:.4f}]")
    tp = sum(1 for r in allrows if r['PI'] < 1 and r['ok'])
    fp = sum(1 for r in allrows if r['PI'] < 1 and not r['ok'])
    fn = sum(1 for r in allrows if r['PI'] >= 1 and r['ok'])
    tn = sum(1 for r in allrows if r['PI'] >= 1 and not r['ok'])
    print(f"\nPI < 1 pooled:  TP {tp}  FP {fp}  FN {fn}  TN {tn}   "
          f"accuracy {(tp+tn)/len(allrows):.4f}")
    print(f"  precision {tp/(tp+fp) if tp+fp else float('nan'):.4f}   "
          f"recall {tp/(tp+fn) if tp+fn else float('nan'):.4f}")
    print(f"Spearman(log NU, log nu_hat) pooled = "
          f"{spearman([math.log(r['NU']) for r in allrows], [math.log(r['nuhat']) for r in allrows]):.4f}")
    print("(Thread 24 W6 measured this ~0 at fixed eff; PI is a product of "
          "two near-independent factors.)")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP Y6: what PI < 1 actually says")
    print("-" * 78)
    print("PI = NU*nu_hat < 1  <=>  NU < 1/nu_hat.  Nearest-plane is")
    print("guaranteed only for NU < 1 (a theorem).  So 1/nu_hat is the factor")
    print("by which Kannan-LLL empirically RELAXES the nearest-plane bound,")
    print("and it is the skewness sqrt(det L2)/lambda_1(L2) of the lambda-block.\n")
    HERMITE = (4.0 / 3.0) ** 0.25
    nh = [r['nuhat'] for r in allrows]
    print(f"nu_hat over all {len(allrows)}: min {min(nh):.4f}  "
          f"median {sorted(nh)[len(nh)//2]:.4f}  max {max(nh):.4f}")
    print(f"Hermite bound for rank 2, (4/3)^(1/4) = {HERMITE:.4f}  "
          f"-> violations: {sum(1 for v in nh if v > HERMITE + 1e-9)}")
    rel = sorted(1.0 / v for v in nh)
    print(f"relaxation factor 1/nu_hat: min {rel[0]:.4f}  "
          f"median {rel[len(rel)//2]:.4f}  90th pct {rel[int(0.9*len(rel))]:.4f}  "
          f"max {rel[-1]:.4f}")

    print("\nThe cases that matter: NU > 1 (no nearest-plane guarantee) but "
          "PI < 1.")
    print(f"{'bits':>5} {'N(NU>1)':>8} {'of those PI<1':>14} "
          f"{'recovered':>11} {'rate':>7} | {'NU>1 & PI>=1 rate':>18}")
    for bits, m in SIZES:
        rows = tables[bits]
        hard = [r for r in rows if r['NU'] > 1.0]
        resc = [r for r in hard if r['PI'] < 1.0]
        rest = [r for r in hard if r['PI'] >= 1.0]
        rr = sum(1 for r in resc if r['ok'])
        r2 = sum(1 for r in rest if r['ok'])
        print(f"{bits:>5} {len(hard):>8} {len(resc):>14} {rr:>11} "
              f"{rr/len(resc) if resc else float('nan'):>7.4f} | "
              f"{r2/len(rest) if rest else float('nan'):>18.4f}")
    print("\nIf the left rate is high and the right rate is low, 1/nu_hat is")
    print("exactly the right amount of relaxation and PI is not merely NU")
    print("with noise added.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
