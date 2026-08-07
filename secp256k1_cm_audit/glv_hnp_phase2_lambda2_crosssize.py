"""
GLV-HNP Phase 2, Thread 25c: does the lambda_2 law survive a change of size?

Thread 25b established, on a single 1800-instance table at 17 bits, that

        lambda_2(L2)  large  ->  Kannan-LLL recovers

with pooled rank-AUC 0.950, beating the exact BDD certificate NU (0.828) and
Thread 20's nu_hat = mu/sqrt(det L2) (0.674).  A one-size regularity is not a
law.  This script re-runs the same measurement at 12, 14, 17 and 20 bits
(dim 16, 20, 24, 28) and reports, per size:

  * AUC(+lambda_2), pooled and held-out by seed parity
  * the same for NU, nu_hat, mu
  * AUC(-n) as the size control — near 0.5 means the effect is geometry
  * the freely fitted exponent alpha in  score = mu / det(L2)^alpha.
    nu_hat asserts alpha = 0.5; the lambda_2 law asserts alpha = 1, because
    Thread 24 W2 gives det(L2) = mu * lambda_2 to within 4.3%.

Run: python3 glv_hnp_phase2_lambda2_crosssize.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import search_curves
from glv_hnp_phase2_projected import run_new
from glv_hnp_phase2_gsprofile import instance, auc
from glv_hnp_phase2_nu_conditioned import SEEDS25, EFFS, logistic_fit

SIZES = ((12, 8), (14, 10), (17, 12), (20, 14))
NSEED = 10


def collect(bits, m, seeds):
    curves = search_curves(1 << (bits - 1), 1 << bits, per_bin=2, nbins=10)
    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in seeds:
                d = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), m, d, k1b, seed, exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), m, d, k1b, seed)
                rows.append({'ok': bool(rk['ok']), 'l2': r['l2'], 'NU': r['NU'],
                             'mu': r['mu'], 'nuhat': r['nuhat'],
                             'det2': r['det2'], 'n': n, 'K1': k1b,
                             'effq': eff, 'seed': seed})
    return curves, rows


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25c — cross-size validation of the lambda_2 law")
    print("=" * 78)
    seeds = SEEDS25[:NSEED]
    train = set(seeds[0::2])
    summary = []
    for bits, m in SIZES:
        t0 = time.time()
        curves, rows = collect(bits, m, seeds)
        pos = [r for r in rows if r['ok']]
        neg = [r for r in rows if not r['ok']]
        print(f"\n=== {bits} bits (m={m}, dim {2*m}), {len(curves)} curves, "
              f"{len(rows)} instances, {len(pos)}/{len(rows)} recovered, "
              f"{time.time()-t0:.1f}s")
        if not pos or not neg:
            print("    degenerate — one class empty")
            continue
        B = [r for r in rows if r['seed'] not in train]
        A = [r for r in rows if r['seed'] in train]
        hp = [r for r in B if r['ok']]
        hn = [r for r in B if not r['ok']]
        a_l2 = auc([-r['l2'] for r in pos], [-r['l2'] for r in neg])
        h_l2 = auc([-r['l2'] for r in hp], [-r['l2'] for r in hn])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_n = auc([r['n'] for r in pos], [r['n'] for r in neg])
        bb, _ = logistic_fit(
            [[math.log(r['NU']), math.log(r['mu']), math.log(r['det2'])]
             for r in A], [1.0 if r['ok'] else 0.0 for r in A])
        alpha = -bb[3] / bb[2] if abs(bb[2]) > 1e-12 else float('nan')
        print(f"  AUC(+lambda_2) pooled {a_l2:.4f}   held-out {h_l2:.4f}")
        print(f"  AUC(-NU)       pooled {a_nu:.4f}")
        print(f"  AUC(-nu_hat)   pooled {a_nh:.4f}")
        print(f"  AUC(-mu)       pooled {a_mu:.4f}")
        print(f"  AUC(-n)        pooled {a_n:.4f}   [size control]")
        print(f"  free-fit alpha (train) = {alpha:.4f}   "
              f"[nu_hat: 0.5 | lambda_2 law: 1.0]")
        for e in EFFS:
            s = [r for r in rows if r['effq'] == e]
            p2 = [r for r in s if r['ok']]
            n2 = [r for r in s if not r['ok']]
            if p2 and n2:
                print(f"    eff={e:.2f} {len(p2):>3}/{len(s):<4} "
                      f"AUC(+l2)={auc([-r['l2'] for r in p2], [-r['l2'] for r in n2]):.4f}")
        summary.append((bits, a_l2, h_l2, a_nu, a_nh, a_mu, a_n, alpha))

    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(f"{'bits':>5} {'AUC +l2':>9} {'held-out':>9} {'AUC -NU':>9} "
          f"{'AUC -nuhat':>11} {'AUC -mu':>9} {'AUC -n':>8} {'alpha':>7}")
    for b, a, h, nu, nh, mu, nn, al in summary:
        print(f"{b:>5} {a:>9.4f} {h:>9.4f} {nu:>9.4f} {nh:>11.4f} "
              f"{mu:>9.4f} {nn:>8.4f} {al:>7.4f}")
    print("\ndone")
