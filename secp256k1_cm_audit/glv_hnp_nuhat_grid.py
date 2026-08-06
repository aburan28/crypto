"""
GLV-HNP Phase 2 -- Thread 23 secondary: does the nu_hat cut generalise?

Thread 20d (2026-07-29) established nu_hat as a Phase 2 separator against the
June C1/C2 classes at a SINGLE operating point: K1 = 72, m = 12, eff = 0.0993.
Every C2 curve in that sample had nu_hat <= 0.645, and no curve with
nu_hat > 0.645 was C2.  The open question, in the 2026-07-29 words:

    "20d used one m and one K1.  Check whether the nu_hat cut generalises by
     re-running 20d at K1 in {36, 72, 144} and m in {10, 12, 14} -- if the C2
     ceiling stays near nu_hat ~ 0.65 across all nine cells, the cut is a
     property of the lattice family, not of the sample."

This script runs that 3x3 grid on one shared curve sample, so the cells differ
only in (K1, m) and not in which curves were drawn.

nu_hat is computed with the Thread 23 exact-integer CF closed form
(glv_hnp_nuhat_cf_closed_form.nu_hat_exact), which E2 of that script verified
is bitwise identical to the float Lagrange-Gauss incumbent on 300/300 curves.

Class convention (matches June Exp S / 20d, and is the opposite of the naive
reading): **C2 is the EASY class**.  20d's decile table settles it -- decile 1
(lowest nu_hat) has "C2 rate 1.00, mean wins/6 = 6.00", i.e. full recovery.  So
here a curve is C2 in a cell iff it wins ALL SEEDS_N trials, C1 otherwise.  This
also matches the direction of the effect: low nu_hat = easier attack.

Sanity anchor: at K1=72, m=12 this reproduces 20d's "74/100 C1" (we get 73).

Run: python3 glv_hnp_nuhat_grid.py
"""

import importlib.util
import math
import os
import random
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))

def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_base = _load("_base", "glv_hnp_nuhat_base.py")
_cf = _load("_cf", "glv_hnp_nuhat_cf_closed_form.py")

eisenstein_decompose = _base.eisenstein_decompose
j0_traces = _base.j0_traces
glv_eigenvalues = _base.glv_eigenvalues
mu_of = _base.mu_of
identify_twist = _base.identify_twist
find_generator = _base.find_generator
run_experiment = _base.run_experiment
nu_hat_exact = _cf.nu_hat_exact
global_max_a = _cf.global_max_a

import sympy

K1_GRID = [36, 72, 144]
M_GRID = [10, 12, 14]
SEEDS_N = 6
N_CURVES = 100
BITS_LO, BITS_HI = 2 ** 19, 2 ** 20

def harvest(want):
    """Full j=0 CM curves (need b and G to actually run signatures)."""
    out, seen = [], set()
    p = int(sympy.nextprime(BITS_LO))
    while p < BITS_HI and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n in seen or n % 3 != 1 or not sympy.isprime(n):
                        continue
                    roots = glv_eigenvalues(n)
                    if roots is None:
                        continue
                    b = identify_twist(p, n)
                    if b is None:
                        break
                    G = find_generator(p, b, n)
                    if G is None:
                        break
                    seen.add(n)
                    out.append({'p': p, 'b': b, 'n': n, 'lam': roots[0], 'G': G})
                    break
        p = int(sympy.nextprime(p))
    return out

def auc(scores, labels):
    """AUC of `scores` separating labels (True = C2 = easy = low nu_hat)."""
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float('nan')
    wins = sum((a < b) + 0.5 * (a == b) for a in pos for b in neg)
    return wins / (len(pos) * len(neg))

def best_cut(scores, labels):
    """Best single threshold accuracy (predict C2 when score <= cut)."""
    best = (0.0, None)
    for c in sorted(set(scores)):
        acc = sum((s <= c) == l for s, l in zip(scores, labels)) / len(scores)
        if acc > best[0]:
            best = (acc, c)
    return best

def main():
    t0 = time.time()
    print("=" * 78)
    print("GLV-HNP Thread 23 secondary: nu_hat cut across K1 x m")
    print("=" * 78)
    print(f"  {N_CURVES} shared 20-bit curves, {SEEDS_N} seeds/cell-curve, "
          f"K1 in {K1_GRID}, m in {M_GRID}")
    print("  nu_hat via Thread 23 exact-integer CF closed form")

    curves = harvest(N_CURVES)
    print(f"  harvested {len(curves)} curves in {time.time()-t0:.1f}s")

    print(f"\n{'K1':>5} {'m':>3} {'eff':>7} {'C2':>5} {'base':>6} "
          f"{'AUC':>6} {'acc':>6} {'cut':>6} {'C2 ceil':>8} {'C1 min':>7}")
    print("-" * 72)

    summary = []
    for k1 in K1_GRID:
        for m in M_GRID:
            nus, labels = [], []
            for c in curves:
                n, lam = c['n'], c['lam']
                k2 = math.isqrt(n) + 1
                nh, _, _, _ = nu_hat_exact(n, lam, k1, k2)
                wins = 0
                for s in range(SEEDS_N):
                    d = random.Random(s + 7777 + 1000 * m).randint(1, n - 1)
                    wins += run_experiment(c['p'], n, lam, c['G'], m, d,
                                           k1, k2, s + 100 * m)
                nus.append(nh)
                labels.append(wins == SEEDS_N)     # True = C2 (easy: full recovery)
            n_c2 = sum(labels)
            base = max(n_c2, len(labels) - n_c2) / len(labels)
            a = auc(nus, labels)
            acc, cut = best_cut(nus, labels)
            c2_ceil = max((s for s, l in zip(nus, labels) if l), default=float('nan'))
            c1_min = min((s for s, l in zip(nus, labels) if not l), default=float('nan'))
            eff = k1 * (math.isqrt(curves[0]['n']) + 1) / curves[0]['n']
            print(f"{k1:>5} {m:>3} {eff:>7.4f} {n_c2:>5} {base:>6.2f} "
                  f"{a:>6.3f} {acc:>6.3f} {cut:>6.3f} {c2_ceil:>8.3f} {c1_min:>7.3f}")
            summary.append((k1, m, n_c2, base, a, acc, cut, c2_ceil))

    print("\n  C2 = curve winning ALL %d/%d seeds (easy). 'cut' = best single threshold" % (SEEDS_N, SEEDS_N))
    print("  (predict C2 when nu_hat <= cut). 'C2 ceil' = largest nu_hat among C2 curves.")

    valid = [s for s in summary if not math.isnan(s[7]) and 0 < s[2] < N_CURVES]
    if valid:
        ceils = [s[7] for s in valid]
        aucs = [s[4] for s in valid]
        lifts = [s[5] - s[3] for s in valid]
        print(f"\n  cells with both classes present: {len(valid)}/9")
        print(f"  C2 ceiling: min={min(ceils):.3f} max={max(ceils):.3f} "
              f"mean={sum(ceils)/len(ceils):.3f}   (20d single-cell value: 0.645)")
        print(f"  AUC:        min={min(aucs):.3f} max={max(aucs):.3f} "
              f"mean={sum(aucs)/len(aucs):.3f}   (20d: 0.935)")
        print(f"  acc - baseline: min={min(lifts):+.3f} max={max(lifts):+.3f} "
              f"mean={sum(lifts)/len(lifts):+.3f}")

    print(f"\nTotal wall time {time.time()-t0:.1f}s")
    print("Done.")

if __name__ == '__main__':
    sys.exit(main())
