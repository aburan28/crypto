"""
GLV-HNP — Thread 23 (T6): is the EXPLANATORY invariant also PREDICTIVE?

Thread 23 (glv_hnp_nuhat_cf_closed_form.py) proved that nu_hat is a continued
fraction quantity of lam/n, and that it is small exactly when lam/n has a large
partial quotient a_loc located at the scale q* = sqrt(n*K2/K1).  On random
instances spearman(nu_hat, a_loc) = -0.77 while the scale-free max_a manages
only -0.30.

That is a statement about nu_hat, not about the attack.  This script closes the
loop by replaying the 2026-06-29 Exp S protocol exactly as Thread 20d did
(K1=72, m=12, 6 seeds, 100 fresh 20-bit j=0 CM curves; C2 = recovers d on all
6 seeds) and scoring FOUR invariants head to head against the C1/C2 labels:

    nu_hat  the Thread 20 predictor          (AUC 0.935 on 2026-07-29)
    a_loc   scale-localised partial quotient (Thread 23 -- new)
    max_a   scale-free max partial quotient  (falsified Exp S 2026-06-29)
    mu      lam/n symmetrised               (falsified twice)

Pre-registered prediction: a_loc lands near nu_hat and far above max_a.  If
a_loc scores at the trivial baseline, the scale-localisation story explains
nu_hat but has no independent predictive content, and should be reported as
explanatory only.
"""

import math
import random
import sys
import time
import importlib.util

import sympy

_spec = importlib.util.spec_from_file_location(
    "_nuhat", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_lib.py")
_nuhat = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_nuhat)

_spec2 = importlib.util.spec_from_file_location(
    "_t23", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_cf_closed_form.py")
_t23 = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(_t23)

eisenstein_decompose = _nuhat.eisenstein_decompose
j0_traces = _nuhat.j0_traces
glv_eigenvalues = _nuhat.glv_eigenvalues
mu_of = _nuhat.mu_of
identify_twist = _nuhat.identify_twist
find_generator = _nuhat.find_generator
rival_sublattice_nu = _nuhat.rival_sublattice_nu
run_experiment = _nuhat.run_experiment
scales = _nuhat.scales

max_partial_quotient = _t23.max_partial_quotient
scale_local_partial_quotient = _t23.scale_local_partial_quotient
nu_hat_cf = _t23.nu_hat_cf

# Exp S protocol (2026-06-29), identical to Thread 20d
K1_FIXED = 72
M_FIXED = 12
SEEDS_N = 6
N_CURVES = 100
BITS_LO, BITS_HI = 2 ** 19, 2 ** 20


def harvest(n_curves):
    out, seen = [], set()
    p = int(sympy.nextprime(BITS_LO))
    while p < BITS_HI and len(out) < n_curves:
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
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float('nan')
    tot = 0.0
    for a in pos:
        for b in neg:
            tot += 1.0 if a > b else 0.5 if a == b else 0.0
    return tot / (len(pos) * len(neg))


def best_cut_accuracy(scores, labels):
    best, cut = 0.0, None
    for c in sorted(set(scores)):
        for sign in (1, -1):
            acc = sum(1 for s, l in zip(scores, labels)
                      if ((s > c) if sign > 0 else (s <= c)) == l) / len(labels)
            if acc > best:
                best, cut = acc, (c, sign)
    return best, cut


def main():
    print("=" * 78)
    print("GLV-HNP Thread 23 (T6): a_loc vs nu_hat vs the falsified invariants")
    print("=" * 78)
    print(f"Exp S protocol: K1={K1_FIXED}, m={M_FIXED}, {SEEDS_N} seeds, "
          f"{N_CURVES} fresh 20-bit j=0 CM curves")
    print("C2 = recovers d on all seeds; C1 = does not (Exp S convention)")

    t0 = time.time()
    curves = harvest(N_CURVES)
    print(f"\n  harvested {len(curves)} curves in {time.time()-t0:.1f}s")

    t1 = time.time()
    for c in curves:
        n = c['n']
        k2 = math.isqrt(n) + 1
        wins = 0
        for s in range(SEEDS_N):
            d = random.Random(s + 7777).randint(1, n - 1)
            wins += run_experiment(c['p'], n, c['lam'], c['G'], M_FIXED, d,
                                   K1_FIXED, k2, s)
        c['wins'] = wins
        c['C2'] = (wins == SEEDS_N)
        _, _, c['nu_hat'] = rival_sublattice_nu(n, c['lam'], K1_FIXED, k2)
        # cross-check the closed form on every curve of the sample
        nh_cf, _, _, _ = nu_hat_cf(n, c['lam'], K1_FIXED, k2)
        c['cf_err'] = abs(nh_cf - c['nu_hat'])
        A, _, B, _ = scales(n, K1_FIXED, k2)
        a_loc, _, _ = scale_local_partial_quotient(n, c['lam'], A, B)
        c['a_loc'] = a_loc
        c['mu'] = mu_of(c['lam'], n)
        c['max_a'] = max_partial_quotient(c['lam'], n)
    print(f"  {len(curves)*SEEDS_N} LLL trials in {time.time()-t1:.1f}s")
    print(f"  max |nu_hat(closed form) - nu_hat(Gauss)| over the sample = "
          f"{max(c['cf_err'] for c in curves):.3e}   (Theorem T23 re-check)")

    c2 = [c for c in curves if c['C2']]
    c1 = [c for c in curves if not c['C2']]
    base = max(len(c1), len(c2)) / len(curves)
    print(f"\n  C1 (fails): {len(c1)}/{len(curves)}   C2 (6/6): {len(c2)}/{len(curves)}")
    print(f"  trivial baseline (always predict majority) = {base*100:.1f}%")

    print("\n--- Head to head against the C1/C2 labels ---")
    print(f"{'invariant':>10} {'C1 range':>20} {'C2 range':>20} {'AUC':>7} {'best acc':>9}")
    labels = [c['C2'] for c in curves]
    for key in ('nu_hat', 'a_loc', 'max_a', 'mu'):
        sc = [c[key] for c in curves]
        a = auc(sc, labels)
        # report the orientation-free AUC (invariants point in opposite ways)
        a_or = max(a, 1 - a)
        acc, _cut = best_cut_accuracy(sc, labels)
        r1 = [c[key] for c in c1]
        r2 = [c[key] for c in c2]
        def rng(v):
            return f"[{min(v):.3f},{max(v):.3f}]" if v else "n/a"
        print(f"{key:>10} {rng(r1):>20} {rng(r2):>20} {a_or:>7.3f} {acc*100:>8.1f}%")
    print(f"{'':>10} {'':>20} {'trivial baseline':>20} {'':>7} {base*100:>8.1f}%")

    print("\n--- a_loc response (mean wins/6 by bucket) ---")
    print(f"  {'a_loc':>8} {'count':>6} {'C2 rate':>9} {'mean wins/6':>12} {'mean nu_hat':>12}")
    for lo, hi in ((1, 1), (2, 2), (3, 4), (5, 8), (9, 16), (17, 10 ** 9)):
        grp = [c for c in curves if lo <= c['a_loc'] <= hi]
        if not grp:
            continue
        lab = f"{lo}" if lo == hi else f"{lo}-{hi if hi < 10**9 else '+'}"
        print(f"  {lab:>8} {len(grp):>6} "
              f"{sum(1 for c in grp if c['C2'])/len(grp):>9.2f} "
              f"{sum(c['wins'] for c in grp)/len(grp):>12.2f} "
              f"{sum(c['nu_hat'] for c in grp)/len(grp):>12.3f}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
