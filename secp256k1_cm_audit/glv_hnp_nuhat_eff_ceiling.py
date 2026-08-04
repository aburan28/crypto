"""
GLV-HNP — Thread 21b: is the nu_hat "C2 ceiling" universal, or per-eff?

Why
---
Thread 20d (2026-07-29) measured, at a single operating point (K1=72 fixed,
sample mean eff = K1*K2/n = 0.0993, m=12, 6 seeds, 20-bit):

    AUC(nu_hat) = 0.935, and every C2 curve had nu_hat <= 0.645,

and concluded from a 20/24-bit extrapolation that secp256k1 (nu_hat = 0.664)
"is not in the low-nu_hat easy tail at any eff tested".

Thread 21 replaced the extrapolation with an exact closed form, nu_hat =
F(tau, theta) with tau = S_K1/S_K2 ~ 1/eff and theta the Eisenstein CM angle,
and computed secp256k1 exactly:

    eff     0.0200  0.0500  0.0993  0.1500  0.2500  0.5000
    nu_hat  0.8782  0.8709  0.6639  0.5990  0.5852  0.6851

So secp256k1 sits BELOW 0.645 at eff = 0.15 and 0.25.  Taken at face value that
contradicts the 20d conclusion.  But 0.645 was calibrated at eff ~ 0.0993 only,
and Thread 20a already showed the attack itself weakens fast with eff (at
eff=0.15 only 3/20 curves recovered; at eff=0.25, 0/20).  A fixed numeric
ceiling cannot be transported across eff.

H-ceiling: the ceiling max{nu_hat : curve is C2} is eff-dependent, and at large
eff the C2 class is empty, making "nu_hat <= 0.645" vacuous rather than
favourable.

Falsifier: if C2 curves still exist at eff = 0.15/0.25 with a ceiling near
0.645, then the ceiling IS transportable and secp256k1 really is in the easy
tail there.

Run: python3 glv_hnp_nuhat_eff_ceiling.py
"""

import math
import random
import sys
import time

import sympy

import importlib.util


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, __file__.rsplit("/", 1)[0] + "/" + fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_t20c = _load("_t20c", "glv_hnp_thread20c_lib.py")
_t21 = _load("_t21", "glv_hnp_nuhat_closed_form.py")

eisenstein_decompose = _t20c.eisenstein_decompose
j0_traces = _t20c.j0_traces
glv_eigenvalues = _t20c.glv_eigenvalues
identify_twist = _t20c.identify_twist
find_generator = _t20c.find_generator
rival_sublattice_nu = _t20c.rival_sublattice_nu
run_experiment = _t20c.run_experiment

nu_hat_from_cm = _t21.nu_hat_from_cm
SECP256K1_N = _t21.SECP256K1_N

BITS_LO, BITS_HI = 2 ** 19, 2 ** 20
N_CURVES = 120
SEEDS_N = 6
M_LIST = [12, 20]
EFF_LIST = [0.02, 0.05, 0.0993, 0.15, 0.25]


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


def main():
    t0 = time.time()
    print("=" * 78)
    print("GLV-HNP Thread 21b: eff-dependence of the nu_hat C2 ceiling")
    print("=" * 78)
    curves = harvest(N_CURVES)
    print(f"  harvested {len(curves)} fresh 20-bit j=0 curves in "
          f"{time.time()-t0:.1f}s\n")

    trials = 0
    for m in M_LIST:
        print("=" * 78)
        print(f"m = {m} signatures, {SEEDS_N} seeds/curve, C2 = all seeds recover")
        print("=" * 78)
        print(f"  {'eff':>7} {'C2 rate':>9} {'AUC':>7} {'best acc':>9} "
              f"{'baseline':>9} {'C2 ceiling':>11} {'secp nu_hat':>12} {'below?':>7}")
        for eff in EFF_LIST:
            labels, scores = [], []
            for c in curves:
                n = c['n']
                k2 = math.isqrt(n) + 1
                k1 = max(2, int(eff * n / k2))
                wins = 0
                for s in range(SEEDS_N):
                    d = random.Random(s + 7777).randint(1, n - 1)
                    wins += run_experiment(c['p'], n, c['lam'], c['G'], m, d,
                                           k1, k2, s)
                    trials += 1
                labels.append(wins == SEEDS_N)
                scores.append(rival_sublattice_nu(n, c['lam'], k1, k2)[2])
            c2 = [s for s, l in zip(scores, labels) if l]
            rate = len(c2) / len(labels)
            base = max(sum(labels), len(labels) - sum(labels)) / len(labels)
            a = auc(scores, labels)
            a = max(a, 1 - a) if a == a else a
            best = 0.0
            for cut in sorted(set(scores)):
                for sgn in (1, -1):
                    acc = sum(1 for s, l in zip(scores, labels)
                              if ((s > cut) if sgn > 0 else (s <= cut)) == l)
                    best = max(best, acc / len(labels))
            ceil_txt = f"{max(c2):.4f}" if c2 else "-- (empty)"

            k1s = max(2, int(eff * SECP256K1_N / (math.isqrt(SECP256K1_N) + 1)))
            nh_s = nu_hat_from_cm(SECP256K1_N, glv_eigenvalues(SECP256K1_N)[0],
                                  k1s, math.isqrt(SECP256K1_N) + 1)[0]
            below = ("vacuous" if not c2
                     else ("YES" if nh_s <= max(c2) else "no"))
            auc_txt = f"{a:.3f}" if a == a else "  n/a"
            print(f"  {eff:>7.4f} {rate:>8.1%} {auc_txt:>7} {best:>8.1%} "
                  f"{base:>8.1%} {ceil_txt:>11} {nh_s:>12.4f} {below:>7}")
        print()

    print(f"  total LLL trials: {trials}")
    print("""
Reading
-------
"C2 ceiling" is max nu_hat over curves that recover on all seeds -- the
quantity Thread 20d reported as 0.645.  Where the C2 class is empty the
ceiling is undefined and the comparison to secp256k1 is vacuous: no curve
of any angle is attackable at that eff, so a low nu_hat buys nothing.""")
    print(f"\nTotal wall time {time.time()-t0:.1f}s")
    print("Done.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
