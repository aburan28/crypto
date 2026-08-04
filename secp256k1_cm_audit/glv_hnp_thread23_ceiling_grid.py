"""
GLV-HNP -- Thread 23b: does the nu_hat C2 ceiling generalise across (K1, m)?

Thread 20d (2026-07-29, e845207) measured the C1/C2 separation at ONE cell of
the parameter space: K1 = 72, m = 12.  It found

    every C2 curve had nu_hat <= 0.645; no curve with nu_hat > 0.645 was C2
    AUC 0.935, best single-cut accuracy 89%

The 2026-07-29 next-step proposal asked for the cheap generalisation check:
re-run at K1 in {36, 72, 144} and m in {10, 12, 14}.  If the ceiling stays near
0.65 across all nine cells, it is a property of the lattice family rather than
of that one sample.

Prediction under the Thread 23a closed form: it should NOT stay fixed.  nu_hat
depends on (n, lam, K1, K2) through S1/S2 = K2/K1, so changing K1 moves the
critical scale q_crit = sqrt(n*K2/K1) and rescales nu_hat itself.  m does not
enter L2 at all, but larger m gives the attacker more planted structure, so the
ceiling should rise with m at fixed K1.  Both are tested here.

Run: python3 glv_hnp_thread23_ceiling_grid.py
"""

import importlib.util
import math
import random
import sys
import time

_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_lib.py")
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

_s23 = importlib.util.spec_from_file_location(
    "_t23", __file__.rsplit("/", 1)[0] + "/glv_hnp_thread23_cf_closedform.py")
_t23 = importlib.util.module_from_spec(_s23)
_s23.loader.exec_module(_t23)

rival_sublattice_nu = _t20a.rival_sublattice_nu
run_experiment = _t20a.run_experiment
mu_of = _t20a.mu_of
auc = None  # defined below

K1_GRID = (36, 72, 144)
M_GRID = (10, 12, 14)
SEEDS_N = 6
N_CURVES = 100
BITS_LO, BITS_HI = 2 ** 19, 2 ** 20


def harvest(n_curves):
    import sympy
    out, seen = [], set()
    p = int(sympy.nextprime(BITS_LO))
    while p < BITS_HI and len(out) < n_curves:
        if p % 3 == 1:
            eis = _t20a.eisenstein_decompose(p)
            if eis is not None:
                for t in _t20a.j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n in seen or n % 3 != 1 or not sympy.isprime(n):
                        continue
                    if _t20a.glv_eigenvalues(n) is None:
                        continue
                    b = _t20a.identify_twist(p, n)
                    if b is None:
                        break
                    G = _t20a.find_generator(p, b, n)
                    if G is None:
                        break
                    seen.add(n)
                    out.append({'p': p, 'b': b, 'n': n,
                                'lam': _t20a.glv_eigenvalues(n)[0], 'G': G})
                    break
        p = int(sympy.nextprime(p))
    return out


def auc_of(scores, labels):
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
    print("=" * 90)
    print("GLV-HNP Thread 23b -- nu_hat C2 ceiling across the (K1, m) grid")
    print("=" * 90)
    print(f"{N_CURVES} fresh 20-bit j=0 CM curves, {SEEDS_N} seeds/cell; "
          f"C2 = recovers d on all {SEEDS_N} seeds")
    print("anchor cell (K1=72, m=12) reproduced 2026-07-29: ceiling 0.645, AUC 0.935\n")

    t0 = time.time()
    curves = harvest(N_CURVES)
    print(f"harvested {len(curves)} curves in {time.time()-t0:.1f}s\n")

    print(f"{'K1':>5} {'m':>4} {'C2/N':>8} {'C2 ceiling':>11} {'C1 floor':>9} "
          f"{'AUC':>7} {'best acc':>9} {'baseline':>9} {'gain':>6} "
          f"{'nu_hat med':>11} {'sep?':>5}")
    rows = []
    for k1 in K1_GRID:
        for m in M_GRID:
            for c in curves:
                n = c['n']
                k2 = math.isqrt(n) + 1
                wins = 0
                for s in range(SEEDS_N):
                    d = random.Random(s + 7777).randint(1, n - 1)
                    wins += bool(run_experiment(c['p'], n, c['lam'], c['G'],
                                                m, d, k1, k2, s))
                c['wins'] = wins
                c['C2'] = (wins == SEEDS_N)
                _, _, c['nu_hat'] = rival_sublattice_nu(n, c['lam'], k1, k2)
            c2 = [c for c in curves if c['C2']]
            c1 = [c for c in curves if not c['C2']]
            nh = sorted(c['nu_hat'] for c in curves)
            med = nh[len(nh) // 2]
            ceil_ = max((c['nu_hat'] for c in c2), default=float('nan'))
            floor_ = min((c['nu_hat'] for c in c1), default=float('nan'))
            a = auc_of([c['nu_hat'] for c in curves], [c['C2'] for c in curves])
            a_or = max(a, 1 - a) if a == a else float('nan')
            best = 0.0
            for cut in nh:
                for sign in (1, -1):
                    acc = sum(1 for c in curves
                              if ((c['nu_hat'] > cut) if sign > 0
                                  else (c['nu_hat'] <= cut)) == c['C2']) / len(curves)
                    best = max(best, acc)
            sep = "yes" if (floor_ == floor_ and ceil_ == ceil_ and floor_ > ceil_) else "no"
            base = max(len(c1), len(c2)) / len(curves)
            print(f"{k1:>5} {m:>4} {len(c2):>4}/{len(curves):<3} {ceil_:>11.3f} "
                  f"{floor_:>9.3f} {a_or:>7.3f} {best*100:>8.1f}% {base*100:>8.1f}% "
                  f"{(best-base)*100:>+5.1f} {med:>11.3f} {sep:>5}")
            rows.append((k1, m, len(c2), ceil_, a_or, best, med))

    print("\n--- ceiling normalised by the sample median nu_hat at that K1 ---")
    print(f"{'K1':>5} {'m':>4} {'ceiling':>9} {'median':>9} {'ceiling/median':>15}")
    for (k1, m, nc2, ceil_, a_or, best, med) in rows:
        r = ceil_ / med if med else float('nan')
        print(f"{k1:>5} {m:>4} {ceil_:>9.3f} {med:>9.3f} {r:>15.3f}")

    print("\nDone.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
