"""
Thread 23 (secondary) — is the nu_hat ~ 0.65 C2 ceiling a property of the
lattice family, or an artefact of the single (K1, m) cell used in 20d?

2026-07-29 run #2 measured the C2 ceiling at ONE cell: K1 = 72, m = 12
(`glv_hnp_nuhat_vs_c1c2.py`), finding "every C2 curve has nu_hat <= 0.645".
Its own next-step proposal asked for K1 in {36, 72, 144} x m in {10, 12, 14}:

    "if the C2 ceiling stays near nu_hat ~ 0.65 across all nine cells, the cut
     is a property of the lattice family, not of the sample."

This script runs that 3x3 grid on ONE shared curve sample (so the cells differ
only in (K1, m), not in which curves were drawn), and reports per cell:

    C2 rate, AUC of nu_hat against the C2 label, and the C2 ceiling
    (max nu_hat over C2 curves).

nu_hat is computed with the Thread 23 closed form
(`thread23_nuhat_cf_closed_form.nu_hat_cf`), which E1 of that script verified
equals the Lagrange-Gauss value exactly on 640/640 instances.  The grid also
records a_at_t, the partial quotient of lam/n straddling t = sqrt(n*S_K1/S_K2),
which the closed form identifies as the mechanism behind nu_hat.  Note a_at_t
depends on K1 (through S_K1) but not on m.

Run: python3 secp256k1_cm_audit/thread23_nuhat_k1m_grid.py
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

_t20a = _load("_t20a", "glv_hnp_thread20a_nuhat_base.py")
_t23 = _load("_t23", "thread23_nuhat_cf_closed_form.py")

eisenstein_decompose = _t20a.eisenstein_decompose
j0_traces = _t20a.j0_traces
glv_eigenvalues = _t20a.glv_eigenvalues
identify_twist = _t20a.identify_twist
find_generator = _t20a.find_generator
run_experiment = _t20a.run_experiment
rival_sublattice_nu = _t20a.rival_sublattice_nu

nu_hat_cf = _t23.nu_hat_cf

# ---------------------------------------------------------------------------

K1_GRID = [36, 72, 144]
M_GRID = [10, 12, 14]
SEEDS_N = 4
N_CURVES = 40
BITS_LO, BITS_HI = 2 ** 19, 2 ** 20
CEILING_REF = 0.645          # the 20d C2 ceiling to be tested

def harvest(n_curves):
    """20-bit j=0 CM curves; identical procedure to glv_hnp_nuhat_vs_c1c2.harvest."""
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
    print("=" * 78)
    print("Thread 23 (secondary): nu_hat C2 ceiling across the K1 x m grid")
    print("=" * 78)
    print(f"K1 in {K1_GRID}, m in {M_GRID}, {SEEDS_N} seeds, "
          f"{N_CURVES} shared 20-bit j=0 CM curves")
    print(f"C2 = LLL recovers d on all {SEEDS_N} seeds (Exp S convention)")
    print(f"nu_hat via the Thread 23 CF closed form; reference ceiling "
          f"{CEILING_REF} from 20d (K1=72, m=12)\n")

    t0 = time.time()
    curves = harvest(N_CURVES)
    print(f"harvested {len(curves)} curves in {time.time()-t0:.1f}s")

    rng = random.Random(20260806)
    results = {}
    t1 = time.time()
    for k1 in K1_GRID:
        for c in curves:
            n = c['n']
            k2 = math.isqrt(n) + 1
            cf = nu_hat_cf(n, c['lam'], k1, k2)
            c[f'nu_{k1}'] = cf['nu_hat']
            c[f'a_{k1}'] = cf['a_at_t']
            # sanity: closed form must agree with the shipped Gauss reduction
            _, _, ng = rival_sublattice_nu(n, c['lam'], k1, k2)
            assert abs(ng - cf['nu_hat']) < 1e-9, (n, k1, ng, cf['nu_hat'])
        for m in M_GRID:
            wins_l, c2_l, nu_l = [], [], []
            for c in curves:
                n = c['n']
                k2 = math.isqrt(n) + 1
                wins = 0
                for s in range(SEEDS_N):
                    d = rng.randrange(1, n)
                    wins += run_experiment(c['p'], n, c['lam'], c['G'], m, d,
                                           k1, k2, s)
                wins_l.append(wins)
                c2_l.append(wins == SEEDS_N)
                nu_l.append(c[f'nu_{k1}'])
            results[(k1, m)] = (wins_l, c2_l, nu_l)
            print(f"  cell K1={k1:<4} m={m:<3} done "
                  f"({sum(c2_l)}/{len(curves)} C2)  t={time.time()-t1:.0f}s")

    print(f"\n{len(K1_GRID)*len(M_GRID)*len(curves)*SEEDS_N} LLL trials in "
          f"{time.time()-t1:.1f}s\n")

    print(f"{'K1':>5} {'m':>4} {'C2 rate':>9} {'AUC(nu_hat)':>12} "
          f"{'C2 ceiling':>11} {'C1 min nu':>10} {'q(C2rate)':>10} {'mean a_at_t':>12}")
    print("-" * 78)
    ceilings = []
    for k1 in K1_GRID:
        for m in M_GRID:
            wins_l, c2_l, nu_l = results[(k1, m)]
            n_c2 = sum(c2_l)
            a = auc(nu_l, c2_l)
            a = max(a, 1 - a) if a == a else a
            ceil = max((v for v, l in zip(nu_l, c2_l) if l), default=float('nan'))
            c1min = min((v for v, l in zip(nu_l, c2_l) if not l), default=float('nan'))
            mean_a = sum(c[f'a_{k1}'] for c in curves) / len(curves)
            # If nu_hat ranked perfectly, the ceiling would be exactly the
            # C2-rate quantile of the sample's own nu_hat distribution.
            srt = sorted(nu_l)
            qc = srt[min(len(srt) - 1, max(0, n_c2 - 1))] if n_c2 else float('nan')
            if n_c2:
                ceilings.append(ceil)
            print(f"{k1:>5} {m:>4} {n_c2:>4}/{len(curves):<4} {a:>12.3f} "
                  f"{ceil:>11.3f} {c1min:>10.3f} {qc:>10.3f} {mean_a:>12.2f}")

    print("-" * 78)
    if ceilings:
        print(f"C2 ceiling across {len(ceilings)} non-degenerate cells: "
              f"min={min(ceilings):.3f} max={max(ceilings):.3f} "
              f"mean={sum(ceilings)/len(ceilings):.3f}")
        print(f"20d reference (K1=72, m=12, 100 curves, 6 seeds): {CEILING_REF}")
        within = sum(1 for c in ceilings if abs(c - CEILING_REF) <= 0.15)
        print(f"cells with ceiling within +-0.15 of {CEILING_REF}: "
              f"{within}/{len(ceilings)}")

    print("\nSAMPLE CAVEAT: harvest() walks primes upward from 2^19 and takes the")
    print("first N, exactly as glv_hnp_nuhat_vs_c1c2.harvest does. With N=40 vs")
    print("that script's N=100, this sample is a PREFIX SUBSET of the 20d sample,")
    print("not an independent draw. The K1 x m comparison is still valid (all")
    print("cells share one curve set, which is the point), but agreement with the")
    print("20d numbers at K1=72/m=12 is a consistency check, not a replication.")

    print("\nPer-K1 a_at_t (scale-dependent, m-independent by construction):")
    for k1 in K1_GRID:
        vals = sorted(c[f'a_{k1}'] for c in curves)
        print(f"  K1={k1:<4} median a_at_t={vals[len(vals)//2]:>4} "
              f"range [{vals[0]}, {vals[-1]}]")
    return 0

if __name__ == '__main__':
    sys.exit(main())
