"""
GLV-HNP — Thread 20d: does nu_hat separate the June C1/C2 classes?

Why this matters
----------------
`glv_hnp_delta_threshold.py:224 build_lattice()` and the Phase 2 builder in
`glv_hnp_phase2_lambda_threshold.py` are the SAME (2m+2)-dimensional
column-scaled lattice -- identical rows, identical S_K1/S_K2/S_KANNAN, identical
recovery test.  The only difference is that the June runs fixed K1=72 in
absolute terms while Phase 2 sets K1 from a target eff = K1*K2/n.

So the June "C1/C2" classification (C2 = LLL recovers d on all seeds at K1=72,
m=12; C1 = it does not) is a classification of *the same lattice*, and the
Thread 20c predictor nu_hat can be tested directly against it.

That matters because the 2026-06-30 log declared Thread 5 a DEAD END after six
consecutive invariants failed to separate C1 from C2:

    delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n

and the 2026-06-29 entry closed with the (never executed) conjecture that the
separator "requires a lattice-geometric computation: specifically, whether the
BV lattice for (n, lam) admits a short non-planted vector".  nu_hat is exactly
that computation.

This script replicates the Exp S protocol (K1=72, m=12, 6 seeds, fresh 20-bit
j=0 CM curves) and puts nu_hat head to head with max_a on the same sample.

Run: python3 glv_hnp_nuhat_vs_c1c2.py
"""

import math
import random
import sys
import time

import sympy

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_thread20c_lib.py")
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

eisenstein_decompose = _t20a.eisenstein_decompose
j0_traces = _t20a.j0_traces
glv_eigenvalues = _t20a.glv_eigenvalues
mu_of = _t20a.mu_of
identify_twist = _t20a.identify_twist
find_generator = _t20a.find_generator
rival_sublattice_nu = _t20a.rival_sublattice_nu
run_experiment = _t20a.run_experiment

# Exp S protocol (2026-06-29)
K1_FIXED = 72
M_FIXED = 12
SEEDS_N = 6
N_CURVES = 100
BITS_LO, BITS_HI = 2 ** 19, 2 ** 20

def max_partial_quotient(a, b):
    """Largest partial quotient in the continued fraction of a/b (the
    'max_a' invariant falsified by Exp S, 2026-06-29)."""
    best = 0
    while b:
        q, r = divmod(a, b)
        best = max(best, q)
        a, b = b, r
    return best

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
    """Probability a random positive outranks a random negative (ties=0.5)."""
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
    print("GLV-HNP Thread 20d: nu_hat vs the June C1/C2 classes")
    print("=" * 78)
    print(f"Exp S protocol replication: K1={K1_FIXED}, m={M_FIXED}, "
          f"{SEEDS_N} seeds, {N_CURVES} fresh 20-bit j=0 CM curves")
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
        c['mu'] = mu_of(c['lam'], n)
        c['max_a'] = max_partial_quotient(c['lam'], n)
    print(f"  {len(curves)*SEEDS_N} LLL trials in {time.time()-t1:.1f}s")

    c2 = [c for c in curves if c['C2']]
    c1 = [c for c in curves if not c['C2']]
    base = max(len(c1), len(c2)) / len(curves)
    print(f"\n  C1 (fails): {len(c1)}/{len(curves)}   C2 (6/6): {len(c2)}/{len(curves)}")
    print(f"  trivial baseline (always predict majority) = {base*100:.1f}%")
    print("  (Exp S 2026-06-29 saw 40/50 C1, 10/50 C2 -- same regime)")

    print("\n--- Head to head: nu_hat vs the falsified invariants ---")
    print(f"{'invariant':>10} {'C1 range':>18} {'C2 range':>18} {'AUC':>7} {'best acc':>9}")
    for key in ('nu_hat', 'mu', 'max_a'):
        s1 = [c[key] for c in c1]
        s2 = [c[key] for c in c2]
        # AUC oriented so that "higher score => C2" ; report max(auc, 1-auc)
        a = auc([c[key] for c in curves], [c['C2'] for c in curves])
        a_or = max(a, 1 - a)
        acc, cut = best_cut_accuracy([c[key] for c in curves],
                                     [c['C2'] for c in curves])
        r1 = f"[{min(s1):.3f},{max(s1):.3f}]" if s1 else "-"
        r2 = f"[{min(s2):.3f},{max(s2):.3f}]" if s2 else "-"
        print(f"{key:>10} {r1:>18} {r2:>18} {a_or:>7.3f} {acc*100:>8.1f}%")

    print("\n--- nu_hat decile response (C2 rate and mean wins/6) ---")
    order = sorted(curves, key=lambda c: c['nu_hat'])
    k = max(1, len(order) // 10)
    print(f"{'decile':>7} {'nu_hat mid':>11} {'C2 rate':>9} {'mean wins/6':>12}")
    for i in range(0, len(order), k):
        grp = order[i:i + k]
        if len(grp) < 2:
            continue
        print(f"{i//k+1:>7} {sum(c['nu_hat'] for c in grp)/len(grp):>11.3f} "
              f"{sum(c['C2'] for c in grp)/len(grp):>9.2f} "
              f"{sum(c['wins'] for c in grp)/len(grp):>12.2f}")

    print("\n--- secp256k1 placement at the same eff ---")
    n_s = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    lam_s = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
    assert (lam_s * lam_s + lam_s + 1) % n_s == 0
    k2_s = math.isqrt(n_s) + 1
    eff_here = K1_FIXED * (math.isqrt(curves[0]['n']) + 1) / curves[0]['n']
    k1_s = max(2, int(eff_here * math.isqrt(n_s)))
    _, _, nh_s = rival_sublattice_nu(n_s, lam_s, k1_s, k2_s)
    below = sum(1 for c in curves if c['nu_hat'] < nh_s) / len(curves)
    print(f"  sample eff = {eff_here:.4f}; secp256k1 at that eff: "
          f"nu_hat={nh_s:.4f} (percentile {below*100:.0f}% of sample)")

    print("\nDone.")

if __name__ == '__main__':
    sys.exit(main())
