"""
GLV-HNP Phase 2 — Thread 20b: mu-response curve, well powered.

Motivation
----------
`glv_hnp_phase2_lambda_threshold.py` (Thread 20a) falsified the 2026-07-26
conjecture that small mu = min(lam, n-lam)/n obstructs the Phase 2 attack:
at 20 bits with eff = K1*K2/n = 0.05, 27/28 curves reach 3/3, and the
smallest mu in the sample (mu = 0.0200) was among the *fastest* (3/3 at m=7).

But the graded response -- the smallest m reaching 3/3 -- was not flat.  Binned
by eye on the 20a sample:

    mu < 0.04            first-3/3 m ~ 7.5   (2 curves)
    0.04 <= mu <= 0.19   first-3/3 m ~ 11.7  (10 curves)
    mu > 0.20            first-3/3 m ~ 8.6   (15 curves)

i.e. a non-monotone *hard band* at intermediate mu, not an obstruction at
mu -> 0.  With 3 seeds per curve that is far too noisy to believe.

This script tests it properly.

Design
------
- Harvest CURVES_PER_BIN curves in each of BINS uniform mu-bins over (0, 0.5].
- Unit of analysis is the CURVE: for each curve estimate the success
  probability p_hat at fixed m in M_FIXED over SEEDS_N independent seeds.
  Fixed m avoids the right-censoring that biases "first 3/3 m".
- Compare bands with a curve-level permutation test (curves are exchangeable
  under H0; seeds within a curve are not independent observations of a
  band-level effect).
- Also correlate p_hat against nu_hat (the lattice-geometric invariant
  falsified as a binary separator in 20a) to see if it explains the graded
  response instead.

H-band  : p_hat at fixed m is lower for mu in [0.04, 0.20] than outside it.
Falsifier: permutation p-value > 0.05, or the effect vanishes at a second m.

Run: python3 glv_hnp_phase2_mu_response.py
"""

import math
import random
import sys
import time

import sympy

# Reuse the verified Thread 20a core rather than re-implementing it.
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_core.py")
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

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BITS_LO, BITS_HI = 2 ** 19, 2 ** 20
EFF_TARGET = 0.05
BINS = 20                       # uniform bins over (0, 0.5]
CURVES_PER_BIN = 5
SEEDS_N = 24                    # seeds per (curve, m)
M_FIXED = [8, 9, 10, 11]
MAX_PRIMES_SCANNED = 40000
PERM_ITERS = 20000
BAND_LO, BAND_HI = 0.04, 0.20   # the conjectured hard band

# ---------------------------------------------------------------------------

def harvest(n_bins, per_bin):
    edges = [(i / n_bins * 0.5, (i + 1) / n_bins * 0.5) for i in range(n_bins)]
    buckets = [[] for _ in edges]
    seen_n = set()
    scanned = 0
    p = int(sympy.nextprime(BITS_LO))
    t0 = time.time()
    while p < BITS_HI and scanned < MAX_PRIMES_SCANNED:
        scanned += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n in seen_n or n % 3 != 1 or not sympy.isprime(n):
                        continue
                    roots = glv_eigenvalues(n)
                    if roots is None:
                        continue
                    lam = roots[0]
                    mu = mu_of(lam, n)
                    for bi, (lo, hi) in enumerate(edges):
                        if lo < mu <= hi and len(buckets[bi]) < per_bin:
                            b = identify_twist(p, n)
                            if b is None:
                                break
                            G = find_generator(p, b, n)
                            if G is None:
                                break
                            seen_n.add(n)
                            buckets[bi].append({'p': p, 'b': b, 'n': n,
                                                'lam': lam, 'mu': mu, 'bin': bi, 'G': G})
                            break
        if all(len(v) >= per_bin for v in buckets):
            break
        p = int(sympy.nextprime(p))
    total = sum(len(v) for v in buckets)
    print(f"  harvested {total} curves from {scanned} primes in {time.time()-t0:.1f}s")
    for bi, (lo, hi) in enumerate(edges):
        if len(buckets[bi]) < per_bin:
            print(f"    bin {lo:.3f}-{hi:.3f}: only {len(buckets[bi])}/{per_bin}")
    return [c for v in buckets for c in v], edges

def measure(curve, m, seeds_n):
    n = curve['n']
    k2 = math.isqrt(n) + 1
    k1 = max(2, int(EFF_TARGET * math.sqrt(n)))
    wins = 0
    for s in range(seeds_n):
        seed = 1000 * m + s
        d = random.Random(seed + 7777).randint(1, n - 1)
        wins += run_experiment(curve['p'], n, curve['lam'], curve['G'],
                               m, d, k1, k2, seed)
    return wins / seeds_n

def perm_test(values, labels, iters, rng):
    """Two-sided permutation test on difference of means between label groups."""
    a = [v for v, l in zip(values, labels) if l]
    b = [v for v, l in zip(values, labels) if not l]
    if not a or not b:
        return float('nan'), float('nan')
    obs = sum(a) / len(a) - sum(b) / len(b)
    pool = list(values)
    na = len(a)
    hits = 0
    for _ in range(iters):
        rng.shuffle(pool)
        d = sum(pool[:na]) / na - sum(pool[na:]) / (len(pool) - na)
        if abs(d) >= abs(obs) - 1e-15:
            hits += 1
    return obs, (hits + 1) / (iters + 1)

def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float('nan')

def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))

def main():
    print("=" * 78)
    print("GLV-HNP Phase 2 — Thread 20b: mu-response curve (well powered)")
    print("=" * 78)
    print(f"20-bit j=0 CM curves, eff={EFF_TARGET}, {SEEDS_N} seeds/(curve,m), "
          f"m in {M_FIXED}, {BINS} bins x {CURVES_PER_BIN} curves")

    print("\n--- Part A: harvest ---")
    curves, edges = harvest(BINS, CURVES_PER_BIN)

    print("\n--- Part B: measure p_hat at fixed m ---")
    t0 = time.time()
    for c in curves:
        n, k2 = c['n'], math.isqrt(c['n']) + 1
        k1 = max(2, int(EFF_TARGET * math.sqrt(n)))
        _, _, c['nu_hat'] = rival_sublattice_nu(n, c['lam'], k1, k2)
        c['p_hat'] = {m: measure(c, m, SEEDS_N) for m in M_FIXED}
    trials = len(curves) * len(M_FIXED) * SEEDS_N
    print(f"  {trials} LLL trials in {time.time()-t0:.1f}s")

    print("\n--- Part C: mu-bin response ---")
    print(f"{'mu bin':>14} {'n':>3} " + " ".join(f"{'m='+str(m):>7}" for m in M_FIXED))
    for bi, (lo, hi) in enumerate(edges):
        grp = [c for c in curves if c['bin'] == bi]
        if not grp:
            continue
        cells = " ".join(f"{sum(c['p_hat'][m] for c in grp)/len(grp):>7.3f}"
                         for m in M_FIXED)
        print(f"{lo:.3f}-{hi:.3f} {len(grp):>3} {cells}")

    print("\n--- Part D: H-band permutation test ---")
    print(f"  band = [{BAND_LO}, {BAND_HI}] (conjectured HARD); "
          f"comparison = mu outside the band")
    rng = random.Random(20260729)
    in_band = [BAND_LO <= c['mu'] <= BAND_HI for c in curves]
    print(f"  {sum(in_band)} curves in band, {len(curves)-sum(in_band)} outside")
    for m in M_FIXED:
        vals = [c['p_hat'][m] for c in curves]
        obs, pv = perm_test(vals, in_band, PERM_ITERS, rng)
        inb = [v for v, l in zip(vals, in_band) if l]
        out = [v for v, l in zip(vals, in_band) if not l]
        print(f"  m={m:<3} p_hat in-band={sum(inb)/len(inb):.3f}  "
              f"out-band={sum(out)/len(out):.3f}  diff={obs:+.3f}  perm p={pv:.4f}")

    print("\n--- Part E: monotone correlations (mu, nu_hat) vs p_hat ---")
    for m in M_FIXED:
        vals = [c['p_hat'][m] for c in curves]
        print(f"  m={m:<3} spearman(mu, p_hat)={spearman([c['mu'] for c in curves], vals):+.3f}"
              f"   spearman(nu_hat, p_hat)="
              f"{spearman([c['nu_hat'] for c in curves], vals):+.3f}")

    print("\n--- Part F: per-curve variance check ---")
    m0 = M_FIXED[len(M_FIXED) // 2]
    vals = sorted(c['p_hat'][m0] for c in curves)
    mean = sum(vals) / len(vals)
    var_obs = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
    var_binom = mean * (1 - mean) / SEEDS_N
    print(f"  at m={m0}: mean p_hat={mean:.3f}")
    print(f"    observed between-curve variance = {var_obs:.5f}")
    print(f"    variance if all curves identical = {var_binom:.5f}  "
          f"(binomial, {SEEDS_N} seeds)")
    print(f"    overdispersion ratio = {var_obs/var_binom:.2f}"
          f"  ({'real per-curve variation' if var_obs/var_binom > 1.5 else 'consistent with pure seed noise'})")
    print(f"    p_hat spread: min={vals[0]:.3f} median={vals[len(vals)//2]:.3f} "
          f"max={vals[-1]:.3f}")

    print("\nDone.")

if __name__ == '__main__':
    sys.exit(main())
