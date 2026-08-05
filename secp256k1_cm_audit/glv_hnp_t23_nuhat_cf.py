"""
GLV-HNP -- Thread 23: closed form for lambda_1(L2), and why the June continued
fraction invariants failed.

Background
----------
2026-07-29 (commit e845207) found the first predictor to separate the June
C1/C2 classes: nu_hat = lambda_1(L2)/sqrt(det L2), AUC 0.935, where

    L2 = < (n*S1, 0), (-lam*S1, S2) >,     S1 = n//K1,  S2 = max(1, n//K2)

is the non-planted 2D sublattice sitting on coordinate pair (i, m+1+i) of the
Phase 2 lattice (`glv_hnp_phase2_20bit.py:263 build_glv_lattice`).  det L2 =
n*S1*S2 is independent of lam, so nu_hat isolates the lam-dependence.

That entry closed by proposing Thread 23: nu_hat is *computed*, not
*understood*.  Show that lambda_1(L2) is a continued-fraction quantity, and that
the relevant approximation scale is set by S1/S2 -- which would retroactively
explain why the six scale-free June invariants (delta/n, kappa(M), q_cf,
max_q_cf, max_a, a_corn/n) all failed to separate C1 from C2.

The proposal asked for a >= 0.9 correlation between CF-predicted and computed
nu_hat.  The claim proved and tested here is stronger: the identity is *exact*.

Theorem (Thread 23)
-------------------
Every vector of L2 is  w(a,b) = ((a*n - b*lam)*S1, b*S2),  a,b in Z.  For fixed
b the first coordinate is minimised over a at |b*lam|_n := min_p |b*lam - p*n|,
the centred residue.  Hence

    lambda_1(L2)^2 = min_{b}  [ (|b*lam|_n * S1)^2 + (b*S2)^2 ],      (*)

the b = 0 term contributing (n*S1)^2.  A minimiser b of (*) must be a *record*:
if 1 <= b' < b had |b'*lam|_n <= |b*lam|_n then w(.,b') would be strictly
shorter.  The records of b |-> |b*lam|_n are exactly the best approximations of
the second kind to lam/n, i.e. the convergent denominators q_j of lam/n.
Therefore

    lambda_1(L2)^2 = min_j [ (|q_j*lam - p_j*n| * S1)^2 + (q_j*S2)^2 ]   (**)

over the O(log n) convergents, with (q,r) = (0,n) included for b = 0.  The
extended Euclidean algorithm on (n, lam) emits exactly the pairs
(|s_j|, r_j) = (q_j, |q_j*lam - p_j*n|).

Scale decomposition
-------------------
Write, for each convergent,

    D_j = q_j * |q_j*lam|_n / n           (scale-free; ~ 1/a_{j+1})
    R_j = q_j * S2 / (|q_j*lam|_n * S1)   (carries the scale S2/S1 = K1/K2)

Then the two terms of (**) have product D_j and ratio R_j, so

    nu_hat^2 = min_j  D_j * (R_j + 1/R_j).                            (***)

R_j = 1 at the *balanced* index, where q_j*q_{j+1} ~ n*S1/S2 = n*K2/K1, i.e.
q_j ~ sqrt(n*K2/K1).  So nu_hat is small iff the CF expansion of lam/n has a
large partial quotient *at that particular scale*.  Every June invariant
aggregated over the whole expansion (max over all j, or the corner index), which
is why they carried no signal: the minimum in (***) is dominated by one index
j*, and j* is a function of K1/K2, not of the curve.

Falsifiable predictions tested here
-----------------------------------
  A. (**) reproduces exact Lagrange-Gauss lambda_1 with zero error, at every
     bit size.  [Falsifier: any mismatch.]
  B. j* is at or adjacent to the balanced index j_bal defined by
     q_j ~ sqrt(n*S1/S2).  [Falsifier: |j*-j_bal| frequently > 1.]
  C. (***) holds numerically, and D_j ~ 1/a_{j+1} to within the Legendre bound.
  D. The *local* partial quotient a_{j*+1} separates the June C1/C2 classes,
     where the *global* max_a (falsified 2026-06-30) does not.
     [Falsifier: AUC(a_local) at the trivial baseline like max_a.]
  E. Exact-integer nu_hat for secp256k1 agrees with the float value reported
     2026-07-29 (`gauss_reduce_2d` divides in f64; at 256 bits the squared
     norms reach ~2^780 and q may round wrong).

Run: python3 glv_hnp_t23_nuhat_cf.py
Deps: sympy, fpylll (Part D only).
"""

import importlib.util
import math
import os
import random
import sys
import time

import sympy

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "_t20a", os.path.join(_HERE, "glv_hnp_t20a_lib.py"))
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

eisenstein_decompose = _t20a.eisenstein_decompose
j0_traces = _t20a.j0_traces
glv_eigenvalues = _t20a.glv_eigenvalues
mu_of = _t20a.mu_of
identify_twist = _t20a.identify_twist
find_generator = _t20a.find_generator
scales = _t20a.scales
run_experiment = _t20a.run_experiment
gauss_reduce_2d_float = _t20a.gauss_reduce_2d

SECP_P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

# ---------------------------------------------------------------------------
# Exact integer Lagrange-Gauss (reference implementation)
# ---------------------------------------------------------------------------

def _rdiv(a, b):
    """Nearest integer to a/b, exact, b > 0. Ties round up."""
    q, r = divmod(a, b)
    return q + 1 if 2 * r >= b else q

def gauss_reduce_2d_exact(u, v):
    """Lagrange-Gauss reduction in exact integer arithmetic.

    Returns (lambda_1_squared, shortest_vector).  Unlike the f64 version in
    glv_hnp_t20a_lib.py this never loses precision, so it is the ground truth
    for Part A.
    """
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    if n2(u) < n2(v):
        u, v = v, u
    while True:
        nv = n2(v)
        if nv == 0:
            return n2(u), u
        q = _rdiv(u[0] * v[0] + u[1] * v[1], nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if n2(u) <= n2(v):
            return n2(u), u

# ---------------------------------------------------------------------------
# Continued fraction frontier of lam/n
# ---------------------------------------------------------------------------

def cf_frontier(n, lam):
    """The Pareto frontier of b |-> (b, |b*lam|_n), via extended Euclid.

    Returns a list of records {q, r, a} where q = |s_j| is a convergent
    denominator of lam/n, r = |q*lam - p*n| = |q*lam|_n, and a is the partial
    quotient produced at that step (a_{j+1} in the usual indexing, i.e. the one
    that *follows* this convergent).  The list starts with (q, r) = (0, n),
    which is the b = 0 vector, and ends with r = 0, q = n.
    """
    r_prev, r = n, lam % n
    s_prev, s = 0, 1
    out = [{'q': 0, 'r': r_prev, 'a': None}]
    while r != 0:
        a = r_prev // r
        out.append({'q': abs(s), 'r': r, 'a': a})
        r_prev, r = r, r_prev - a * r
        s_prev, s = s, s_prev - a * s
    out.append({'q': abs(s), 'r': 0, 'a': None})
    return out

def nu_from_cf(n, lam, k1_bound, k2_bound):
    """lambda_1(L2)^2 and diagnostics from the CF frontier -- equation (**).

    Returns dict with lam1_sq, nu_hat, j_star, q_star, r_star, a_local (the
    partial quotient a_{j*+1}), D_star, R_star, j_bal, max_a (global).
    """
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    fr = cf_frontier(n, lam)
    det = n * s1 * s2

    best = None
    for j, e in enumerate(fr):
        # b = 0 contributes (n*S1, 0); b = q_j contributes (r_j*S1, q_j*S2)
        first = e['r'] * s1 if e['q'] != 0 else n * s1
        norm2 = first * first + (e['q'] * s2) ** 2
        if best is None or norm2 < best[0]:
            best = (norm2, j, e)
    lam1_sq, j_star, e_star = best

    # balanced index: q_j*S2 closest to r_j*S1  (equivalently R_j closest to 1)
    j_bal, bal_gap = None, None
    for j, e in enumerate(fr):
        if e['q'] == 0 or e['r'] == 0:
            continue
        gap = abs(math.log((e['q'] * s2) / (e['r'] * s1)))
        if bal_gap is None or gap < bal_gap:
            j_bal, bal_gap = j, gap

    q, r = e_star['q'], e_star['r']
    D = (q * r) / n if q and r else None
    R = (q * s2) / (r * s1) if q and r else None
    a_local = e_star['a']
    if a_local is None and j_star + 1 < len(fr):
        a_local = fr[j_star + 1]['a']
    partials = [e['a'] for e in fr if e['a'] is not None]
    return {
        'lam1_sq': lam1_sq, 'det': det,
        'nu_hat': math.sqrt(lam1_sq / det) if det else float('nan'),
        'j_star': j_star, 'q_star': q, 'r_star': r,
        'a_local': a_local, 'D_star': D, 'R_star': R,
        'j_bal': j_bal, 'n_conv': len(fr),
        'max_a': max(partials) if partials else 0,
        'sqrt_scale': math.sqrt(n * s1 / s2) if s2 else float('nan'),
    }

def nu_exact_gauss(n, lam, k1_bound, k2_bound):
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    l2, _ = gauss_reduce_2d_exact((n * s1, 0), (-lam * s1, s2))
    return l2, n * s1 * s2

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def auc(scores, labels):
    """AUC of `scores` predicting label == 1, by rank (ties averaged)."""
    pairs = sorted(zip(scores, labels))
    n1 = sum(labels)
    n0 = len(labels) - n1
    if n1 == 0 or n0 == 0:
        return float('nan')
    ranks, i = [0.0] * len(pairs), 0
    while i < len(pairs):
        j = i
        while j + 1 < len(pairs) and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1
    s1 = sum(rk for rk, (_, lb) in zip(ranks, pairs) if lb == 1)
    return (s1 - n1 * (n1 + 1) / 2.0) / (n1 * n0)

def best_cut_acc(scores, labels):
    """Accuracy of the best single threshold (either polarity)."""
    best = 0.0
    for t in sorted(set(scores)):
        for pol in (1, -1):
            acc = sum(1 for s, l in zip(scores, labels)
                      if ((s <= t) if pol == 1 else (s > t)) == (l == 1))
            best = max(best, acc / len(labels))
    return best

def rand_lambda(n, rng):
    return rng.randrange(1, n)

def k1k2_for_eff(n, eff):
    k2 = math.isqrt(n) + 1
    k1 = max(2, int(eff * n / k2))
    return k1, k2

# ---------------------------------------------------------------------------
# Part A -- exactness of the CF closed form
# ---------------------------------------------------------------------------

def part_a(rng):
    print("\n" + "=" * 78)
    print("PART A -- CF closed form (**) vs exact Lagrange-Gauss")
    print("=" * 78)
    print("  Prediction: identical, every sample, every bit size.")
    print(f"  {'bits':>5} {'samples':>8} {'CF==exact':>10} {'f64==exact':>11} "
          f"{'max f64 rel err':>16}")

    all_ok = True
    f64_bad_at = []
    for bits in (16, 20, 24, 32, 64, 96, 128, 192, 256):
        lo = 1 << (bits - 1)
        n = int(sympy.nextprime(lo + rng.randrange(lo // 2)))
        agree = f64_agree = 0
        worst = 0.0
        trials = 200 if bits <= 64 else 60
        for _ in range(trials):
            lam = rand_lambda(n, rng)
            eff = rng.choice([0.02, 0.05, 0.1, 0.25, 0.5])
            k1, k2 = k1k2_for_eff(n, eff)
            l2_exact, det = nu_exact_gauss(n, lam, k1, k2)
            cf = nu_from_cf(n, lam, k1, k2)
            if cf['lam1_sq'] == l2_exact:
                agree += 1
            try:
                f = gauss_reduce_2d_float((n * (n // k1), 0),
                                          (-lam * (n // k1), max(1, n // k2)))
                rel = abs(f * f - l2_exact) / l2_exact
            except (OverflowError, ValueError, ZeroDivisionError):
                rel = float('inf')
            if rel < 1e-12:
                f64_agree += 1
            worst = max(worst, rel if rel != float('inf') else 1e99)
        if agree != trials:
            all_ok = False
        if f64_agree != trials:
            f64_bad_at.append((bits, trials - f64_agree))
        w = "inf" if worst >= 1e99 else f"{worst:.3e}"
        print(f"  {bits:>5} {trials:>8} {agree:>7}/{trials:<3} "
              f"{f64_agree:>8}/{trials:<3} {w:>16}")

    print(f"\n  CF closed form exact everywhere: {all_ok}")
    if f64_bad_at:
        print("  f64 gauss_reduce_2d WRONG at: " +
              ", ".join(f"{b} bits ({c} of them)" for b, c in f64_bad_at))
    else:
        print("  f64 gauss_reduce_2d agreed everywhere tested")
    return all_ok

# ---------------------------------------------------------------------------
# Part B -- where the minimum sits
# ---------------------------------------------------------------------------

def part_b(rng):
    print("\n" + "=" * 78)
    print("PART B -- j* vs the balanced index; q_{j*} vs sqrt(n*S1/S2)")
    print("=" * 78)
    print("  Prediction: |j* - j_bal| <= 1 almost always; q_{j*} within a small")
    print("  factor of sqrt(n*S1/S2) = sqrt(n*K2/K1).")
    print(f"\n  {'bits':>5} {'eff':>6} {'N':>5} {'j*==j_bal':>10} "
          f"{'|d|<=1':>8} {'med q*/sqrt(nS1/S2)':>21} {'med j*':>7} {'med #conv':>10}")

    rows = []
    for bits in (20, 24, 64, 256):
        lo = 1 << (bits - 1)
        n = int(sympy.nextprime(lo + rng.randrange(lo // 2)))
        for eff in (0.05, 0.25):
            k1, k2 = k1k2_for_eff(n, eff)
            eq = near = 0
            ratios, jstars, nconv = [], [], []
            N = 400 if bits <= 64 else 150
            for _ in range(N):
                lam = rand_lambda(n, rng)
                cf = nu_from_cf(n, lam, k1, k2)
                if cf['j_bal'] is not None:
                    d = abs(cf['j_star'] - cf['j_bal'])
                    eq += (d == 0)
                    near += (d <= 1)
                if cf['q_star']:
                    ratios.append(cf['q_star'] / cf['sqrt_scale'])
                jstars.append(cf['j_star'])
                nconv.append(cf['n_conv'])
            ratios.sort(); jstars.sort(); nconv.sort()
            med = ratios[len(ratios) // 2] if ratios else float('nan')
            print(f"  {bits:>5} {eff:>6.2f} {N:>5} {eq / N:>9.1%} {near / N:>7.1%} "
                  f"{med:>21.3f} {jstars[len(jstars)//2]:>7} "
                  f"{nconv[len(nconv)//2]:>10}")
            rows.append((bits, eff, eq / N, near / N, med))
    return rows

# ---------------------------------------------------------------------------
# Part C -- the scale decomposition (***)
# ---------------------------------------------------------------------------

def part_c(rng):
    print("\n" + "=" * 78)
    print("PART C -- nu_hat^2 = min_j D_j*(R_j + 1/R_j);  D_j vs 1/a_{j+1}")
    print("=" * 78)

    bits = 24
    lo = 1 << (bits - 1)
    n = int(sympy.nextprime(lo + rng.randrange(lo // 2)))
    k1, k2 = k1k2_for_eff(n, 0.05)
    print(f"  n = {n} ({bits} bits), K1 = {k1}, K2 = {k2}, "
          f"S1/S2 = {(n//k1)/(max(1,n//k2)):.1f}")
    print(f"\n  {'lam':>10} {'nu_hat':>8} {'D*':>7} {'R*':>8} "
          f"{'sqrt(D*(R*+1/R*))':>18} {'a_local':>8} {'1/a_loc':>8} {'max_a':>7}")

    ok = True
    for _ in range(12):
        lam = rand_lambda(n, rng)
        cf = nu_from_cf(n, lam, k1, k2)
        if cf['D_star'] is None:
            continue
        pred = math.sqrt(cf['D_star'] * (cf['R_star'] + 1.0 / cf['R_star']))
        if abs(pred - cf['nu_hat']) > 1e-9 * max(1.0, cf['nu_hat']):
            ok = False
        al = cf['a_local']
        print(f"  {lam:>10} {cf['nu_hat']:>8.4f} {cf['D_star']:>7.4f} "
              f"{cf['R_star']:>8.4f} {pred:>18.4f} {al if al else 0:>8} "
              f"{1.0/al if al else 0:>8.4f} {cf['max_a']:>7}")

    print(f"\n  identity (***) holds to 1e-9 on all rows: {ok}")

    # D_j ~ 1/a_{j+1}:  Legendre gives 1/(a+2) < D < 1/a
    viol = tot = 0
    for _ in range(2000):
        lam = rand_lambda(n, rng)
        fr = cf_frontier(n, lam)
        for j in range(1, len(fr) - 1):
            e, a = fr[j], fr[j]['a']
            if not (e['q'] and e['r'] and a):
                continue
            D = e['q'] * e['r'] / n
            tot += 1
            if not (1.0 / (a + 2) <= D <= 1.0 / a + 1e-12):
                viol += 1
    print(f"  D_j in [1/(a+2), 1/a] on {tot - viol}/{tot} convergents "
          f"({(tot-viol)/tot:.4%})")
    return ok

# ---------------------------------------------------------------------------
# Part D -- does the LOCAL partial quotient separate C1/C2?
# ---------------------------------------------------------------------------

EXPS_K1 = 72
EXPS_M = 12
EXPS_SEEDS = [11, 42, 1234, 5150, 9999, 31337]

def harvest_plain(bits_lo, bits_hi, want, max_primes=200000):
    """First `want` 20-bit j=0 CM curves, no mu bucketing (20d protocol)."""
    out, seen = [], set()
    p = int(sympy.nextprime(bits_lo))
    scanned = 0
    while p < bits_hi and len(out) < want and scanned < max_primes:
        scanned += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n in seen or not sympy.isprime(n) or n % 3 != 1:
                        continue
                    roots = glv_eigenvalues(n)
                    if roots is None:
                        continue
                    b = identify_twist(p, n)
                    if b is None:
                        continue
                    G = find_generator(p, b, n)
                    if G is None:
                        continue
                    seen.add(n)
                    out.append({'p': p, 'b': b, 'n': n, 'lam': roots[0], 'G': G,
                                'mu': mu_of(roots[0], n)})
                    if len(out) >= want:
                        break
        p = int(sympy.nextprime(p))
    return out

def part_d(rng, n_curves=400):
    print("\n" + "=" * 78)
    print("PART D -- local a_{j*+1} vs global max_a vs nu_hat on the C1/C2 split")
    print("=" * 78)
    print(f"  Exp S protocol: K1={EXPS_K1}, m={EXPS_M}, {len(EXPS_SEEDS)} seeds, "
          f"{n_curves} fresh 20-bit curves.")
    print("  C2 := LLL recovers d on ALL seeds.  Prediction: a_local separates,")
    print("  max_a (falsified 2026-06-30) does not.")

    t0 = time.time()
    curves = harvest_plain(2 ** 19, 2 ** 20, n_curves)
    print(f"  harvested {len(curves)} curves in {time.time()-t0:.1f}s", flush=True)

    recs = []
    t0 = time.time()
    for idx, c in enumerate(curves):
        n, lam, p, b, G = c['n'], c['lam'], c['p'], c['b'], c['G']
        k2 = math.isqrt(n) + 1
        cf = nu_from_cf(n, lam, EXPS_K1, k2)
        wins = 0
        for seed in EXPS_SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            wins += run_experiment(p, n, lam, G, EXPS_M, d_trial,
                                   EXPS_K1, k2, seed)
        recs.append({'n': n, 'lam': lam, 'mu': c['mu'], 'wins': wins,
                     'c2': 1 if wins == len(EXPS_SEEDS) else 0,
                     'nu_hat': cf['nu_hat'], 'a_local': cf['a_local'] or 1,
                     'max_a': cf['max_a'], 'q_star': cf['q_star'],
                     'j_star': cf['j_star'], 'D_star': cf['D_star'] or 1.0})
        if (idx + 1) % 10 == 0:
            print(f"    {idx+1}/{len(curves)} curves "
                  f"({time.time()-t0:.0f}s)", flush=True)

    labels = [r['c2'] for r in recs]
    n_c2 = sum(labels)
    base = max(n_c2, len(labels) - n_c2) / len(labels)
    print(f"\n  C2 = {n_c2}/{len(recs)}   trivial baseline = {base:.1%}")

    print(f"\n  {'invariant':>12} {'C1 range':>22} {'C2 range':>22} "
          f"{'AUC':>7} {'best acc':>9}")
    summary = {}
    for name, key, flip in (("nu_hat", 'nu_hat', True),
                            ("a_local", 'a_local', False),
                            ("D_star", 'D_star', True),
                            ("max_a", 'max_a', False),
                            ("mu", 'mu', False)):
        vals = [r[key] for r in recs]
        # orient so that AUC >= 0.5 means "predicts C2"
        sc = [(-v if flip else v) for v in vals]
        a = auc(sc, labels)
        acc = best_cut_acc(vals, labels)
        c1v = [v for v, l in zip(vals, labels) if l == 0]
        c2v = [v for v, l in zip(vals, labels) if l == 1]
        rng_s = lambda xs: (f"[{min(xs):.3f}, {max(xs):.3f}]" if xs else "-")
        print(f"  {name:>12} {rng_s(c1v):>22} {rng_s(c2v):>22} "
              f"{a:>7.3f} {acc:>8.1%}")
        summary[name] = (a, acc)

    # a_local response table
    print(f"\n  a_local response ({'a_local':>8} -> C2 rate, mean wins/6):")
    buckets = [(1, 1), (2, 2), (3, 4), (5, 8), (9, 10 ** 9)]
    for lo, hi in buckets:
        grp = [r for r in recs if lo <= r['a_local'] <= hi]
        if not grp:
            continue
        lbl = f"{lo}" if lo == hi else (f">={lo}" if hi > 10 ** 8 else f"{lo}-{hi}")
        print(f"    a_local {lbl:>6}  N={len(grp):>3}  "
              f"C2 rate={sum(r['c2'] for r in grp)/len(grp):>5.2f}  "
              f"mean wins={sum(r['wins'] for r in grp)/len(grp):>4.2f}")

    # --- significance: bootstrap CIs, and paired bootstrap on AUC gaps -----
    print("\n  Bootstrap (2000 resamples of curves) and label-permutation tests:")
    print(f"  {'invariant':>12} {'AUC':>7} {'95% CI':>16} {'perm p':>9}")
    ORIENT = {'nu_hat': -1, 'a_local': 1, 'D_star': -1, 'max_a': 1, 'mu': 1}
    sc = {k: [ORIENT[k] * r[k] for r in recs] for k in ORIENT}
    labels_l = labels
    B = 2000
    brng = random.Random(4242)
    idx_all = range(len(recs))
    boots = {k: [] for k in ORIENT}
    for _ in range(B):
        idx = [brng.randrange(len(recs)) for _ in idx_all]
        lb = [labels_l[i] for i in idx]
        if sum(lb) in (0, len(lb)):
            continue
        for k in ORIENT:
            boots[k].append(auc([sc[k][i] for i in idx], lb))
    for k in ORIENT:
        bs = sorted(boots[k])
        lo, hi = bs[int(0.025 * len(bs))], bs[int(0.975 * len(bs))]
        obs = auc(sc[k], labels_l)
        # permutation test: AUC under shuffled labels
        prng = random.Random(99)
        ge = 0
        for _ in range(B):
            sh = labels_l[:]
            prng.shuffle(sh)
            if auc(sc[k], sh) >= obs:
                ge += 1
        print(f"  {k:>12} {obs:>7.3f} [{lo:.3f}, {hi:.3f}] {(ge+1)/(B+1):>9.4f}")

    print("\n  Paired bootstrap on AUC gaps (does the theory's ordering hold?):")
    for a_k, b_k in (('nu_hat', 'a_local'), ('a_local', 'max_a'),
                     ('nu_hat', 'max_a'), ('a_local', 'mu')):
        d = [x - y for x, y in zip(boots[a_k], boots[b_k])]
        d.sort()
        lo, hi = d[int(0.025 * len(d))], d[int(0.975 * len(d))]
        frac = sum(1 for x in d if x > 0) / len(d)
        print(f"    AUC({a_k}) - AUC({b_k}) = "
              f"{auc(sc[a_k], labels_l) - auc(sc[b_k], labels_l):+.3f}  "
              f"95% CI [{lo:+.3f}, {hi:+.3f}]  P(gap>0) = {frac:.3f}")
    return summary, base, recs

# ---------------------------------------------------------------------------
# Part E -- secp256k1, exact
# ---------------------------------------------------------------------------

def part_e():
    print("\n" + "=" * 78)
    print("PART E -- secp256k1 exact-integer nu_hat (checks the f64 numbers)")
    print("=" * 78)
    n = SECP_N
    roots = glv_eigenvalues(n)
    assert roots is not None
    lam = roots[0]
    assert (lam * lam + lam + 1) % n == 0
    print(f"  lam = {hex(lam)}")
    print(f"  mu  = {mu_of(lam, n):.6f}")
    print(f"\n  {'eff':>6} {'nu_hat(exact)':>14} {'nu_hat(f64)':>13} "
          f"{'rel err':>10} {'j*':>4} {'a_local':>8} {'max_a':>7} "
          f"{'q*/sqrt(nS1/S2)':>16}")
    out = []
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1, k2 = k1k2_for_eff(n, eff)
        l2, det = nu_exact_gauss(n, lam, k1, k2)
        nu_ex = math.sqrt(l2 / det)
        cf = nu_from_cf(n, lam, k1, k2)
        assert cf['lam1_sq'] == l2, "CF disagrees with exact Gauss on secp256k1"
        try:
            f = gauss_reduce_2d_float((n * (n // k1), 0),
                                      (-lam * (n // k1), max(1, n // k2)))
            nu_f = f / math.sqrt(det)
            rel = abs(nu_f - nu_ex) / nu_ex
        except (OverflowError, ValueError):
            nu_f, rel = float('nan'), float('nan')
        print(f"  {eff:>6.2f} {nu_ex:>14.6f} {nu_f:>13.6f} {rel:>10.2e} "
              f"{cf['j_star']:>4} {cf['a_local'] or 0:>8} {cf['max_a']:>7} "
              f"{cf['q_star']/cf['sqrt_scale']:>16.4f}")
        out.append((eff, nu_ex, nu_f, cf))
    print("\n  Caveat (unchanged from 2026-07-29): this is conditional on a")
    print("  non-standard nonce generator k = k1 + lam*k2 with k1 bounded, and")
    print("  is a 20/24-bit -> 256-bit extrapolation.  NOT an attack on secp256k1.")
    return out

# ---------------------------------------------------------------------------

def main():
    do_d = "--no-partd" not in sys.argv
    rng = random.Random(20260805)
    print("=" * 78)
    print("GLV-HNP Thread 23 -- closed form for lambda_1(L2)")
    print("=" * 78)

    ok_a = part_a(rng)
    part_b(rng)
    ok_c = part_c(rng)
    if do_d:
        part_d(rng)
    part_e()

    print("\n" + "=" * 78)
    print(f"VERDICT: CF closed form exact = {ok_a}; decomposition (***) = {ok_c}")
    print("=" * 78)

if __name__ == "__main__":
    main()
