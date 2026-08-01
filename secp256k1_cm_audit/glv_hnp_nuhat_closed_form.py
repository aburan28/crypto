"""
GLV-HNP -- Thread 23: closed form for nu_hat, and why six scale-free
continued-fraction invariants failed in June.

Background
----------
2026-07-29 (commit e845207) found that

    nu_hat = lambda_1(L2) / sqrt(det L2),
    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >

separates the June C1/C2 classes with AUC 0.935, where six earlier invariants
(delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n) and later lam/n and mu all
sat at the trivial baseline.  That entry closed with:

    "nu_hat is currently computed, not understood.  ... minimising
     |a*lam mod n|*S_K1 against a*S_K2 is best rational approximation to lam/n
     *at the scale* S_K2/S_K1.  That is why the raw CF invariants failed in
     June -- they are scale-free."

    Falsifier: if predicted-from-CF nu_hat correlates <0.9 with computed
    nu_hat, the convergent-scale story is wrong.

This script executes that test, and goes one step further: if the story holds,
the *scale-local* partial quotient a_{j*+1} (j* = the convergent index nearest
the critical scale) should be a genuine algebraic invariant with the predictive
power that scale-free max_a lacks.

Parts
-----
A  Theorem check: the minimiser of the L2 objective is a CF convergent
   denominator of lam/n.  Exhaustive brute force vs CF-restricted search.
B  Closed form: nu_hat_cf (CF convergents only) vs nu_hat (Lagrange-Gauss).
C  Critical index: is the winning j the one with q_j nearest sqrt(n*S_K1/S_K2)?
D  Scale-local a_{j*+1} vs C1/C2, head to head with scale-free max_a.
E  secp256k1 decomposition.

Run: python3 glv_hnp_nuhat_closed_form.py
"""

import importlib.util
import math
import random
import sys
import time

import sympy

_spec = importlib.util.spec_from_file_location(
    "_core", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_core.py")
_core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_core)

scales = _core.scales
rival_sublattice_nu = _core.rival_sublattice_nu
gauss_reduce_2d = _core.gauss_reduce_2d
eisenstein_decompose = _core.eisenstein_decompose
j0_traces = _core.j0_traces
glv_eigenvalues = _core.glv_eigenvalues
mu_of = _core.mu_of
identify_twist = _core.identify_twist
find_generator = _core.find_generator
run_experiment = _core.run_experiment

# Exp S protocol, identical to glv_hnp_nuhat_vs_c1c2.py (deterministic sample)
K1_FIXED = 72
M_FIXED = 12
SEEDS_N = 6
N_CURVES = 100
BITS_LO, BITS_HI = 2 ** 19, 2 ** 20


# ---------------------------------------------------------------------------
# Continued fractions
# ---------------------------------------------------------------------------

def cf_convergents(a, b):
    """Continued fraction of a/b.  Returns (quotients, convergents) where
    convergents is the list of (p_j, q_j) with p_j/q_j -> a/b, j = 0, 1, ..."""
    quots, convs = [], []
    pm1, qm1, pm2, qm2 = 1, 0, 0, 1
    x, y = a, b
    while y:
        q, r = divmod(x, y)
        quots.append(q)
        p_j = q * pm1 + pm2
        q_j = q * qm1 + qm2
        convs.append((p_j, q_j))
        pm2, qm2, pm1, qm1 = pm1, qm1, p_j, q_j
        x, y = y, r
    return quots, convs


def centered_residue(b, lam, n):
    """|b*lam mod (+/-) n| -- the distance from b*lam to the nearest multiple
    of n.  Equals min over a of |a*n - b*lam|."""
    r = (b * lam) % n
    return min(r, n - r)


def objective(b, lam, n, s_k1, s_k2):
    """Squared norm of the shortest L2 vector with second coordinate b*S_K2."""
    r = centered_residue(b, lam, n)
    return s_k1 * s_k1 * r * r + s_k2 * s_k2 * b * b


def nu_hat_from_cf(n, lam, k1_bound, k2_bound):
    """
    Closed-form nu_hat: minimise the L2 objective over CF convergent
    denominators of lam/n only (plus the b=0 vector (n*S_K1, 0)).

    Returns (nu_hat_cf, j_win, q_win, ncand) -- j_win = -1 means b=0 won.
    """
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    best = (n * s_k1) ** 2          # b = 0 branch
    j_win, q_win = -1, 0
    _, convs = cf_convergents(lam, n)
    for j, (_, q_j) in enumerate(convs):
        if q_j == 0:
            continue
        if s_k2 * q_j >= n * s_k1:   # already longer than the b=0 vector
            break
        val = objective(q_j, lam, n, s_k1, s_k2)
        if val < best:
            best, j_win, q_win = val, j, q_j
    det = n * s_k1 * s_k2
    return math.sqrt(best / det), j_win, q_win, len(convs)


def critical_index(n, lam, k1_bound, k2_bound):
    """
    j* = argmin_j |log(q_j / sqrt(n * S_K1 / S_K2))|, the convergent nearest
    the scale at which S_K1*delta and S_K2*q balance, together with the
    scale-local partial quotient a_{j*+1}.
    """
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    target = math.sqrt(n * s_k1 / s_k2)
    quots, convs = cf_convergents(lam, n)
    best_j, best_d = 0, float('inf')
    for j, (_, q_j) in enumerate(convs):
        if q_j <= 0:
            continue
        d = abs(math.log(q_j / target))
        if d < best_d:
            best_j, best_d = j, d
    # quots[j+1] is the partial quotient a_{j+1} that governs delta_j
    a_next = quots[best_j + 1] if best_j + 1 < len(quots) else quots[-1]
    return best_j, convs[best_j][1], a_next, target


def max_partial_quotient(a, b):
    """Scale-free max_a -- falsified by Exp S, 2026-06-29.  Kept for the
    head-to-head."""
    best = 0
    while b:
        q, r = divmod(a, b)
        best = max(best, q)
        a, b = b, r
    return best


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

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
    best = 0.0
    for c in sorted(set(scores)):
        for sign in (1, -1):
            acc = sum(1 for s, l in zip(scores, labels)
                      if ((s > c) if sign > 0 else (s <= c)) == l) / len(labels)
            best = max(best, acc)
    return best


def pearson(xs, ys):
    k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return float('nan')
    return sxy / math.sqrt(sxx * syy)


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


# ---------------------------------------------------------------------------
# Curve harvesting -- identical walk to glv_hnp_nuhat_vs_c1c2.py
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Part A -- the minimiser is a convergent denominator
# ---------------------------------------------------------------------------

def part_a():
    print("=" * 78)
    print("PART A -- is the L2 minimiser always a CF convergent denominator?")
    print("=" * 78)
    print("For each (n, lam, K1): brute force b = 0..B vs CF-restricted search.")
    print("B is the point where S_K2*b exceeds the b=0 vector n*S_K1.\n")

    rng = random.Random(20260801)
    primes = [int(sympy.nextprime(rng.randint(2000, 40000))) for _ in range(40)]
    trials, mismatch, notconv = 0, 0, 0
    worst = []
    for n in primes:
        for k1 in (8, 72, 512):
            k2 = math.isqrt(n) + 1
            s_k1, _, s_k2, _ = scales(n, k1, k2)
            if s_k1 == 0 or s_k2 == 0:
                continue
            for _ in range(3):
                lam = rng.randrange(1, n)
                bmax = min(n - 1, int(n * s_k1 / s_k2) + 2)
                if bmax > 4_000_000:
                    continue
                best_val, best_b = (n * s_k1) ** 2, 0
                for b in range(1, bmax + 1):
                    v = objective(b, lam, n, s_k1, s_k2)
                    if v < best_val:
                        best_val, best_b = v, b
                nu_bf = math.sqrt(best_val / (n * s_k1 * s_k2))
                nu_cf, _, q_win, _ = nu_hat_from_cf(n, lam, k1, k2)
                _, convs = cf_convergents(lam, n)
                qs = {q for _, q in convs}
                trials += 1
                if abs(nu_bf - nu_cf) > 1e-12 * max(1.0, nu_bf):
                    mismatch += 1
                    worst.append((n, lam, k1, nu_bf, nu_cf))
                if best_b != 0 and best_b not in qs:
                    notconv += 1
    print(f"  trials                                   : {trials}")
    print(f"  brute-force minimiser NOT a convergent   : {notconv}")
    print(f"  nu_hat_cf != nu_hat_bruteforce           : {mismatch}")
    if worst:
        for w in worst[:5]:
            print(f"    n={w[0]} lam={w[1]} K1={w[2]} bf={w[3]:.6f} cf={w[4]:.6f}")
    print("  VERDICT:", "CONFIRMED" if (mismatch == 0 and notconv == 0)
          else "REFUTED")
    return mismatch == 0 and notconv == 0


# ---------------------------------------------------------------------------
# Part B/C/D -- the 20d sample
# ---------------------------------------------------------------------------

def part_bcd():
    print()
    print("=" * 78)
    print("PARTS B/C/D -- closed form and scale-local a_{j*+1} on the 20d sample")
    print("=" * 78)
    t0 = time.time()
    curves = harvest(N_CURVES)
    print(f"  harvested {len(curves)} curves in {time.time()-t0:.1f}s "
          f"(deterministic walk from nextprime(2^19))")

    # --- B: closed form vs Lagrange-Gauss, no LLL needed ---
    for c in curves:
        n, lam = c['n'], c['lam']
        k2 = math.isqrt(n) + 1
        _, _, c['nu_hat'] = rival_sublattice_nu(n, lam, K1_FIXED, k2)
        c['nu_cf'], c['j_win'], c['q_win'], c['ncf'] = \
            nu_hat_from_cf(n, lam, K1_FIXED, k2)
        c['j_star'], c['q_star'], c['a_star'], c['target'] = \
            critical_index(n, lam, K1_FIXED, k2)
        c['mu'] = mu_of(lam, n)
        c['max_a'] = max_partial_quotient(lam, n)

    exact = sum(1 for c in curves
                if abs(c['nu_cf'] - c['nu_hat']) <= 1e-9 * max(1.0, c['nu_hat']))
    r_p = pearson([c['nu_cf'] for c in curves], [c['nu_hat'] for c in curves])
    r_s = spearman([c['nu_cf'] for c in curves], [c['nu_hat'] for c in curves])
    maxdev = max(abs(c['nu_cf'] - c['nu_hat']) for c in curves)
    print("\n--- B: nu_hat_cf (convergents only) vs nu_hat (Lagrange-Gauss) ---")
    print(f"  exact agreement (rel 1e-9)   : {exact}/{len(curves)}")
    print(f"  max |nu_cf - nu_hat|         : {maxdev:.3e}")
    print(f"  pearson r                    : {r_p:.6f}")
    print(f"  spearman rho                 : {r_s:.6f}")
    print(f"  falsifier was r < 0.9        : "
          f"{'PASSED (story holds)' if r_p >= 0.9 else 'FAILED'}")
    print(f"  mean CF length               : "
          f"{sum(c['ncf'] for c in curves)/len(curves):.1f} quotients; "
          f"winning index j mean = {sum(c['j_win'] for c in curves)/len(curves):.1f}")

    # --- C: is the winner the convergent nearest the critical scale? ---
    off = [c['j_win'] - c['j_star'] for c in curves]
    hist = {}
    for o in off:
        hist[o] = hist.get(o, 0) + 1
    within1 = sum(1 for o in off if abs(o) <= 1)
    ratio = [c['q_win'] / c['target'] for c in curves if c['q_win'] > 0]
    print("\n--- C: winning convergent vs the critical scale sqrt(n*S_K1/S_K2) ---")
    print(f"  j_win - j_star histogram     : "
          f"{dict(sorted(hist.items()))}")
    print(f"  |j_win - j_star| <= 1        : {within1}/{len(curves)}")
    if ratio:
        ratio.sort()
        print(f"  q_win / sqrt(n*S_K1/S_K2)    : "
              f"min {ratio[0]:.3f}  q25 {ratio[len(ratio)//4]:.3f}  "
              f"med {ratio[len(ratio)//2]:.3f}  "
              f"q75 {ratio[3*len(ratio)//4]:.3f}  max {ratio[-1]:.3f}")

    # --- D: LLL labels, then head-to-head ---
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
    print(f"\n  {len(curves)*SEEDS_N} LLL trials in {time.time()-t1:.1f}s")

    n_c2 = sum(1 for c in curves if c['C2'])
    base = max(n_c2, len(curves) - n_c2) / len(curves)
    print(f"  C1 (fails): {len(curves)-n_c2}/{len(curves)}   "
          f"C2 (6/6): {n_c2}/{len(curves)}   baseline = {base*100:.1f}%")

    print("\n--- D: head to head against the C1/C2 labels ---")
    print(f"{'invariant':>22} {'kind':>12} {'AUC':>7} {'best acc':>9}")
    labels = [c['C2'] for c in curves]
    rows = [
        ('nu_hat (Lagrange-Gauss)', 'lattice', [c['nu_hat'] for c in curves]),
        ('nu_hat_cf (closed form)', 'lattice', [c['nu_cf'] for c in curves]),
        ('a_{j*+1} scale-local', 'algebraic', [float(c['a_star']) for c in curves]),
        ('log a_{j*+1}', 'algebraic',
         [math.log(max(1, c['a_star'])) for c in curves]),
        ('max_a scale-free', 'algebraic', [float(c['max_a']) for c in curves]),
        ('mu', 'algebraic', [c['mu'] for c in curves]),
    ]
    for name, kind, sc in rows:
        a = auc(sc, labels)
        print(f"{name:>22} {kind:>12} {max(a, 1-a):>7.3f} "
              f"{best_cut_accuracy(sc, labels)*100:>8.1f}%")
    print(f"{'(trivial baseline)':>22} {'':>12} {'':>7} {base*100:>8.1f}%")

    print("\n--- D2: a_{j*+1} response (C2 rate by scale-local quotient) ---")
    buckets = [(1, 1), (2, 2), (3, 4), (5, 9), (10, 10 ** 9)]
    print(f"{'a_{j*+1}':>12} {'curves':>7} {'C2 rate':>9} {'mean wins/6':>12} "
          f"{'mean nu_hat':>12}")
    for lo, hi in buckets:
        grp = [c for c in curves if lo <= c['a_star'] <= hi]
        if not grp:
            continue
        tag = f"{lo}" if lo == hi else (f">={lo}" if hi > 10 ** 8 else f"{lo}-{hi}")
        print(f"{tag:>12} {len(grp):>7} "
              f"{sum(c['C2'] for c in grp)/len(grp):>9.2f} "
              f"{sum(c['wins'] for c in grp)/len(grp):>12.2f} "
              f"{sum(c['nu_hat'] for c in grp)/len(grp):>12.3f}")
    return curves


# ---------------------------------------------------------------------------
# Part E -- secp256k1
# ---------------------------------------------------------------------------

def part_e(curves):
    print()
    print("=" * 78)
    print("PART E -- secp256k1 CF decomposition of nu_hat")
    print("=" * 78)
    n_s = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    lam_s = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
    assert (lam_s * lam_s + lam_s + 1) % n_s == 0
    k2_s = math.isqrt(n_s) + 1
    quots, _ = cf_convergents(lam_s, n_s)
    print(f"  CF(lam/n) has {len(quots)} partial quotients; "
          f"max_a = {max(quots)} (scale-free, falsified)")
    print(f"\n{'eff':>7} {'K1':>22} {'nu_hat':>9} {'nu_cf':>9} {'j*':>4} "
          f"{'a_{j*+1}':>9} {'q*/scale':>9}")
    for eff in (0.05, 0.0993, 0.10, 0.25, 0.50):
        k1_s = max(2, int(eff * math.isqrt(n_s)))
        _, _, nh = rival_sublattice_nu(n_s, lam_s, k1_s, k2_s)
        nc, j_win, q_win, _ = nu_hat_from_cf(n_s, lam_s, k1_s, k2_s)
        j_star, q_star, a_star, target = critical_index(n_s, lam_s, k1_s, k2_s)
        print(f"{eff:>7.4f} {k1_s:>22} {nh:>9.4f} {nc:>9.4f} {j_star:>4} "
              f"{a_star:>9} {q_win/target if q_win else 0:>9.3f}")

    print("\n  Reference: every C2 curve in the 20d sample has nu_hat <= "
          f"{max(c['nu_hat'] for c in curves if c['C2']):.3f}")
    print("  (heuristic 20-bit -> 256-bit extrapolation; conditional on a "
          "non-standard\n   k = k1 + lam*k2 nonce generator. NOT an attack on "
          "secp256k1.)")


def main():
    ok = part_a()
    curves = part_bcd()
    part_e(curves)
    print("\nDone.  Part A verdict:", "CONFIRMED" if ok else "REFUTED")
    return 0


if __name__ == '__main__':
    sys.exit(main())
