"""
GLV-HNP -- Thread 23-CF: closed form for nu_hat via continued fractions.

Background
----------
The 2026-07-29 run (commit e845207) found an empirical separator for the
Phase 2 / June-C1C2 lattice attack:

    L2     = < (n*S1, 0), (-lam*S1, S2) >        (a non-planted 2D sublattice)
    nu_hat = lambda_1(L2) / sqrt(det L2),        det L2 = n*S1*S2

with S1 = n//K1, S2 = max(1, n//K2)  (`scales()` in
glv_hnp_phase2_lambda_threshold_20a.py:229).  Empirically small nu_hat =>
attack succeeds; AUC 0.935 against the June C1/C2 classes.  nu_hat was
*computed* (Lagrange-Gauss) but not *understood*, and its next-step proposal
was to derive it in closed form.  This script does that.

Theory (proved here, verified numerically in Part A)
----------------------------------------------------
A general element of L2 is

    a*(-lam*S1, S2) + b*(n*S1, 0) = ( (b*n - a*lam)*S1 , a*S2 ).

For fixed a the best b makes the first coordinate S1*||a*lam||_n, where
||x||_n is the distance from x to the nearest multiple of n.  Hence

    lambda_1(L2)^2 = min_{a >= 0} [ (S1 * ||a*lam||_n)^2 + (S2 * a)^2 ].

Now let p_j/q_j be the continued-fraction convergents of lam/n and

    eta_j = |q_j*lam - p_j*n|   (so eta_j = ||q_j*lam||_n).

The convergent denominators are exactly the best approximations of the second
kind: for 0 < a < q_{j+1} we have ||a*lam||_n >= eta_j where q_j is the largest
convergent denominator <= a.  So (q_j, eta_j) dominates (a, ||a*lam||_n)
coordinate-wise, and since the objective is increasing in both coordinates the
minimum is always attained at a convergent.  Including the index j = -1
(q_{-1} = 0, eta_{-1} = n), which is the a = 0 vector (n*S1, 0):

    ***  lambda_1(L2) = min_j sqrt( (S1*eta_j)^2 + (S2*q_j)^2 )  ***

Normalising by sqrt(det L2) = sqrt(n*S1*S2) and setting

    x_j = q_j   / sqrt(n*S1/S2)          y_j = eta_j / sqrt(n*S2/S1)

gives the scale-free form

    ***  nu_hat^2 = min_j ( x_j^2 + y_j^2 ),      x_j * y_j = q_j*eta_j/n  ***

Three consequences, each tested below:

 (1) The minimising convergent sits at the *critical scale*
     q_j ~ sqrt(n*S1/S2) = sqrt(n*K2/K1).   [Part B]

 (2) Since x*y is fixed per convergent, x^2 + y^2 >= 2*x*y, so

         nu_hat^2 >= 2 * c_{j*},     c_j := q_j*eta_j/n.

     The three-term identity eta_j*q_{j+1} + eta_{j+1}*q_j = n gives

         c_j = 1 / ( q_{j+1}/q_j + eta_{j+1}/eta_j ),

     so c_j is small exactly when q_{j+1}/q_j is large, i.e. when the partial
     quotient a_{j+1} is large.  Therefore:

     *** nu_hat is small  <=>  lam/n has a large partial quotient AT the
         critical scale sqrt(n*K2/K1) ***                          [Part C]

 (3) This retro-explains the June falsifications.  q_cf, max_q_cf and max_a are
     scale-free functionals of the CF of lam/n -- max_a is the largest partial
     quotient *anywhere*.  The predictive quantity is the partial quotient at
     one specific index.  A large a_j far from the critical scale does nothing;
     a moderate a_j at the critical scale does everything.  Part D exhibits
     curve pairs with identical max_a and widely different nu_hat, which is a
     direct proof that no function of max_a alone can reproduce nu_hat.

Part E derives a new *algebraic* predictor a_crit (the partial quotient at the
critical index) and checks how much of nu_hat's separating power it retains.

Run: python3 glv_hnp_nuhat_cf.py
"""

import math
import random
import sys

import sympy

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_lambda_threshold_20a.py")
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

scales = _t20a.scales
rival_sublattice_nu = _t20a.rival_sublattice_nu
gauss_reduce_2d = _t20a.gauss_reduce_2d
eisenstein_decompose = _t20a.eisenstein_decompose
j0_traces = _t20a.j0_traces
glv_eigenvalues = _t20a.glv_eigenvalues

SECP_P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


# ---------------------------------------------------------------------------
# Continued fractions of lam/n
# ---------------------------------------------------------------------------

def cf_convergents_clean(lam, n):
    """
    Convergents p_j/q_j of lam/n as records (j, a_j, q_j, eta_j) with
    eta_j = |q_j*lam - p_j*n|.  Index j = -1 is the pseudo-convergent
    (q, eta) = (0, n), representing the a = 0 lattice vector (n*S1, 0).
    The list ends at eta = 0 (q = n), reached since gcd(lam, n) = 1.
    """
    a_list = []
    x, y = lam, n
    while y:
        a = x // y
        a_list.append(a)
        x, y = y, x - a * y

    # standard seeds: p_{-1}=1, q_{-1}=0 ; p_{-2}=0, q_{-2}=1
    out = [(-1, 0, 0, n)]
    p_prev2, q_prev2 = 0, 1
    p_prev1, q_prev1 = 1, 0
    for j, a in enumerate(a_list):
        p = a * p_prev1 + p_prev2
        q = a * q_prev1 + q_prev2
        eta = abs(q * lam - p * n)
        out.append((j, a, q, eta))
        p_prev2, q_prev2 = p_prev1, q_prev1
        p_prev1, q_prev1 = p, q
    return out


def nu_hat_cf(n, lam, k1_bound, k2_bound, want_detail=False):
    """
    Closed-form nu_hat via continued fractions.  Returns nu_hat, or
    (nu_hat, detail) where detail describes the minimising convergent.
    """
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    det = n * s1 * s2
    conv = cf_convergents_clean(lam % n, n)
    best = None
    for (j, a, q, eta) in conv:
        val = (s1 * eta) ** 2 + (s2 * q) ** 2
        if best is None or val < best[0]:
            best = (val, j, a, q, eta)
    val, j, a, q, eta = best
    nu = math.sqrt(val)
    nh = nu / math.sqrt(det)
    if not want_detail:
        return nh
    crit_scale = math.sqrt(n * s1 / s2)
    # a_crit: the partial quotient immediately AFTER the winning convergent
    a_next = None
    for (jj, aa, qq, ee) in conv:
        if jj == j + 1:
            a_next = aa
    detail = {
        'j': j, 'a_j': a, 'q': q, 'eta': eta,
        'a_next': a_next,
        'c': q * eta / n,
        'x': q / crit_scale, 'y': eta / (n / crit_scale),
        'crit_scale': crit_scale,
        'q_over_crit': (q / crit_scale) if crit_scale else float('nan'),
        'lower_bound': math.sqrt(2.0 * q * eta / n),
        'n_conv': len(conv),
    }
    return nh, detail


def max_partial_quotient(a, b):
    """The falsified June invariant `max_a`: largest partial quotient of a/b."""
    best = 0
    while b:
        q, r = divmod(a, b)
        best = max(best, q)
        a, b = b, r
    return best


# ---------------------------------------------------------------------------
# Sample harvesting (20-bit j=0 CM curves, Exp S regime)
# ---------------------------------------------------------------------------

def harvest_curves(n_curves, lo=2 ** 19, hi=2 ** 20):
    """(p, n, lam) for prime-order j=0 CM curves.  No generator needed here --
    this script only touches (n, lam), not the EC group."""
    out, seen = [], set()
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < n_curves:
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
                    seen.add(n)
                    out.append({'p': p, 'n': n, 'lam': roots[0]})
                    break
        p = int(sympy.nextprime(p))
    return out


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float('nan')


# ---------------------------------------------------------------------------
# Part A -- exactness of the closed form
# ---------------------------------------------------------------------------

def part_a(trials=4000):
    print("=" * 78)
    print("PART A -- closed form vs Lagrange-Gauss (exactness)")
    print("=" * 78)
    print("Claim: lambda_1(L2) = min_j sqrt((S1*eta_j)^2 + (S2*q_j)^2) over the")
    print("continued-fraction convergents of lam/n, exactly (not approximately).")
    print()

    rng = random.Random(20260805)
    worst = 0.0
    worst_case = None
    mismatches = 0
    bitsizes = [16, 20, 24, 32, 48, 64, 128, 256]
    per = trials // len(bitsizes)
    for bits in bitsizes:
        wb = 0.0
        for _ in range(per):
            n = int(sympy.nextprime(rng.randrange(2 ** (bits - 1), 2 ** bits)))
            lam = rng.randrange(1, n)
            k1 = rng.randrange(2, max(3, int(n ** rng.uniform(0.2, 0.8))))
            k2 = rng.randrange(2, max(3, n // max(2, k1)))
            nu_gauss, det, nh_gauss = rival_sublattice_nu(n, lam, k1, k2)
            nh_cf = nu_hat_cf(n, lam, k1, k2)
            rel = abs(nh_cf - nh_gauss) / max(nh_gauss, 1e-300)
            if rel > 1e-9:
                mismatches += 1
                if rel > worst:
                    worst, worst_case = rel, (n, lam, k1, k2, nh_gauss, nh_cf)
            wb = max(wb, rel)
        print(f"  {bits:>4}-bit: {per:>4} trials, max rel.dev = {wb:.3e}")
    print()
    print(f"  total trials   : {per * len(bitsizes)}")
    print(f"  mismatches     : {mismatches} (tolerance 1e-9)")
    if worst_case:
        print(f"  worst case     : {worst_case}")
    print(f"  VERDICT: {'CONFIRMED -- closed form is exact' if mismatches == 0 else 'FALSIFIED'}")
    print()
    return mismatches == 0


# ---------------------------------------------------------------------------
# Part B -- the winning convergent sits at the critical scale
# ---------------------------------------------------------------------------

def part_b(trials=2000):
    print("=" * 78)
    print("PART B -- location of the minimising convergent")
    print("=" * 78)
    print("Claim: the winner q_{j*} straddles the critical scale")
    print("       sqrt(n*S1/S2) = sqrt(n*K2/K1), i.e. q_{j*} <= crit <= q_{j*+1}.")
    print()

    rng = random.Random(4242)
    straddle = 0
    ratios = []
    for _ in range(trials):
        bits = rng.choice([20, 24, 32, 64, 128, 256])
        n = int(sympy.nextprime(rng.randrange(2 ** (bits - 1), 2 ** bits)))
        lam = rng.randrange(1, n)
        k1 = rng.randrange(2, max(3, int(n ** rng.uniform(0.25, 0.75))))
        k2 = rng.randrange(2, max(3, n // max(2, k1)))
        _, d = nu_hat_cf(n, lam, k1, k2, want_detail=True)
        conv = cf_convergents_clean(lam % n, n)
        q_next = None
        for (jj, aa, qq, ee) in conv:
            if jj == d['j'] + 1:
                q_next = qq
        if q_next is not None and d['q'] <= d['crit_scale'] <= q_next:
            straddle += 1
        if d['crit_scale'] > 0 and d['q'] > 0:
            ratios.append(math.log(d['q'] / d['crit_scale']))

    ratios.sort()
    def q(p):
        return ratios[min(len(ratios) - 1, int(p * len(ratios)))]
    print(f"  trials                              : {trials}")
    print(f"  q_{{j*}} <= crit_scale <= q_{{j*+1}} : {straddle}/{trials} "
          f"({100.0 * straddle / trials:.1f}%)")
    print(f"  log(q_{{j*}} / crit_scale) quantiles :")
    print(f"      q05={q(0.05):+.3f}  q25={q(0.25):+.3f}  q50={q(0.50):+.3f}  "
          f"q75={q(0.75):+.3f}  q95={q(0.95):+.3f}")
    print(f"  (a value of 0 means the winner sits exactly at the critical scale)")
    print()
    return straddle / trials


# ---------------------------------------------------------------------------
# Part C -- nu_hat^2 = x^2 + y^2 and the partial-quotient bound
# ---------------------------------------------------------------------------

def part_c(n_curves=200):
    print("=" * 78)
    print("PART C -- nu_hat^2 >= 2c and the role of the critical partial quotient")
    print("=" * 78)
    print("nu_hat^2 = x^2 + y^2  with  x*y = c = q*eta/n,  so nu_hat >= sqrt(2c).")
    print("c = 1/(q_{j+1}/q_j + eta_{j+1}/eta_j), so c is small exactly when the")
    print("NEXT partial quotient a_{j*+1} is large.  Hence: small nu_hat <=> a large")
    print("partial quotient at the critical scale.")
    print()

    curves = harvest_curves(n_curves)
    K1, K2 = 72, None
    rows = []
    for c in curves:
        n, lam = c['n'], c['lam']
        k2 = math.isqrt(n) + 1
        nh, d = nu_hat_cf(n, lam, K1, k2, want_detail=True)
        rows.append((nh, d, n, lam))

    tight = [abs(r[0] - r[1]['lower_bound']) / r[0] for r in rows]
    tight.sort()
    print(f"  sample: {len(rows)} 20-bit j=0 CM curves, K1={K1}, K2=isqrt(n)+1")
    print(f"  relative gap between nu_hat and its bound sqrt(2c):")
    print(f"      median={tight[len(tight)//2]:.4f}   "
          f"q90={tight[int(0.9*len(tight))]:.4f}   max={tight[-1]:.4f}")
    print()

    print("  nu_hat vs a_crit (partial quotient right after the winning convergent):")
    print(f"  {'a_crit':>8} {'count':>6} {'mean nu_hat':>12} {'min':>8} {'max':>8}")
    buckets = {}
    for nh, d, n, lam in rows:
        a = d['a_next'] if d['a_next'] is not None else 0
        key = 1 if a <= 1 else 2 if a == 2 else 3 if a <= 4 else 4 if a <= 8 else 5
        buckets.setdefault(key, []).append(nh)
    labels = {1: '1', 2: '2', 3: '3-4', 4: '5-8', 5: '>8'}
    for k in sorted(buckets):
        v = buckets[k]
        print(f"  {labels[k]:>8} {len(v):>6} {sum(v)/len(v):>12.3f} "
              f"{min(v):>8.3f} {max(v):>8.3f}")
    print()

    nhs = [r[0] for r in rows]
    acs = [(r[1]['a_next'] or 0) for r in rows]
    cs = [r[1]['c'] for r in rows]
    print(f"  spearman(a_crit, nu_hat) = {spearman(acs, nhs):+.3f}")
    print(f"  spearman(c,      nu_hat) = {spearman(cs, nhs):+.3f}")
    print()
    return rows


# ---------------------------------------------------------------------------
# Part D -- why max_a had to fail (constructive)
# ---------------------------------------------------------------------------

def part_d(rows):
    print("=" * 78)
    print("PART D -- why the June scale-free CF invariants had to fail")
    print("=" * 78)
    print("max_a is the largest partial quotient ANYWHERE in the CF of lam/n.")
    print("nu_hat depends on the partial quotient at ONE index (the critical scale).")
    print("If that is right, curves sharing a max_a value must still show a wide")
    print("spread of nu_hat -- so no function of max_a can reproduce nu_hat.")
    print()

    by_maxa = {}
    for nh, d, n, lam in rows:
        ma = max_partial_quotient(lam, n)
        by_maxa.setdefault(ma, []).append((nh, n, lam, d))
    print(f"  {'max_a':>7} {'count':>6} {'nu_hat min':>11} {'nu_hat max':>11} {'spread':>8}")
    shown = 0
    for ma in sorted(by_maxa, key=lambda k: -len(by_maxa[k])):
        v = by_maxa[ma]
        if len(v) < 3:
            continue
        nhs = [t[0] for t in v]
        print(f"  {ma:>7} {len(v):>6} {min(nhs):>11.3f} {max(nhs):>11.3f} "
              f"{max(nhs)-min(nhs):>8.3f}")
        shown += 1
        if shown >= 8:
            break
    print()

    all_nh = [r[0] for r in rows]
    all_ma = [max_partial_quotient(r[3], r[2]) for r in rows]
    print(f"  spearman(max_a, nu_hat) over the whole sample = "
          f"{spearman(all_ma, all_nh):+.3f}")
    print(f"  (compare spearman(a_crit, nu_hat) from Part C)")
    print()

    # explicit witness pair: same max_a, most different nu_hat
    best = None
    for ma, v in by_maxa.items():
        if len(v) < 2:
            continue
        v2 = sorted(v, key=lambda t: t[0])
        gap = v2[-1][0] - v2[0][0]
        if best is None or gap > best[0]:
            best = (gap, ma, v2[0], v2[-1])
    if best:
        gap, ma, lo, hi = best
        print(f"  WITNESS PAIR (identical max_a = {ma}):")
        print(f"    n={lo[1]:>8} lam={lo[2]:>8}  nu_hat={lo[0]:.4f}  "
              f"a_crit={lo[3]['a_next']}  q*={lo[3]['q']}")
        print(f"    n={hi[1]:>8} lam={hi[2]:>8}  nu_hat={hi[0]:.4f}  "
              f"a_crit={hi[3]['a_next']}  q*={hi[3]['q']}")
        print(f"    nu_hat differs by {gap:.4f} at identical max_a.")
    print()


# ---------------------------------------------------------------------------
# Part E -- secp256k1 in closed form
# ---------------------------------------------------------------------------

def part_e():
    print("=" * 78)
    print("PART E -- secp256k1 in closed form")
    print("=" * 78)
    n = SECP_N
    roots = glv_eigenvalues(n)
    assert roots is not None
    lam = roots[0]
    assert (lam * lam + lam + 1) % n == 0
    print(f"  n   = {n}")
    print(f"  lam = {lam}   (lam^2+lam+1 = 0 mod n verified)")
    print()
    print(f"  {'eff':>6} {'K1':>8} {'nu_hat(CF)':>11} {'nu_hat(GG)':>11} "
          f"{'j*':>4} {'a_crit':>7} {'log2 q*':>8} {'log2 crit':>10}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff * n / k2))
        nh, d = nu_hat_cf(n, lam, k1, k2, want_detail=True)
        _, _, nh_g = rival_sublattice_nu(n, lam, k1, k2)
        print(f"  {eff:>6.2f} {k1.bit_length():>7}b {nh:>11.4f} {nh_g:>11.4f} "
              f"{d['j']:>4} {str(d['a_next']):>7} "
              f"{math.log2(max(d['q'],1)):>8.1f} {math.log2(d['crit_scale']):>10.1f}")
    print()
    print("  Partial quotients of lam/n around the critical index (eff=0.10):")
    k2 = math.isqrt(n) + 1
    k1 = max(2, int(0.10 * n / k2))
    nh, d = nu_hat_cf(n, lam, k1, k2, want_detail=True)
    conv = cf_convergents_clean(lam, n)
    js = d['j']
    print(f"  {'j':>4} {'a_j':>8} {'log2 q_j':>9} {'log2 eta_j':>11} {'c_j':>8}")
    for (j, a, q, eta) in conv:
        if js - 4 <= j <= js + 4:
            mark = "  <-- winner" if j == js else ""
            print(f"  {j:>4} {a:>8} {math.log2(max(q,1)):>9.2f} "
                  f"{math.log2(max(eta,1)):>11.2f} {q*eta/n:>8.4f}{mark}")
    print()
    print(f"  total CF length of lam/n: {len(conv)-1} partial quotients")
    print(f"  max_a (largest anywhere): {max_partial_quotient(lam, n)}")
    print()
    print("  CAVEATS (unchanged from 2026-07-29): this is a property of a lattice")
    print("  that only arises under a NON-STANDARD nonce generator k = k1 + lam*k2")
    print("  with k1 bounded.  It is not an attack on secp256k1 with correct nonces.")
    print()


# ---------------------------------------------------------------------------
# Part F -- null distribution of nu_hat, CF prediction vs random 2D lattices
# ---------------------------------------------------------------------------

def part_f(draws=4000):
    print("=" * 78)
    print("PART F -- null distribution of nu_hat over random lam")
    print("=" * 78)
    print("If nu_hat is governed by CF statistics at a scale, its null law should")
    print("match the lambda_1/sqrt(det) law of random 2D lattices (Gauss-Kuzmin /")
    print("equidistribution).  Compare against the 2026-07-29 measured quantiles.")
    print()
    n = SECP_N
    k2 = math.isqrt(n) + 1
    k1 = max(2, int(0.05 * n / k2))
    rng = random.Random(31415)
    vals = []
    for _ in range(draws):
        lam = rng.randrange(1, n)
        vals.append(nu_hat_cf(n, lam, k1, k2))
    vals.sort()
    def qq(p):
        return vals[min(len(vals) - 1, int(p * len(vals)))]
    print(f"  256-bit, eff=0.05, {draws} random lam (CF closed form):")
    print(f"    q05={qq(0.05):.3f} q10={qq(0.10):.3f} q25={qq(0.25):.3f} "
          f"q50={qq(0.50):.3f} q75={qq(0.75):.3f} q90={qq(0.90):.3f} q95={qq(0.95):.3f}")
    print(f"  2026-07-29 measured (Lagrange-Gauss, same setting):")
    print(f"    q05=0.232 q10=0.322 q25=0.496 q50=0.715 q75=0.885 q90=0.974 q95=0.999")
    print()
    below45 = sum(1 for v in vals if v < 0.45) / len(vals)
    above90 = sum(1 for v in vals if v > 0.90) / len(vals)
    print(f"    P(nu_hat < 0.45) = {below45*100:.1f}%   (2026-07-29: 20.5%)")
    print(f"    P(nu_hat > 0.90) = {above90*100:.1f}%   (2026-07-29: 22.4%)")
    print()
    # Siegel small-ball prediction: E#{primitive pts of norm <= t*sqrt(det)} = pi t^2/zeta(2)
    # counting +/- pairs => P(nu_hat <= t) ~ (pi/2) t^2 / zeta(2) = 3 t^2 / pi for small t
    print("  small-t check against the Siegel mean-value prediction")
    print("      P(nu_hat <= t) ~ (pi/2)*t^2/zeta(2) = 3t^2/pi   (t -> 0)")
    print(f"  {'t':>6} {'empirical':>10} {'3t^2/pi':>10}")
    for t in (0.05, 0.10, 0.15, 0.20, 0.30):
        emp = sum(1 for v in vals if v <= t) / len(vals)
        pred = 3.0 * t * t / math.pi
        print(f"  {t:>6.2f} {emp:>10.4f} {pred:>10.4f}")
    print()


def main():
    print()
    print("#" * 78)
    print("# GLV-HNP Thread 23-CF: closed form for nu_hat via continued fractions")
    print("#" * 78)
    print()
    ok = part_a()
    part_b()
    rows = part_c()
    part_d(rows)
    part_e()
    part_f()
    print("=" * 78)
    print("DONE.  Closed form exact:", ok)
    print("=" * 78)


if __name__ == "__main__":
    main()
