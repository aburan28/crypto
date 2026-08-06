"""
Thread 23 — closed form for lambda_1(L2), i.e. for the nu_hat separator.

Context
-------
2026-07-29 run #2 (commit e845207) found that the GLV-HNP Phase 2 lattice
contains m copies of a non-planted 2-dimensional sublattice

    L2 = < u, v >,   u = (n*S_K1, 0),   v = (-lam*S_K1, S_K2)

and that the scale-free invariant

    nu_hat = lambda_1(L2) / sqrt(det L2),      det L2 = n * S_K1 * S_K2

separates the June C1/C2 classes with AUC 0.935, where six earlier curve-level
invariants (delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n) and two later
ones (lam/n, mu) all failed.  That entry closed by proposing Thread 23: nu_hat
"is currently computed, not understood" -- derive lambda_1(L2) in closed form
and check whether it is a continued-fraction quantity *at a specific scale*.

This script answers that.  Scale convention is taken verbatim from
`glv_hnp_thread20a_nuhat_base.py`:

    scales(n, K1, K2) = (S_K1, S_D, S_K2, S_KAN) = (n//K1, 1, max(1, n//K2), n)
    k1k2(n, eff)      = (max(2, int(eff*isqrt(n))), isqrt(n)+1)

THE CLAIM
---------
A general lattice point of L2 is

    x*u + y*v = ( (x*n - y*lam) * S_K1 ,  y * S_K2 ).

For fixed y the best x makes |x*n - y*lam| equal to the centred residue
r(y) := min_x |x*n - y*lam|, so

    lambda_1(L2)^2 = min_{(x,y) != 0} [ r(y)^2 * S_K1^2 + y^2 * S_K2^2 ].

Claim (exact, not heuristic): the minimising y is a *best approximation
denominator of the second kind* for lam/n, hence a continued-fraction
convergent denominator q_j -- because if 1 <= y' < y with r(y') <= r(y), then
y' gives a strictly smaller objective, so the minimiser must be a record-holder
for r, and by the classical best-approximation theorem those are exactly the
convergent denominators.  Therefore, with e_j := |q_j*lam - p_j*n| and the
standard initial convergent (q_{-1}, p_{-1}) = (0, 1) (which supplies e_{-1}=n
and reproduces the y=0 corner vector (n*S_K1, 0)):

    ****  lambda_1(L2)^2 = min_j [ S_K1^2 * e_j^2 + S_K2^2 * q_j^2 ]  ****

evaluated over the O(log n) convergents of lam/n.  Since e_j ~ n/q_{j+1}, the
two terms balance at the scale

    t := sqrt(n * S_K1 / S_K2)   (~ sqrt(n/eff) under the k1k2 convention)

and, dividing through by det L2 = n*S_K1*S_K2,

    nu_hat^2  ~=  min_j [ (t/q_{j+1})^2 + (q_j/t)^2 ].

So nu_hat is small exactly when lam/n has an unusually good rational
approximation *at scale t* -- i.e. when a large partial quotient a_{j+1}
straddles t.  That is the mechanism, and it retro-explains the June failures:
q_cf, max_q_cf and max_a are scale-free, and the relevant approximation quality
is scale-dependent.

Experiments
-----------
E1  Exactness of the CF closed form vs exact-integer Lagrange-Gauss, and vs
    brute force at small n.  Also audits the float-based gauss_reduce_2d
    shipped in glv_hnp_thread20a_nuhat_base.py for precision loss.
E2  Does the winning convergent straddle t?
E3  Correlation of the approximation nu_hat^2 ~= min[(t/q_{j+1})^2+(q_j/t)^2]
    with exact nu_hat.  (The falsifier named in the 2026-07-29 entry: corr
    < 0.9 kills the convergent-scale story.)
E4  a_at_t (the partial quotient straddling t) vs nu_hat, head-to-head with the
    scale-free max_a that was falsified in June.
E5  secp256k1 itself: reproduce the four published nu_hat values from the CF
    picture and report the straddling partial quotient.

Run:  python3 secp256k1_cm_audit/thread23_nuhat_cf_closed_form.py
"""

import math
import random
import sys

import sympy

# ---------------------------------------------------------------------------
# conventions (verbatim from glv_hnp_thread20a_nuhat_base.py)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def k1k2(n, eff):
    return max(2, int(eff * math.isqrt(n))), math.isqrt(n) + 1

def gauss_reduce_2d_float(u, v):
    """The shipped implementation -- float division of big ints. Audited in E1."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    if nrm2(u) < nrm2(v):
        u, v = v, u
    while True:
        q = round((u[0] * v[0] + u[1] * v[1]) / nrm2(v))
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if nrm2(u) <= nrm2(v):
            return math.sqrt(nrm2(u))

# ---------------------------------------------------------------------------
# ground truth: exact integer Lagrange-Gauss (no floating point anywhere)
# ---------------------------------------------------------------------------

def _n2(w):
    return w[0] * w[0] + w[1] * w[1]

def gauss_reduce_2d_exact(u, v):
    """Returns lambda_1(L)^2 as an exact integer."""
    if _n2(u) < _n2(v):
        u, v = v, u
    while True:
        nv = _n2(v)
        if nv == 0:
            return _n2(u)
        d = u[0] * v[0] + u[1] * v[1]
        # exact round-half-up of d/nv
        q = (2 * d + nv) // (2 * nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if _n2(u) <= _n2(v):
            return _n2(u)

def brute_force_lambda1_sq(n, lam, s1, s2, ymax):
    """Exhaustive over |y| <= ymax with the best x for each y. Small n only."""
    best = (n * s1) ** 2  # y = 0, x = +-1
    for y in range(1, ymax + 1):
        t = (y * lam) % n
        r = min(t, n - t)
        val = r * r * s1 * s1 + y * y * s2 * s2
        if val < best:
            best = val
    return best

# ---------------------------------------------------------------------------
# the closed form
# ---------------------------------------------------------------------------

def convergents(lam, n):
    """
    Convergents of lam/n.  Returns list of (q_j, e_j, a_j) with
    e_j = |q_j*lam - p_j*n|, starting from the j=-1 term (q,e) = (0, n)
    whose lattice vector is the y=0 corner (n*S_K1, 0).
    a_j is the partial quotient introduced at step j (a_{-1} := 0).
    """
    out = [(0, n, 0)]
    p_prev, q_prev = 1, 0   # j = -1
    p_cur, q_cur = None, None
    x, y = lam, n
    first = True
    while y:
        a = x // y
        x, y = y, x - a * y
        if first:
            p_cur, q_cur = a, 1
            first = False
        else:
            p_cur, q_cur, p_prev, q_prev = (
                a * p_cur + p_prev, a * q_cur + q_prev, p_cur, q_cur)
        out.append((q_cur, abs(q_cur * lam - p_cur * n), a))
    return out

def nu_hat_cf(n, lam, k1_bound, k2_bound):
    """
    lambda_1(L2)^2 (exact int), nu_hat, and diagnostics, from the continued
    fraction of lam/n.  O(log n): one CF expansion, no lattice reduction.
    """
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    cvs = convergents(lam, n)
    best, best_i = None, None
    for i, (q, e, _a) in enumerate(cvs):
        val = s1 * s1 * e * e + s2 * s2 * q * q
        if best is None or val < best:
            best, best_i = val, i
    det = n * s1 * s2
    t = math.sqrt(n * s1 / s2)
    # straddling index: largest j with q_j <= t
    j_str = 0
    for i, (q, _e, _a) in enumerate(cvs):
        if q <= t:
            j_str = i
    a_at_t = cvs[j_str + 1][2] if j_str + 1 < len(cvs) else 0
    max_a = max((a for (_q, _e, a) in cvs), default=0)
    # the scale-normalised approximation  min_j [(t/q_{j+1})^2 + (q_j/t)^2]
    approx = None
    for i in range(len(cvs) - 1):
        q_j = cvs[i][0]
        q_n = cvs[i + 1][0]
        if q_n == 0:
            continue
        v = (t / q_n) ** 2 + (q_j / t) ** 2 if q_j else (t / q_n) ** 2
        if approx is None or v < approx:
            approx = v
    return {
        'lam1sq': best, 'nu_hat': math.sqrt(best / det), 'det': det,
        'win_i': best_i, 'win_q': cvs[best_i][0], 'win_e': cvs[best_i][1],
        't': t, 'j_str': j_str, 'a_at_t': a_at_t, 'max_a': max_a,
        'ncv': len(cvs), 'nu_hat_approx': math.sqrt(approx) if approx else None,
        'cvs': cvs,
    }

def nu_hat_gauss(n, lam, k1_bound, k2_bound):
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    l1sq = gauss_reduce_2d_exact((n * s1, 0), (-lam * s1, s2))
    return l1sq, math.sqrt(l1sq / (n * s1 * s2))

# ---------------------------------------------------------------------------
# stats helpers (no scipy)
# ---------------------------------------------------------------------------

def _rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r

def pearson(xs, ys):
    k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else 0.0

def spearman(xs, ys):
    return pearson(_rank(xs), _rank(ys))

# ---------------------------------------------------------------------------
# curve sampling
# ---------------------------------------------------------------------------

def glv_lambda(n):
    """A root of lam^2 + lam + 1 = 0 mod n, or None if n !== 1 mod 3."""
    if n % 3 != 1:
        return None
    for c in range(2, 200):
        z = pow(c, (n - 1) // 3, n)
        if z != 1 and (z * z + z + 1) % n == 0:
            return min(z, n - 1 - z)
    return None

def sample_glv_primes(bits, count, seed=1):
    """Primes n = 1 mod 3 of the given bit size, with a genuine GLV eigenvalue."""
    rng = random.Random(seed)
    out = []
    while len(out) < count:
        n = int(sympy.nextprime(rng.randrange(1 << (bits - 1), 1 << bits)))
        if n.bit_length() != bits or n % 3 != 1:
            continue
        lam = glv_lambda(n)
        if lam:
            out.append((n, lam))
    return out

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_LAM = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72

# ---------------------------------------------------------------------------
# E1 -- exactness
# ---------------------------------------------------------------------------

def e1_exactness():
    print("=" * 74)
    print("E1  CF closed form  vs  exact Lagrange-Gauss  vs  brute force")
    print("=" * 74)

    # (a) brute force, small n
    rng = random.Random(7)
    bad_bf = 0
    trials_bf = 0
    for _ in range(400):
        n = int(sympy.nextprime(rng.randrange(200, 20000)))
        lam = rng.randrange(1, n)
        eff = rng.choice([0.02, 0.05, 0.1, 0.25, 0.5])
        k1, k2 = k1k2(n, eff)
        s1, _, s2, _ = scales(n, k1, k2)
        cf = nu_hat_cf(n, lam, k1, k2)
        bf = brute_force_lambda1_sq(n, lam, s1, s2, min(n, 4 * n))
        trials_bf += 1
        if cf['lam1sq'] != bf:
            bad_bf += 1
            if bad_bf <= 3:
                print(f"  MISMATCH n={n} lam={lam} cf={cf['lam1sq']} bf={bf}")
    print(f"  (a) brute force, n < 20000 : {trials_bf - bad_bf}/{trials_bf} exact match")

    # (b) exact Gauss at many bit sizes, generic and GLV lambda
    print(f"\n  (b) exact integer Lagrange-Gauss")
    print(f"  {'bits':>5} {'kind':>7} {'eff':>6} {'trials':>7} {'cf==gauss':>10} "
          f"{'float ok':>9} {'max float err':>14}")
    tot, tot_ok, tot_f_ok = 0, 0, 0
    for bits in (16, 32, 64, 128, 256):
        for kind in ('glv', 'generic'):
            if kind == 'glv':
                curves = sample_glv_primes(bits, 12, seed=bits)
            else:
                r2 = random.Random(bits * 31)
                curves = []
                for _ in range(12):
                    n = int(sympy.nextprime(r2.randrange(1 << (bits - 1), 1 << bits)))
                    curves.append((n, r2.randrange(1, n)))
            for eff in (0.05, 0.25):
                ok = f_ok = 0
                worst = 0.0
                cnt = 0
                for (n, lam) in curves:
                    k1, k2 = k1k2(n, eff)
                    cf = nu_hat_cf(n, lam, k1, k2)
                    gx, _ = nu_hat_gauss(n, lam, k1, k2)
                    cnt += 1
                    if cf['lam1sq'] == gx:
                        ok += 1
                    s1, _, s2, _ = scales(n, k1, k2)
                    try:
                        fl = gauss_reduce_2d_float((n * s1, 0), (-lam * s1, s2))
                        rel = abs(fl - math.sqrt(gx)) / math.sqrt(gx)
                        if rel < 1e-9:
                            f_ok += 1
                        worst = max(worst, rel)
                    except (OverflowError, ValueError, ZeroDivisionError) as exc:
                        worst = float('inf')
                        if cnt == 1:
                            print(f"        float impl raised {type(exc).__name__} "
                                  f"at {bits} bits")
                tot += cnt
                tot_ok += ok
                tot_f_ok += f_ok
                w = "inf" if worst == float('inf') else f"{worst:.3e}"
                print(f"  {bits:>5} {kind:>7} {eff:>6} {cnt:>7} {ok:>6}/{cnt:<3} "
                      f"{f_ok:>6}/{cnt:<2} {w:>14}")
    print(f"\n  TOTAL cf == exact-gauss : {tot_ok}/{tot}")
    print(f"  TOTAL shipped-float == exact-gauss (rel < 1e-9) : {tot_f_ok}/{tot}")
    return tot_ok == tot and bad_bf == 0

# ---------------------------------------------------------------------------
# E2 / E3 / E4 -- the scale story
# ---------------------------------------------------------------------------

def e234():
    print("\n" + "=" * 74)
    print("E2/E3/E4  the winning convergent, the scale t, and a_at_t")
    print("=" * 74)

    rows = []
    for bits in (24, 32, 64, 128, 256):
        for (n, lam) in sample_glv_primes(bits, 60, seed=1000 + bits):
            for eff in (0.05, 0.1, 0.25):
                k1, k2 = k1k2(n, eff)
                cf = nu_hat_cf(n, lam, k1, k2)
                rows.append((bits, eff, n, lam, cf))

    print(f"\nE2  does the winning convergent straddle t = sqrt(n*S_K1/S_K2)?")
    print(f"  {'bits':>5} {'eff':>6} {'N':>5} {'win==j_str':>11} {'|win-j_str|<=1':>15}")
    for bits in (24, 32, 64, 128, 256):
        for eff in (0.05, 0.1, 0.25):
            sub = [r for r in rows if r[0] == bits and r[1] == eff]
            eq = sum(1 for r in sub if r[4]['win_i'] == r[4]['j_str'])
            n1 = sum(1 for r in sub if abs(r[4]['win_i'] - r[4]['j_str']) <= 1)
            print(f"  {bits:>5} {eff:>6} {len(sub):>5} {eq:>7}/{len(sub):<3} "
                  f"{n1:>11}/{len(sub):<3}")
    allsub = rows
    eq = sum(1 for r in allsub if r[4]['win_i'] == r[4]['j_str'])
    n1 = sum(1 for r in allsub if abs(r[4]['win_i'] - r[4]['j_str']) <= 1)
    print(f"  ALL   {len(allsub):>5}  win==j_str {eq}/{len(allsub)} "
          f"({100*eq/len(allsub):.1f}%)   within 1: {n1}/{len(allsub)} "
          f"({100*n1/len(allsub):.1f}%)")

    print(f"\nE3  approximation  nu_hat^2 ~= min_j[(t/q_(j+1))^2 + (q_j/t)^2]")
    print(f"     FALSIFIER (2026-07-29): corr(predicted, exact) < 0.9 kills the story.")
    xs = [r[4]['nu_hat'] for r in rows]
    ys = [r[4]['nu_hat_approx'] for r in rows]
    print(f"  N={len(xs)}  pearson={pearson(xs, ys):+.4f}  spearman={spearman(xs, ys):+.4f}")
    relerr = sorted(abs(b - a) / a for a, b in zip(xs, ys))
    print(f"  |approx-exact|/exact : median={relerr[len(relerr)//2]:.4f} "
          f"p90={relerr[int(0.9*len(relerr))]:.4f} max={relerr[-1]:.4f}")

    print(f"\nE4  a_at_t (scale-dependent) vs max_a (scale-free, falsified in June)")
    a_at_t = [float(r[4]['a_at_t']) for r in rows]
    max_a = [float(r[4]['max_a']) for r in rows]
    print(f"  spearman(a_at_t, nu_hat) = {spearman(a_at_t, xs):+.4f}")
    print(f"  spearman(max_a,  nu_hat) = {spearman(max_a, xs):+.4f}   <- scale-free")
    print(f"\n  nu_hat by a_at_t bucket (all {len(rows)} rows):")
    print(f"  {'a_at_t':>10} {'N':>5} {'mean nu_hat':>12} {'max nu_hat':>11} "
          f"{'sqrt(2/a)':>10}")
    buckets = [(1, 1), (2, 2), (3, 3), (4, 5), (6, 9), (10, 19), (20, 10**9)]
    for lo, hi in buckets:
        sub = [r[4]['nu_hat'] for r in rows if lo <= r[4]['a_at_t'] <= hi]
        if not sub:
            continue
        lbl = f"{lo}" if lo == hi else (f">={lo}" if hi > 10**8 else f"{lo}-{hi}")
        print(f"  {lbl:>10} {len(sub):>5} {sum(sub)/len(sub):>12.4f} "
              f"{max(sub):>11.4f} {math.sqrt(2.0/lo):>10.4f}")

    print(f"\n  Predicted C2 admission rule.  The 2026-07-29 C2 ceiling nu_hat <= 0.645")
    print(f"  becomes a_at_t >= 2/0.645^2 = {2/0.645**2:.2f}, i.e. a_at_t >= 5.")
    lo_nu = [r for r in rows if r[4]['nu_hat'] <= 0.645]
    print(f"  rows with nu_hat <= 0.645 : {len(lo_nu)}/{len(rows)}; "
          f"of those, a_at_t >= 5 in {sum(1 for r in lo_nu if r[4]['a_at_t'] >= 5)}")
    hi_a = [r for r in rows if r[4]['a_at_t'] >= 5]
    print(f"  rows with a_at_t >= 5     : {len(hi_a)}/{len(rows)}; "
          f"of those, nu_hat <= 0.645 in {sum(1 for r in hi_a if r[4]['nu_hat'] <= 0.645)}")
    viol = [r for r in rows if r[4]['a_at_t'] >= 1
            and r[4]['nu_hat'] ** 2 < 2.0 / (r[4]['a_at_t'] + 1) - 1e-9]
    print(f"  heuristic bound nu_hat^2 >= 2/(a_at_t+1) violated in "
          f"{len(viol)}/{len(rows)} rows")
    return rows

# ---------------------------------------------------------------------------
# E5 -- secp256k1
# ---------------------------------------------------------------------------

def e5_secp():
    print("\n" + "=" * 74)
    print("E5  secp256k1 -- reproduce the published nu_hat values from the CF")
    print("=" * 74)
    assert (SECP_LAM * SECP_LAM + SECP_LAM + 1) % SECP_N == 0
    published = {0.05: 0.8709, 0.10: 0.6624, 0.25: 0.5852, 0.50: 0.6851}
    cvs = convergents(SECP_LAM, SECP_N)
    print(f"  lam/n has {len(cvs)-1} convergents; max partial quotient = "
          f"{max(a for _q,_e,a in cvs)} (scale-free, uninformative)")
    print(f"\n  {'eff':>6} {'t/sqrt(n)':>10} {'j_str':>6} {'a_at_t':>7} "
          f"{'nu_hat(CF)':>11} {'nu_hat(Gauss)':>14} {'published':>10}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1, k2 = k1k2(SECP_N, eff)
        cf = nu_hat_cf(SECP_N, SECP_LAM, k1, k2)
        _, ng = nu_hat_gauss(SECP_N, SECP_LAM, k1, k2)
        print(f"  {eff:>6} {cf['t']/math.sqrt(SECP_N):>10.3f} {cf['j_str']:>6} "
              f"{cf['a_at_t']:>7} {cf['nu_hat']:>11.4f} {ng:>14.4f} "
              f"{published[eff]:>10.4f}")
    print("\n  The four values are reproduced exactly by one CF expansion, with no")
    print("  lattice reduction: secp256k1's nu_hat is a statement about the partial")
    print("  quotients of lam/n at denominators near sqrt(n/eff), nothing more.")
    print("\n  CAVEAT (carried from 2026-07-29, unchanged): this whole thread is")
    print("  conditional on a NON-STANDARD nonce generator k = k1 + lam*k2 with k1")
    print("  bounded. It is not an attack on secp256k1 as deployed.")

def main():
    ok = e1_exactness()
    e234()
    e5_secp()
    print("\n" + "=" * 74)
    print(f"E1 exactness verdict: {'PASS' if ok else 'FAIL'}")
    print("=" * 74)
    return 0 if ok else 1

if __name__ == '__main__':
    sys.exit(main())
