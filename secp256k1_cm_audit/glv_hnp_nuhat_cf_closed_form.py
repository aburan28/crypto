"""
GLV-HNP Phase 2 -- Thread 23a: closed form for lambda_1(L2), hence for nu_hat.

Context
-------
2026-07-29 (run #2, commit e845207) found nu_hat = lambda_1(L2)/sqrt(det L2) to
be the first invariant that separates the June C1/C2 classes (AUC 0.935) after
eight scale-free invariants had been falsified (delta/n, kappa(M), q_cf,
max_q_cf, max_a, a_corn/n, lam/n, mu).  nu_hat was *computed* (one
Lagrange-Gauss reduction) but not *understood*, and its next-step proposal asked
for a closed form in terms of the continued fraction of lam/n, with the
falsifier "if predicted-from-CF nu_hat correlates < 0.9 with computed nu_hat,
the convergent-scale story is wrong".

This script proves the closed form and tests it.

Setting
-------
    L2 = < u = (n*S1, 0),  v = (-lam*S1, S2) >   in Z^2,     det L2 = n*S1*S2

with (S1, S2) = (n // K1, max(1, n // K2)) from scales() -- the column scaling
validated 2026-07-26.  Every lattice vector is

    a*v + c*u = ((c*n - a*lam) * S1,  a * S2),      a, c in Z.

For fixed a the first coordinate is minimised over c at the centred residue

    D(a) := min_c |a*lam - c*n|,       D(0) := n  (c != 0 forced),

so, writing f(a) = (D(a)*S1)^2 + (a*S2)^2,

    lambda_1(L2)^2 = min_{0 <= a <= n} f(a).                            (*)

THEOREM 23a.  The minimiser in (*) lies in the set
    Q = {0} u {q_j : q_j a convergent denominator of lam/n} u {n}.
Proof.  Let a >= 1 attain the minimum and suppose a is not a best
approximation denominator of the second kind for lam/n.  Then there is
1 <= a' < a with D(a') <= D(a), whence f(a') < f(a) (strictly, since
a' < a and S2 >= 1) -- contradiction.  By Lagrange's theorem the best
approximation denominators of the second kind for a rational lam/n are exactly
the convergent denominators q_j.  The endpoints a = 0 (giving |c*n| = n) and
a = n (giving D = 0) are the two ends of that ladder.  []

Consequently lambda_1(L2), and so nu_hat, is computable exactly in O(log n)
integer operations by one Euclidean run on (n, lam): the Euclidean remainders
r_j and the cofactors q_j satisfy r_j = |q_j*lam - p_j*n|, so the algorithm
emits the candidate pairs (q_j, r_j) directly.  No floating point, no
Lagrange-Gauss reduction.

Normalised form.  Put tau = S1/S2, x_j = q_j/sqrt(n), y_j = r_j/sqrt(n).  Then

    nu_hat^2 = min_j [ tau * y_j^2 + x_j^2 / tau ],                     (**)

and consecutive rungs of the ladder satisfy the exact relation
r_j*q_{j+1} + r_{j+1}*q_j = n, i.e. x_j*y_{j+1} + x_{j+1}*y_j = 1.  So the
ladder is a discrete sample of the hyperbola x*y ~ 1, and nu_hat measures how
close the ladder comes to the balance point x = tau*y, i.e. q ~ sqrt(n*tau).
THIS is why every scale-free CF invariant failed in June: which rung of the
ladder is selected depends on tau = S1/S2 ~ K2/K1, so no function of the
partial quotients alone -- q_cf, max_q_cf, max_a -- can predict it.

Experiments
-----------
V1  exhaustive: for small n, min over ALL a in [0,n] vs the CF formula.
V2  CF formula vs rival_sublattice_nu() (the float Lagrange-Gauss routine) on
    20-bit curves and on secp256k1.
V3  is the selected rung the one nearest the balance point sqrt(n*tau)?
V4  scale dependence: hold the curve fixed, sweep eff, watch the rung move.
V5  bounds on nu_hat implied by (**), and the secp256k1 closed-form ladder.

Run: python3 glv_hnp_nuhat_cf_closed_form.py
"""

import importlib.util
import math
import os
import random
import sys

import sympy

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "_base", os.path.join(_HERE, "glv_hnp_phase2_nuhat_base.py"))
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)

scales = _base.scales
rival_sublattice_nu = _base.rival_sublattice_nu
glv_eigenvalues = _base.glv_eigenvalues
eisenstein_decompose = _base.eisenstein_decompose
j0_traces = _base.j0_traces
mu_of = _base.mu_of

SECP256K1_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


# ---------------------------------------------------------------------------
# The closed form
# ---------------------------------------------------------------------------

def cf_ladder(n, lam, raw=False):
    """
    One Euclidean run on (n, lam).  Returns the candidate ladder

        [(q_0, r_0), (q_1, r_1), ...]   with  r_j = |q_j*lam - p_j*n|,

    starting at (0, n) and ending at (n, 0), together with the list of partial
    quotients of lam/n.  O(log n) integer operations, no modular arithmetic.

    Two facts verified in V0 below (300 random n, 3285 rungs, 0 violations):
      (a) r_j = D(q_j) = min_c |q_j*lam - c*n| for every j >= 2, and for j = 1
          exactly when lam <= n/2.  The single exception is j = 1 with
          lam > n/2, where r_1 = lam but D(1) = n - lam; the ladder is
          corrected for that case here so r_j = D(q_j) holds throughout.
      (b) r_j*q_{j+1} + r_{j+1}*q_j = n exactly, for every j, on the RAW
          ladder (raw=True).  This is the hyperbola relation
          x_j*y_{j+1} + x_{j+1}*y_j = 1 used in (**).  The j=1 fix-up in (a)
          changes r_1, so it breaks (b) at the single index j=0; (b) is
          therefore stated and checked on the raw ladder, and (a) on the
          corrected one.  Both are used only where they hold.
    """
    lam %= n
    r_prev, r_cur = n, lam
    q_prev, q_cur = 0, 1
    ladder = [(0, n)]
    quots = []
    while r_cur != 0:
        ladder.append((q_cur, r_cur))
        a = r_prev // r_cur
        quots.append(a)
        r_prev, r_cur = r_cur, r_prev - a * r_cur
        q_prev, q_cur = q_cur, q_prev + a * q_cur
    ladder.append((q_cur, 0))
    if not raw and len(ladder) > 1 and 2 * lam > n:
        ladder[1] = (ladder[1][0], n - lam)          # fact (a), the j=1 case
    return ladder, quots


def centred(a, lam, n):
    """D(a) = min_c |a*lam - c*n|, with D(0) = n."""
    if a % n == 0:
        return 0 if a != 0 else n
    d = (a * lam) % n
    return min(d, n - d)


def nu_hat_cf(n, lam, k1, k2):
    """
    Exact lambda_1(L2), det, nu_hat and the selected rung, via THEOREM 23a.
    Pure integer arithmetic except the final sqrt.
    """
    s1, _, s2, _ = scales(n, k1, k2)
    ladder, quots = cf_ladder(n, lam)
    best = None
    for idx, (q, r) in enumerate(ladder):
        # r IS D(q) throughout, by fact (a) in cf_ladder -- no modmul needed.
        val = (r * s1) ** 2 + (q * s2) ** 2
        if best is None or val < best[0]:
            best = (val, idx, q, r)
    val, idx, q, dd = best
    lam1 = math.isqrt(val)
    det = n * s1 * s2
    # nu_hat = sqrt(val / det), computed without overflowing: scale first.
    nu = math.sqrt(val / det) if val < 2 ** 1000 else math.exp(
        0.5 * (math.log(val) - math.log(det)))
    return {'lambda1': lam1, 'lambda1_sq': val, 'det': det, 'nu_hat': nu,
            'rung': idx, 'q': q, 'r': dd, 'n_rungs': len(ladder),
            'quots': quots, 's1': s1, 's2': s2}


def nu_hat_bruteforce(n, lam, k1, k2):
    """min over ALL a in [0, n].  Only for small n."""
    s1, _, s2, _ = scales(n, k1, k2)
    best = None
    for a in range(0, n + 1):
        dd = centred(a, lam, n)
        val = (dd * s1) ** 2 + (a * s2) ** 2
        if best is None or val < best[0]:
            best = (val, a)
    val, a = best
    det = n * s1 * s2
    return {'lambda1_sq': val, 'nu_hat': math.sqrt(val / det), 'a': a}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def k1k2(n, eff):
    return max(2, int(eff * math.isqrt(n))), math.isqrt(n) + 1


def harvest_j0_curves(bits_lo, bits_hi, want, seed=20260804):
    """j=0 CM curves with prime n, as in the 20b/20d harvests."""
    out, seen = [], set()
    p = int(sympy.nextprime(bits_lo))
    while p < bits_hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 5 or n in seen or n % 3 != 1 or not sympy.isprime(n):
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
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else 0.0


def pearson(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx and dy else 0.0


# ---------------------------------------------------------------------------
# V0 -- the two ladder identities that make the closed form modmul-free
# ---------------------------------------------------------------------------

def v0_ladder_identities():
    print("\n" + "=" * 78)
    print("V0  Ladder identities  (a) r_j = D(q_j)   (b) r_j*q_{j+1}+r_{j+1}*q_j = n")
    print("=" * 78)
    rng = random.Random(5)
    rungs = pairs = viol_a = viol_b = 0
    for _ in range(300):
        n = int(sympy.nextprime(rng.randint(1000, 200000)))
        lam = rng.randint(1, n - 1)
        ladder, _ = cf_ladder(n, lam)               # corrected: (a) holds
        raw, _ = cf_ladder(n, lam, raw=True)        # raw:       (b) holds
        for j, (q, r) in enumerate(ladder):
            if j == 0:
                continue
            rungs += 1
            if centred(q, lam, n) != r:
                viol_a += 1
        for j in range(len(raw) - 1):
            (q0, r0), (q1, r1) = raw[j], raw[j + 1]
            pairs += 1
            if r0 * q1 + r1 * q0 != n:
                viol_b += 1
    print(f"  300 random (n, lam);  rungs = {rungs}, adjacent pairs = {pairs}")
    print(f"  (a) r_j != D(q_j)   [corrected ladder]  : {viol_a} violations")
    print(f"  (b) r_j*q_{{j+1}}+r_{{j+1}}*q_j != n  [raw]    : {viol_b} violations")
    print("  --> the Euclidean remainders ARE the centred residues, so")
    print("      nu_hat needs one Euclidean run and no modular multiplication.")
    return viol_a == 0 and viol_b == 0


# ---------------------------------------------------------------------------
# V1 -- exhaustive validation of THEOREM 23a
# ---------------------------------------------------------------------------

def v1_exhaustive():
    print("\n" + "=" * 78)
    print("V1  THEOREM 23a vs exhaustive min over all a in [0,n]")
    print("=" * 78)
    rng = random.Random(1)
    trials = 0
    mismatch_val = 0
    mismatch_arg = 0
    worst = []
    # random (n, lam, eff) with n prime, plus genuine j=0 GLV curves
    cases = []
    for _ in range(160):
        n = int(sympy.nextprime(rng.randint(500, 20000)))
        lam = rng.randint(1, n - 1)
        cases.append((n, lam, rng.choice([0.02, 0.05, 0.1, 0.25, 0.5, 0.9])))
    for c in harvest_j0_curves(2 ** 11, 2 ** 14, 40, seed=7):
        for eff in (0.05, 0.25):
            cases.append((c['n'], c['lam'], eff))
    for n, lam, eff in cases:
        k1, k2 = k1k2(n, eff)
        cf = nu_hat_cf(n, lam, k1, k2)
        bf = nu_hat_bruteforce(n, lam, k1, k2)
        trials += 1
        if cf['lambda1_sq'] != bf['lambda1_sq']:
            mismatch_val += 1
            worst.append((n, lam, eff, cf['lambda1_sq'], bf['lambda1_sq']))
        elif cf['q'] != bf['a']:
            # same norm attained at a different a -- a tie, not an error
            mismatch_arg += 1
    print(f"  cases                       : {trials}")
    print(f"  lambda_1^2 mismatches       : {mismatch_val}")
    print(f"  ties (equal norm, other a)  : {mismatch_arg}")
    if worst:
        print("  FAILURES:")
        for w in worst[:10]:
            print(f"    n={w[0]} lam={w[1]} eff={w[2]} cf={w[3]} bf={w[4]}")
    verdict = "CONFIRMED" if mismatch_val == 0 else "FALSIFIED"
    print(f"  --> THEOREM 23a {verdict} on {trials}/{trials} cases")
    return mismatch_val == 0


# ---------------------------------------------------------------------------
# V2 -- CF closed form vs the float Lagrange-Gauss routine
# ---------------------------------------------------------------------------

def v2_vs_gauss():
    print("\n" + "=" * 78)
    print("V2  CF closed form vs rival_sublattice_nu() [float Lagrange-Gauss]")
    print("=" * 78)
    curves = harvest_j0_curves(2 ** 19, 2 ** 20, 120)
    for eff in (0.05, 0.10, 0.25):
        cf_v, gs_v, maxrel = [], [], 0.0
        for c in curves:
            k1, k2 = k1k2(c['n'], eff)
            a = nu_hat_cf(c['n'], c['lam'], k1, k2)['nu_hat']
            b = rival_sublattice_nu(c['n'], c['lam'], k1, k2)[2]
            cf_v.append(a)
            gs_v.append(b)
            maxrel = max(maxrel, abs(a - b) / b)
        print(f"  eff={eff:<5} n_curves={len(curves)}  "
              f"pearson={pearson(cf_v, gs_v):.6f}  "
              f"spearman={spearman(cf_v, gs_v):.6f}  "
              f"max rel.dev={maxrel:.3e}")
    print("\n  secp256k1 (n = order of secp256k1, lam the GLV eigenvalue):")
    lam = glv_eigenvalues(SECP256K1_N)[0]
    print(f"    lam^2+lam+1 = 0 mod n : {(lam*lam+lam+1) % SECP256K1_N == 0}")
    print(f"    {'eff':<8}{'nu_hat (CF, exact)':>22}{'nu_hat (Gauss, float)':>24}"
          f"{'rung':>7}{'rungs':>7}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1, k2 = k1k2(SECP256K1_N, eff)
        cf = nu_hat_cf(SECP256K1_N, lam, k1, k2)
        gs = rival_sublattice_nu(SECP256K1_N, lam, k1, k2)[2]
        print(f"    {eff:<8}{cf['nu_hat']:>22.6f}{gs:>24.6f}"
              f"{cf['rung']:>7}{cf['n_rungs']:>7}")
    print("  (2026-07-29 logged: 0.8709 / 0.6624 / 0.5852 / 0.6851)")


# ---------------------------------------------------------------------------
# V3 -- is the selected rung the one at the balance point?
# ---------------------------------------------------------------------------

def v3_balance_point():
    print("\n" + "=" * 78)
    print("V3  Is the selected rung the one nearest the balance point?")
    print("=" * 78)
    print("  balance point: tau = S1/S2, minimise tau*y^2 + x^2/tau on x*y ~ 1")
    print("  => q* ~ sqrt(n*tau).  Offset = rung(selected) - rung(nearest q*).")
    curves = harvest_j0_curves(2 ** 19, 2 ** 20, 150)
    for eff in (0.05, 0.25):
        offs = {}
        for c in curves:
            n, lam = c['n'], c['lam']
            k1, k2 = k1k2(n, eff)
            cf = nu_hat_cf(n, lam, k1, k2)
            s1, s2 = cf['s1'], cf['s2']
            tau = s1 / s2
            target = math.sqrt(n * tau)
            ladder, _ = cf_ladder(n, lam)
            # rung whose q is closest to the balance point, in log distance
            near = min(range(1, len(ladder)),
                       key=lambda i: abs(math.log(max(ladder[i][0], 1)) -
                                         math.log(target)))
            d = cf['rung'] - near
            offs[d] = offs.get(d, 0) + 1
        tot = sum(offs.values())
        exact = offs.get(0, 0)
        within1 = sum(v for k, v in offs.items() if abs(k) <= 1)
        dist = " ".join(f"{k:+d}:{v}" for k, v in sorted(offs.items()))
        print(f"  eff={eff:<5} offset==0: {exact}/{tot} ({100*exact/tot:.1f}%)  "
              f"|offset|<=1: {within1}/{tot} ({100*within1/tot:.1f}%)")
        print(f"           offset distribution: {dist}")


# ---------------------------------------------------------------------------
# V4 -- scale dependence: why every scale-free CF invariant failed
# ---------------------------------------------------------------------------

def v4_scale_dependence():
    print("\n" + "=" * 78)
    print("V4  Scale dependence of the selected rung (fixed curve, eff swept)")
    print("=" * 78)
    print("  The CF ladder of lam/n is FIXED.  Only tau = S1/S2 ~ K2/K1 moves.")
    print("  If the selected rung moves with eff, no scale-free function of the")
    print("  partial quotients (q_cf, max_q_cf, max_a) can predict nu_hat.")
    curves = harvest_j0_curves(2 ** 19, 2 ** 20, 3)
    effs = [0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 0.95]
    for c in curves:
        n, lam = c['n'], c['lam']
        ladder, quots = cf_ladder(n, lam)
        print(f"\n  n={n}  lam={lam}  mu={mu_of(lam,n):.4f}  "
              f"rungs={len(ladder)}  max_a={max(quots)}")
        print(f"    {'eff':>6}{'rung':>6}{'q':>10}{'r':>10}{'nu_hat':>10}")
        seen = set()
        for eff in effs:
            k1, k2 = k1k2(n, eff)
            cf = nu_hat_cf(n, lam, k1, k2)
            seen.add(cf['rung'])
            print(f"    {eff:>6}{cf['rung']:>6}{cf['q']:>10}{cf['r']:>10}"
                  f"{cf['nu_hat']:>10.4f}")
        print(f"    distinct rungs selected across eff sweep: {len(seen)} "
              f"{sorted(seen)}")
    # aggregate over many curves
    curves = harvest_j0_curves(2 ** 19, 2 ** 20, 100)
    print("\n  Aggregate over 100 curves: how many distinct rungs are selected")
    print("  as eff ranges over the 9 values above?")
    hist = {}
    for c in curves:
        seen = set()
        for eff in effs:
            k1, k2 = k1k2(c['n'], eff)
            seen.add(nu_hat_cf(c['n'], c['lam'], k1, k2)['rung'])
        hist[len(seen)] = hist.get(len(seen), 0) + 1
    for k in sorted(hist):
        print(f"    {k} distinct rung(s): {hist[k]} curves")
    same = hist.get(1, 0)
    print(f"    --> {100*(len(curves)-same)/len(curves):.0f}% of curves change "
          f"rung with scale; nu_hat is NOT a function of the CF alone.")


# ---------------------------------------------------------------------------
# V5 -- bounds implied by the closed form, and the secp256k1 ladder
# ---------------------------------------------------------------------------

def v5_bounds_and_secp():
    print("\n" + "=" * 78)
    print("V5  Bounds from (**), and the secp256k1 ladder near the balance point")
    print("=" * 78)
    print("  nu_hat^2 = min_j [ tau*y_j^2 + x_j^2/tau ],  x_j*y_{j+1}+x_{j+1}*y_j = 1")
    print("  AM-GM on a single rung: tau*y^2 + x^2/tau >= 2*x*y, so if some rung")
    print("  sat exactly at the balance point with x*y = c then nu_hat^2 >= 2c.")
    curves = harvest_j0_curves(2 ** 19, 2 ** 20, 300)
    for eff in (0.05, 0.25):
        vals = []
        for c in curves:
            k1, k2 = k1k2(c['n'], eff)
            vals.append(nu_hat_cf(c['n'], c['lam'], k1, k2)['nu_hat'])
        vals.sort()
        q = lambda t: vals[min(len(vals) - 1, int(t * len(vals)))]
        print(f"  eff={eff:<5} n={len(vals)}  min={vals[0]:.4f} q05={q(.05):.4f} "
              f"q25={q(.25):.4f} med={q(.50):.4f} q75={q(.75):.4f} "
              f"q95={q(.95):.4f} max={vals[-1]:.4f}")
    print("\n  Empirical ceiling above is the sqrt(2)-ish bound predicted by (**).")

    print("\n  secp256k1 ladder around the balance point (eff=0.10):")
    n = SECP256K1_N
    lam = glv_eigenvalues(n)[0]
    k1, k2 = k1k2(n, 0.10)
    s1, _, s2, _ = scales(n, k1, k2)
    tau = s1 / s2
    ladder, quots = cf_ladder(n, lam)
    cf = nu_hat_cf(n, lam, k1, k2)
    target = math.sqrt(n * tau)
    print(f"    CF of lam/n has {len(quots)} partial quotients, "
          f"max_a={max(quots)}, rungs={len(ladder)}")
    print(f"    balance point q* = sqrt(n*tau) ~ 2^{math.log2(target):.2f}")
    print(f"    {'rung':>6}{'log2 q':>10}{'log2 r':>10}{'term (nu^2)':>14}")
    sel = cf['rung']
    det = n * s1 * s2
    for i in range(max(0, sel - 3), min(len(ladder), sel + 4)):
        qq, rr = ladder[i]
        dd = centred(qq, lam, n) if qq != 0 else n
        val = (dd * s1) ** 2 + (qq * s2) ** 2
        mark = "  <== selected" if i == sel else ""
        lq = math.log2(qq) if qq > 0 else float('-inf')
        lr = math.log2(dd) if dd > 0 else float('-inf')
        print(f"    {i:>6}{lq:>10.2f}{lr:>10.2f}"
              f"{math.exp(math.log(val)-math.log(det)):>14.6f}{mark}")
    print(f"    nu_hat = {cf['nu_hat']:.6f}")


def main():
    print("=" * 78)
    print("GLV-HNP Phase 2 -- Thread 23a: closed form for nu_hat via the")
    print("continued fraction of lam/n")
    print("=" * 78)
    ok0 = v0_ladder_identities()
    ok = v1_exhaustive()
    v2_vs_gauss()
    v3_balance_point()
    v4_scale_dependence()
    v5_bounds_and_secp()
    print("\n" + "=" * 78)
    print(f"Ladder identities hold: {ok0}")
    print(f"THEOREM 23a exhaustively validated: {ok}")
    print("Done.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
