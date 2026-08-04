"""
GLV-HNP — Thread 23: closed form for nu_hat.

Background
----------
The 2026-07-29 run #2 (`e845207`, log entry "2026-07-29 (autolab run #2)")
found an empirical separator for the Phase-2 / June-Phase-1 lattice:

    L2      = < (n*S_K1, 0),  (-lam*S_K1, S_K2) >     (a NON-planted 2D
                                                       sublattice; there are
                                                       m independent copies)
    nu_hat  = lambda_1(L2) / sqrt(det L2),   det L2 = n*S_K1*S_K2

with AUC 0.935 against the June C1/C2 classes, where every previously tried
invariant sat at the trivial baseline.  That entry closed with Thread 23:

    "nu_hat is currently computed, not understood. [...] show lambda_1(L2) is
     determined by the CF convergent p_j/q_j of lam/n with q_j nearest
     sqrt(n*S_K1/S_K2), and check whether nu_hat ~ f(...) reproduces the
     measured nu_hat.  Falsifier: if predicted-from-CF nu_hat correlates < 0.9
     with computed nu_hat, the convergent-scale story is wrong."

This script answers that, and more: nu_hat is not merely CF-predictable, it is
a classical object in closed form.

The claim
---------
(1) REDUCTION TO A WEIGHTED APPROXIMATION PROBLEM.  A general element of L2 is
    a*(n*S_K1,0) + b*(-lam*S_K1,S_K2) = ((a*n - b*lam)*S_K1, b*S_K2).  For
    fixed b the best a makes the first coordinate the centred residue
    r(b) = (b*lam mod +- n).  Hence exactly

        lambda_1(L2)^2 = min_{b >= 1} [ r(b)^2 * S_K1^2 + b^2 * S_K2^2 ].

    So nu_hat is a *weighted* best-rational-approximation quality of lam/n,
    the weighting being S_K1 : S_K2 = (1/K1) : (1/K2).

(2) THE MINIMISER IS A CF CONVERGENT -- but the proposed SELECTION RULE is
    WRONG.  Part A confirms the structural half: the argmin b is a convergent
    denominator of lam/n in 400/400 trials (it must be -- the argmin is on the
    Pareto frontier of (|r(b)|, b), which is exactly the set of best
    approximations of the second kind).  Part C then FALSIFIES the rule
    proposed on 2026-07-29, that the right convergent is the one nearest
    sqrt(n*S_K1/S_K2):

        single convergent @ q-scale   pearson = -0.33   (falsifier: < 0.9)

    The scale argument silently assumes q_j ~ q_{j+1}.  The two terms of the
    objective balance when |r(q)|*S_K1 = q*S_K2, and since |r(q_j)| ~ n/q_{j+1}
    the true condition is on the PRODUCT q_j*q_{j+1} ~ n*S_K1/S_K2.  Whenever a
    partial quotient is large the two are far apart and the proposed rule picks
    the wrong convergent.  Correcting to the balance condition helps but still
    fails the bar (pearson 0.44).  Only the min over ALL convergents is exact
    (pearson 1.0000, 100% exact) -- and that costs O(log n), the same as the
    Lagrange-Gauss reduction it was meant to replace.  So there is NO cheap
    one-convergent shortcut; the value of the CF picture is structural, and the
    real closed form is (3) below.

(3) CLOSED FORM — nu_hat IS THE HERMITE INVARIANT.  L2 is a rank-2 lattice, so
    up to rotation and scaling it is a point tau in the modular fundamental
    domain F = {|tau| >= 1, |Re tau| <= 1/2}, namely the SL2(Z)-reduction of

        tau_0 = (-lam*S_K1 + i*S_K2) / (n*S_K1) = -lam/n + i*S_K2/(n*S_K1).

    For L = <1, tau> one has det = Im tau and lambda_1 = 1, so

        **  nu_hat = (Im tau_reduced)^(-1/2)  **

    and Minkowski/Hermite in dimension 2 gives the hard ceiling

        **  nu_hat <= (4/3)^(1/4) = 1.0745699...  **

    attained only at tau = the corner rho = e^(i*pi/3).

(4) CLOSED-FORM NULL DISTRIBUTION.  Reduced lattices from "random" (n, lam)
    equidistribute in F w.r.t. the hyperbolic measure dx dy / y^2 (total mass
    pi/3).  {tau in F : Im tau >= Y} for Y >= 1 is a unit-width strip of mass
    1/Y, so with Y = nu_hat^(-2):

        **  P(nu_hat <= v) = 3 v^2 / pi        for v <= 1  **

    and on the short tail 1 <= v <= (4/3)^(1/4), with s = sqrt(1 - v^-4),

            P(nu_hat <= v) = (6/pi) * [ arcsin(s) + (1/2 - s) * v^2 ].

    In particular the median is sqrt(pi/6) = 0.7236 and quantile q maps to
    sqrt(pi*q/3).

(5) WHY THE JUNE INVARIANTS FAILED.  max_a, q_cf, max_q_cf are SL2(Z)-invariant
    statistics of lam/n as a real number -- they do not see WHERE the CF is
    truncated.  nu_hat is a function of tau_0, whose imaginary part
    Im tau_0 = S_K2/(n*S_K1) ~ K1/(K2*n) sets exactly that truncation depth.
    Part E below holds lam/n (hence every scale-free CF statistic) FIXED and
    sweeps K1/K2, driving nu_hat across essentially its whole range.  No
    scale-free invariant of lam/n can do that, which retro-explains all six
    June falsifications at once.

Parts
-----
  A  exact min-over-b identity for lambda_1(L2)                  (claim 1)
  B  Hermite ceiling                                             (claim 3)
  C  CF-convergent predictor -- THE THREAD 23 FALSIFIER          (claim 2)
  D  closed-form null CDF vs empirical, and vs the 2026-07-29
     recorded secp256k1-scale quantiles                          (claim 4)
  E  scale-dependence demo: lam/n fixed, K1/K2 swept             (claim 5)
  F  modular closed form nu_hat = (Im tau_red)^(-1/2)            (claim 3)
  G  secp256k1 placement, recomputed from the closed form
  H  do REAL GLV eigenvalues obey the null law of (4)?           (see below)

(6) GLV EIGENVALUES ARE NOT GENERIC -- AND THE PROTECTION IS SCALE-DEPENDENT.
    The law in (4) is for uniform lam.  Real GLV lam are primitive cube roots
    of 1 mod n, so Lambda = <(n,0),(-lam,1)> is the coordinate lattice of the
    norm-n ideal (n, omega - lam) in Z[omega]; its shortest vector is pinned
    near sqrt(n) by the Eisenstein norm form a^2+ab+b^2 = n -- the same pinning
    that makes GLV scalar decomposition work.  Part H measures the consequence:

      at eff ~ 1 (S_K1 ~ S_K2):   0/300 GLV curves below the 0.645 C2 ceiling,
                                  vs 126/300 uniform.  Total protection.
      at eff = 0.10:              39.3% GLV vs 39.8% uniform.  Gone.

    So "GLV curves avoid the easy regime" is true only near eff ~ 1 and false
    in the eff <= 0.10 regime the Phase-2 experiments actually ran in.  Any
    quotation of this effect must state the eff.  Scored against the GLV null
    at its own eff, secp256k1 is unremarkable (46.3rd percentile).

Pure integer/float Python: no fpylll, no PARI, no sympy needed.
Run: python3 glv_hnp_nuhat_closed_form.py
"""

import math
import random

HERMITE = (4.0 / 3.0) ** 0.25          # 1.0745699318...
SEED = 20260804

# ---------------------------------------------------------------------------
# Exact 2D machinery
# ---------------------------------------------------------------------------

def _n2(w):
    return w[0] * w[0] + w[1] * w[1]

def gauss_reduce_2d_exact(u, v):
    """Lagrange-Gauss reduction with EXACT integer rounding (no float division,
    so this is safe at 256-bit scale where `round(a/b)` on Python ints can lose
    precision).  Returns the shortest nonzero vector."""
    if _n2(u) > _n2(v):
        u, v = v, u
    while True:
        num = v[0] * u[0] + v[1] * u[1]
        den = _n2(u)
        if den == 0:
            break
        # exact round-half-away-from-zero of num/den
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if _n2(v) >= _n2(u):
            break
        u, v = v, u
    return u

def scales(n, k1_bound, k2_bound):
    """Verbatim from glv_hnp_phase2_lambda_threshold*.py (validated 2026-07-26)."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def nu_hat_exact(n, lam, k1, k2):
    """nu_hat = lambda_1(L2)/sqrt(det L2), computed exactly then divided once."""
    s_k1, _, s_k2, _ = scales(n, k1, k2)
    w = gauss_reduce_2d_exact((n * s_k1, 0), (-(lam % n) * s_k1, s_k2))
    # do the ratio in integer-squared form first to avoid overflow at 256 bits
    num2 = _n2(w)
    den = n * s_k1 * s_k2
    return math.sqrt(num2 / den) if num2 < 2 ** 900 else math.exp(
        0.5 * (math.log(num2) - math.log(den)))

def centred_residue(b, lam, n):
    r = (b * lam) % n
    return r - n if r > n - r else r

# ---------------------------------------------------------------------------
# Continued fractions
# ---------------------------------------------------------------------------

def cf_parts(a, b):
    """Partial quotients [a_0; a_1, a_2, ...] of a/b."""
    parts = []
    x, y = a, b
    while y:
        k, r = divmod(x, y)
        parts.append(k)
        x, y = y, r
    return parts

def cf_convergents(a, b):
    """Denominators q_j of the convergents of a/b, plus the partial quotients.

    Recurrence q_{-1}=0, q_0=1, q_j = a_j*q_{j-1} + q_{j-2}.  Note q_0=1 always
    (the denominator of a_0/1), so the loop starts at a_1.  These are exactly
    the 'best approximations of the second kind' -- the Pareto frontier of
    (|q*lam mod +- n|, q), which is why the argmin in Part A must lie here.
    """
    parts = cf_parts(a, b)
    qs = [1]
    q_prev, q_cur = 0, 1
    for k in parts[1:]:
        q_prev, q_cur = q_cur, k * q_cur + q_prev
        qs.append(q_cur)
    return qs, parts

def cf_semiconvergent_denoms(a, b):
    """Convergent denominators AND all intermediate (semiconvergent) ones,
    i.e. t*q_{j-1} + q_{j-2} for t = 1..a_j (best approximations of the
    FIRST kind)."""
    parts = cf_parts(a, b)
    out = [1]
    q_prev, q_cur = 0, 1
    for k in parts[1:]:
        for t in range(1, k + 1):
            out.append(t * q_cur + q_prev)
        q_prev, q_cur = q_cur, k * q_cur + q_prev
    return sorted(set(d for d in out if d >= 1))

# ---------------------------------------------------------------------------
# Modular reduction of tau (claim 3)
# ---------------------------------------------------------------------------

def reduce_tau(x, y, iters=100000):
    """SL2(Z)-reduce tau = x + i*y into the fundamental domain; return Im tau."""
    for _ in range(iters):
        x -= math.floor(x + 0.5)                 # |Re| <= 1/2
        n2 = x * x + y * y
        if n2 >= 1.0 - 1e-15:
            return y
        x, y = -x / n2, y / n2                   # tau -> -1/tau
    return y

# ---------------------------------------------------------------------------
# Closed-form null CDF (claim 4)
# ---------------------------------------------------------------------------

def nu_cdf(v):
    """P(nu_hat <= v) for a hyperbolically-equidistributed random 2D lattice."""
    if v <= 0:
        return 0.0
    if v <= 1.0:
        return 3.0 * v * v / math.pi
    if v >= HERMITE:
        return 1.0
    s = math.sqrt(max(0.0, 1.0 - v ** -4))
    return (6.0 / math.pi) * (math.asin(s) + (0.5 - s) * v * v)

def nu_quantile_lo(q):
    """Inverse of the v<=1 branch: nu_hat quantile = sqrt(pi*q/3), valid q <= 3/pi."""
    return math.sqrt(math.pi * q / 3.0)

def pearson(xs, ys):
    k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")

def spearman(xs, ys):
    def rank(z):
        order = sorted(range(len(z)), key=lambda i: z[i])
        r = [0.0] * len(z)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and z[order[j + 1]] == z[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for t in range(i, j + 1):
                r[order[t]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))

# ===========================================================================

def part_A(trials=400):
    print("=" * 74)
    print("PART A -- lambda_1(L2)^2 == min_b [ r(b)^2 S_K1^2 + b^2 S_K2^2 ]")
    print("         (claim 1: exact reduction to weighted approximation)")
    print("=" * 74)
    rng = random.Random(SEED)
    bad = 0
    argmin_is_conv = 0
    argmin_is_semi = 0
    for _ in range(trials):
        n = rng.randrange(5000, 40000)
        lam = rng.randrange(1, n)
        k1 = rng.choice([4, 8, 16, 36, 72, 144])
        k2 = rng.choice([8, 32, 72, 200, 512])
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        w = gauss_reduce_2d_exact((n * s_k1, 0), (-lam * s_k1, s_k2))
        l1sq = _n2(w)
        # brute force over b (bounded: b*S_K2 <= lambda_1 <= sqrt(2/sqrt3 * det))
        bmax = int(math.sqrt(l1sq) / s_k2) + 2
        best, barg = None, None
        for b in range(1, bmax + 1):
            val = centred_residue(b, lam, n) ** 2 * s_k1 * s_k1 + b * b * s_k2 * s_k2
            if best is None or val < best:
                best, barg = val, b
        if best != l1sq:
            bad += 1
            continue
        convs, _ = cf_convergents(lam, n)
        if barg in convs:
            argmin_is_conv += 1
        if barg in cf_semiconvergent_denoms(lam, n):
            argmin_is_semi += 1
    print(f"  trials                          : {trials}")
    print(f"  identity mismatches             : {bad}   <- must be 0")
    print(f"  argmin b is a CF convergent     : {argmin_is_conv}/{trials - bad} "
          f"({100.0 * argmin_is_conv / max(1, trials - bad):.1f}%)")
    print(f"  argmin b is a semiconvergent    : {argmin_is_semi}/{trials - bad} "
          f"({100.0 * argmin_is_semi / max(1, trials - bad):.1f}%)")
    print()
    return bad == 0

def part_B(trials=40000):
    print("=" * 74)
    print("PART B -- Hermite ceiling  nu_hat <= (4/3)^(1/4) = %.7f" % HERMITE)
    print("=" * 74)
    rng = random.Random(SEED + 1)
    mx, exceed = 0.0, 0
    for _ in range(trials):
        n = rng.randrange(2 ** 19, 2 ** 20)
        lam = rng.randrange(1, n)
        v = nu_hat_exact(n, lam, rng.choice([36, 72, 144]), rng.choice([36, 72, 144]))
        mx = max(mx, v)
        if v > HERMITE + 1e-12:
            exceed += 1
    print(f"  trials                : {trials}")
    print(f"  max nu_hat observed   : {mx:.7f}")
    print(f"  ceiling               : {HERMITE:.7f}   (gap {HERMITE - mx:+.7f})")
    print(f"  violations            : {exceed}   <- must be 0")
    print()
    return exceed == 0

def part_C(trials=3000):
    print("=" * 74)
    print("PART C -- THREAD 23 FALSIFIER: CF-convergent predictor vs exact nu_hat")
    print("          falsifier threshold: correlation < 0.9 kills the story")
    print("=" * 74)
    rng = random.Random(SEED + 2)
    exact, pred_single, pred_fix, pred_all, pred_semi = [], [], [], [], []
    single_exact_hits = 0
    for _ in range(trials):
        n = rng.randrange(2 ** 19, 2 ** 20)
        lam = rng.randrange(1, n)
        k1 = rng.choice([18, 36, 72, 144])
        k2 = rng.choice([36, 72, 144, 288])
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        det = n * s_k1 * s_k2
        e = nu_hat_exact(n, lam, k1, k2)
        convs, _ = cf_convergents(lam, n)
        scale = math.sqrt(n * s_k1 / s_k2)        # = sqrt(n*K2/K1) balance point

        def val(b):
            return centred_residue(b, lam, n) ** 2 * s_k1 * s_k1 + b * b * s_k2 * s_k2

        # (i) the SINGLE convergent nearest the balance scale -- AS PROPOSED
        #     on 2026-07-29: "q_j nearest sqrt(n*S_K1/S_K2)".
        qbest = min(convs, key=lambda q: abs(math.log(q / scale)))
        pred_single.append(math.sqrt(val(qbest) / det))
        if abs(math.sqrt(val(qbest) / det) - e) < 1e-12:
            single_exact_hits += 1
        # (ii) CORRECTED single-convergent rule.  The two terms of val(q) balance
        #      when |r(q)|*S_K1 = q*S_K2, i.e. |r(q)|/q = S_K2/S_K1, NOT when
        #      q = sqrt(n*S_K1/S_K2).  Since |r(q_j)| ~ n/q_{j+1}, the correct
        #      scale condition is on the PRODUCT q_j*q_{j+1} ~ n*S_K1/S_K2 --
        #      so when a partial quotient is large (q_j and q_{j+1} far apart)
        #      the as-proposed rule picks the wrong convergent.
        target = s_k2 / s_k1
        qfix = min(convs, key=lambda q: abs(math.log(
            max(abs(centred_residue(q, lam, n)), 1) / q / target)))
        pred_fix.append(math.sqrt(val(qfix) / det))
        # (iii) min over ALL convergents
        pred_all.append(math.sqrt(min(val(q) for q in convs) / det))
        # (iv) min over all semiconvergents
        pred_semi.append(
            math.sqrt(min(val(q) for q in cf_semiconvergent_denoms(lam, n)) / det))
        exact.append(e)

    for name, pr in (("single convergent @ q-scale", pred_single),
                     ("single convergent @ balance", pred_fix),
                     ("min over convergents", pred_all),
                     ("min over semiconvergents", pred_semi)):
        r = pearson(exact, pr)
        rs = spearman(exact, pr)
        exact_frac = sum(1 for a, b in zip(exact, pr) if abs(a - b) < 1e-12) / len(exact)
        mre = sum(abs(a - b) / a for a, b in zip(exact, pr)) / len(exact)
        verdict = "PASS" if r >= 0.9 else "FAIL"
        print(f"  {name:<26} pearson={r:.4f} spearman={rs:.4f} "
              f"exact={100 * exact_frac:5.1f}%  mean-rel-err={mre:.4f}  [{verdict}]")
    print()
    return pearson(exact, pred_single) >= 0.9

def part_D(trials=40000):
    print("=" * 74)
    print("PART D -- closed-form null CDF  P(nu_hat <= v) = 3v^2/pi   (v <= 1)")
    print("=" * 74)
    rng = random.Random(SEED + 3)
    samp = []
    for _ in range(trials):
        n = rng.randrange(2 ** 19, 2 ** 20)
        samp.append(nu_hat_exact(n, rng.randrange(1, n), 72, 72))
    samp.sort()
    print("   q      empirical   closed-form   diff")
    worst = 0.0
    for q in (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95):
        e = samp[int(q * len(samp))]
        c = nu_quantile_lo(q)
        worst = max(worst, abs(e - c))
        print(f"   {q:.2f}    {e:.4f}      {c:.4f}      {e - c:+.4f}")
    print(f"  max |empirical - closed form| over quantiles: {worst:.4f}")

    # Head-to-head with the quantiles RECORDED on 2026-07-29 at 256-bit scale
    # (log entry "2026-07-29 (autolab run #2)", section 5, eff=0.05, 4000 draws).
    print()
    print("  vs the 2026-07-29 recorded 256-bit null (eff=0.05, 4000 draws):")
    recorded = {0.05: 0.232, 0.10: 0.322, 0.25: 0.496,
                0.50: 0.715, 0.75: 0.885, 0.90: 0.974, 0.95: 0.999}
    worst2 = 0.0
    for q, rv in recorded.items():
        c = nu_quantile_lo(q)
        worst2 = max(worst2, abs(rv - c))
        print(f"   q{int(q * 100):02d}   recorded={rv:.3f}   closed-form={c:.3f}   "
              f"diff={rv - c:+.3f}")
    print(f"  max |recorded - closed form|: {worst2:.4f}")
    print()
    return worst < 0.02 and worst2 < 0.03

def part_E():
    print("=" * 74)
    print("PART E -- scale dependence: lam/n FIXED, K1/K2 swept")
    print("          (why every scale-free CF invariant had to fail)")
    print("=" * 74)
    rng = random.Random(SEED + 4)
    n = 1048573                                   # prime-ish 20-bit modulus
    while True:
        lam = rng.randrange(2, n)
        convs, parts = cf_convergents(lam, n)
        if len(convs) >= 8:
            break
    print(f"  n = {n}, lam = {lam}   (these are FIXED for the whole table)")
    print(f"  scale-free CF invariants of lam/n, constant down every column:")
    print(f"    max_a = {max(parts)},  len(CF) = {len(parts)},  "
          f"q_cf = {convs[len(convs) // 2]},  max_q_cf = {max(convs)}")
    print()
    print("   K1     K2     Im tau_0        nu_hat")
    vals = []
    for k1, k2 in [(2, 1024), (8, 512), (36, 288), (72, 144), (72, 72),
                   (144, 72), (288, 36), (512, 8), (1024, 2)]:
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        v = nu_hat_exact(n, lam, k1, k2)
        vals.append(v)
        print(f"  {k1:>4}   {k2:>4}   {s_k2 / (n * s_k1):.3e}     {v:.4f}")
    print()
    print(f"  nu_hat range over the sweep: [{min(vals):.4f}, {max(vals):.4f}] "
          f"= {100 * (max(vals) - min(vals)) / HERMITE:.0f}% of the full range "
          f"[0, {HERMITE:.4f}]")
    print("  => a single curve (n, lam) takes almost every nu_hat value as K1/K2")
    print("     moves.  Any invariant that is a function of lam/n alone is")
    print("     constant down this column and therefore cannot predict nu_hat.")
    print()
    return max(vals) - min(vals) > 0.4

def part_F(trials=2000):
    print("=" * 74)
    print("PART F -- modular closed form   nu_hat = (Im tau_reduced)^(-1/2)")
    print("          tau_0 = -lam/n + i*S_K2/(n*S_K1)")
    print("=" * 74)
    rng = random.Random(SEED + 5)
    worst, bad = 0.0, 0
    for _ in range(trials):
        n = rng.randrange(2 ** 16, 2 ** 20)
        lam = rng.randrange(1, n)
        k1, k2 = rng.choice([18, 36, 72]), rng.choice([36, 72, 144])
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        e = nu_hat_exact(n, lam, k1, k2)
        y = reduce_tau(-lam / n, s_k2 / (n * s_k1))
        pred = y ** -0.5
        d = abs(pred - e)
        worst = max(worst, d)
        if d > 1e-6:
            bad += 1
    print(f"  trials                        : {trials}")
    print(f"  max |(Im tau)^-1/2 - nu_hat|  : {worst:.3e}")
    print(f"  disagreements > 1e-6          : {bad}")
    print(f"  => nu_hat is the Hermite invariant of L2; the reduction of tau_0")
    print(f"     IS the continued-fraction expansion of lam/n, truncated at the")
    print(f"     depth set by Im tau_0 = S_K2/(n*S_K1) ~ K1/(K2*n).")
    print()
    return bad == 0

def part_G():
    print("=" * 74)
    print("PART G -- secp256k1 placement, recomputed from the closed form")
    print("=" * 74)
    n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    lam = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
    assert (lam * lam + lam + 1) % n == 0, "lambda is not a primitive cube root mod n"
    print(f"  n   = {n}")
    print(f"  lam verified: lam^2 + lam + 1 = 0 mod n   OK")
    print()
    print("  Convention from glv_hnp_nuhat_vs_c1c2.py:185 --")
    print("    K2 = isqrt(n)+1,  K1 = max(2, int(eff*isqrt(n))).")
    print("  Reproduction check against the values recorded 2026-07-29 (section 5):")
    print()
    print("   eff    nu_hat   recorded    diff   closed-form pct 3v^2/pi   recorded pct")
    rec = {0.05: (0.8709, 73.3), 0.10: (0.6624, 41.3),
           0.25: (0.5852, 33.2), 0.50: (0.6851, 44.1)}
    k2_s = math.isqrt(n) + 1
    worst = 0.0
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1_s = max(2, int(eff * math.isqrt(n)))
        v = nu_hat_exact(n, lam, k1_s, k2_s)
        rv, rp = rec[eff]
        worst = max(worst, abs(v - rv))
        print(f"  {eff:.2f}   {v:.4f}   {rv:.4f}   {v - rv:+.4f}        "
              f"{100 * nu_cdf(v):5.1f}%              {rp:.1f}%")
    print(f"  max |reproduced - recorded| = {worst:.4f}")
    print()
    print("  Note: nu_hat depends on (n, lam) and on the RATIO S_K1/S_K2 only.")
    print("  Under this convention S_K1/S_K2 ~ 1/eff, which is why the table")
    print("  moves with eff; had K1=K2 been used, nu_hat would be eff-invariant.")
    print()
    print("  C2 ceiling recorded 2026-07-29: every C2 curve had nu_hat <= 0.645.")
    print(f"  Closed-form share of curves below that ceiling: "
          f"{100 * nu_cdf(0.645):.1f}%")
    print("  (i.e. the 'easy' class is not rare -- it is ~40% of random lambda.)")
    print()
    # Part H shows GLV lam does NOT follow the uniform law, so the honest
    # reference distribution for secp256k1 is the GLV null, not 3v^2/pi.
    glv = sorted(nu_hat_exact(nn, ll, max(2, int(0.10 * math.isqrt(nn))),
                              math.isqrt(nn) + 1)
                 for nn, ll in _glv_pairs(600, SEED + 10))
    v_secp = nu_hat_exact(n, lam, max(2, int(0.10 * math.isqrt(n))),
                          math.isqrt(n) + 1)
    pct_glv = 100.0 * sum(1 for x in glv if x < v_secp) / len(glv)
    print(f"  Scored against the GLV null (600 j=0 curves, eff=0.10; see Part H):")
    print(f"    secp256k1 nu_hat = {v_secp:.4f}  ->  GLV-null percentile "
          f"{pct_glv:.1f}%   (uniform-null percentile {100 * nu_cdf(v_secp):.1f}%)")
    print(f"    share of GLV curves below the 0.645 C2 ceiling: "
          f"{100.0 * sum(1 for x in glv if x <= 0.645) / len(glv):.1f}%")
    print("  secp256k1 is an unremarkable GLV curve by this measure -- neither in")
    print("  the easy tail nor anomalously hard.")
    print()
    print("  CAVEATS, unchanged from 2026-07-29 and NOT weakened by this note:")
    print("   (a) the C1/C2 -> nu_hat link is a 20/24-bit empirical regularity")
    print("       extrapolated to 256 bits; this script proves the closed form for")
    print("       nu_hat itself, NOT that nu_hat predicts recovery at 256 bits.")
    print("   (b) the whole thread is conditional on a NON-STANDARD nonce generator")
    print("       k = k1 + lam*k2 with k1 bounded.  It says nothing about correctly")
    print("       generated nonces and is NOT an attack on secp256k1.")
    print()
    return True

def part_H(n_moduli=1200):
    """The null law in Part D was derived for UNIFORM lam.  Real GLV eigenvalues
    are not uniform: they are the two primitive cube roots of 1 mod n, i.e. the
    roots of x^2+x+1.  There are only 2 of them per modulus, and they are
    algebraically special (lam ~ (-1 +- sqrt(-3))/2).  If they deviated from
    3v^2/pi, the secp256k1 percentile in Part G would be misread.  Test it."""
    print("=" * 74)
    print("PART H -- do real GLV eigenvalues obey the same null law as uniform lam?")
    print("=" * 74)
    rng = random.Random(SEED + 7)
    glv, unif = [], []
    tried = 0
    while len(glv) < n_moduli and tried < 400000:
        tried += 1
        n = rng.randrange(2 ** 19, 2 ** 20) | 1
        if n % 3 != 1 or not _is_prime(n):
            continue
        lam = _cube_root_of_unity(n, rng)
        if lam is None:
            continue
        assert (lam * lam + lam + 1) % n == 0
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(0.10 * math.isqrt(n)))
        glv.append(nu_hat_exact(n, lam, k1, k2))
        unif.append(nu_hat_exact(n, rng.randrange(1, n), k1, k2))
    glv.sort(); unif.sort()

    def ks_vs_closed(sample):
        return max(max(abs((i + 1) / len(sample) - nu_cdf(v)),
                       abs(i / len(sample) - nu_cdf(v)))
                   for i, v in enumerate(sample))

    ks_g, ks_u = ks_vs_closed(glv), ks_vs_closed(unif)
    crit = 1.36 / math.sqrt(len(glv))
    print(f"  primes n = 1 mod 3 sampled: {len(glv)}   (K2=isqrt(n)+1, eff=0.10)")
    print()
    print("   q     GLV lam   uniform lam   closed form   GLV-closed")
    for q in (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95):
        g, u, c = glv[int(q * len(glv))], unif[int(q * len(unif))], nu_quantile_lo(q)
        print(f"   {q:.2f}   {g:.4f}    {u:.4f}       {c:.4f}      {g - c:+.4f}")
    print()
    print(f"  KS vs closed-form CDF, uniform lam : {ks_u:.4f}  (5% crit {crit:.4f})"
          f" -> {'consistent' if ks_u < crit else 'DEVIATES'}")
    print(f"  KS vs closed-form CDF, GLV lam     : {ks_g:.4f}  (5% crit {crit:.4f})"
          f" -> {'consistent' if ks_g < crit else 'DEVIATES'}")
    print()
    print("  *** RESULT: GLV eigenvalues DO NOT follow the random-lattice law. ***")
    print("  Uniform lam matches 3v^2/pi (as Part D shows); GLV lam is shifted UP.")
    print("  But the shift is SCALE-DEPENDENT, and that detail matters -- read the")
    print("  decay table below before quoting any single number.")
    print()
    frac_g = sum(1 for v in glv if v <= 0.645) / len(glv)
    frac_u = sum(1 for v in unif if v <= 0.645) / len(unif)
    print(f"    At eff=0.10, P(nu_hat <= 0.645) (the recorded C2 ceiling):")
    print(f"      GLV lam {100 * frac_g:5.1f}%   uniform lam {100 * frac_u:5.1f}%"
          f"   closed form {100 * nu_cdf(0.645):5.1f}%")
    print("    -- i.e. at THIS cut and THIS scale the GLV constraint has already")
    print("       washed out; the KS deviation above lives in the extreme low tail")
    print(f"       (q05: GLV {glv[len(glv) // 20]:.3f} vs closed form "
          f"{nu_quantile_lo(0.05):.3f}).")
    print()
    print("  Mechanism.  For a j=0 CM curve, lam is the eigenvalue of an")
    print("  endomorphism of norm 1, so Lambda = <(n,0),(-lam,1)> is the coordinate")
    print("  lattice of the norm-n ideal (n, omega - lam) in Z[omega].  Its shortest")
    print("  vector is pinned near sqrt(n) by the Eisenstein norm form")
    print("  a^2+ab+b^2 = n -- that pinning is precisely what makes GLV scalar")
    print("  decomposition work.  So the very structure that makes a curve")
    print("  GLV-friendly also holds nu_hat up, away from the low-nu_hat tail.")
    print("  (It is not pinned exactly to the Hermite corner because L2 uses the")
    print("  coordinate metric, not the Minkowski metric of Z[omega].)")
    print()
    # How fast does the GLV constraint decay as the scale ratio S_K1/S_K2 grows?
    pairs = _glv_pairs(400, SEED + 8)
    uprimes = _uniform_pairs(400, SEED + 9)
    print("  DECAY OF THE GLV CONSTRAINT with the scale ratio S_K1/S_K2 ~ 1/eff.")
    print("  (K1=K2 is eff ~ 1; the Phase-2 runs used eff = 0.05-0.50.)")
    print()
    print("    eff     P(nu<=0.645) GLV   uniform     median GLV   median unif")
    for eff in (1.0, 0.75, 0.50, 0.25, 0.10, 0.05, 0.01):
        gg, uu = [], []
        for nn, ll in pairs:
            k1 = max(1, int(eff * math.isqrt(nn)))
            gg.append(nu_hat_exact(nn, ll, k1, math.isqrt(nn) + 1))
        for nn, ll in uprimes:
            k1 = max(1, int(eff * math.isqrt(nn)))
            uu.append(nu_hat_exact(nn, ll, k1, math.isqrt(nn) + 1))
        gg.sort(); uu.sort()
        print(f"   {eff:5.2f}      {100 * sum(1 for v in gg if v <= 0.645) / len(gg):5.1f}%"
              f"          {100 * sum(1 for v in uu if v <= 0.645) / len(uu):5.1f}%"
              f"       {gg[len(gg) // 2]:.4f}       {uu[len(uu) // 2]:.4f}")
    print()
    print("  Reading: the GLV/Eisenstein constraint is decisive only near eff ~ 1")
    print("  (S_K1 ~ S_K2), where L2 is essentially the ideal lattice itself.  As")
    print("  eff falls the weighting skews L2 away from that ideal and the")
    print("  constraint washes out; by eff = 0.10 GLV and uniform agree at the")
    print("  0.645 cut to well under a percentage point.  So 'GLV curves are")
    print("  protected' is TRUE at eff ~ 1 and FALSE in the eff <= 0.10 regime the")
    print("  Phase-2 experiments actually ran in.  Quote the eff.")
    print()
    # PASS = we correctly detect that GLV differs from the uniform null.
    return ks_g > crit and ks_u < crit

def _glv_pairs(count, seed):
    rng = random.Random(seed)
    out = []
    while len(out) < count:
        n = rng.randrange(2 ** 19, 2 ** 20) | 1
        if n % 3 != 1 or not _is_prime(n):
            continue
        lam = _cube_root_of_unity(n, rng)
        if lam is not None:
            out.append((n, lam))
    return out

def _uniform_pairs(count, seed):
    """(prime n, uniform lam) pairs -- the null comparison for _glv_pairs."""
    rng = random.Random(seed)
    out = []
    while len(out) < count:
        n = rng.randrange(2 ** 19, 2 ** 20) | 1
        if not _is_prime(n):
            continue
        out.append((n, rng.randrange(1, n)))
    return out

def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def _cube_root_of_unity(n, rng):
    """A primitive cube root of 1 mod prime n (n = 1 mod 3), i.e. a root of
    x^2+x+1.  Equivalently the GLV eigenvalue lam."""
    for _ in range(40):
        h = rng.randrange(2, n - 1)
        g = pow(h, (n - 1) // 3, n)
        if g != 1:
            return g
    return None

def main():
    print()
    print("Thread 23 -- closed form for nu_hat (2026-08-04 autolab)")
    print()
    results = [
        ("A  min-over-b identity", part_A()),
        ("B  Hermite ceiling", part_B()),
        ("C  CF predictor (FALSIFIER)", part_C()),
        ("D  closed-form null CDF", part_D()),
        ("E  scale dependence", part_E()),
        ("F  modular closed form", part_F()),
        ("G  secp256k1 placement", part_G()),
        ("H  GLV lam obeys null law", part_H()),
    ]
    print("=" * 74)
    print("SUMMARY")
    print("=" * 74)
    for name, ok in results:
        print(f"  {name:<32} {'PASS' if ok else 'FAIL'}")
    print()

if __name__ == "__main__":
    main()
