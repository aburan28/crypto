"""
GLV-HNP Phase 2 -- Thread 23: closed form for lambda_1(L2), and why eight
scale-free invariants had to fail.

Background
----------
Thread 20 (2026-07-29, commit e845207) established nu_hat as the first working
Phase 2 separator (AUC 0.935 against the June C1/C2 classes):

    L2     = < (n*S_K1, 0), (-lam*S_K1, S_K2) >      (the non-planted 2D sublattice)
    det L2 = n * S_K1 * S_K2                          (independent of lam)
    nu_hat = lambda_1(L2) / sqrt(det L2)

but nu_hat was *computed* (Lagrange-Gauss), not *understood*.  The 2026-07-29
next-step proposal: show lambda_1(L2) is a continued-fraction quantity of lam/n
evaluated at the scale sqrt(n*K2/K1), and check the prediction reproduces the
measured nu_hat.  Falsifier: correlation < 0.9 between CF-predicted and computed
nu_hat.

The derivation being tested
---------------------------
Every vector of L2 is  b*u + a*v = ( S_K1*(b*n - a*lam),  a*S_K2 ).  For fixed a
the optimal b makes the first coordinate S_K1*r(a) with

    r(a) = min_b |a*lam - b*n|      (centered residue)

so lambda_1(L2)^2 = min_{a >= 0} [ S_K1^2 * r(a)^2 + S_K2^2 * a^2 ].         (*)

The pairs (a, r(a)) that are Pareto-minimal are exactly the classical best
approximations to lam/n, i.e. the CF convergent denominators q_j (possibly with
semiconvergents).  Writing c = S_K2/S_K1 and normalising

    x = a * sqrt(c/n),   y = r(a) / sqrt(n*c)   =>   nu_hat^2 = x^2 + y^2,

the balance point x = y is a* = sqrt(n/c) = sqrt(n*K2/K1), and on a convergent
q_j the classical r(q_j) ~ n/q_{j+1} gives

    x*y = a*r(a)/n ~ q_j/q_{j+1} ~ 1/a_{j+1}.

Hence  nu_hat^2 >= 2*x*y ~ 2/a_{j+1}:  nu_hat is small exactly when lam/n has a
LARGE partial quotient at the convergent nearest the scale a*.

That is the retro-explanation for June.  q_cf, max_q_cf and max_a are scale-free
-- they ask whether a large partial quotient exists ANYWHERE in the expansion.
The quantity that matters is the partial quotient WHERE the expansion crosses
a* = sqrt(n*K2/K1), which moves with K1/K2.  A global max_a can be huge at a
convergent far from a* and buy nothing.

Experiments
-----------
E1  Structure: brute-force argmin of (*) over all a < A_MAX on small n; check the
    minimiser is a convergent (or semiconvergent) denominator of lam/n.
E2  Closed form: CF-only predictor vs exact lambda_1 on a 20-bit curve sample.
    Exact-match rate + Pearson/Spearman.  This is the stated falsifier.
E3  Single-convergent rule: use only the one convergent nearest a*.
E4  Mechanism: nu_hat vs LOCAL partial quotient at the scale, head-to-head with
    the June-falsified GLOBAL max_a; plus a same-curve scale sweep showing
    nu_hat move as eff changes while the curve (and max_a) is fixed.
E5  secp256k1: exact integer nu_hat, the achieving convergent, and a check of
    the float Lagrange-Gauss value reported 2026-07-29.

All arithmetic here is exact integer arithmetic; the incumbent
`rival_sublattice_nu` uses float Lagrange-Gauss, which E5 checks at 256 bits.

Run: python3 glv_hnp_nuhat_cf_closed_form.py
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
    "_t20a", os.path.join(_HERE, "glv_hnp_nuhat_base.py"))
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

eisenstein_decompose = _t20a.eisenstein_decompose
j0_traces = _t20a.j0_traces
glv_eigenvalues = _t20a.glv_eigenvalues
mu_of = _t20a.mu_of
scales = _t20a.scales
rival_sublattice_nu = _t20a.rival_sublattice_nu

# ---------------------------------------------------------------------------
# Continued fractions of lam/n
# ---------------------------------------------------------------------------

def cf_convergents(lam, n):
    """
    Convergents of lam/n, as (q_j, r_j, a_j) with

        r_j = |q_j*lam - p_j*n|   (exact int, the centered-residue sequence)

    via the classical Euclidean recurrence seeded at
    (q_{-1}, r_{-1}) = (0, n) and (q_0, r_0) = (1, lam):

        a_j     = r_{j-1} // r_j
        r_{j+1} = r_{j-1} - a_j*r_j
        q_{j+1} = q_{j-1} + a_j*q_j

    a_j is the partial quotient consumed to REACH entry j+1.  The j = -1 term
    (q=0, r=n) is included so the a=0 lattice vector (n*S_K1, 0) is covered.
    """
    out = [(0, n, 0)]
    q_prev, r_prev = 0, n
    q_cur, r_cur = 1, lam
    while r_cur:
        a = r_prev // r_cur
        out.append((q_cur, r_cur, a))
        q_prev, q_cur = q_cur, q_prev + a * q_cur
        r_prev, r_cur = r_cur, r_prev - a * r_cur
    return out

def cf_semiconvergents(lam, n):
    """
    Convergents plus the intermediate fractions
        q = q_{j-1} + t*q_j,  r = r_{j-1} - t*r_j,   1 <= t < a_j.
    Both coordinates stay positive, so these are the remaining vertices of the
    Klein polygon and the only other candidates for the minimum of (*).
    """
    conv = cf_convergents(lam, n)
    out = list(conv)
    for j in range(1, len(conv)):
        q_pp, r_pp, _ = conv[j - 1]
        q_p, r_p, a_j = conv[j]
        for t in range(1, a_j):
            out.append((q_pp + t * q_p, r_pp - t * r_p, a_j))
    return out

# ---------------------------------------------------------------------------
# lambda_1(L2): exact, via the quadratic form (*)
# ---------------------------------------------------------------------------

def centered(q, lam, n):
    """r(q) = min_b |q*lam - b*n|, the exact centered residue."""
    t = (q * lam) % n
    return min(t, n - t)

def form_val(s_k1, s_k2, a, r):
    """Exact squared norm S_K1^2*r^2 + S_K2^2*a^2."""
    return s_k1 * s_k1 * r * r + s_k2 * s_k2 * a * a

def lambda1_exact_bruteforce(n, lam, s_k1, s_k2, a_max):
    """Minimise (*) over every a in [0, a_max].  Reference implementation."""
    best, best_a = None, None
    for a in range(a_max + 1):
        t = (a * lam) % n
        r = min(t, n - t)
        val = form_val(s_k1, s_k2, a, r)
        if a == 0:
            val = form_val(s_k1, s_k2, 0, n)
        if best is None or val < best:
            best, best_a = val, a
    return best, best_a

def lambda1_from_cf(n, lam, s_k1, s_k2, use_semi=False):
    """
    Minimise (*) restricted to CF convergent denominators (optionally including
    semiconvergents).  r is recomputed as the exact centered residue rather than
    taken from the recurrence: the recurrence's r_j equals |q_j*lam - p_j*n| for
    the convergent's own p_j, which for j=0 (q=1, r=lam) is NOT the minimum over
    b when lam > n/2.
    """
    cands = cf_semiconvergents(lam, n) if use_semi else cf_convergents(lam, n)
    best, best_a, best_r, best_pq = None, None, None, None
    for (q, _r_rec, a_pq) in cands:
        r = n if q == 0 else centered(q, lam, n)
        val = form_val(s_k1, s_k2, q, r)
        if best is None or val < best:
            best, best_a, best_r, best_pq = val, q, r, a_pq
    return best, best_a, best_r, best_pq

def nu_hat_exact(n, lam, k1, k2, use_semi=True):
    """Exact-integer nu_hat via the CF closed form."""
    s_k1, _, s_k2, _ = scales(n, k1, k2)
    val, a, r, pq = lambda1_from_cf(n, lam, s_k1, s_k2, use_semi)
    det = n * s_k1 * s_k2
    return math.sqrt(val / det), a, r, pq

def scale_point(n, k1, k2):
    """Balance scale a* = sqrt(n*S_K1/S_K2)."""
    s_k1, _, s_k2, _ = scales(n, k1, k2)
    return math.sqrt(n * s_k1 / s_k2)

def local_partial_quotient(lam, n, k1, k2):
    """
    The partial quotient a_{j+1} at the convergent q_j nearest (below) the scale
    a*.  This is the LOCAL, scale-dependent analogue of the June max_a.
    """
    conv = cf_convergents(lam, n)
    a_star = scale_point(n, k1, k2)
    best_j, best_d = 0, None
    for j, (q, r, _) in enumerate(conv):
        d = abs(math.log((q if q else 0.5) / a_star))
        if best_d is None or d < best_d:
            best_j, best_d = j, d
    # entry j stores a_j = r_{j-1}//r_j, the quotient that governs how far r
    # drops at this convergent -- i.e. the one attached to q_j itself.
    return conv[best_j][2], best_j, len(conv)

def global_max_a(lam, n):
    """The June scale-free invariant: largest partial quotient anywhere."""
    return max((a for (_, _, a) in cf_convergents(lam, n)), default=0)

# ---------------------------------------------------------------------------
# Curve sampling
# ---------------------------------------------------------------------------

def harvest_curves(bits_lo, bits_hi, want, max_primes=400000):
    """j=0 CM curves with prime n = 1 mod 3 (so a GLV eigenvalue exists)."""
    out, seen = [], set()
    p = int(sympy.nextprime(bits_lo))
    scanned = 0
    while p < bits_hi and scanned < max_primes and len(out) < want:
        scanned += 1
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

def k1k2(n, eff):
    return max(2, int(eff * math.sqrt(n))), math.isqrt(n) + 1

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def pearson(xs, ys):
    k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else 0.0

def rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for t in range(i, j + 1):
            r[order[t]] = avg
        i = j + 1
    return r

def spearman(xs, ys):
    return pearson(rank(xs), rank(ys))

# ---------------------------------------------------------------------------
# E1 -- structure: is the minimiser always a convergent?
# ---------------------------------------------------------------------------

def exp1(rng):
    print("\n" + "=" * 78)
    print("E1  STRUCTURE -- brute-force argmin of (*) vs CF convergents")
    print("=" * 78)
    print("  For each curve, minimise S_K1^2*r(a)^2 + S_K2^2*a^2 over EVERY")
    print("  a in [0, A_MAX] and check the minimiser is a CF denominator.")

    curves = harvest_curves(2 ** 13, 2 ** 15, 40)
    print(f"  {len(curves)} curves, 14-15 bit (brute force is O(n) per curve)")

    for eff in (0.05, 0.10, 0.25):
        conv_hit = semi_hit = total = 0
        val_match_c = val_match_s = 0
        for c in curves:
            n, lam = c['n'], c['lam']
            k1, k2 = k1k2(n, eff)
            s_k1, _, s_k2, _ = scales(n, k1, k2)
            # a beyond sqrt(det)/S_K2 can never win: S_K2^2 a^2 > det.
            a_max = min(n - 1, int(math.sqrt(n * s_k1 / s_k2)) * 8 + 50)
            bf_val, bf_a = lambda1_exact_bruteforce(n, lam, s_k1, s_k2, a_max)
            qs_c = {q for (q, _, _) in cf_convergents(lam, n)}
            qs_s = {q for (q, _, _) in cf_semiconvergents(lam, n)}
            cval, _, _, _ = lambda1_from_cf(n, lam, s_k1, s_k2, False)
            sval, _, _, _ = lambda1_from_cf(n, lam, s_k1, s_k2, True)
            total += 1
            conv_hit += bf_a in qs_c
            semi_hit += bf_a in qs_s
            val_match_c += cval == bf_val
            val_match_s += sval == bf_val
        print(f"\n  eff={eff}")
        print(f"    argmin is a CF convergent      : {conv_hit}/{total}")
        print(f"    argmin is a semiconvergent     : {semi_hit}/{total}")
        print(f"    convergents reproduce min VALUE: {val_match_c}/{total}")
        print(f"    +semiconvergents reproduce it  : {val_match_s}/{total}")

# ---------------------------------------------------------------------------
# E2 -- the stated falsifier
# ---------------------------------------------------------------------------

def exp2(rng, curves):
    print("\n" + "=" * 78)
    print("E2  CLOSED FORM -- CF predictor vs computed nu_hat  [FALSIFIER]")
    print("=" * 78)
    print("  Falsifier (2026-07-29): corr(CF-predicted, computed) < 0.9.")

    for eff in (0.05, 0.10, 0.25, 0.50):
        exact, incumbent, conv_only = [], [], []
        agree_s = agree_c = 0
        for c in curves:
            n, lam = c['n'], c['lam']
            k1, k2 = k1k2(n, eff)
            s_k1, _, s_k2, _ = scales(n, k1, k2)
            det = n * s_k1 * s_k2
            vs, _, _, _ = lambda1_from_cf(n, lam, s_k1, s_k2, True)
            vc, _, _, _ = lambda1_from_cf(n, lam, s_k1, s_k2, False)
            _, _, nh_inc = rival_sublattice_nu(n, lam, k1, k2)
            exact.append(math.sqrt(vs / det))
            conv_only.append(math.sqrt(vc / det))
            incumbent.append(nh_inc)
            agree_s += abs(math.sqrt(vs / det) - nh_inc) < 1e-9
            agree_c += vc == vs
        r_p = pearson(exact, incumbent)
        r_s = spearman(exact, incumbent)
        print(f"\n  eff={eff}  ({len(curves)} curves)")
        print(f"    pearson (CF+semi, Lagrange-Gauss) = {r_p:.6f}")
        print(f"    spearman                          = {r_s:.6f}")
        print(f"    bitwise-identical nu_hat          = {agree_s}/{len(curves)}")
        print(f"    convergents alone suffice         = {agree_c}/{len(curves)}")
        print(f"    verdict: {'CONFIRMED' if r_p >= 0.9 else 'FALSIFIED'}")

# ---------------------------------------------------------------------------
# E3 -- single nearest convergent
# ---------------------------------------------------------------------------

def exp3(curves):
    print("\n" + "=" * 78)
    print("E3  SINGLE-CONVERGENT RULE -- only the convergent nearest a*")
    print("=" * 78)
    print("  a* = sqrt(n*S_K1/S_K2) = sqrt(n*K2/K1).  If one convergent suffices,")
    print("  nu_hat is a closed form, not a minimisation.")

    for eff in (0.05, 0.10, 0.25):
        rs_one, rs_full, hit = [], [], 0
        for c in curves:
            n, lam = c['n'], c['lam']
            k1, k2 = k1k2(n, eff)
            s_k1, _, s_k2, _ = scales(n, k1, k2)
            det = n * s_k1 * s_k2
            conv = cf_convergents(lam, n)
            a_star = scale_point(n, k1, k2)
            j_near = min(range(len(conv)),
                         key=lambda j: abs(math.log((conv[j][0] or 0.5) / a_star)))
            q = conv[j_near][0]
            r = n if q == 0 else centered(q, lam, n)
            v_one = form_val(s_k1, s_k2, q, r)
            v_full, _, _, _ = lambda1_from_cf(n, lam, s_k1, s_k2, True)
            rs_one.append(math.sqrt(v_one / det))
            rs_full.append(math.sqrt(v_full / det))
            hit += v_one == v_full
        print(f"\n  eff={eff}")
        print(f"    nearest convergent is the minimiser: {hit}/{len(curves)}"
              f" ({100.0*hit/len(curves):.1f}%)")
        ratios = [a / b for a, b in zip(rs_one, rs_full) if b > 0]
        print(f"    pearson(one-convergent, exact)     = {pearson(rs_one, rs_full):.4f}")
        print(f"    mean ratio one/exact               = "
              f"{sum(ratios)/len(ratios):.4f}" if ratios else "    (undefined)")

# ---------------------------------------------------------------------------
# E4 -- mechanism: local vs global partial quotient
# ---------------------------------------------------------------------------

def exp4(curves):
    print("\n" + "=" * 78)
    print("E4  MECHANISM -- local (scale-dependent) vs global (June) partial quotient")
    print("=" * 78)
    print("  Prediction: nu_hat^2 >~ 2/a_local.  The June invariant max_a is the")
    print("  GLOBAL max and should carry far less information.")

    for eff in (0.05, 0.10, 0.25):
        nus, loc, glob = [], [], []
        for c in curves:
            n, lam = c['n'], c['lam']
            k1, k2 = k1k2(n, eff)
            nh, _, _, _ = nu_hat_exact(n, lam, k1, k2)
            a_loc, _, _ = local_partial_quotient(lam, n, k1, k2)
            nus.append(nh)
            loc.append(1.0 / math.sqrt(a_loc) if a_loc > 0 else 1.0)
            glob.append(1.0 / math.sqrt(global_max_a(lam, n)))
        print(f"\n  eff={eff}")
        print(f"    spearman(nu_hat, 1/sqrt(a_LOCAL))  = {spearman(nus, loc):+.3f}")
        print(f"    spearman(nu_hat, 1/sqrt(max_a))    = {spearman(nus, glob):+.3f}"
              "   <- June invariant")

    print("\n  --- Same-curve scale sweep: curve FIXED, only eff varies ---")
    print("  max_a is a property of lam/n alone, so it cannot move.  If nu_hat")
    print("  moves, every scale-free invariant is provably dead.")
    for c in curves[:6]:
        n, lam = c['n'], c['lam']
        vals = []
        for eff in (0.02, 0.05, 0.10, 0.25, 0.50, 1.00):
            k1, k2 = k1k2(n, eff)
            nh, _, _, _ = nu_hat_exact(n, lam, k1, k2)
            vals.append(nh)
        print(f"    n={n:<9} max_a={global_max_a(lam, n):<6} mu={mu_of(lam, n):.3f}  "
              "nu_hat: " + " ".join(f"{v:.3f}" for v in vals))
    print("    (eff = 0.02 0.05 0.10 0.25 0.50 1.00)")

    print("\n  --- E4b: partial quotient at the ACHIEVING convergent ---")
    print("  The mechanism predicts nu_hat^2 ~ x^2+y^2 with x*y = q*r/n at the")
    print("  convergent that ACHIEVES the min -- not the one nearest a*.")
    for eff in (0.05, 0.10, 0.25):
        nus, ach, near, tight = [], [], [], []
        for c in curves:
            n, lam = c['n'], c['lam']
            k1, k2 = k1k2(n, eff)
            s_k1, _, s_k2, _ = scales(n, k1, k2)
            det = n * s_k1 * s_k2
            v, q, r, pq = lambda1_from_cf(n, lam, s_k1, s_k2, True)
            nh2 = v / det
            nus.append(math.sqrt(nh2))
            ach.append(1.0 / math.sqrt(pq) if pq > 0 else 1.0)
            a_near, _, _ = local_partial_quotient(lam, n, k1, k2)
            near.append(1.0 / math.sqrt(a_near) if a_near > 0 else 1.0)
            xy = 2.0 * (s_k1 * r) * (s_k2 * q) / det
            if xy > 0:
                tight.append(nh2 / xy)
        print(f"    eff={eff}  spearman(nu_hat, 1/sqrt(a_ACHIEVING)) = "
              f"{spearman(nus, ach):+.3f}   "
              f"[nearest-to-a* was {spearman(nus, near):+.3f}]")
        if tight:
            tight.sort()
            print(f"              AM-GM tightness nu_hat^2/(2xy): "
                  f"median={tight[len(tight)//2]:.3f} "
                  f"q10={tight[len(tight)//10]:.3f} "
                  f"q90={tight[9*len(tight)//10]:.3f}")

    print("\n  --- Lower bound check: nu_hat^2 >= 2*x*y ---")
    viol = 0
    tot = 0
    for c in curves:
        n, lam = c['n'], c['lam']
        for eff in (0.05, 0.25):
            k1, k2 = k1k2(n, eff)
            s_k1, _, s_k2, _ = scales(n, k1, k2)
            det = n * s_k1 * s_k2
            v, a, r, _ = lambda1_from_cf(n, lam, s_k1, s_k2, True)
            nh2 = v / det
            xy = 2.0 * (s_k1 * r) * (s_k2 * a) / det
            tot += 1
            viol += nh2 + 1e-12 < xy
    print(f"    AM-GM bound nu_hat^2 >= 2*x*y violated in {viol}/{tot} cases"
          " (0 expected)")

# ---------------------------------------------------------------------------
# E5 -- secp256k1
# ---------------------------------------------------------------------------

def exp5():
    print("\n" + "=" * 78)
    print("E5  secp256k1 -- exact-integer nu_hat and the achieving convergent")
    print("=" * 78)
    n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    roots = glv_eigenvalues(n)
    lam = roots[0]
    assert (lam * lam + lam + 1) % n == 0, "lambda is not a GLV eigenvalue"
    print(f"  n   = {n}")
    print(f"  lam = {lam}   (lam^2+lam+1 = 0 mod n verified)")
    print(f"  mu  = {mu_of(lam, n):.6f}")
    conv = cf_convergents(lam, n)
    print(f"  CF of lam/n: {len(conv)-1} convergents, "
          f"max partial quotient = {global_max_a(lam, n)}")

    print("\n  eff    a*(scale)      nu_hat(exact)  nu_hat(float LG)   idx  a_local")
    for eff in (0.05, 0.0993, 0.10, 0.25, 0.50):
        k1, k2 = k1k2(n, eff)
        nh, a, r, pq = nu_hat_exact(n, lam, k1, k2)
        _, _, nh_f = rival_sublattice_nu(n, lam, k1, k2)
        a_loc, j, _ = local_partial_quotient(lam, n, k1, k2)
        idx = next((i for i, (q, _, _) in enumerate(conv) if q == a), -1)
        print(f"  {eff:<6} {scale_point(n,k1,k2):.4e}   {nh:.6f}       "
              f"{nh_f:.6f}      {idx:<4} {a_loc}")
    print("\n  (2026-07-29 reported nu_hat = 0.8709 at eff=0.05, 0.6624 at eff=0.10,")
    print("   0.5852 at eff=0.25, 0.6851 at eff=0.50, via float Lagrange-Gauss.)")

# ---------------------------------------------------------------------------

def main():
    t0 = time.time()
    print("=" * 78)
    print("GLV-HNP Phase 2 -- Thread 23: closed form for lambda_1(L2)")
    print("=" * 78)
    rng = random.Random(20260806)

    exp1(rng)

    print("\n  harvesting 20-bit sample for E2-E4 ...")
    curves = harvest_curves(2 ** 19, 2 ** 20, 300)
    print(f"  {len(curves)} curves")

    exp2(rng, curves)
    exp3(curves)
    exp4(curves)
    exp5()

    print(f"\nTotal wall time {time.time()-t0:.1f}s")
    print("Done.")

if __name__ == '__main__':
    sys.exit(main())
