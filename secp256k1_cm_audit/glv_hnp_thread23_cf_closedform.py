"""
GLV-HNP -- Thread 23: closed form for lambda_1(L2), and why six June
invariants failed.

Background
----------
Thread 20c/20d (2026-07-29, commit e845207) found the first working separator
for the Phase 2 / June-Exp-S lattice after eight falsified invariants:

    L2   = < (n*S_K1, 0), (-lam*S_K1, S_K2) >        (non-planted 2D sublattice)
    nu   = lambda_1(L2)                              (Lagrange-Gauss, O(log n))
    nu_hat = nu / sqrt(det L2),   det L2 = n*S_K1*S_K2

AUC 0.935 against the June C1/C2 classes vs 0.523 (mu) and 0.748 (max_a).
But nu_hat was *computed*, not *understood*.  The 2026-07-29 next-step proposal
asked for a closed form and set the falsifier at "predicted-from-CF nu_hat
correlates < 0.9 with computed nu_hat".

The claim tested here is much stronger than that falsifier: not a correlation,
an identity.

Theory
------
A general element of L2 is

    a*(n*S1, 0) + b*(-lam*S1, S2) = ( (a*n - b*lam)*S1,  b*S2 ),   a,b in Z.

For fixed b the inner coordinate is minimised by the centred residue

    r(b) = |b*lam mod^{+-} n|  in [0, n/2],

so, writing S1 = S_K1 and S2 = S_K2,

    lambda_1(L2)^2 = min{ (n*S1)^2 ,  min_{b>=1} [ (r(b)*S1)^2 + (b*S2)^2 ] }.

f(b) = (r(b)*S1)^2 + (b*S2)^2 is strictly increasing in each of b and r(b), so
its minimiser must be Pareto-optimal for (b, r(b)): no smaller b may have a
residue that is no larger.  Those Pareto points are exactly the *best
approximations of the second kind* to lam/n, and by the classical theorem those
are exactly the continued-fraction convergent denominators q_j of lam/n.  Hence

    lambda_1(L2)^2 = min{ (n*S1)^2, min_j [ (r_j*S1)^2 + (q_j*S2)^2 ] },
    r_j = |q_j*lam - p_j*n|,                                            (*)

an O(log n) formula over the ~log n convergents -- no lattice reduction.

Localisation.  Since r_j * q_{j+1} ~ n, the two terms of (*) balance near

    q_j ~ sqrt(n * S1 / S2) = sqrt(n * K2 / K1)   =:  q_crit,

so a single convergent -- the one at scale q_crit -- decides nu_hat.

*This is the retroactive explanation of the June falsifications.*  q_cf, max_q_cf
and max_a are scale-FREE: max_a is the largest partial quotient anywhere in the
expansion.  A large partial quotient only makes the attack easier if it sits at
the critical index.  The scale-LOCAL partial quotient a_{j*+1} should therefore
predict what max_a could not.

Experiments
-----------
A  exactness      -- (*) vs Lagrange-Gauss over many (n, lam, K1, K2), to 256 bits
A2 Pareto/brute   -- convergents really are the only candidates (small-n brute force)
B  localisation   -- how close is argmin j* to the q_crit index?
C  scale-local    -- a_local = a_{j*+1} vs max_a as a predictor of nu_hat
D  secp256k1      -- critical index and a_local at the eff values of 20d

Run: python3 glv_hnp_thread23_cf_closedform.py
"""

import importlib.util
import math
import random
import sys
import time

_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_lib.py")
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

gauss_reduce_2d = _t20a.gauss_reduce_2d
scales = _t20a.scales
rival_sublattice_nu = _t20a.rival_sublattice_nu

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_LAM = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72

# ---------------------------------------------------------------------------
# Continued fractions
# ---------------------------------------------------------------------------

def cf_convergents(lam, n):
    """
    Convergents p_j/q_j of lam/n.

    Returns a list of records (j, a_j, q_j, r_j) where a_j is the j-th partial
    quotient, q_j the convergent denominator and r_j = |q_j*lam - p_j*n| the
    (uncentred) approximation defect.  j starts at 0 with a_0 = lam // n.
    """
    out = []
    x, y = lam, n
    p_prev, p_cur = 0, 1      # p_{-2}, p_{-1}
    q_prev, q_cur = 1, 0      # q_{-2}, q_{-1}
    j = 0
    while y:
        a = x // y
        x, y = y, x - a * y
        p_prev, p_cur = p_cur, a * p_cur + p_prev
        q_prev, q_cur = q_cur, a * q_cur + q_prev
        out.append((j, a, q_cur, abs(q_cur * lam - p_cur * n)))
        j += 1
    return out

def lambda1_cf(n, lam, s1, s2):
    """
    lambda_1(L2) by the closed form (*), plus the winning convergent.

    Returns (lambda_1, j_star, q_star, a_next) where j_star is the index of the
    minimising convergent (-1 for the b=0 vector (n*S1, 0)) and a_next is the
    partial quotient immediately after it -- the 'scale-local' quotient.
    """
    best = (n * s1) ** 2
    j_star, q_star, a_next = -1, 0, 0
    cvs = cf_convergents(lam, n)
    for idx, (j, _a, q, r) in enumerate(cvs):
        cand = (r * s1) ** 2 + (q * s2) ** 2
        if cand < best:
            best = cand
            j_star, q_star = j, q
            a_next = cvs[idx + 1][1] if idx + 1 < len(cvs) else 0
    return math.sqrt(best), j_star, q_star, a_next

def lambda1_brute(n, lam, s1, s2):
    """Reference implementation: scan every b in [0, n)."""
    best = (n * s1) ** 2
    arg = 0
    for b in range(1, n):
        r = (b * lam) % n
        r = min(r, n - r)
        cand = (r * s1) ** 2 + (b * s2) ** 2
        if cand < best:
            best, arg = cand, b
    return math.sqrt(best), arg

def max_partial_quotient(a, b):
    """The scale-free 'max_a' invariant falsified by Exp S (2026-06-29)."""
    best = 0
    while b:
        q, r = divmod(a, b)
        best = max(best, q)
        a, b = b, r
    return best

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def rankdata(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks

def pearson(xs, ys):
    nn = len(xs)
    mx, my = sum(xs) / nn, sum(ys) / nn
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    if sxx <= 0 or syy <= 0:
        return float('nan')
    return sxy / math.sqrt(sxx * syy)

def spearman(xs, ys):
    return pearson(rankdata(xs), rankdata(ys))

def perm_p(xs, ys, stat=spearman, trials=2000, seed=99):
    obs = stat(xs, ys)
    rng = random.Random(seed)
    ys2 = list(ys)
    hits = 0
    for _ in range(trials):
        rng.shuffle(ys2)
        if abs(stat(xs, ys2)) >= abs(obs) - 1e-15:
            hits += 1
    return obs, (hits + 1) / (trials + 1)

# ---------------------------------------------------------------------------
# Instance generation
# ---------------------------------------------------------------------------

def random_instances(count, bits, seed, eff=0.10):
    """
    Random (n prime, lam, K1, K2) with K2 = isqrt(n)+1 and K1 set from eff.

    lam is drawn uniformly -- Thread 20c Arm 5 established that the lattice does
    not require lam to be a GLV eigenvalue for the HNP instance to be well posed,
    so uniform lam is the right null model for a geometry question.
    """
    import sympy
    rng = random.Random(seed)
    out = []
    while len(out) < count:
        n = int(sympy.nextprime(rng.randrange(2 ** (bits - 1), 2 ** bits)))
        lam = rng.randrange(1, n)
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff * n / k2))
        out.append((n, lam, k1, k2))
    return out

# ---------------------------------------------------------------------------
# Experiments
# ---------------------------------------------------------------------------

def exp_a_exactness():
    print("=" * 78)
    print("EXP A -- is the CF closed form (*) EXACT against Lagrange-Gauss?")
    print("=" * 78)
    print("falsifier from 2026-07-29 was corr < 0.9; the claim tested is equality\n")
    print(f"{'bits':>6} {'eff':>6} {'cases':>6} {'exact':>7} {'max rel err':>13} {'ms/CF':>8} {'ms/LG':>8}")
    allx, ally = [], []
    total, exact = 0, 0
    for bits in (20, 24, 32, 64, 128, 256):
        for eff in (0.05, 0.10, 0.25, 0.50):
            inst = random_instances(40, bits, seed=hash((bits, eff)) & 0xffff, eff=eff)
            worst = 0.0
            t_cf = t_lg = 0.0
            ok = 0
            for (n, lam, k1, k2) in inst:
                s1, _, s2, _ = scales(n, k1, k2)
                t0 = time.time()
                l_cf, _, _, _ = lambda1_cf(n, lam, s1, s2)
                t_cf += time.time() - t0
                t0 = time.time()
                l_lg = gauss_reduce_2d((n * s1, 0), (-lam * s1, s2))
                t_lg += time.time() - t0
                rel = abs(l_cf - l_lg) / l_lg if l_lg else 0.0
                worst = max(worst, rel)
                if rel < 1e-9:
                    ok += 1
                det = n * s1 * s2
                allx.append(l_cf / math.sqrt(det))
                ally.append(l_lg / math.sqrt(det))
            total += len(inst)
            exact += ok
            print(f"{bits:>6} {eff:>6.2f} {len(inst):>6} {ok:>3}/{len(inst):<3} "
                  f"{worst:>13.2e} {1000*t_cf/len(inst):>8.3f} {1000*t_lg/len(inst):>8.3f}")
    print(f"\n  TOTAL exact: {exact}/{total}")
    print(f"  pearson(nu_hat_cf, nu_hat_gauss) = {pearson(allx, ally):.6f}")
    return exact == total

def exp_a2_pareto():
    print("\n" + "=" * 78)
    print("EXP A2 -- brute force: is the minimiser always a CF convergent?")
    print("=" * 78)
    import sympy
    rng = random.Random(4242)
    bad = 0
    tested = 0
    off_conv = 0
    for _ in range(300):
        n = int(sympy.nextprime(rng.randrange(500, 4000)))
        lam = rng.randrange(1, n)
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(rng.choice([0.05, 0.1, 0.25, 0.5]) * n / k2))
        s1, _, s2, _ = scales(n, k1, k2)
        l_cf, _, q_star, _ = lambda1_cf(n, lam, s1, s2)
        l_bf, b_bf = lambda1_brute(n, lam, s1, s2)
        tested += 1
        if abs(l_cf - l_bf) > 1e-9 * max(1.0, l_bf):
            bad += 1
        qs = {q for (_j, _a, q, _r) in cf_convergents(lam, n)}
        if b_bf != 0 and b_bf not in qs:
            off_conv += 1
    print(f"  instances tested                : {tested}")
    print(f"  CF value != brute-force value   : {bad}")
    print(f"  brute-force argmin NOT a q_j    : {off_conv}")
    print("  (both zero => the Pareto/best-approximation argument is sound)")
    return bad == 0 and off_conv == 0

def exp_b_localisation():
    print("\n" + "=" * 78)
    print("EXP B -- localisation: is the winner the convergent at q_crit?")
    print("=" * 78)
    print("  q_crit = sqrt(n*S1/S2) = sqrt(n*K2/K1)\n")
    print(f"{'bits':>6} {'eff':>6} {'N':>4} {'dj=0':>6} {'|dj|<=1':>8} "
          f"{'median q*/q_crit':>18} {'prod dj=0':>10}")
    for bits in (24, 64, 256):
        for eff in (0.05, 0.10, 0.25, 0.50):
            inst = random_instances(60, bits, seed=hash((bits, eff, 'B')) & 0xffff, eff=eff)
            d0 = d1 = d0p = 0
            ratios = []
            for (n, lam, k1, k2) in inst:
                s1, _, s2, _ = scales(n, k1, k2)
                _, j_star, q_star, _ = lambda1_cf(n, lam, s1, s2)
                cvs = cf_convergents(lam, n)
                q_crit = math.sqrt(n * s1 / s2)
                # naive predictor: convergent whose q is nearest q_crit (log scale)
                j_pred = min(cvs, key=lambda t: abs(math.log(max(t[2], 1)) - math.log(q_crit)))[0]
                # refined predictor: the two terms of (*) balance when
                #   r_j*S1 = q_j*S2  and  r_j ~ n/q_{j+1},  i.e.
                #   q_j * q_{j+1} ~ n*S1/S2 = q_crit^2  -- a condition on the
                # PRODUCT of consecutive denominators, not on q_j alone.  That
                # is why q*/q_crit sits systematically below 1 (by ~sqrt(a)).
                tgt = math.log(n * s1 / s2)
                j_pred2 = min(
                    range(len(cvs) - 1),
                    key=lambda i: abs(math.log(max(cvs[i][2], 1))
                                      + math.log(max(cvs[i + 1][2], 1)) - tgt),
                    default=0)
                if j_star >= 0:
                    if j_star == j_pred:
                        d0 += 1
                    if abs(j_star - j_pred) <= 1:
                        d1 += 1
                    if j_star == cvs[j_pred2][0]:
                        d0p += 1
                    ratios.append(q_star / q_crit)
            ratios.sort()
            med = ratios[len(ratios) // 2] if ratios else float('nan')
            print(f"{bits:>6} {eff:>6.2f} {len(inst):>4} {d0:>6} {d1:>8} "
                  f"{med:>18.3f} {d0p:>10}")

def exp_c_scale_local():
    print("\n" + "=" * 78)
    print("EXP C -- scale-LOCAL vs scale-FREE partial quotient as nu_hat predictor")
    print("=" * 78)
    print("  a_local = a_{j*+1}, the partial quotient just past the winning")
    print("  convergent;  max_a = largest partial quotient anywhere (June, falsified)\n")
    print(f"{'bits':>6} {'eff':>6} {'N':>4} {'sp(a_local,nu_hat)':>20} {'p':>8} "
          f"{'sp(max_a,nu_hat)':>18} {'p':>8}")
    for bits in (24, 64, 256):
        for eff in (0.05, 0.10, 0.25):
            inst = random_instances(200, bits, seed=hash((bits, eff, 'C')) & 0xffff, eff=eff)
            nh, al, ma = [], [], []
            for (n, lam, k1, k2) in inst:
                s1, _, s2, _ = scales(n, k1, k2)
                l1, _, _, a_next = lambda1_cf(n, lam, s1, s2)
                nh.append(l1 / math.sqrt(n * s1 * s2))
                al.append(float(a_next))
                ma.append(float(max_partial_quotient(lam, n)))
            r1, p1 = perm_p(al, nh, trials=1000, seed=7)
            r2, p2 = perm_p(ma, nh, trials=1000, seed=7)
            print(f"{bits:>6} {eff:>6.2f} {len(inst):>4} {r1:>20.3f} {p1:>8.4f} "
                  f"{r2:>18.3f} {p2:>8.4f}")
    print("\n  Interpretation: a large partial quotient helps the attacker only")
    print("  when it sits AT the critical scale.  max_a cannot see the scale,")
    print("  which is why it (and q_cf, max_q_cf) failed Exp S in June.")

def exp_d_secp():
    print("\n" + "=" * 78)
    print("EXP D -- secp256k1: where is the critical index?")
    print("=" * 78)
    assert (SECP_LAM * SECP_LAM + SECP_LAM + 1) % SECP_N == 0
    cvs = cf_convergents(SECP_LAM, SECP_N)
    print(f"  CF of lam/n has {len(cvs)} convergents;  "
          f"max_a = {max_partial_quotient(SECP_LAM, SECP_N)}")
    k2 = math.isqrt(SECP_N) + 1
    print(f"\n{'eff':>6} {'nu_hat':>9} {'j*':>4} {'a_local':>9} {'q*/q_crit':>11} {'agrees LG':>10}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1 = max(2, int(eff * SECP_N / k2))
        s1, _, s2, _ = scales(SECP_N, k1, k2)
        l1, j_star, q_star, a_next = lambda1_cf(SECP_N, SECP_LAM, s1, s2)
        det = SECP_N * s1 * s2
        _, _, nh_ref = rival_sublattice_nu(SECP_N, SECP_LAM, k1, k2)
        q_crit = math.sqrt(SECP_N * s1 / s2)
        agree = abs(l1 / math.sqrt(det) - nh_ref) < 1e-9
        print(f"{eff:>6.2f} {l1/math.sqrt(det):>9.4f} {j_star:>4} {a_next:>9} "
              f"{(q_star/q_crit if q_star else float('nan')):>11.3f} {str(agree):>10}")
    print("\n  Caveats unchanged from 2026-07-29: the whole thread is conditional")
    print("  on a non-standard nonce generator k = k1 + lam*k2 with k1 bounded.")
    print("  It is not an attack on secp256k1.")

def main():
    print("GLV-HNP Thread 23 -- closed form for lambda_1(L2)")
    print("continuation of Thread 20d (e845207, 2026-07-29)\n")
    ok_a = exp_a_exactness()
    ok_a2 = exp_a2_pareto()
    exp_b_localisation()
    exp_c_scale_local()
    exp_d_secp()
    print("\n" + "=" * 78)
    print(f"VERDICT: closed form exact on all random instances : {ok_a}")
    print(f"         convergent-only argument verified by brute: {ok_a2}")
    print("=" * 78)
    return 0 if (ok_a and ok_a2) else 1

if __name__ == '__main__':
    sys.exit(main())
