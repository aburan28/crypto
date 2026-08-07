"""
GLV-HNP — Thread 23: closed form for nu_hat via the continued fraction of lam/n.

Background
----------
2026-07-29 (session 2, commit e845207) found nu_hat, the first predictor of GLV-HNP
Phase 2 success to beat the trivial baseline (AUC 0.935 against the June C1/C2
classes) after eight falsified invariants.  It is defined on the NON-PLANTED
2-dimensional sublattice of the (2m+2)-dim Phase 2 lattice,

    L2 = < (n*A, 0), (-lam*A, B) >,   A = S_K1 = n//K1,  B = S_K2 = max(1, n//K2),
    det L2 = n*A*B  (independent of lam),
    nu_hat = lambda_1(L2) / sqrt(det L2).

That entry's next-step proposal (Thread 23) was: nu_hat is *computed*, not
*understood*.  Derive lambda_1(L2) in closed form from the continued fraction of
lam/n and pin the threshold.  Stated falsifier there: "if predicted-from-CF nu_hat
correlates <0.9 with computed nu_hat, the convergent-scale story is wrong."

What this script tests
----------------------
A general vector of L2 is

    q*(-lam*A, B) + b*(n*A, 0) = ( A*(b*n - q*lam),  q*B ),

so for each integer q >= 0 the best b gives first coordinate A*delta(q), where
delta(q) = dist(q*lam, n*Z) is the centered residue.  Hence

    lambda_1(L2)^2 = min_{q >= 0} [ (A*delta(q))^2 + (q*B)^2 ],   (q,b) != (0,0)

with q=0 contributing the vector (n*A, 0).  Minimising a positive-definite
quadratic in delta(q) and q is *best rational approximation to lam/n at a scale*,
so the minimiser should be a denominator of a semiconvergent of lam/n.  Define

    Q = sqrt(n*A/B) = sqrt(n*K2/K1)      (the balance scale)

Normalising by sqrt(det) = sqrt(n*A*B) gives the scale-free form tested here:

    CLAIM 1 (exact).  nu_hat^2 = min_{q in D} [ (Q*delta(q)/n)^2 + (q/Q)^2 ],
      D = {0} union {denominators of semiconvergents of lam/n with q <= 1.08*Q}.

    CLAIM 2 (asymptotic).  With u_j = q_j/Q and delta(q_j) ~ n/q_{j+1},
      nu_hat^2 ~ min_j [ u_j^2 + 1/u_{j+1}^2 ].
      This is small exactly when some j has q_j << Q << q_{j+1}, i.e. when a LARGE
      partial quotient STRADDLES the scale Q.  Predictor:
          a_straddle = a_{j*+1} where q_{j*} <= Q < q_{j*+1}.

    CLAIM 3 (explanatory).  a_straddle is scale-DEPENDENT (it moves with K1/K2),
      whereas q_cf / max_q_cf / max_a -- falsified in June 2026-06-27..06-29 --
      are scale-FREE functions of lam/n alone.  That is why they had to fail.

Experiments
-----------
  A  exactness of CLAIM 1: exact integer Gauss reduction vs CF enumeration,
     convergents-only and semiconvergents, over random and GLV lam.
  B  brute-force control at small n: enumerate every q <= 1.08*Q.
  C  CLAIM 2: accuracy of the u_j^2 + 1/u_{j+1}^2 approximation.
  D  CLAIM 2 predictor: spearman(a_straddle, nu_hat) vs spearman(max_a, nu_hat).
  E  CLAIM 3: fix the curve, sweep eff -- does j* move while max_a stays fixed?
  F  secp256k1: reproduce the four nu_hat values logged 2026-07-29 and name the
     convergent index that attains each.
  G  translate the empirical C2 ceiling nu_hat <= 0.645 into an a_straddle cut.

No fpylll / no PARI needed: everything is exact integer arithmetic plus one
Lagrange-Gauss reduction per instance, O(log n).

Run: python3 glv_hnp_nuhat_cf.py
"""

import math
import random
import sys

import sympy

# ---------------------------------------------------------------------------
# exact 2D lattice machinery
# ---------------------------------------------------------------------------

def gauss_reduce_2d_exact(u, v):
    """Lagrange/Gauss reduction over Z.  Returns the exact shortest nonzero
    vector as a tuple of ints (norms compared as exact integers)."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        # exact round-half-away integer division
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def scales(n, k1_bound, k2_bound):
    """Column-diagonal scaling validated 2026-07-26 (S_K1, S_D, S_K2, S_KANNAN)."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def lambda1_sq_exact(n, lam, A, B):
    """Exact lambda_1(L2)^2 as a Python int."""
    w = gauss_reduce_2d_exact((n * A, 0), (-(lam % n) * A, B))
    return w[0] * w[0] + w[1] * w[1]

def nu_hat_exact(n, lam, k1_bound, k2_bound):
    A, _, B, _ = scales(n, k1_bound, k2_bound)
    l1sq = lambda1_sq_exact(n, lam, A, B)
    det = n * A * B
    return math.sqrt(l1sq / det), l1sq, A, B, det

# ---------------------------------------------------------------------------
# continued fraction of lam/n
# ---------------------------------------------------------------------------

def cf_expansion(num, den):
    """Partial quotients of num/den (num < den), plus convergent denominators.
    Returns (a_list, q_list) where q_list[j] is the denominator of the j-th
    convergent, q_list[0] = 1 corresponds to a_list[0] (= a_1 in [0;a1,a2,...])."""
    a, q = [], []
    q_prev, q_cur = 0, 1     # q_{-1}, q_0
    x, y = num, den
    # lam/n < 1 so a_0 = 0; start the loop at a_1
    while y:
        ai = x // y
        x, y = y, x - ai * y
        if not a and ai == 0:
            # this is a_0 = 0 for lam/n < 1; skip, do not emit a convergent
            continue
        a.append(ai)
        q_prev, q_cur = q_cur, ai * q_cur + q_prev
        q.append(q_cur)
    return a, q

def centered_residue(q, lam, n):
    """delta(q) = dist(q*lam, n*Z), exact."""
    r = (q * lam) % n
    return min(r, n - r)

def semiconvergent_denoms(a, q, cap):
    """All semiconvergent (intermediate-fraction) denominators <= cap.
    q_{j-1}*t + q_{j-2} for t = 1..a_j; t = a_j reproduces the convergent q_j."""
    out = set()
    q_prev2, q_prev1 = 0, 1    # q_{-1}, q_0 ... but q_0 = 1 only after a_0 = 0
    # rebuild recurrence alongside the emitted convergents
    qm2, qm1 = 0, 1
    for j, aj in enumerate(a):
        for t in range(1, aj + 1):
            cand = t * qm1 + qm2
            if cand <= cap:
                out.add(cand)
            else:
                break
        qm2, qm1 = qm1, aj * qm1 + qm2
        if qm1 > cap:
            # later convergents only grow; intermediates already covered above
            break
    out.add(1)
    return sorted(out)

def nu_hat_from_cf(n, lam, A, B, mode='semi'):
    """CF-predicted lambda_1^2 (exact int) and the minimising q.
    mode='conv'  -> convergent denominators only
    mode='semi'  -> all semiconvergent denominators"""
    det = n * A * B
    # lambda_1 <= (4/3)^(1/4) sqrt(det)  =>  q*B <= 1.0746 sqrt(det)
    cap = int(1.08 * math.isqrt(det) // B) + 2
    a, q = cf_expansion(lam % n, n)
    cands = q if mode == 'conv' else semiconvergent_denoms(a, q, cap)
    best, bestq = (n * A) ** 2, 0          # the q = 0 vector (n*A, 0)
    for qq in cands:
        if qq == 0 or qq > cap:
            continue
        d = centered_residue(qq, lam, n)
        val = (A * d) ** 2 + (qq * B) ** 2
        if val < best:
            best, bestq = val, qq
    return best, bestq

def straddle_data(n, lam, A, B):
    """Q, the straddling index j*, a_straddle = a_{j*+1}, max_a, and u_j."""
    Q = math.sqrt(n * A / B)
    a, q = cf_expansion(lam % n, n)
    jstar = -1
    for j, qj in enumerate(q):
        if qj <= Q:
            jstar = j
        else:
            break
    a_str = a[jstar + 1] if 0 <= jstar + 1 < len(a) else (a[0] if a else 1)
    return Q, jstar, a_str, (max(a) if a else 1), a, q

# ---------------------------------------------------------------------------
# statistics helpers
# ---------------------------------------------------------------------------

def rank(xs):
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

def spearman(xs, ys):
    rx, ry = rank(xs), rank(ys)
    nn = len(xs)
    mx, my = sum(rx) / nn, sum(ry) / nn
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(nn))
    dx = math.sqrt(sum((v - mx) ** 2 for v in rx))
    dy = math.sqrt(sum((v - my) ** 2 for v in ry))
    return num / (dx * dy) if dx * dy else 0.0

def pearson(xs, ys):
    nn = len(xs)
    mx, my = sum(xs) / nn, sum(ys) / nn
    num = sum((xs[i] - mx) * (ys[i] - my) for i in range(nn))
    dx = math.sqrt(sum((v - mx) ** 2 for v in xs))
    dy = math.sqrt(sum((v - my) ** 2 for v in ys))
    return num / (dx * dy) if dx * dy else 0.0

# ---------------------------------------------------------------------------
# instance generation
# ---------------------------------------------------------------------------

def glv_eigenvalue(n):
    """Smaller root of x^2 + x + 1 = 0 mod n, or None."""
    if n % 3 != 1:
        return None
    # need a cube root of unity: g^((n-1)/3)
    for g in range(2, 200):
        w = pow(g, (n - 1) // 3, n)
        if w != 1:
            return min(w, n - 1 - w) if (w * w + w + 1) % n == 0 else None
    return None

def k1k2(n, eff):
    """K1 = eff*sqrt(n) (rounded), K2 = sqrt(n); K1*K2/n = eff."""
    return max(2, int(eff * math.isqrt(n))), math.isqrt(n) + 1

def random_instances(bits, count, rng, glv=False):
    """(n, lam) pairs.  glv=True restricts to n = 1 mod 3 with the GLV eigenvalue."""
    out = []
    while len(out) < count:
        n = int(sympy.nextprime(rng.randrange(1 << (bits - 1), 1 << bits)))
        if glv:
            if n % 3 != 1:
                continue
            lam = glv_eigenvalue(n)
            if lam is None:
                continue
        else:
            lam = rng.randrange(1, n)
        out.append((n, lam))
    return out

# ---------------------------------------------------------------------------
# EXPERIMENTS
# ---------------------------------------------------------------------------

def exp_A(rng):
    print("=" * 78)
    print("EXP A -- CLAIM 1: exact CF characterisation of lambda_1(L2)")
    print("=" * 78)
    print("Exact integer Gauss reduction vs CF enumeration.  'match' = exact")
    print("equality of the integer lambda_1^2, not a correlation.\n")
    print(f"{'bits':>5} {'lam':>6} {'eff':>6} {'N':>5} {'conv match':>11} {'semi match':>11}"
          f" {'conv relerr(max)':>17}")
    rows = []
    for bits in (20, 24, 32, 64, 128, 256):
        for glv in (False, True):
            insts = random_instances(bits, 60, rng, glv=glv)
            for eff in (0.05, 0.10, 0.25):
                nconv = nsemi = 0
                worst = 0.0
                for (n, lam) in insts:
                    k1, k2 = k1k2(n, eff)
                    A, _, B, _ = scales(n, k1, k2)
                    if A == 0 or B == 0:
                        continue
                    exact = lambda1_sq_exact(n, lam, A, B)
                    c, _ = nu_hat_from_cf(n, lam, A, B, 'conv')
                    s, _ = nu_hat_from_cf(n, lam, A, B, 'semi')
                    nconv += (c == exact)
                    nsemi += (s == exact)
                    worst = max(worst, math.sqrt(c / exact) - 1.0)
                    rows.append((n, lam, eff, exact, s))
                N = len(insts)
                print(f"{bits:>5} {'GLV' if glv else 'rand':>6} {eff:>6.2f} {N:>5}"
                      f" {nconv:>6}/{N:<4} {nsemi:>6}/{N:<4} {worst:>17.3e}")
    # global tally
    tot = len(rows)
    ok = sum(1 for (_, _, _, e, s) in rows if e == s)
    print(f"\nCLAIM 1 (semiconvergents): {ok}/{tot} exact "
          f"({100.0*ok/tot:.2f}%)")
    return rows

def exp_B(rng):
    print()
    print("=" * 78)
    print("EXP B -- brute-force control (small n): is the minimiser always a")
    print("         semiconvergent denominator?")
    print("=" * 78)
    bad = 0
    tested = 0
    for bits in (16, 18, 20):
        for (n, lam) in random_instances(bits, 40, rng):
            for eff in (0.02, 0.05, 0.10, 0.25, 0.5):
                k1, k2 = k1k2(n, eff)
                A, _, B, _ = scales(n, k1, k2)
                if A == 0 or B == 0:
                    continue
                det = n * A * B
                cap = int(1.08 * math.isqrt(det) // B) + 2
                if cap > 400000:
                    continue
                best = (n * A) ** 2
                bq = 0
                for q in range(1, cap + 1):
                    d = centered_residue(q, lam, n)
                    v = (A * d) ** 2 + (q * B) ** 2
                    if v < best:
                        best, bq = v, q
                s, sq = nu_hat_from_cf(n, lam, A, B, 'semi')
                tested += 1
                if s != best:
                    bad += 1
                    if bad <= 5:
                        print(f"  MISMATCH n={n} lam={lam} eff={eff} "
                              f"brute q={bq} cf q={sq}")
    print(f"  brute-force cases: {tested}, mismatches: {bad}")
    print(f"  => the minimiser is a semiconvergent denominator in "
          f"{tested-bad}/{tested} cases")

def exp_C(rng):
    print()
    print("=" * 78)
    print("EXP C -- CLAIM 2: accuracy of  nu_hat^2 ~ min_j [ u_j^2 + 1/u_{j+1}^2 ]")
    print("=" * 78)
    print(f"{'bits':>5} {'eff':>6} {'N':>4} {'median relerr':>14} {'p90 relerr':>11}"
          f" {'pearson':>9} {'spearman':>9}")
    for bits in (24, 64, 256):
        for eff in (0.05, 0.10, 0.25):
            ex, ap = [], []
            for (n, lam) in random_instances(bits, 80, rng):
                k1, k2 = k1k2(n, eff)
                A, _, B, _ = scales(n, k1, k2)
                if A == 0 or B == 0:
                    continue
                l1sq = lambda1_sq_exact(n, lam, A, B)
                det = n * A * B
                nu = math.sqrt(l1sq / det)
                Q, jstar, a_str, max_a, a, q = straddle_data(n, lam, A, B)
                # min over j of u_j^2 + 1/u_{j+1}^2  with u_j = q_j/Q
                best = None
                us = [qq / Q for qq in q]
                for j in range(len(us) - 1):
                    v = us[j] ** 2 + 1.0 / (us[j + 1] ** 2)
                    best = v if best is None else min(best, v)
                if best is None:
                    continue
                ex.append(nu)
                ap.append(math.sqrt(best))
            rel = sorted(abs(ap[i] - ex[i]) / ex[i] for i in range(len(ex)))
            print(f"{bits:>5} {eff:>6.2f} {len(ex):>4} {rel[len(rel)//2]:>14.4f}"
                  f" {rel[int(0.9*len(rel))]:>11.4f}"
                  f" {pearson(ex, ap):>9.4f} {spearman(ex, ap):>9.4f}")

def exp_D(rng):
    print()
    print("=" * 78)
    print("EXP D -- CLAIM 2 predictor: a_straddle vs the falsified scale-free max_a")
    print("=" * 78)
    print("a_straddle = a_{j*+1} where q_{j*} <= Q < q_{j*+1}.  Small nu_hat should")
    print("mean a LARGE straddling partial quotient, so the spearman is negative.\n")
    print(f"{'bits':>5} {'eff':>6} {'N':>4} {'sp(a_str,nu)':>13} {'sp(max_a,nu)':>13}"
          f" {'sp(q_cf,nu)':>12}")
    for bits in (24, 64, 256):
        for eff in (0.05, 0.10, 0.25):
            nus, astr, maxa, qcf = [], [], [], []
            for (n, lam) in random_instances(bits, 200, rng):
                k1, k2 = k1k2(n, eff)
                A, _, B, _ = scales(n, k1, k2)
                if A == 0 or B == 0:
                    continue
                l1sq = lambda1_sq_exact(n, lam, A, B)
                nus.append(math.sqrt(l1sq / (n * A * B)))
                Q, jstar, a_s, m_a, a, q = straddle_data(n, lam, A, B)
                astr.append(float(a_s))
                maxa.append(float(m_a))
                qcf.append(float(len(a)))       # CF length, the June 'q_cf'
            print(f"{bits:>5} {eff:>6.2f} {len(nus):>4} {spearman(astr, nus):>13.4f}"
                  f" {spearman(maxa, nus):>13.4f} {spearman(qcf, nus):>12.4f}")

def exp_E(rng):
    print()
    print("=" * 78)
    print("EXP E -- CLAIM 3: nu_hat is scale-DEPENDENT, so no scale-free CF")
    print("         invariant of lam/n can predict it")
    print("=" * 78)
    print("One FIXED (n, lam); only eff = K1*K2/n varies.  max_a and the CF itself")
    print("are constant by construction; j*, a_straddle and nu_hat are not.\n")
    for (n, lam) in random_instances(64, 3, rng):
        _, _, _, max_a, a, q = straddle_data(n, lam, n // 2, n // 2)
        print(f"n = {n}  lam = {lam}")
        print(f"  CF length = {len(a)}, max_a = {max_a}  (both eff-independent)")
        print(f"  {'eff':>7} {'Q':>16} {'j*':>4} {'a_straddle':>11} {'nu_hat':>9}"
              f" {'argmin q':>16}")
        for eff in (0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 0.90):
            k1, k2 = k1k2(n, eff)
            A, _, B, _ = scales(n, k1, k2)
            if A == 0 or B == 0:
                continue
            l1sq = lambda1_sq_exact(n, lam, A, B)
            nu = math.sqrt(l1sq / (n * A * B))
            Q, jstar, a_s, _, _, _ = straddle_data(n, lam, A, B)
            _, argq = nu_hat_from_cf(n, lam, A, B, 'semi')
            print(f"  {eff:>7.2f} {Q:>16.4e} {jstar:>4} {a_s:>11} {nu:>9.4f}"
                  f" {argq:>16}")
        print()

SECP256K1_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP256K1_LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72

def exp_F():
    print()
    print("=" * 78)
    print("EXP F -- secp256k1: reproduce the four nu_hat values logged 2026-07-29")
    print("=" * 78)
    n, lam = SECP256K1_N, SECP256K1_LAM
    assert (lam * lam + lam + 1) % n == 0, "lam is not a GLV eigenvalue"
    print(f"  lam^2 + lam + 1 = 0 mod n  VERIFIED")
    logged = {0.05: 0.8709, 0.10: 0.6624, 0.25: 0.5852, 0.50: 0.6851}
    _, _, _, max_a, a_all, _ = straddle_data(n, lam, n // 2, n // 2)
    print(f"  CF of lam/n: length {len(a_all)}, max partial quotient {max_a}\n")
    print(f"  {'eff':>6} {'nu_hat (this run)':>18} {'logged 2026-07-29':>18}"
          f" {'j*':>4} {'a_straddle':>11} {'argmin q bits':>14}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1, k2 = k1k2(n, eff)
        A, _, B, _ = scales(n, k1, k2)
        l1sq = lambda1_sq_exact(n, lam, A, B)
        nu = math.sqrt(l1sq / (n * A * B))
        Q, jstar, a_s, _, _, _ = straddle_data(n, lam, A, B)
        _, argq = nu_hat_from_cf(n, lam, A, B, 'semi')
        cf_sq, _ = nu_hat_from_cf(n, lam, A, B, 'semi')
        assert cf_sq == l1sq, "CLAIM 1 fails on secp256k1"
        print(f"  {eff:>6.2f} {nu:>18.4f} {logged[eff]:>18.4f} {jstar:>4}"
              f" {a_s:>11} {argq.bit_length():>14}")
    print("\n  (CLAIM 1 verified exactly at all four eff on the real secp256k1 lattice)")
    print("  Reminder: conditional on a NON-STANDARD nonce k = k1 + lam*k2 with k1")
    print("  bounded.  This is not an attack on secp256k1.")

def exp_G(rng):
    print()
    print("=" * 78)
    print("EXP G -- translate the empirical C2 ceiling nu_hat <= 0.645 into an")
    print("         a_straddle cut")
    print("=" * 78)
    eff = 0.05
    for bits in (20, 64, 256):
        rows = []
        for (n, lam) in random_instances(bits, 400, rng):
            k1, k2 = k1k2(n, eff)
            A, _, B, _ = scales(n, k1, k2)
            if A == 0 or B == 0:
                continue
            nu = math.sqrt(lambda1_sq_exact(n, lam, A, B) / (n * A * B))
            _, _, a_s, _, _, _ = straddle_data(n, lam, A, B)
            rows.append((nu, a_s))
        lo = [r for r in rows if r[0] <= 0.645]
        hi = [r for r in rows if r[0] > 0.645]
        def med(v):
            v = sorted(v)
            return v[len(v) // 2] if v else float('nan')
        # smallest a_straddle cut that captures most of the nu<=0.645 set
        best_t, best_acc = None, 0.0
        for t in range(2, 60):
            acc = sum(1 for (nu, a) in rows if (a >= t) == (nu <= 0.645)) / len(rows)
            if acc > best_acc:
                best_acc, best_t = acc, t
        print(f"  {bits:>3} bits, eff={eff}:  n={len(rows)}  "
              f"P(nu<=0.645)={len(lo)/len(rows):.3f}")
        print(f"      median a_straddle | nu<=0.645 : {med([a for _,a in lo])}")
        print(f"      median a_straddle | nu >0.645 : {med([a for _,a in hi])}")
        print(f"      best cut 'a_straddle >= {best_t}' reproduces the nu cut with "
              f"accuracy {best_acc:.3f}")

def exp_H(rng):
    print()
    print("=" * 78)
    print("EXP H -- consequence of CLAIM 1: the eff-optimised nu_hat, and why max_a")
    print("         *nearly* worked in June")
    print("=" * 78)
    print("For a FIXED winning convergent (q_j, delta_j), CLAIM 1 gives")
    print("    nu_hat^2 = (Q*delta_j/n)^2 + (q_j/Q)^2,")
    print("a function of Q alone.  Minimising over Q:")
    print("    Q_opt^2 = n*q_j/delta_j,   nu_hat_min^2 = 2*q_j*delta_j/n.")
    print("With K2 = sqrt(n) and K1 = eff*sqrt(n) we have Q^2 = n/eff, so")
    print("    eff_opt = delta_j / q_j     (feasible iff delta_j <= q_j),")
    print("    nu_hat_min = sqrt(2 * min_j q_j*delta_j / n).")
    print("Since q_j*delta_j/n ~ q_j/q_{j+1} ~ 1/a_{j+1},")
    print("    nu_hat_min ~ sqrt(2 / max_a).")
    print("So max_a DOES govern nu_hat -- but only after optimising eff.  At the")
    print("FIXED eff of every June experiment, the relevant quotient is a_straddle,")
    print("not max_a.  That is the precise sense in which the June invariant was")
    print("the right quantity measured at the wrong scale.\n")

    def analyse(n, lam, label):
        a, q = cf_expansion(lam % n, n)
        best = None
        for j, qq in enumerate(q):
            d = centered_residue(qq, lam, n)
            if d == 0 or d > qq:
                continue                      # eff_opt = d/q must be <= 1
            val = 2.0 * qq * d / n
            if best is None or val < best[0]:
                best = (val, j, qq, d)
        if best is None:
            print(f"  {label}: no feasible convergent")
            return
        val, j, qq, d = best
        nu_min = math.sqrt(val)
        eff_opt = d / qq
        max_a = max(a)
        k1_implied = eff_opt * math.isqrt(n)
        print(f"  {label}")
        print(f"    max_a = {max_a}   sqrt(2/max_a) = {math.sqrt(2.0/max_a):.4f}")
        print(f"    argmin j = {j} (q has {qq.bit_length()} bits), a_{{j+1}} = "
              f"{a[j+1] if j+1 < len(a) else '-'}")
        print(f"    UNCONSTRAINED nu_hat_min = {nu_min:.4f}  at eff_opt = "
              f"{eff_opt:.4e}  (implied K1 = {k1_implied:.3e})")
        if k1_implied < 2:
            print(f"    -> INFEASIBLE: K1 < 2, i.e. the k1 nonce part carries <1 bit.")
        # feasible regime: eff in [EFF_MIN, 1]  <=>  Q in [sqrt(n), sqrt(n/EFF_MIN)]
        EFF_MIN = 0.01
        bestf = None
        for jj, q2 in enumerate(q):
            d2 = centered_residue(q2, lam, n)
            if d2 == 0:
                continue
            for eff in (EFF_MIN, min(1.0, max(EFF_MIN, d2 / q2)), 1.0):
                Q2 = math.sqrt(n / eff)
                v = (Q2 * d2 / n) ** 2 + (q2 / Q2) ** 2
                if bestf is None or v < bestf[0]:
                    bestf = (v, jj, eff)
        vf, jf, efff = bestf
        print(f"    FEASIBLE (eff in [{EFF_MIN}, 1]) nu_hat_min = "
              f"{math.sqrt(vf):.4f}  at eff = {efff:.4f}, j = {jf}")
        k1 = max(2, int(efff * math.isqrt(n)))
        k2 = math.isqrt(n) + 1
        A, _, B, _ = scales(n, k1, k2)
        if A > 0 and B > 0:
            nu = math.sqrt(lambda1_sq_exact(n, lam, A, B) / (n * A * B))
            print(f"    exact nu_hat at that eff = {nu:.4f}  (check, "
                  f"integer-K1 rounding aside)")

    analyse(SECP256K1_N, SECP256K1_LAM, "secp256k1 (real n, lam)")
    print()
    for i, (n, lam) in enumerate(random_instances(64, 3, rng)):
        analyse(n, lam, f"random 64-bit #{i+1}  n={n}")
        print()

    # how well does sqrt(2/max_a) predict the eff-optimised nu_hat, unconstrained
    # vs restricted to the operationally relevant eff range?
    xs, yu, yf = [], [], []
    EFF_MIN = 0.01
    for (n, lam) in random_instances(64, 200, rng):
        a, q = cf_expansion(lam % n, n)
        bu = bf = None
        for qq in q:
            d = centered_residue(qq, lam, n)
            if d == 0:
                continue
            if d <= qq:
                v = 2.0 * qq * d / n
                bu = v if bu is None else min(bu, v)
            for eff in (EFF_MIN, min(1.0, max(EFF_MIN, d / qq)), 1.0):
                Q2 = math.sqrt(n / eff)
                v2 = (Q2 * d / n) ** 2 + (qq / Q2) ** 2
                bf = v2 if bf is None else min(bf, v2)
        if bu is None or bf is None:
            continue
        xs.append(math.sqrt(2.0 / max(a)))
        yu.append(math.sqrt(bu))
        yf.append(math.sqrt(bf))
    print(f"  64-bit, N={len(xs)}:")
    print(f"    spearman(sqrt(2/max_a), nu_hat_min UNCONSTRAINED) = "
          f"{spearman(xs, yu):>7.4f}   pearson {pearson(xs, yu):>7.4f}")
    print(f"    spearman(sqrt(2/max_a), nu_hat_min FEASIBLE eff)  = "
          f"{spearman(xs, yf):>7.4f}   pearson {pearson(xs, yf):>7.4f}")
    print("  Contrast with EXP D, where sp(max_a, nu_hat) at fixed eff was "
          "-0.01 .. -0.32.")

def main():
    seed = 20260807
    print("Thread 23 -- closed form for nu_hat from the continued fraction of lam/n")
    print(f"seed = {seed}\n")
    exp_A(random.Random(seed))
    exp_B(random.Random(seed + 1))
    exp_C(random.Random(seed + 2))
    exp_D(random.Random(seed + 3))
    exp_E(random.Random(seed + 4))
    exp_F()
    exp_G(random.Random(seed + 5))
    exp_H(random.Random(seed + 6))

if __name__ == '__main__':
    main()
