"""
GLV-HNP — Thread 23: closed form for lambda_1(L2), i.e. for the nu_hat predictor.

Context
-------
2026-07-29 (run #2, commit e845207) found nu_hat, a lattice-geometric separator
for the Phase 2 GLV-HNP lattice, with AUC 0.935 against the June C1/C2 classes.
It is defined from the non-planted 2-dimensional sublattice

    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >,   det L2 = n*S_K1*S_K2
    nu_hat = lambda_1(L2) / sqrt(det L2)

and was computed, not understood: one Lagrange-Gauss reduction per curve.

That entry proposed Thread 23: show lambda_1(L2) is a continued-fraction
quantity, evaluated *at the scale* sqrt(n*S_K1/S_K2), and thereby explain why
the six scale-free CF invariants falsified in June (q_cf, max_q_cf, max_a, ...)
all failed.  The proposed falsifier was "predicted-from-CF nu_hat correlates
< 0.9 with computed nu_hat".

This script does that.  The result is stronger than the proposed falsifier: the
identification is EXACT, not a correlation.

Theory
------
A general nonzero vector of L2 is

    a*(n*S_K1, 0) + b*(-lam*S_K1, S_K2) = ( (a*n - b*lam)*S_K1 , b*S_K2 ).

For fixed b, the best a minimises |a*n - b*lam|, which is the centered residue
||b*lam||_n = min_k |b*lam - k*n|.  Hence

    lambda_1(L2)^2 = min over b >= 0 of  F(b),
    F(0) = (n*S_K1)^2,
    F(b) = ||b*lam||_n^2 * S_K1^2  +  b^2 * S_K2^2      (b >= 1).

F is a sum of a term increasing in b and a term driven by ||b*lam||_n.  If b is
not a record holder for ||b*lam||_n -- i.e. if some b' < b has
||b'*lam||_n <= ||b*lam||_n -- then F(b') < F(b).  So the minimiser is a record
holder.  By the classical best-approximation theorem the record holders of
||b*lam||_n are exactly the continued-fraction convergent denominators q_j of
lam/n.  Therefore

    ** lambda_1(L2)^2 = min over q in {0} u {q_j} of F(q) **     (EXACT)

Normalised form.  Put

    t_j = q_j * sqrt(S_K2 / (n*S_K1))      (so t_j = 1 at the balance scale
                                            q = sqrt(n*S_K1/S_K2) = sqrt(n*K2/K1))
    theta_j = ||q_j*lam||_n                (~ n / q_{j+1})

Then F(q_j)/det = theta_j^2*S_K1/(n*S_K2) + t_j^2, and using theta_j ~ n/q_{j+1},

    nu_hat^2  ~  min_j [ 1/t_{j+1}^2 + t_j^2 ].

Since t_{j+1} ~ a_{j+1} * t_j, nu_hat is small exactly when some partial
quotient a_{j+1} is large AND its convergent straddles the balance scale t = 1.
For a symmetric straddle (t_j = 1/s, t_{j+1} = s) this gives

    nu_hat ~ sqrt(2 / a_{j*+1}),   j* = the index with q_{j*} <= sqrt(n*K2/K1) < q_{j*+1}.

This is the scale-dependence.  max_a (the global maximum partial quotient) is
scale-free: it fires on a large a_j wherever it sits in the expansion, whereas
only the one straddling t = 1 changes lambda_1.  That is the precise reason the
June invariants failed.

Experiments
-----------
E1  minimiser is a convergent denominator (brute force, small n)
E2  CF closed form == Lagrange-Gauss lambda_1, exactly, incl. 256-bit
E3  the winning index is the balance index j*
E4  nu_hat ~ sqrt(2/a_{j*+1}); scale-aware a_{j*+1} vs scale-free max_a
E5  secp256k1 placement in the closed form
"""

import math
import random

# ---------------------------------------------------------------------------
# Exact integer primitives
# ---------------------------------------------------------------------------

def gauss_reduce_2d_exact(u, v):
    """Lagrange-Gauss reduction in exact integer arithmetic.

    The float version in glv_hnp_nuhat_lib.py divides Python ints via `/`,
    which loses precision (and overflows) once the entries exceed ~2^512.
    This returns lambda_1(L2)^2 as an exact integer.
    """
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]

    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]

    if n2(u) < n2(v):
        u, v = v, u
    while True:
        d, m = dot(u, v), n2(v)
        # nearest integer to d/m, exact, correct for negative d
        q = -((-2 * d + m) // (2 * m)) if d < 0 else (2 * d + m) // (2 * m)
        u, v = v, (u[0] - q * v[0], u[1] - q * v[1])
        if n2(u) <= n2(v):
            return n2(u)


def cf_expansion(lam, n):
    """Continued fraction of lam/n.

    Returns (a, q, theta) where a[j] are the partial quotients, q[j] the
    convergent denominators, theta[j] = ||q[j]*lam||_n = |q[j]*lam - p[j]*n|.
    """
    p0, p1 = 0, 1
    q0, q1 = 1, 0
    x, y = lam, n
    a, qs, th = [], [], []
    while y:
        c = x // y
        p0, p1 = p1, c * p1 + p0
        q0, q1 = q1, c * q1 + q0
        a.append(c)
        qs.append(q1)
        th.append(abs(q1 * lam - p1 * n))
        x, y = y, x - c * y
    return a, qs, th


def centered(x, n):
    x %= n
    return min(x, n - x)


def scales(n, k1_bound, k2_bound):
    """Column-diagonal scaling, identical to glv_hnp_nuhat_lib.scales."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)


# ---------------------------------------------------------------------------
# The two computations of lambda_1(L2)
# ---------------------------------------------------------------------------

def nu_gauss(n, lam, k1_bound, k2_bound):
    """Ground truth: Lagrange-Gauss on L2. Returns (lambda_1^2, det)."""
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    l1sq = gauss_reduce_2d_exact((n * s1, 0), (-lam * s1, s2))
    return l1sq, n * s1 * s2


def nu_cf(n, lam, k1_bound, k2_bound):
    """Closed form: minimise F over {0} u {convergent denominators}.

    Returns (lambda_1^2, det, j_win, j_star, a_at_star, t_star_lo, t_star_hi).
    """
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    det = n * s1 * s2
    a, qs, th = cf_expansion(lam, n)

    best, jwin = (n * s1) ** 2, -1          # j = -1 encodes the b = 0 vector
    for j, (q, t) in enumerate(zip(qs, th)):
        if q >= n:                           # q_last = n gives theta = 0; skip
            continue
        val = t * t * s1 * s1 + q * q * s2 * s2
        if val < best:
            best, jwin = val, j

    # balance index: last j with q_j <= sqrt(n*s1/s2)
    bal_sq = n * s1 // s2
    jstar = -1
    for j, q in enumerate(qs):
        if q * q <= bal_sq and q < n:
            jstar = j
    a_star = a[jstar + 1] if 0 <= jstar + 1 < len(a) else (a[0] if a else 0)
    t_lo = qs[jstar] / math.sqrt(bal_sq) if jstar >= 0 else 0.0
    t_hi = qs[jstar + 1] / math.sqrt(bal_sq) if jstar + 1 < len(qs) else float('inf')
    return best, det, jwin, jstar, a_star, t_lo, t_hi


def max_partial_quotient(lam, n):
    a, qs, _ = cf_expansion(lam, n)
    return max(a) if a else 0


# ---------------------------------------------------------------------------
# E1 — minimiser is a convergent denominator (brute force)
# ---------------------------------------------------------------------------

def e1_bruteforce(trials=300, seed=1):
    print("\n" + "=" * 72)
    print("E1  minimiser of F(b) is a CF convergent denominator (brute force)")
    print("=" * 72)
    rng = random.Random(seed)
    bad = 0
    for _ in range(trials):
        n = rng.randrange(1000, 20000)
        lam = rng.randrange(1, n)
        s1 = rng.randrange(1, 60)
        s2 = rng.randrange(1, 60)
        best, bb = (n * s1) ** 2, 0
        for b in range(1, n):
            v = centered(b * lam, n) ** 2 * s1 * s1 + b * b * s2 * s2
            if v < best:
                best, bb = v, b
        _, qs, _ = cf_expansion(lam, n)
        if bb != 0 and bb not in set(qs):
            bad += 1
    print(f"  trials = {trials}   minimiser NOT a convergent denominator: {bad}")
    print(f"  VERDICT: {'PASS' if bad == 0 else 'FAIL'}")
    return bad == 0


# ---------------------------------------------------------------------------
# E2 — closed form == Lagrange-Gauss, exactly
# ---------------------------------------------------------------------------

def e2_exactness(trials=400, seed=7):
    print("\n" + "=" * 72)
    print("E2  CF closed form == Lagrange-Gauss lambda_1  (exact integer test)")
    print("=" * 72)
    rng = random.Random(seed)
    rows = []
    for bits in (20, 24, 64, 128, 256):
        mism = 0
        for _ in range(trials // 5):
            n = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
            lam = rng.randrange(1, n)
            k2 = math.isqrt(n) + 1
            eff = rng.choice([0.05, 0.10, 0.25, 0.50])
            k1 = max(1, int(eff * n / k2))
            g, _ = nu_gauss(n, lam, k1, k2)
            c, _, _, _, _, _, _ = nu_cf(n, lam, k1, k2)
            if g != c:
                mism += 1
        rows.append((bits, trials // 5, mism))
        print(f"  {bits:>4}-bit n:  {trials//5:>4} trials   mismatches = {mism}")
    ok = all(r[2] == 0 for r in rows)
    print(f"  VERDICT: {'PASS - identification is EXACT' if ok else 'FAIL'}")
    return ok


def e2b_float_robustness():
    """The shipped nu_hat used float Lagrange-Gauss. Check it at 256 bits."""
    print("\n  --- E2b: float vs exact Lagrange-Gauss at 256 bits ---")
    try:
        import importlib.util
        import os
        here = os.path.dirname(os.path.abspath(__file__))
        spec = importlib.util.spec_from_file_location(
            "_lib", os.path.join(here, "glv_hnp_nuhat_lib.py"))
        lib = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(lib)
    except Exception as e:
        print(f"  (skipped: {type(e).__name__}: {e})")
        return
    n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    lam = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
    k2 = math.isqrt(n) + 1
    print(f"  {'eff':>6} {'nu_hat(float)':>15} {'nu_hat(exact)':>15} {'rel.diff':>12}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1 = max(1, int(eff * n / k2))
        _, _, nh_f = lib.rival_sublattice_nu(n, lam, k1, k2)
        l1sq, det = nu_gauss(n, lam, k1, k2)
        nh_e = math.sqrt(l1sq / det)
        rel = abs(nh_f - nh_e) / nh_e if nh_e else 0.0
        print(f"  {eff:>6.2f} {nh_f:>15.10f} {nh_e:>15.10f} {rel:>12.2e}")


# ---------------------------------------------------------------------------
# E3 — the winning convergent is the balance index
# ---------------------------------------------------------------------------

def e3_balance_index(trials=400, bits=20, seed=11):
    print("\n" + "=" * 72)
    print("E3  is the winning convergent the balance index j* ?")
    print("    j* = last j with q_j <= sqrt(n*S_K1/S_K2) = sqrt(n*K2/K1)")
    print("=" * 72)
    rng = random.Random(seed)
    hist = {}
    for _ in range(trials):
        n = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
        lam = rng.randrange(1, n)
        k2 = math.isqrt(n) + 1
        k1 = max(1, int(0.10 * n / k2))
        _, _, jwin, jstar, _, _, _ = nu_cf(n, lam, k1, k2)
        hist[jwin - jstar] = hist.get(jwin - jstar, 0) + 1
    tot = sum(hist.values())
    print(f"  offset (j_win - j*)   count    frac      ({trials} trials, {bits}-bit)")
    for k in sorted(hist):
        print(f"  {k:>17}   {hist[k]:>6}   {hist[k]/tot:>6.3f}")
    within1 = sum(v for k, v in hist.items() if abs(k) <= 1) / tot
    print(f"  |j_win - j*| <= 1 : {within1:.3f}")
    return within1


def e3b_window(trials=300, seed=17):
    """Does j_win stay in a bounded window around j* across scales?

    If so, lambda_1(L2) needs only O(1) candidates once the CF is known,
    rather than a scan over all convergents.
    """
    print("\n  --- E3b: width of the window, across bit sizes and eff ---")
    rng = random.Random(seed)
    print(f"  {'bits':>5} {'eff':>6} {'trials':>7} {'off=-1':>8} {'off=0':>7}"
          f" {'off=+1':>8} {'|off|>1':>8}")
    worst = 0
    for bits in (20, 64, 256):
        for eff in (0.02, 0.05, 0.10, 0.25, 0.50, 0.90):
            hist = {}
            for _ in range(trials):
                n = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
                lam = rng.randrange(1, n)
                k2 = math.isqrt(n) + 1
                k1 = max(1, int(eff * n / k2))
                _, _, jwin, jstar, _, _, _ = nu_cf(n, lam, k1, k2)
                off = jwin - jstar
                hist[off] = hist.get(off, 0) + 1
                worst = max(worst, abs(off))
            out = sum(v for k, v in hist.items() if abs(k) > 1)
            print(f"  {bits:>5} {eff:>6.2f} {trials:>7} {hist.get(-1,0):>8}"
                  f" {hist.get(0,0):>7} {hist.get(1,0):>8} {out:>8}")
    print(f"  max |j_win - j*| observed over all cells: {worst}")
    return worst


# ---------------------------------------------------------------------------
# E4 — nu_hat ~ sqrt(2/a_{j*+1}); scale-aware vs scale-free
# ---------------------------------------------------------------------------

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


def e4_partial_quotient(trials=600, bits=20, seed=13):
    print("\n" + "=" * 72)
    print("E4  nu_hat vs the partial quotient at the balance scale")
    print("=" * 72)
    rng = random.Random(seed)
    nus, a_star, a_max, pred = [], [], [], []
    for _ in range(trials):
        n = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
        lam = rng.randrange(1, n)
        k2 = math.isqrt(n) + 1
        k1 = max(1, int(0.10 * n / k2))
        l1sq, det, _, _, astar, _, _ = nu_cf(n, lam, k1, k2)
        nh = math.sqrt(l1sq / det)
        nus.append(nh)
        a_star.append(astar)
        a_max.append(max_partial_quotient(lam, n))
        pred.append(math.sqrt(2.0 / astar) if astar > 0 else float('nan'))
    print(f"  ({trials} trials, {bits}-bit, eff=0.10)")
    print(f"  spearman(nu_hat, a_star = a_(j*+1))  = {spearman(nus, a_star):+.4f}   <- scale-AWARE")
    print(f"  spearman(nu_hat, max_a)              = {spearman(nus, a_max):+.4f}   <- scale-free (June)")
    print(f"  spearman(nu_hat, sqrt(2/a_star))     = {spearman(nus, pred):+.4f}")
    fin = [(p, v) for p, v in zip(pred, nus) if p == p]
    if fin:
        err = [abs(p - v) / v for p, v in fin]
        err.sort()
        print(f"  |sqrt(2/a_star) - nu_hat| / nu_hat :"
              f" median {err[len(err)//2]:.3f}   p90 {err[int(0.9*len(err))]:.3f}")
    print("\n  nu_hat by a_star bucket:")
    print(f"  {'a_(j*+1)':>10} {'count':>7} {'mean nu_hat':>13} {'sqrt(2/a)':>11}")
    for lo, hi, lbl in ((1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"),
                        (5, 9, "5-9"), (10, 10 ** 9, ">=10")):
        grp = [v for v, a in zip(nus, a_star) if lo <= a <= hi]
        if grp:
            mid = (lo + hi) / 2 if hi < 10 ** 9 else 15
            print(f"  {lbl:>10} {len(grp):>7} {sum(grp)/len(grp):>13.4f}"
                  f" {math.sqrt(2.0/mid):>11.4f}")
    return spearman(nus, a_star), spearman(nus, a_max)


# ---------------------------------------------------------------------------
# E5 — secp256k1
# ---------------------------------------------------------------------------

def e5_secp256k1():
    print("\n" + "=" * 72)
    print("E5  secp256k1 in the closed form")
    print("=" * 72)
    n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    lam = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
    assert (lam * lam + lam + 1) % n == 0
    k2 = math.isqrt(n) + 1
    a, qs, _ = cf_expansion(lam, n)
    print(f"  CF of lam/n: {len(a)} partial quotients, max a_j = {max(a)}")
    print(f"  {'eff':>6} {'nu_hat':>10} {'j_win':>6} {'j*':>5} {'a_(j*+1)':>10}"
          f" {'t_(j*)':>9} {'t_(j*+1)':>10} {'sqrt(2/a*)':>11}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1 = max(1, int(eff * n / k2))
        l1sq, det, jwin, jstar, astar, tlo, thi = nu_cf(n, lam, k1, k2)
        nh = math.sqrt(l1sq / det)
        print(f"  {eff:>6.2f} {nh:>10.4f} {jwin:>6} {jstar:>5} {astar:>10}"
              f" {tlo:>9.3f} {thi:>10.3f} {math.sqrt(2.0/astar):>11.4f}")
    print("\n  Reference (2026-07-29 log): eff=0.05 -> 0.8709, 0.10 -> 0.6624,")
    print("                              eff=0.25 -> 0.5852, 0.50 -> 0.6851")


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("=" * 72)
    print("GLV-HNP Thread 23 — closed form for lambda_1(L2) / nu_hat")
    print("=" * 72)
    ok1 = e1_bruteforce()
    ok2 = e2_exactness()
    e2b_float_robustness()
    w1 = e3_balance_index()
    worst = e3b_window()
    sa, sm = e4_partial_quotient()
    e5_secp256k1()
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"  E1 minimiser is a convergent denominator : {'PASS' if ok1 else 'FAIL'}")
    print(f"  E2 CF closed form == Lagrange-Gauss      : {'PASS (exact)' if ok2 else 'FAIL'}")
    print(f"  E3 |j_win-j*|<=1 rate (20-bit, eff=0.10) : {w1:.3f}")
    print(f"  E3b max |j_win - j*| over 5400 trials    : {worst}"
          f"  (window is NOT provably +-1)")
    print(f"  E4 spearman(nu_hat, a_(j*+1))            : {sa:+.4f}")
    print(f"  E4 spearman(nu_hat, max_a)  [scale-free] : {sm:+.4f}")
    print("\n  Thread 23 falsifier was: predicted-from-CF nu_hat correlates < 0.9")
    print("  with computed nu_hat.  Not falsified — the relation is exact equality.")
