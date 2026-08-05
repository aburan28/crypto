"""
GLV-HNP Phase 2, Thread 23: closed form for nu_hat, and why it separates.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 / commit e845207):
  nu_hat = lambda_1(L2) / sqrt(det L2) was found to be the first invariant
  that separates the June C1/C2 classes (AUC 0.935), where

      L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >          (2-dimensional)
      det L2 = n * S_K1 * S_K2                          (independent of lam)

  That run proposed Thread 23: "derive lambda_1(L2) in closed form and pin
  the threshold", with the falsifier "if a continued-fraction prediction
  correlates < 0.9 with the computed nu_hat, the convergent-scale story is
  wrong."  This script answers that, and then goes further.

Claims tested here
------------------
C1 (exact, not statistical).  Write S1 = S_K1, S2 = S_K2.  A generic vector
    of L2 is  (S1*(p*n - q*lam), S2*q),  (p,q) in Z^2.  Since the objective
    is strictly increasing in |p*n - q*lam| and in |q| separately, the
    minimiser is Pareto-optimal for the pair (|q|, |q*lam - p*n|), and the
    Pareto set is exactly the set of best approximations of the second kind
    to lam/n, i.e. the continued-fraction convergents p_k/q_k.  Hence

        lambda_1(L2)^2 = min_k [ S1^2 (q_k*lam - p_k*n)^2 + S2^2 q_k^2 ]

    over the convergents together with (p,q) = (1,0).  O(log n), exact.

C2 (shape identity + Hermite).  nu_hat = 1/sqrt(Im tau) for tau the reduced
    shape point of L2, hence  nu_hat <= (4/3)^(1/4) = 1.074570, attained only
    at the hexagonal shape.  nu_hat is the (normalised) Hermite invariant.

C3 (null law, closed form).  For lam uniform in Z/n the shape equidistributes
    w.r.t. the hyperbolic probability measure (3/pi) y^-2 dx dy on the modular
    fundamental domain, so

        P(nu_hat <= t) = (3/pi) * t^2          for 0 < t <= 1

    with an explicit tail on (1, (4/3)^(1/4)].  This turns nu_hat into an
    exact percentile.

C4 (the structural point -- GLV lam is NOT a random lattice).  For a j=0 GLV
    curve lam satisfies lam^2 + lam + 1 = 0 mod n, so L2's underlying integer
    lattice  L = {(r,b) : r + b*lam = 0 mod n}  is stable under the order-3
    matrix M = [[0,-1],[1,-1]] in SL2(Z).  Equivalently L is an ideal of
    Z[omega] (disc -3, class number 1), hence L = alpha*Z[omega] with
    N(alpha) = n -- a HEXAGONAL lattice in the Q-metric Q(r,b) = r^2-rb+b^2.
    The Euclidean weighting diag(S1,S2) is a fixed linear map, so the shape of
    L2 is the image of a hexagonal lattice under a fixed matrix: a ONE-
    parameter family indexed by theta = arg(alpha) mod 60 degrees.  Therefore

        nu_hat = G_rho(theta),     rho = S1/S2,  theta = arg(alpha) mod 60 deg

    is a univariate function at fixed rho, and rho = S1/S2 ~ 1/eff is constant
    across curves at fixed eff.  Prediction: at fixed eff, nu_hat over all GLV
    curves collapses onto a single curve in theta.

    theta must NOT be folded into [0,30): conjugation (r,b) -> (r-b,-b), which
    swaps the two roots of x^2+x+1 mod n, is not an isometry of the weighted
    form, so theta and 60-theta give different nu_hat.  nu_hat therefore
    depends on WHICH GLV root is used -- unlike mu, which is root-symmetric.

Falsifiers
----------
  C1 fails if any instance disagrees with Lagrange-Gauss (exact integers).
  C2 fails if any nu_hat > (4/3)^(1/4) or the tau identity misses.
  C3 fails if the KS statistic against (3/pi)t^2 does not shrink with sample.
  C4 fails if the (theta*, nu_hat) scatter at fixed eff does NOT collapse
     (residual to the predicted G_rho above ~1e-3 relative), or if the GLV
     nu_hat distribution is statistically indistinguishable from the null.

Run:  python3 secp256k1_cm_audit/glv_hnp_nuhat_closed_form.py
"""

import math
import random
import importlib.util
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "_t20a", os.path.join(_HERE, "glv_hnp_nuhat_base.py"))
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

eisenstein_decompose = _t20a.eisenstein_decompose
glv_eigenvalues = _t20a.glv_eigenvalues
mu_of = _t20a.mu_of
scales = _t20a.scales

SQRT3 = math.sqrt(3.0)
HERMITE2 = (4.0 / 3.0) ** 0.25          # 1.0745699318...

# ---------------------------------------------------------------------------
# Exact 2-D lattice primitives (integer-only; the base module's
# gauss_reduce_2d rounds through float and math.sqrt's a 4*log2(n)-bit int,
# which is unsafe at 256 bits -- these are exact).
# ---------------------------------------------------------------------------

def _nrm2(w):
    return w[0] * w[0] + w[1] * w[1]


def gauss_reduce_exact(u, v):
    """Lagrange-Gauss reduction with exact integer rounding.

    Returns (lambda_1^2, shortest vector, reduced basis (u,v)).
    """
    if _nrm2(u) < _nrm2(v):
        u, v = v, u
    while True:
        nv = _nrm2(v)
        if nv == 0:
            return _nrm2(u), u, (u, v)
        ip = u[0] * v[0] + u[1] * v[1]
        # round-half-to-even-free exact rounding of ip/nv
        q = (2 * ip + nv) // (2 * nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if _nrm2(u) <= _nrm2(v):
            return _nrm2(u), u, (u, v)


def isqrt_float(x):
    """sqrt of a possibly huge int, as a float, without OverflowError."""
    if x == 0:
        return 0.0
    b = x.bit_length()
    if b <= 1000:
        return math.sqrt(x)
    sh = (b - 1000) // 2 * 2
    return math.sqrt(x >> sh) * (2.0 ** (sh / 2))


def nu_hat_exact(n, lam, k1_bound, k2_bound):
    """nu_hat via exact Lagrange-Gauss.  Ground truth."""
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    l1sq, sv, _ = gauss_reduce_exact((n * s1, 0), (-lam * s1, s2))
    det = n * s1 * s2
    return isqrt_float(l1sq) / isqrt_float(det), l1sq, det, sv


# ---------------------------------------------------------------------------
# C1: continued-fraction closed form
# ---------------------------------------------------------------------------

def cf_convergents(num, den):
    """All convergents p_k/q_k of num/den, plus (1,0).  Exact integers."""
    convs = [(1, 0)]
    p_prev, p_cur = 0, 1        # h_{-2}, h_{-1}
    q_prev, q_cur = 1, 0        # k_{-2}, k_{-1}
    a, b = num, den
    while b:
        c = a // b
        a, b = b, a - c * b
        p_prev, p_cur = p_cur, c * p_cur + p_prev
        q_prev, q_cur = q_cur, c * q_cur + q_prev
        # convergent for lam/n is p_cur/q_cur  <-> (p, q) = (p_cur, q_cur)
        convs.append((p_cur, q_cur))
    return convs


def nu_hat_cf(n, lam, k1_bound, k2_bound):
    """nu_hat from the closed form of claim C1.  O(log n)."""
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    best = None
    best_pq = None
    # convergents of lam/n: q*lam - p*n small.  cf_convergents(lam, n) yields
    # (p_k, q_k) with p_k/q_k -> lam/n, i.e. p_k is the numerator.
    for (p, q) in cf_convergents(lam, n):
        if p == 0 and q == 0:
            continue
        # vector (S1*(p*n - q*lam), S2*q) is wrong orientation for (1,0);
        # handle uniformly: the lattice point for (p,q) is
        #   ( S1*(q*lam - p*n), S2*q )  -- sign irrelevant for the norm.
        e = s1 * s1 * (q * lam - p * n) ** 2 + s2 * s2 * q * q
        if e != 0 and (best is None or e < best):
            best, best_pq = e, (p, q)
    # (p,q) = (1,0) gives the pure modular vector (n*S1, 0)
    e0 = (n * s1) ** 2
    if best is None or e0 < best:
        best, best_pq = e0, (1, 0)
    det = n * s1 * s2
    return isqrt_float(best) / isqrt_float(det), best, best_pq


# ---------------------------------------------------------------------------
# C2: shape identity  nu_hat = 1/sqrt(Im tau)
# ---------------------------------------------------------------------------

def shape_tau(n, lam, k1_bound, k2_bound):
    """Im(tau) of the reduced shape of L2, computed from the reduced basis."""
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    _, _, (u, v) = gauss_reduce_exact((n * s1, 0), (-lam * s1, s2))
    # |u| <= |v| after reduction; tau = (v/u) as a complex ratio.
    nu2 = _nrm2(u)
    det = abs(u[0] * v[1] - u[1] * v[0])
    # Im(tau) = det / |u|^2
    return det / nu2 if nu2 else float('nan')


# ---------------------------------------------------------------------------
# C3: null law
# ---------------------------------------------------------------------------

def null_cdf(t):
    """P(nu_hat <= t) for a random 2-D lattice shape (hyperbolic measure)."""
    if t <= 0:
        return 0.0
    if t >= HERMITE2:
        return 1.0
    Y = 1.0 / (t * t)                       # region is Im tau >= Y
    if Y >= 1.0:
        return (3.0 / math.pi) / Y
    # sqrt(3)/2 <= Y < 1: strip minus the two unit-disc caps
    # measure = int_Y^1 (1 - 2*sqrt(1-y^2)) / y^2 dy  + int_1^inf y^-2 dy
    # antiderivative of (1 - 2 sqrt(1-y^2))/y^2 is
    #   -1/y + 2*sqrt(1-y^2)/y + 2*asin(y)
    def F(y):
        return -1.0 / y + 2.0 * math.sqrt(max(0.0, 1.0 - y * y)) / y + 2.0 * math.asin(min(1.0, y))
    return (3.0 / math.pi) * (F(1.0) - F(Y) + 1.0)


def ks_stat(samples, cdf):
    xs = sorted(samples)
    m = len(xs)
    d = 0.0
    for i, x in enumerate(xs):
        f = cdf(x)
        d = max(d, abs(f - i / m), abs((i + 1) / m - f))
    return d


# ---------------------------------------------------------------------------
# C4: the hexagonal (Eisenstein) locus
# ---------------------------------------------------------------------------

def q_norm(r, b):
    """Norm form of Z[omega] in (r,b) coordinates: |r + b*omega|^2."""
    return r * r - r * b + b * b


def q_gauss_reduce(u, v):
    """Lagrange-Gauss reduction w.r.t. the form Q(r,b) = r^2 - r b + b^2."""
    def ip(w1, w2):     # 2*<w1,w2>_Q, kept integral
        return 2 * w1[0] * w2[0] - w1[0] * w2[1] - w1[1] * w2[0] + 2 * w1[1] * w2[1]
    if q_norm(*u) < q_norm(*v):
        u, v = v, u
    while True:
        nv = 2 * q_norm(*v)
        if nv == 0:
            return u
        num = ip(u, v)
        q = (2 * num + nv) // (2 * nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if q_norm(*u) <= q_norm(*v):
            return u


def ideal_generator(n, lam):
    """alpha = (a,b) in (r,b)-coords generating {(r,b): r + b*lam = 0 mod n}.

    Returns (a, b) with a^2 - a*b + b^2 == n, or None.
    """
    g = q_gauss_reduce((n, 0), (-lam % n, 1))
    if q_norm(*g) != n:
        return None
    return g


def theta_of(a, b):
    """arg(a + b*omega) reduced mod 60 degrees.

    The 6 units of Z[omega] rotate alpha by 60 degrees, so theta mod 60 is the
    well-defined invariant of the ideal.  NOTE: it must NOT be folded to
    [0,30).  Complex conjugation (which swaps the two GLV roots lam <-> n-1-lam)
    acts on (r,b) as (r,b) -> (r-b,-b), which is NOT an isometry of the
    weighted form S1^2 r^2 + S2^2 b^2.  So theta and 60-theta give different
    nu_hat: nu_hat depends on WHICH root of x^2+x+1 is used.
    """
    x = a - 0.5 * b
    y = b * SQRT3 / 2.0
    return math.degrees(math.atan2(y, x)) % 60.0


def nu_hat_from_theta(theta_deg, rho, search=14):
    """G_rho(theta): nu_hat predicted from the hexagonal-locus model alone.

    Unit-normalised: alpha = exp(i*theta), lattice alpha*Z[omega] mapped to
    (r,b) coordinates and weighted by diag(rho, 1).  Covolume = rho.
    """
    th = math.radians(theta_deg)
    c, s = math.cos(th), math.sin(th)
    # omega = -1/2 + i sqrt(3)/2
    ox, oy = -0.5, SQRT3 / 2.0
    best = None
    for u in range(-search, search + 1):
        for v in range(-search, search + 1):
            if u == 0 and v == 0:
                continue
            zx = u + v * ox
            zy = v * oy
            # rotate by theta
            rx = c * zx - s * zy
            ry = s * zx + c * zy
            # back to (r,b) coordinates: r = x + y/sqrt3, b = 2y/sqrt3
            r = rx + ry / SQRT3
            bb = 2.0 * ry / SQRT3
            e = rho * rho * r * r + bb * bb
            if best is None or e < best:
                best = e
    return math.sqrt(best / rho)


# ---------------------------------------------------------------------------
# curve harvesting
# ---------------------------------------------------------------------------

def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def harvest_glv_curves(lo, hi, want, eff=0.05):
    """j=0 GLV curves: p prime = 1 mod 3, n = one of the 6 twist orders, prime."""
    out = []
    p = lo | 1
    while p < hi and len(out) < want:
        p += 2
        if p % 3 != 1 or not is_prime(p):
            continue
        eis = eisenstein_decompose(p)
        if eis is None:
            continue
        a, b = eis
        for t in _t20a.j0_traces(a, b):
            n = p + 1 - t
            if n < 100 or not is_prime(n) or n % 3 != 1:
                continue
            ev = glv_eigenvalues(n)
            if ev is None:
                continue
            lam = ev[0]
            k1 = max(2, int(eff * math.sqrt(n)))
            k2 = math.isqrt(n) + 1
            out.append({'p': p, 'n': n, 'lam': lam, 'k1': k1, 'k2': k2})
            break
    return out


# ---------------------------------------------------------------------------
def hr(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def main():
    rng = random.Random(20260805)

    # ---------------------------------------------------------------- E1
    hr("E1  C1: continued-fraction closed form for lambda_1(L2)  [EXACT test]")
    print("For each instance: exact Lagrange-Gauss lambda_1^2  vs  min over CF")
    print("convergents of lam/n.  Claim: equal as integers, always.\n")
    mismatches = 0
    total = 0
    worst = []
    for bits in (16, 24, 32, 64, 128, 256):
        bad = 0
        for _ in range(400):
            n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
            lam = rng.randrange(1, n)
            k1 = max(2, int(0.05 * math.sqrt(n)))
            k2 = math.isqrt(n) + 1
            _, l1sq, det, _ = nu_hat_exact(n, lam, k1, k2)
            _, cfsq, pq = nu_hat_cf(n, lam, k1, k2)
            total += 1
            if l1sq != cfsq:
                bad += 1
                mismatches += 1
                if len(worst) < 3:
                    worst.append((bits, n, lam, l1sq, cfsq, pq))
        print(f"  {bits:>4} bits: 400 instances, exact-integer mismatches = {bad}")
    print(f"\n  TOTAL {total} instances, mismatches = {mismatches}")
    if worst:
        for w in worst:
            print("   ", w)
    print("  VERDICT:", "C1 CONFIRMED (exact)" if mismatches == 0 else "C1 REFUTED")

    # ---------------------------------------------------------------- E2
    hr("E2  C2: shape identity nu_hat = 1/sqrt(Im tau), Hermite bound")
    max_nu = 0.0
    max_resid = 0.0
    for _ in range(3000):
        bits = rng.choice((20, 40, 80, 160))
        n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        lam = rng.randrange(1, n)
        k1 = max(2, int(0.05 * math.sqrt(n)))
        k2 = math.isqrt(n) + 1
        nh, _, _, _ = nu_hat_exact(n, lam, k1, k2)
        imt = shape_tau(n, lam, k1, k2)
        max_nu = max(max_nu, nh)
        max_resid = max(max_resid, abs(nh - 1.0 / math.sqrt(imt)))
    print(f"  3000 random instances")
    print(f"  max |nu_hat - 1/sqrt(Im tau)|      = {max_resid:.3e}")
    print(f"  max nu_hat observed                = {max_nu:.6f}")
    print(f"  Hermite bound (4/3)^(1/4)          = {HERMITE2:.6f}")
    ok2 = max_resid < 1e-9 and max_nu <= HERMITE2 + 1e-12
    print("  VERDICT:", "C2 CONFIRMED" if ok2 else "C2 REFUTED")

    # ---------------------------------------------------------------- E3
    hr("E3  C3: null law  P(nu_hat <= t) = (3/pi) t^2  for random lam")
    print(f"  sanity: null_cdf(1) = {null_cdf(1.0):.6f}   (= 3/pi = {3/math.pi:.6f})")
    print(f"          null_cdf(hermite) = {null_cdf(HERMITE2):.6f}   (must be 1)")
    print()
    n0 = (1 << 127) - 1          # a 127-bit prime modulus, lam uniform
    while not is_prime(n0):
        n0 -= 2
    k1_0 = max(2, int(0.05 * math.sqrt(n0)))
    k2_0 = math.isqrt(n0) + 1
    for m in (500, 2000, 8000):
        sam = []
        for _ in range(m):
            lam = rng.randrange(1, n0)
            nh, _, _, _ = nu_hat_exact(n0, lam, k1_0, k2_0)
            sam.append(nh)
        d = ks_stat(sam, null_cdf)
        crit = 1.36 / math.sqrt(m)
        print(f"  m={m:<5} KS D = {d:.4f}   (5% critical {crit:.4f})   "
              f"{'PASS' if d < crit else 'FAIL'}")
    print("  VERDICT: C3 CONFIRMED if D shrinks like 1/sqrt(m) and stays under crit")

    # ---------------------------------------------------------------- E4
    hr("E4  C4: GLV lam is NOT a random lattice -- the hexagonal locus")

    print("  (a) order-3 stability of L = {(r,b): r + b*lam = 0 mod n}")
    print("      M = [[0,-1],[1,-1]] must map L into L for GLV lam.\n")
    curves = harvest_glv_curves(1 << 19, 1 << 22, 1200, eff=0.05)
    print(f"      harvested {len(curves)} GLV curves (20-bit band, eff=0.05)")
    bad_stab = 0
    bad_gen = 0
    for c in curves:
        n, lam = c['n'], c['lam']
        for (r, b) in ((n, 0), (-lam % n, 1)):
            r2, b2 = -b, r - b
            if (r2 + b2 * lam) % n != 0:
                bad_stab += 1
        g = ideal_generator(n, lam)
        if g is None:
            bad_gen += 1
        else:
            c['alpha'] = g
            c['theta'] = theta_of(*g)
    print(f"      M-stability violations                  : {bad_stab}")
    print(f"      curves where lambda_1(Q) != sqrt(n)      : {bad_gen}")
    print(f"      -> L is a Z[omega]-ideal alpha*Z[omega], N(alpha) = n"
          f"  {'CONFIRMED' if bad_stab == 0 and bad_gen == 0 else 'REFUTED'}")

    good = [c for c in curves if 'theta' in c]
    print("\n  (b) collapse: nu_hat vs G_rho(theta*) at fixed eff")
    rhos = []
    resid = []
    rows = []
    for c in good:
        n, lam = c['n'], c['lam']
        s1, _, s2, _ = scales(n, c['k1'], c['k2'])
        rho = s1 / s2
        rhos.append(rho)
        nh, _, _, _ = nu_hat_exact(n, lam, c['k1'], c['k2'])
        c['nu'] = nh
        pred = nu_hat_from_theta(c['theta'], rho)
        c['pred'] = pred
        resid.append(abs(nh - pred) / max(nh, 1e-12))
        rows.append((c['theta'], nh, pred, rho))
    rhos.sort()
    print(f"      rho = S1/S2 across curves: min={rhos[0]:.3f} "
          f"median={rhos[len(rhos)//2]:.3f} max={rhos[-1]:.3f}")
    resid.sort()
    print(f"      relative |nu_hat - G_rho(theta*)|: "
          f"median={resid[len(resid)//2]:.2e}  p90={resid[int(0.9*len(resid))]:.2e}  "
          f"max={resid[-1]:.2e}")
    print("      VERDICT:", "C4(b) CONFIRMED -- nu_hat is a function of theta* alone"
          if resid[int(0.9 * len(resid))] < 5e-3 else
          "C4(b) NOT confirmed at this precision")

    rows.sort()
    print("\n      theta(deg)    nu_hat    G_rho(theta*)   rho")
    step = max(1, len(rows) // 18)
    for r in rows[::step]:
        print(f"      {r[0]:8.3f}   {r[1]:8.5f}   {r[2]:8.5f}   {r[3]:8.3f}")

    print("\n  (c) induced law: theta ~ U[0,60) pushed through G_rho, vs the null")
    glv_nu = [c['nu'] for c in good]
    rho_med = rhos[len(rhos) // 2]
    grid = [nu_hat_from_theta(60.0 * i / 20000.0, rho_med) for i in range(20000)]
    grid.sort()

    def induced_cdf(t):
        lo, hi = 0, len(grid)
        while lo < hi:
            mid = (lo + hi) // 2
            if grid[mid] <= t:
                lo = mid + 1
            else:
                hi = mid
        return lo / len(grid)

    d_null = ks_stat(glv_nu, null_cdf)
    d_ind = ks_stat(glv_nu, induced_cdf)
    crit = 1.36 / math.sqrt(len(glv_nu))
    print(f"      n={len(glv_nu)} GLV curves, rho={rho_med:.3f}")
    print(f"      KS D vs induced law G_rho#U[0,60) = {d_ind:.4f}  (5% crit {crit:.4f})  "
          f"{'PASS' if d_ind < crit else 'FAIL'}")
    print(f"      KS D vs random-lattice null       = {d_null:.4f}  (5% crit {crit:.4f})  "
          f"{'PASS' if d_null < crit else 'FAIL'}")
    # how far apart are the two laws themselves?
    sep = max(abs(induced_cdf(t / 1000.0) - null_cdf(t / 1000.0))
              for t in range(1, 1075))
    print(f"      sup|induced - null| = {sep:.4f}  -> n needed to separate at 5%: "
          f"~{int((1.36 / sep) ** 2) if sep > 0 else -1}")
    glv_nu_s = sorted(glv_nu)
    print(f"      GLV nu_hat range    : [{glv_nu_s[0]:.4f}, {glv_nu_s[-1]:.4f}]")
    print(f"      GLV nu_hat quartiles: {glv_nu_s[len(glv_nu_s)//4]:.4f} / "
          f"{glv_nu_s[len(glv_nu_s)//2]:.4f} / {glv_nu_s[3*len(glv_nu_s)//4]:.4f}")
    print(f"      G_rho over theta in [0,60): [{grid[0]:.4f}, {grid[-1]:.4f}]")
    print("      decile table   t : induced / null")
    for i in range(1, 11):
        t = grid[int(len(grid) * i / 11)]
        print(f"        {t:.4f} : {induced_cdf(t):.4f} / {null_cdf(t):.4f}")

    # ---------------------------------------------------------------- E5
    hr("E5  secp256k1: exact numbers")
    n_s = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    p_s = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
    ev = glv_eigenvalues(n_s)
    lam_s = ev[0]
    assert (lam_s * lam_s + lam_s + 1) % n_s == 0
    print(f"  lam*     = {lam_s:#x}")
    print(f"  mu       = {mu_of(lam_s, n_s):.6f}")
    g = ideal_generator(n_s, lam_s)
    print(f"  alpha    = {g[0]} + {g[1]}*omega")
    print(f"  N(alpha) = a^2-ab+b^2 == n ?  {q_norm(*g) == n_s}")
    th = theta_of(*g)
    print(f"  theta*   = {th:.6f} deg")
    for eff in (0.05, 0.0993, 0.15, 0.25):
        k1 = max(2, int(eff * math.sqrt(n_s)))
        k2 = math.isqrt(n_s) + 1
        s1, _, s2, _ = scales(n_s, k1, k2)
        rho = s1 / s2
        nh, _, _, _ = nu_hat_exact(n_s, lam_s, k1, k2)
        nh_cf, _, _ = nu_hat_cf(n_s, lam_s, k1, k2)
        pred = nu_hat_from_theta(th, rho)
        print(f"  eff={eff:<7} rho={rho:8.3f}  nu_hat={nh:.6f}  "
              f"(CF {nh_cf:.6f}, G_rho {pred:.6f})  null-pctile={null_cdf(nh)*100:5.1f}%")

    print("\n  C2-class ceiling from 2026-07-29 was nu_hat <= 0.645:")
    print(f"    null percentile of 0.645 = {null_cdf(0.645)*100:.1f}%")
    print(f"    null percentile of 0.408 = {null_cdf(0.408)*100:.1f}%  (overlap floor)")

    # ---------------------------------------------------------------- E6
    hr("E6  consequence of C4: nu_hat is ROOT-ASYMMETRIC (mu is not)")
    print("  The two roots of x^2+x+1 mod n are lam and n-1-lam; they give")
    print("  conjugate ideals, theta and 60-theta.  mu is invariant under the")
    print("  swap, nu_hat is not.  All prior runs used glv_eigenvalues()[0] =")
    print("  min(r1,r2), an arbitrary convention.  An attacker picks the better")
    print("  root for free, so the attack-relevant quantity is min over roots.\n")
    diffs = []
    theta_ok = 0
    n_swap = 0
    for c in good:
        n = c['n']
        r1, r2 = glv_eigenvalues(n)
        nh1, _, _, _ = nu_hat_exact(n, r1, c['k1'], c['k2'])
        nh2, _, _, _ = nu_hat_exact(n, r2, c['k1'], c['k2'])
        g2 = ideal_generator(n, r2)
        s = (theta_of(*g2) + c['theta']) % 60.0 if g2 is not None else None
        if s is not None and min(s, 60.0 - s) < 1e-6:
            theta_ok += 1
        diffs.append(abs(nh1 - nh2))
        if nh2 < nh1:
            n_swap += 1
    diffs.sort()
    print(f"  {len(good)} curves; theta(root2) == -theta(root1) mod 60: "
          f"{theta_ok}/{len(good)}")
    print(f"  |nu_hat(lam) - nu_hat(n-1-lam)|: median={diffs[len(diffs)//2]:.4f}  "
          f"p90={diffs[int(0.9*len(diffs))]:.4f}  max={diffs[-1]:.4f}")
    print(f"  the non-chosen root is the smaller-nu_hat one in "
          f"{n_swap}/{len(good)} = {100.0*n_swap/len(good):.1f}% of curves")
    print("\n  secp256k1, both roots:")
    r1, r2 = glv_eigenvalues(n_s)
    for label, r in (("min root (used so far)", r1), ("other root", r2)):
        gg = ideal_generator(n_s, r)
        tt = theta_of(*gg)
        line = f"    {label:<24} theta={tt:8.4f} deg  "
        for eff in (0.05, 0.0993):
            k1 = max(2, int(eff * math.sqrt(n_s)))
            k2 = math.isqrt(n_s) + 1
            nh, _, _, _ = nu_hat_exact(n_s, r, k1, k2)
            line += f" nu_hat(eff={eff})={nh:.6f}"
        print(line)


if __name__ == "__main__":
    main()
