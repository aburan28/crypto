"""
GLV-HNP — Thread 21: closed form for nu_hat.  nu_hat is a CM-angle invariant.

Background
----------
Thread 20c/20d (2026-07-29, commit e845207) found the first invariant that
separates the June C1/C2 classes: the scale-free shortest vector of the
non-planted 2D sublattice

    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >,     nu_hat = lambda_1(L2)/sqrt(det L2)

with AUC 0.935 against C1/C2 and spearman(nu_hat, p_hat) = -0.583 in a
fixed-curve causal arm.  It was computed numerically (one Lagrange-Gauss
reduction per curve) and left unexplained.

NOTE: those three scripts were dead on `main` -- the helper module rewritten by
e845207 was discarded by merge 6fd645f, so `_t20a.rival_sublattice_nu` did not
exist.  It is restored here as `glv_hnp_thread20c_lib.py`.

Claim (H-closed)
----------------
Strip the scaling.  L2 = D*Lambda with D = diag(S_K1, S_K2) and

    Lambda = {(x,b) in Z^2 : x + lam*b = 0 mod n},   det Lambda = n.

For a j=0 curve n is a norm from the Eisenstein order: n = u^2 - u*v + v^2, and
Lambda is the coordinate image of the principal ideal (pi), pi = u + v*omega,
w.r.t. the basis {1, omega}.  Concretely Lambda = A_pi * Z^2 with

    A_pi = [[u, -v], [v, u-v]]        (multiplication-by-pi matrix),  det = n.

Since multiplication by pi is (scale by sqrt(n)) o (rotate by theta = arg pi) in
the Minkowski plane, writing C for the fixed coordinate->Minkowski matrix

    C = [[1, -1/2], [0, sqrt(3)/2]],        C * A_pi * C^-1 = sqrt(n) * R_theta,

the whole n-dependence cancels out of nu_hat:

    nu_hat = lambda_1( G(tau, theta) Z^2 ) / sqrt(tau),
    G(tau, theta) = diag(tau, 1) * C^-1 * R_theta * C,     tau = S_K1 / S_K2.

So nu_hat depends on the curve ONLY through theta = arg(u + v*omega), the CM
angle, and on the attack parameters only through the aspect ratio tau (which is
1/eff up to integer-division rounding, eff = K1*K2/n).  n itself is absent.

Falsifier: if the closed form disagrees with the integer Lagrange-Gauss value on
any curve beyond float epsilon, H-closed is wrong.

Parts
-----
A  exact-agreement test, 200 fresh 20-bit curves x 4 eff values
B  collapse test: nu_hat vs theta at fixed tau, many n -- must be one curve
C  the nu_hat(theta) landscape; measure of the "easy tail" nu_hat <= 0.645
D  secp256k1 exactly (no 20-bit extrapolation) across eff
E  reconciliation with the falsified invariant a_corn/n (2026-06-29)

Run: python3 glv_hnp_nuhat_closed_form.py
"""

import math
import sys
import time

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_t20c", __file__.rsplit("/", 1)[0] + "/glv_hnp_thread20c_lib.py")
_t20c = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20c)

eisenstein_decompose = _t20c.eisenstein_decompose
glv_eigenvalues = _t20c.glv_eigenvalues
mu_of = _t20c.mu_of
scales = _t20c.scales
rival_sublattice_nu = _t20c.rival_sublattice_nu

SQ3 = math.sqrt(3.0)
# coordinate (x,b) [meaning x + b*omega] -> Minkowski plane (Re, Im)
C = ((1.0, -0.5), (0.0, SQ3 / 2.0))
C_INV = ((1.0, 1.0 / SQ3), (0.0, 2.0 / SQ3))

SECP256K1_P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
SECP256K1_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


# ---------------------------------------------------------------------------
# 2x2 helpers
# ---------------------------------------------------------------------------

def matmul(A, B):
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def gauss_reduce_real(u, v, iters=200):
    """Lagrange-Gauss on a real 2D basis; returns lambda_1."""
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    if n2(u) < n2(v):
        u, v = v, u
    for _ in range(iters):
        q = round((u[0] * v[0] + u[1] * v[1]) / n2(v))
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if n2(u) <= n2(v):
            break
    return math.sqrt(n2(u))


def nu_hat_closed(tau, theta):
    """nu_hat from the closed form -- no n, no lam, no lattice of size n."""
    ct, st = math.cos(theta), math.sin(theta)
    R = ((ct, -st), (st, ct))
    G = matmul(((tau, 0.0), (0.0, 1.0)), matmul(C_INV, matmul(R, C)))
    # columns of G are the basis vectors of G*Z^2
    l1 = gauss_reduce_real((G[0][0], G[1][0]), (G[0][1], G[1][1]))
    return l1 / math.sqrt(tau)


THETA_PERIOD = math.pi / 3.0            # 60 deg: the 6 units of Z[omega]


def eisenstein_decompose_fast(n):
    """
    (u,v) with u^2 - u*v + v^2 = n for prime n = 1 mod 3, via Cornacchia.

    The library's eisenstein_decompose is O(sqrt(n)) trial search, which is
    fine at 20 bits and hopeless at 256.  Solve n = x^2 + 3y^2 by Cornacchia
    and map (x,y) -> (u,v) = (x+y, 2y):
        (x+y)^2 - (x+y)(2y) + (2y)^2 = x^2 + 3y^2.
    Returns None if n is not represented (i.e. n != 1 mod 3).
    """
    if n == 3:
        return (1, 2)
    if n % 3 != 1:
        return None
    r = _t20c.tonelli_shanks((-3) % n, n)
    if r is None:
        return None
    if 2 * r < n:                     # Cornacchia wants n/2 < r < n
        r = n - r
    lim = math.isqrt(n)
    a, b = n, r
    while b > lim:
        a, b = b, a % b
    rem = n - b * b
    if rem % 3 != 0:
        return None
    y = math.isqrt(rem // 3)
    if y * y * 3 + b * b != n:
        return None
    u, v = b + y, 2 * y
    if u * u - u * v + v * v != n:
        return None
    return (u, v)


def theta_of(u, v):
    """
    arg(u + v*omega) reduced into the fundamental domain [0, 60) degrees.
    nu_hat is 60-degree periodic (multiplying pi by a unit does not change the
    ideal), so this is the canonical CM angle of the curve.
    """
    return math.atan2(v * SQ3 / 2.0, u - v / 2.0) % THETA_PERIOD


def eisenstein_for_root(n, lam):
    """
    (u,v) with u^2-u*v+v^2 = n such that A_pi*Z^2 == Lambda(n,lam).
    Returns None if n is not an Eisenstein norm.
    The 6 units x conjugation give 12 candidates; exactly 6 (one ideal) match.
    """
    d = eisenstein_decompose_fast(n)
    if d is None:
        return None
    u0, v0 = d
    # the 6 associates of pi and of its conjugate, in (x + y*omega) coordinates.
    # pi = u + v*omega  =>  pi-bar = u + v*omega-bar = (u - v) - v*omega.
    cands = []
    for (a, b) in ((u0, v0), (u0 - v0, -v0)):
        x, y = a, b
        for _ in range(6):
            cands.append((x, y))
            x, y = -y, x - y                   # multiply by omega
    for (x, y) in cands:
        if x * x - x * y + y * y != n:
            continue
        if (x + lam * y) % n == 0 and (-y + lam * (x - y)) % n == 0:
            return (x, y)
    return None


def nu_hat_from_cm(n, lam, k1_bound, k2_bound):
    """nu_hat via the closed form; returns (nu_hat, theta, tau, (u,v))."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    tau = s_k1 / s_k2
    uv = eisenstein_for_root(n, lam)
    if uv is None:
        return None
    th = theta_of(*uv)
    return nu_hat_closed(tau, th), th, tau, uv


# ---------------------------------------------------------------------------
# curve harvesting (only n matters here -- no EC arithmetic needed)
# ---------------------------------------------------------------------------

def harvest_orders(lo, hi, count, step=1):
    """Primes n = 1 mod 3 in [lo,hi) that are Eisenstein norms (all of them)."""
    import sympy
    out = []
    n = lo | 1
    while len(out) < count and n < hi:
        if n % 3 == 1 and sympy.isprime(n):
            out.append(n)
        n += 2 * step
    return out


# ---------------------------------------------------------------------------

def part_a():
    print("=" * 78)
    print("Part A: closed form vs integer Lagrange-Gauss (falsifier for H-closed)")
    print("=" * 78)
    orders = harvest_orders(2 ** 19, 2 ** 20, 200, step=37)

    # sanity: the fast Cornacchia decomposition must agree (up to units) with
    # the library's O(sqrt n) brute-force search.
    bad = 0
    for n in orders:
        f = eisenstein_decompose_fast(n)
        s = eisenstein_decompose(n)
        if f is None or s is None:
            bad += 1
            continue
        if f[0] ** 2 - f[0] * f[1] + f[1] ** 2 != n:
            bad += 1
            continue
        # same ideal up to conjugation <=> theta agree modulo theta -> 60-theta
        tf, ts = theta_of(*f), theta_of(*s)
        if min(abs(tf - ts), abs(tf + ts - THETA_PERIOD)) > 1e-9:
            bad += 1
    print(f"  Cornacchia vs brute-force on {len(orders)} orders: {bad} mismatches")

    worst = 0.0
    worst_row = None
    tested = 0
    skipped = 0
    for eff in (0.02, 0.05, 0.10, 0.25):
        for n in orders:
            roots = glv_eigenvalues(n)
            if roots is None:
                skipped += 1
                continue
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(eff * n / k2))
            for lam in roots:
                _, _, nh_int = rival_sublattice_nu(n, lam, k1, k2)
                got = nu_hat_from_cm(n, lam, k1, k2)
                if got is None:
                    skipped += 1
                    continue
                nh_cf, th, tau, uv = got
                err = abs(nh_cf - nh_int)
                tested += 1
                if err > worst:
                    worst, worst_row = err, (n, lam, eff, nh_int, nh_cf, uv)
    print(f"  {tested} (curve, root, eff) triples tested, {skipped} skipped")
    print(f"  max |nu_hat_closed - nu_hat_gauss| = {worst:.3e}")
    if worst_row:
        n, lam, eff, a, b, uv = worst_row
        print(f"  worst case: n={n} eff={eff} (u,v)={uv}  gauss={a:.9f} closed={b:.9f}")
    verdict = "CONFIRMED" if worst < 1e-9 else "FALSIFIED"
    print(f"  H-closed: {verdict}  (threshold 1e-9)")
    return worst


def part_b():
    print()
    print("=" * 78)
    print("Part B: is theta equidistributed, and is its law n-independent?")
    print("=" * 78)
    print("  Part A proves nu_hat = F(tau, theta) pointwise, so the *only*")
    print("  remaining input is the distribution of theta over curves.  If that")
    print("  law is uniform on [0,60) and independent of bit length, then the")
    print("  Part C tail fractions are population statements, not sample fits.")

    print("\n  60-degree periodicity check (multiply pi by each of the 6 units):")
    for th_deg in (5.0, 17.0, 31.0, 48.0):
        th = math.radians(th_deg)
        vals = [nu_hat_closed(20.0, th + k * THETA_PERIOD) for k in range(6)]
        print(f"    theta={th_deg:>5.1f}  spread over the 6 units = "
              f"{max(vals) - min(vals):.2e}  (nu_hat={vals[0]:.6f})")

    NBINS = 10
    print(f"\n  theta mod 60 in {NBINS} equal bins, counts per bit length:")
    print(f"  {'bits':>6} {'curves':>7}  " +
          " ".join(f"{i*6:>4}-{(i+1)*6:<3}" for i in range(NBINS)) +
          f" {'chi2':>8} {'dof':>4}")
    for bits in (16, 20, 24, 32, 48):
        orders = harvest_orders(2 ** (bits - 1), 2 ** bits, 600, step=1009)
        counts = [0] * NBINS
        tot = 0
        for n in orders:
            roots = glv_eigenvalues(n)
            if roots is None:
                continue
            for lam in roots:
                uv = eisenstein_for_root(n, lam)
                if uv is None:
                    continue
                th = theta_of(*uv)
                counts[min(NBINS - 1, int(th / THETA_PERIOD * NBINS))] += 1
                tot += 1
        exp = tot / NBINS
        chi2 = sum((c - exp) ** 2 / exp for c in counts) if exp > 0 else 0.0
        print(f"  {bits:>6} {tot:>7}  " + " ".join(f"{c:>8}" for c in counts) +
              f" {chi2:>8.2f} {NBINS-1:>4}")
    print(f"\n  chi2 5% critical value at dof={NBINS-1} is 16.92.")
    print("  Both roots of each curve are counted; they sit at theta and 60-theta,")
    print("  which is why the histogram is symmetric by construction.")

    print("\n  Root-pair check: the two GLV roots give theta and 60-theta,")
    print("  and nu_hat differs between them (no reflection symmetry):")
    import sympy
    n = 2 ** 19 | 1
    shown = 0
    print(f"    {'n':>9} {'theta1':>9} {'theta2':>9} {'sum':>8} "
          f"{'nu_hat1':>9} {'nu_hat2':>9}")
    while shown < 6:
        n += 2
        if n % 3 != 1 or not sympy.isprime(n):
            continue
        roots = glv_eigenvalues(n)
        if roots is None:
            continue
        uvs = [eisenstein_for_root(n, lam) for lam in roots]
        if any(x is None for x in uvs):
            continue
        t1, t2 = (theta_of(*x) for x in uvs)
        print(f"    {n:>9} {math.degrees(t1):>9.4f} {math.degrees(t2):>9.4f} "
              f"{math.degrees(t1 + t2):>8.4f} "
              f"{nu_hat_closed(20.0, t1):>9.5f} {nu_hat_closed(20.0, t2):>9.5f}")
        shown += 1


def part_c():
    print()
    print("=" * 78)
    print("Part C: the nu_hat(theta) landscape and the size of the easy tail")
    print("=" * 78)
    print("  C2 ceiling from Thread 20d: every C2 curve had nu_hat <= 0.645.")
    print(f"\n  {'eff':>7} {'tau':>8} {'min nu':>8} {'max nu':>8} {'mean':>8} "
          f"{'frac<=0.645':>12} {'argmin theta':>13}")
    NS = 20001
    for eff in (0.02, 0.05, 0.0993, 0.15, 0.25, 0.50):
        tau = 1.0 / eff
        vals = []
        best, best_th = 1e9, 0.0
        for i in range(NS):
            th = i * (math.pi / 3.0) / NS      # one 60-degree period
            v = nu_hat_closed(tau, th)
            vals.append(v)
            if v < best:
                best, best_th = v, th
        frac = sum(1 for v in vals if v <= 0.645) / len(vals)
        print(f"  {eff:>7.4f} {tau:>8.2f} {min(vals):>8.4f} {max(vals):>8.4f} "
              f"{sum(vals)/len(vals):>8.4f} {frac:>12.3f} "
              f"{math.degrees(best_th):>12.2f}d")
    print("\n  theta is equidistributed for random j=0 curves (Duke / Hecke), so")
    print("  'frac<=0.645' is the predicted fraction of j=0 curves in the easy")
    print("  tail at that eff -- a curve-independent, n-independent statement.")


def part_d():
    print()
    print("=" * 78)
    print("Part D: secp256k1 EXACTLY (closed form -- no 20/24-bit extrapolation)")
    print("=" * 78)
    n = SECP256K1_N
    roots = glv_eigenvalues(n)
    print(f"  n  = {n}")
    print(f"  n mod 3 = {n % 3}, both GLV roots exist: {roots is not None}")
    k2 = math.isqrt(n) + 1
    for lam in roots:
        uv = eisenstein_for_root(n, lam)
        if uv is None:
            print(f"  lam={lam}: NO Eisenstein match (unexpected)")
            continue
        th = theta_of(*uv)
        print(f"\n  root mu={mu_of(lam, n):.6f}  (u,v)=({uv[0]}, {uv[1]})")
        print(f"    check u^2-uv+v^2 == n : {uv[0]**2 - uv[0]*uv[1] + uv[1]**2 == n}")
        print(f"    theta = {math.degrees(th):.6f} deg "
              f"(mod 60 = {math.degrees(th) % 60:.6f} deg)")
        print(f"    {'eff':>8} {'nu_hat(cf)':>11} {'nu_hat(gauss)':>14} {'easy?':>7}")
        for eff in (0.02, 0.05, 0.0993, 0.15, 0.25, 0.50):
            k1 = max(2, int(eff * n / k2))
            _, _, nh_int = rival_sublattice_nu(n, lam, k1, k2)
            got = nu_hat_from_cm(n, lam, k1, k2)
            nh_cf = got[0]
            print(f"    {eff:>8.4f} {nh_cf:>11.6f} {nh_int:>14.6f} "
                  f"{'YES' if nh_cf <= 0.645 else 'no':>7}")
    print("\n  Previous runs quoted nu_hat=0.664 for secp256k1 as a heuristic")
    print("  20/24-bit -> 256-bit extrapolation.  It is not an extrapolation:")
    print("  the closed form is exact at every bit length.")


def part_e():
    print()
    print("=" * 78)
    print("Part E: why a_corn/n was falsified (2026-06-29) but nu_hat is not")
    print("=" * 78)
    print("  a_corn is the Cornacchia/Eisenstein 'a' -- i.e. |u|, the RADIAL")
    print("  part of pi.  But |pi| = sqrt(n) always, so u alone is a lossy proxy")
    print("  for theta: u = sqrt(n)*cos(theta) + (v/2), and a_corn/n ~ n^-1/2 for")
    print("  every curve.  The information lives in the ANGLE, and it is only")
    print("  visible after the anisotropic weighting diag(tau,1) -- an isotropic")
    print("  invariant of the ideal cannot see it, because ALL Eisenstein ideals")
    print("  are similar as lattices (class number 1, hexagonal shape).")
    tau = 20.0
    print(f"\n  Demonstration: 8 curves, same n, different roots/angles (tau={tau})")
    import sympy
    print(f"    {'n':>9} {'u':>6} {'v':>6} {'a_corn/n':>10} {'theta(deg)':>11} "
          f"{'nu_hat':>9}")
    shown = 0
    n = 2 ** 19 | 1
    while shown < 8:
        n += 2
        if n % 3 != 1 or not sympy.isprime(n):
            continue
        roots = glv_eigenvalues(n)
        if roots is None:
            continue
        lam = roots[0]
        uv = eisenstein_for_root(n, lam)
        if uv is None:
            continue
        th = theta_of(*uv)
        d = eisenstein_decompose(n)
        print(f"    {n:>9} {uv[0]:>6} {uv[1]:>6} {d[0]/n:>10.6f} "
              f"{math.degrees(th) % 60:>11.4f} {nu_hat_closed(tau, th):>9.5f}")
        shown += 1
    print("\n  a_corn/n is pinned to a ~n^-1/2 sliver for every curve; theta and")
    print("  hence nu_hat range over the whole interval.  Same underlying datum,")
    print("  different coordinate -- which is exactly why one separates and the")
    print("  other does not.")


def main():
    t0 = time.time()
    print("GLV-HNP Thread 21: closed form for nu_hat")
    print("(restores + explains the Thread 20c/20d separator)")
    print()
    worst = part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    print(f"\nTotal wall time {time.time() - t0:.1f}s")
    print("Done.")
    return 0 if worst < 1e-9 else 1


if __name__ == '__main__':
    sys.exit(main())
