"""
Pure-Python finite-field elliptic-curve + low-degree isogeny primitives.

No Sage/PARI/sympy dependency: everything here is plain-int modular arithmetic,
which is enough to run Phase 0 (invariants/sanity) and a Phase 1 small-degree
isogeny crawl at full P-256 scale.

Conventions
-----------
Short Weierstrass curves  E: y^2 = x^3 + a*x + b  over F_p (p > 3 prime).
Affine points are (x, y) tuples; the point at infinity is the Python object
``INF`` (a unique sentinel).

Primitives provided
--------------------
* modular sqrt (Tonelli-Shanks), inverse
* EC point add / double / scalar-mul
* j-invariant, discriminant
* univariate polynomial root finding mod p  (x^p-x gcd split + Cantor-Zassenhaus)
* 2-isogeny via Velu from a rational 2-torsion point
* quadratic twist + same-order twist selection by a random-point order test

The 2-isogeny path is exact (codomain built by Velu), so an l-isogenous
neighbour automatically has #E'(F_p) = #E(F_p); no point counting needed.
"""

from __future__ import annotations
import random

INF = ("INF",)  # unique sentinel for the point at infinity


# --------------------------------------------------------------------------
# modular arithmetic
# --------------------------------------------------------------------------

def inv_mod(a: int, p: int) -> int:
    return pow(a % p, -1, p)


def sqrt_mod(a: int, p: int):
    """Tonelli-Shanks. Returns a square root of a mod p, or None if non-residue."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    # general Tonelli-Shanks
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = (t2 * t2) % p
            i += 1
            if i == m:
                return None
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, (b * b) % p, (t * b * b) % p, (r * b) % p
    return r


def is_square(a: int, p: int) -> bool:
    a %= p
    return a == 0 or pow(a, (p - 1) // 2, p) == 1


# --------------------------------------------------------------------------
# elliptic-curve group law  E: y^2 = x^3 + a x + b  over F_p
# --------------------------------------------------------------------------

class Curve:
    def __init__(self, a: int, b: int, p: int):
        self.a = a % p
        self.b = b % p
        self.p = p
        if (4 * self.a ** 3 + 27 * self.b ** 2) % p == 0:
            raise ValueError("singular curve (discriminant 0)")

    # discriminant and j-invariant
    def disc(self) -> int:
        p = self.p
        return (-16 * (4 * self.a ** 3 + 27 * self.b ** 2)) % p

    def j_invariant(self) -> int:
        p, a, b = self.p, self.a, self.b
        num = (1728 * 4 * a ** 3) % p
        den = (4 * a ** 3 + 27 * b ** 2) % p
        return (num * inv_mod(den, p)) % p

    def is_on_curve(self, P) -> bool:
        if P is INF:
            return True
        x, y = P
        p = self.p
        return (y * y - (x * x * x + self.a * x + self.b)) % p == 0

    def neg(self, P):
        if P is INF:
            return INF
        x, y = P
        return (x, (-y) % self.p)

    def add(self, P, Q):
        if P is INF:
            return Q
        if Q is INF:
            return P
        p = self.p
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2 and (y1 + y2) % p == 0:
            return INF
        if P == Q:
            lam = (3 * x1 * x1 + self.a) * inv_mod(2 * y1, p) % p
        else:
            lam = (y2 - y1) * inv_mod((x2 - x1) % p, p) % p
        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def mul(self, k: int, P):
        if k < 0:
            return self.mul(-k, self.neg(P))
        R = INF
        Q = P
        while k:
            if k & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

    def random_point(self, rng: random.Random):
        p = self.p
        while True:
            x = rng.randrange(p)
            rhs = (x * x * x + self.a * x + self.b) % p
            y = sqrt_mod(rhs, p)
            if y is not None:
                return (x, y if rng.getrandbits(1) else (-y) % p)

    # quadratic twist by a fixed non-residue d:  y^2 = x^3 + a d^2 x + b d^3
    def quadratic_twist(self):
        p = self.p
        d = 2
        while is_square(d, p):
            d += 1
        return Curve(self.a * d * d % p, self.b * d * d * d % p, p)

    def __repr__(self):
        return f"Curve(a={self.a}, b={self.b}, p=...{self.p % 10**6})"


# --------------------------------------------------------------------------
# univariate polynomials mod p  (coeff lists, index = degree, low to high)
# --------------------------------------------------------------------------

def _poly_trim(c):
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c


def poly_mulmod(a, b, mod, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i + j] = (r[i + j] + ai * bj) % p
    return poly_divmod(r, mod, p)[1]


def poly_divmod(a, b, p):
    a = [x % p for x in a]
    b = _poly_trim([x % p for x in b])
    db = len(b) - 1
    inv_lead = inv_mod(b[-1], p)
    q = [0] * max(0, len(a) - db)
    a = a[:]
    while len(a) - 1 >= db and _poly_trim(a) != [0]:
        da = len(a) - 1
        if da < db:
            break
        coef = (a[-1] * inv_lead) % p
        shift = da - db
        q[shift] = coef
        for i in range(len(b)):
            a[shift + i] = (a[shift + i] - coef * b[i]) % p
        _poly_trim(a)
    return _poly_trim(q), _poly_trim(a)


def poly_gcd(a, b, p):
    a = _poly_trim([x % p for x in a])
    b = _poly_trim([x % p for x in b])
    while b != [0] and not (len(b) == 1 and b[0] == 0):
        a, b = b, poly_divmod(a, b, p)[1]
    # make monic
    if a != [0]:
        inv = inv_mod(a[-1], p)
        a = [(x * inv) % p for x in a]
    return _poly_trim(a)


def poly_powmod(base, e, mod, p):
    """base^e mod (mod) in F_p[x]."""
    result = [1]
    base = poly_divmod(base, mod, p)[1]
    while e:
        if e & 1:
            result = poly_mulmod(result, base, mod, p)
        base = poly_mulmod(base, base, mod, p)
        e >>= 1
    return result


def poly_roots_mod_p(coeffs, p, rng=None):
    """All distinct roots in F_p of the polynomial given by coeffs (low->high)."""
    if rng is None:
        rng = random.Random(0xC0FFEE)
    f = _poly_trim([c % p for c in coeffs])
    if len(f) <= 1:
        return []
    roots = []
    # strip x=0 factor
    if f[0] == 0:
        roots.append(0)
        while len(f) > 1 and f[0] == 0:
            f = f[1:]
        _poly_trim(f)
        if len(f) <= 1:
            return sorted(set(roots))
    # g = gcd(f, x^p - x) = product of (x - r) over roots r != 0
    xp = poly_powmod([0, 1], p, f, p)          # x^p mod f
    xp_minus_x = xp[:]
    while len(xp_minus_x) < 2:
        xp_minus_x.append(0)
    xp_minus_x[1] = (xp_minus_x[1] - 1) % p     # x^p - x
    g = poly_gcd(f, xp_minus_x, p)
    if len(g) <= 1:
        return sorted(set(roots))
    roots.extend(_cz_split(g, p, rng))
    return sorted(set(roots))


def _cz_split(g, p, rng):
    """Cantor-Zassenhaus equal-degree split of a squarefree product of linears."""
    g = _poly_trim(g[:])
    deg = len(g) - 1
    if deg == 0:
        return []
    if deg == 1:                                # x - r  ->  root = -const/lead
        return [(-g[0] * inv_mod(g[1], p)) % p]
    while True:
        # random degree-1 poly  (x + c)
        c = rng.randrange(p)
        h = [c, 1]
        hp = poly_powmod(h, (p - 1) // 2, g, p)
        hp = hp[:]
        if not hp:
            hp = [0]
        hp[0] = (hp[0] - 1) % p                  # h^((p-1)/2) - 1
        d = poly_gcd(g, hp, p)
        dd = len(d) - 1
        if 0 < dd < deg:
            return _cz_split(d, p, rng) + _cz_split(poly_divmod(g, d, p)[0], p, rng)


# --------------------------------------------------------------------------
# 2-isogenies via Velu from a rational 2-torsion point
# --------------------------------------------------------------------------

def two_torsion_x(E: Curve):
    """x-coordinates of rational 2-torsion: roots of x^3 + a x + b over F_p."""
    return poly_roots_mod_p([E.b, E.a, 0, 1], E.p)


def velu_2_isogeny(E: Curve, x0: int) -> Curve:
    """
    Codomain of the 2-isogeny with kernel {O, (x0, 0)} (x0 a root of x^3+ax+b).

    Velu, single 2-torsion kernel point Q=(x0,0):
        v = 3*x0^2 + a ;  u = x0*v
        a' = a - 5*v   ;  b' = b - 7*u
    """
    p = E.p
    v = (3 * x0 * x0 + E.a) % p
    u = (x0 * v) % p
    return Curve((E.a - 5 * v) % p, (E.b - 7 * u) % p, p)


def velu_from_kernel_xs(E: Curve, kernel_xs, order2_xs=()) -> Curve:
    """
    Velu codomain for a rational kernel given by the x-coordinates of one
    representative per {Q, -Q} pair (none of order 2) in ``kernel_xs`` plus any
    order-2 kernel x-coordinates in ``order2_xs``.

    For a non-2-torsion representative Q=(x,y):
        g^x_Q = 3 x^2 + a ,   (2 y)^2 = 4 (x^3 + a x + b)
        v_Q   = 2 g^x_Q ,     u_Q = (2 y)^2
    For a 2-torsion point (x,0):
        v_Q = g^x_Q ,         u_Q = x * g^x_Q
    Codomain:  a' = a - 5 sum(v_Q) ,  b' = b - 7 sum(u_Q).

    Only x-coordinates are needed (u_Q uses y^2 = x^3+a x+b), so the codomain is
    defined over F_p even when the kernel points themselves are not rational.
    """
    p, a, b = E.p, E.a, E.b
    t = w = 0
    for x in kernel_xs:                      # non-2-torsion, one per +-pair
        gx = (3 * x * x + a) % p
        tq = (2 * gx) % p
        uq = (4 * (x * x * x + a * x + b)) % p   # (2y)^2
        t = (t + tq) % p
        w = (w + uq + x * tq) % p            # w_Q = u_Q + x_Q t_Q
    for x in order2_xs:                      # 2-torsion: u_Q = 0
        gx = (3 * x * x + a) % p
        t = (t + gx) % p
        w = (w + x * gx) % p
    return Curve((a - 5 * t) % p, (b - 7 * w) % p, p)


def division_poly_3(E: Curve):
    """psi_3 = 3 x^4 + 6 a x^2 + 12 b x - a^2  (low->high coeffs).

    Its roots in F_p are the x-coordinates of 3-torsion points; each rational
    root is the kernel of a rational 3-isogeny (the subgroup {O, Q, -Q} is
    Galois-stable even if Q itself lives over F_{p^2}).
    """
    p, a, b = E.p, E.a, E.b
    return [(-a * a) % p, (12 * b) % p, (6 * a) % p, 0, 3 % p]


def three_isogenous_neighbors(E: Curve):
    """All rational 3-isogenous codomains of E (via psi_3 roots + Velu)."""
    p = E.p
    roots = poly_roots_mod_p(division_poly_3(E), p)
    return [velu_from_kernel_xs(E, [r]) for r in roots], roots


# --------------------------------------------------------------------------
# order verification by random-point test (cofactor / twist disambiguation)
# --------------------------------------------------------------------------

def order_is(E: Curve, n: int, p: int, trials: int = 8, seed: int = 1) -> bool:
    """
    Probabilistic check that #E(F_p) == n.  Over F_p the two relevant orders
    are n and 2p+2-n (the twist order); a handful of random points with
    [n]P == O essentially pins it down when n != twist order.
    """
    rng = random.Random(seed)
    for _ in range(trials):
        P = E.random_point(rng)
        if E.mul(n, P) is not INF:
            return False
    return True
