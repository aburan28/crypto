"""Deterministic toy Koblitz instances for the self-checks.

Small enough to enumerate the whole group, so a test can check a recovered
discrete log against a planted secret the way ``gpu/ecc2k`` does on its k23
and k41 toy curves.  The fields are the binary Koblitz family
``E_0: y^2 + xy = x^3 + a x^2 + 1`` used throughout this repository
(``AGENTS.md`` 8b); the toy *degrees* here are exploratory sizes only and
carry none of the structural fidelity of m=83/131.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from .reference import INF, Curve, Field, clmul, gf2_reduce

# Irreducible low-term sets (exponents < n, constant term implied) giving a
# sparse modulus t^n + ... + 1.  Verified irreducible by ``is_irreducible``.
_LOW_TERMS: Dict[int, List[int]] = {
    7: [1],
    9: [1],
    11: [2],
    13: [4, 3, 1],
    15: [1],
    17: [3],
    19: [5, 2, 1],
    23: [5],
}


def _poly_mulmod(a: int, b: int, irr: int, n: int) -> int:
    return gf2_reduce(clmul(a, b), irr, n)


def _poly_gcd(a: int, b: int) -> int:
    while b:
        # polynomial remainder a mod b
        db = b.bit_length() - 1
        while a.bit_length() - 1 >= db and a:
            a ^= b << ((a.bit_length() - 1) - db)
        a, b = b, a
    return a


def is_irreducible(irr: int, n: int) -> bool:
    """Irreducible over GF(2): x^(2^n) == x (mod irr) and, for each prime p
    dividing n, gcd(x^(2^(n/p)) - x, irr) == 1."""
    if n == 1:
        return True
    # x^(2^k) mod irr by repeated squaring of the x-power.
    def x_pow_2k(k: int) -> int:
        r = 2  # x
        for _ in range(k):
            r = _poly_mulmod(r, r, irr, n)
        return r

    if x_pow_2k(n) != 2:
        return False
    primes = set()
    m = n
    d = 2
    while d * d <= m:
        while m % d == 0:
            primes.add(d)
            m //= d
        d += 1
    if m > 1:
        primes.add(m)
    for p in primes:
        if _poly_gcd(x_pow_2k(n // p) ^ 2, irr) != 1:
            return False
    return True


def find_field(n: int) -> Field:
    if n in _LOW_TERMS:
        f = Field.from_low_terms(n, _LOW_TERMS[n])
        assert is_irreducible(f.irr, n), f"stored modulus for n={n} is reducible"
        return f
    # search trinomials then pentanomials
    for k in range(1, n):
        irr = (1 << n) | (1 << k) | 1
        if is_irreducible(irr, n):
            return Field(n=n, irr=irr)
    for k1 in range(1, n):
        for k2 in range(1, k1):
            for k3 in range(1, k2):
                irr = (1 << n) | (1 << k1) | (1 << k2) | (1 << k3) | 1
                if is_irreducible(irr, n):
                    return Field(n=n, irr=irr)
    raise ValueError(f"no sparse irreducible found for n={n}")


def make_curve(n: int, a: int = 0) -> Curve:
    assert a in (0, 1)
    return Curve(field=find_field(n), a=a, b=1)


def points_with_x(curve: Curve, x: int) -> List[Tuple[int, int, bool]]:
    """Affine points with abscissa ``x`` (0, 1 or 2 of them)."""
    f = curve.field
    pts = []
    rhs = f.mul(f.sqr(x), x) ^ f.mul(curve.a, f.sqr(x)) ^ curve.b  # x^3+ax^2+b
    for y in range(1 << n_of(f)):
        if f.sqr(y) ^ f.mul(x, y) == rhs:
            pts.append((x, y, False))
    return pts


def n_of(f: Field) -> int:
    return f.n


def enumerate_group(curve: Curve) -> List[Tuple[int, int, bool]]:
    """All points (affine + infinity).  Only call for small degrees."""
    pts = [INF]
    for x in range(1 << curve.field.n):
        pts.extend(points_with_x(curve, x))
    return pts


def group_order(curve: Curve) -> int:
    return len(enumerate_group(curve))


def point_order(curve: Curve, p, group_n: int) -> int:
    """Order of ``p``; naive divisor scan of the group order."""
    # factor group_n
    facs = []
    m = group_n
    d = 2
    while d * d <= m:
        while m % d == 0:
            facs.append(d)
            m //= d
        d += 1
    if m > 1:
        facs.append(m)
    order = group_n
    for prime in set(facs):
        while order % prime == 0 and curve.mul_scalar(p, order // prime) == INF:
            order //= prime
    return order


def find_generator(curve: Curve) -> Tuple[Tuple[int, int, bool], int]:
    """A point of maximal order and that order.

    Binary Koblitz curves carry a small cofactor, so the full group is
    generally not cyclic and no point reaches the whole group order; the
    largest-order point generates the biggest cyclic subgroup, which is
    all a probe sequence ``R = [a]G`` needs.  Returns ``(point, order)``
    with ``order`` the order of that point.
    """
    group_n = group_order(curve)
    best = None
    best_order = 0
    for x in range(1, 1 << curve.field.n):
        for p in points_with_x(curve, x):
            o = point_order(curve, p, group_n)
            if o > best_order:
                best_order = o
                best = p
                if best_order == group_n:
                    return best, best_order
    if best is None:
        raise ValueError("curve has no affine points")
    return best, best_order


def factor_base(curve: Curve, size: int, seed: int = 1) -> List[Tuple[int, int, bool]]:
    """A deterministic set of ``size`` distinct affine base points, taken by
    walking abscissae from a seed-dependent start.  Closed under nothing in
    particular -- it is a plain factor base, which is all the pair-table
    stage needs."""
    out: List[Tuple[int, int, bool]] = []
    seen = set()
    x = (seed * 0x9E3779B1) % (1 << curve.field.n)
    guard = 0
    while len(out) < size and guard < (1 << (curve.field.n + 1)):
        for p in points_with_x(curve, x):
            key = (p[0], p[1])
            if key not in seen:
                seen.add(key)
                out.append(p)
                if len(out) >= size:
                    break
        x = (x + 1) % (1 << curve.field.n)
        guard += 1
    return out
