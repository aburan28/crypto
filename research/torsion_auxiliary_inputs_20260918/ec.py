"""Toy prime-field elliptic curves with exact group-operation accounting.

Every curve addition and doubling passes through ``Counter`` so that the
experiments in this directory can report *operations*, never wall time,
as the repository's AGENTS.md requires.  Affine short Weierstrass
arithmetic is used throughout; an addition and a doubling each cost one
field inversion and are counted as one group operation apiece, the same
convention the repository's rho harnesses use.

Nothing here is optimised.  The point is that the counts are exact and
that every scalar multiplication is visible in them.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

Point = Optional[Tuple[int, int]]  # None is the point at infinity


# ---------------------------------------------------------------------------
# Number theory helpers
# ---------------------------------------------------------------------------


def is_probable_prime(n: int, rounds: int = 24) -> bool:
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    rng = random.Random(0xC0FFEE ^ n)
    for _ in range(rounds):
        a = rng.randrange(2, n - 1)
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


def factorise(n: int) -> Dict[int, int]:
    """Trial division plus Pollard rho; fine for the sizes used here."""
    out: Dict[int, int] = {}
    for q in (2, 3, 5, 7, 11, 13):
        while n % q == 0:
            out[q] = out.get(q, 0) + 1
            n //= q
    if n == 1:
        return out
    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_probable_prime(m):
            out[m] = out.get(m, 0) + 1
            continue
        f = _pollard_rho_factor(m)
        stack.append(f)
        stack.append(m // f)
    return out


def _pollard_rho_factor(n: int) -> int:
    if n % 2 == 0:
        return 2
    rng = random.Random(n)
    while True:
        c = rng.randrange(1, n)
        x = y = rng.randrange(0, n)
        g = 1
        while g == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            g = math.gcd(abs(x - y), n)
        if g != n:
            return g


def sqrt_mod(a: int, p: int) -> Optional[int]:
    """Tonelli–Shanks.  Returns None when ``a`` is a non-residue."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
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
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def primitive_root(p: int, fac: Optional[Dict[int, int]] = None) -> int:
    fac = fac or factorise(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g
    raise ValueError("no primitive root")


def divisors(fac: Dict[int, int]) -> List[int]:
    ds = [1]
    for q, e in fac.items():
        ds = [d * q**k for d in ds for k in range(e + 1)]
    return sorted(ds)


# ---------------------------------------------------------------------------
# Operation counter
# ---------------------------------------------------------------------------


@dataclass
class Counter:
    """Exact counts.  ``ops`` is what the notes report; the rest is detail."""

    add: int = 0
    dbl: int = 0
    smul: int = 0  # number of scalar multiplications requested
    phases: Dict[str, int] = field(default_factory=dict)
    _phase: Optional[str] = None

    @property
    def ops(self) -> int:
        return self.add + self.dbl

    def begin(self, name: str) -> None:
        self._flush()
        self._phase = name
        self._mark = self.ops

    def _flush(self) -> None:
        if self._phase is not None:
            self.phases[self._phase] = self.phases.get(self._phase, 0) + (
                self.ops - self._mark
            )
            self._phase = None

    def end(self) -> None:
        self._flush()

    def snapshot(self) -> Dict[str, int]:
        self._flush()
        return {"ops": self.ops, "add": self.add, "dbl": self.dbl, "smul": self.smul,
                "phases": dict(self.phases)}


# ---------------------------------------------------------------------------
# Curve arithmetic
# ---------------------------------------------------------------------------


class Curve:
    """y² = x³ + a x + b over F_q, affine, every op counted."""

    def __init__(self, q: int, a: int, b: int, counter: Optional[Counter] = None):
        self.q, self.a, self.b = q, a % q, b % q
        self.ctr = counter or Counter()
        self.bits = q.bit_length()

    # -- basic group law -----------------------------------------------------

    def on_curve(self, P: Point) -> bool:
        if P is None:
            return True
        x, y = P
        return (y * y - (x * x * x + self.a * x + self.b)) % self.q == 0

    def neg(self, P: Point) -> Point:
        return None if P is None else (P[0], (-P[1]) % self.q)

    def add(self, P: Point, Q: Point) -> Point:
        if P is None:
            return Q
        if Q is None:
            return P
        q = self.q
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % q == 0:
                return None
            return self.dbl(P)
        self.ctr.add += 1
        lam = (y2 - y1) * pow(x2 - x1, -1, q) % q
        x3 = (lam * lam - x1 - x2) % q
        y3 = (lam * (x1 - x3) - y1) % q
        return (x3, y3)

    def dbl(self, P: Point) -> Point:
        if P is None:
            return None
        q = self.q
        x1, y1 = P
        if y1 == 0:
            return None
        self.ctr.dbl += 1
        lam = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, q) % q
        x3 = (lam * lam - 2 * x1) % q
        y3 = (lam * (x1 - x3) - y1) % q
        return (x3, y3)

    def mul(self, k: int, P: Point) -> Point:
        """Plain left-to-right double-and-add.  Counted in full."""
        self.ctr.smul += 1
        if k < 0:
            k, P = -k, self.neg(P)
        R: Point = None
        for bit in bin(k)[2:]:
            R = self.dbl(R)
            if bit == "1":
                R = self.add(R, P)
        return R

    # -- fixed-base comb ----------------------------------------------------

    def comb_table(self, P: Point, nbits: int, w: int) -> "FixedBase":
        return FixedBase(self, P, nbits, w)

    # -- point utilities ------------------------------------------------------

    def random_point(self, rng: random.Random) -> Point:
        while True:
            x = rng.randrange(self.q)
            rhs = (x * x * x + self.a * x + self.b) % self.q
            y = sqrt_mod(rhs, self.q)
            if y is not None:
                if rng.random() < 0.5:
                    y = (-y) % self.q
                return (x, y)

    def lift_x(self, x: int) -> Point:
        rhs = (x * x * x + self.a * x + self.b) % self.q
        y = sqrt_mod(rhs, self.q)
        return None if y is None else (x % self.q, y)


class FixedBase:
    """Fixed-base comb: [k]P costs ⌈nbits/w⌉ − 1 additions after a table
    of (2^w − 1)·⌈nbits/w⌉ points is built.  Both are counted through the
    curve's counter, the table under whatever phase is open when it is
    built (call ``begin('precompute')`` first)."""

    def __init__(self, curve: Curve, P: Point, nbits: int, w: int):
        self.E = curve
        self.w = w
        self.chunks = -(-nbits // w)
        self.table: List[List[Point]] = []
        base = P
        for _ in range(self.chunks):
            row: List[Point] = [None]
            for v in range(1, 1 << w):
                row.append(curve.add(row[-1], base))
            self.table.append(row)
            for _ in range(w):
                base = curve.dbl(base)
        self.size = self.chunks * ((1 << w) - 1)

    def mul(self, k: int) -> Point:
        self.E.ctr.smul += 1
        R: Point = None
        mask = (1 << self.w) - 1
        for i in range(self.chunks):
            v = (k >> (self.w * i)) & mask
            if v:
                R = self.E.add(R, self.table[i][v])
        return R


# ---------------------------------------------------------------------------
# Group order and curve search
# ---------------------------------------------------------------------------


def point_order_candidates(E: Curve, P: Point) -> List[int]:
    """All k in the Hasse interval with [k]P = O, by baby-step giant-step.
    The search is not charged to any experiment: it is curve *generation*."""
    q = E.q
    w = math.isqrt(q) + 1
    lo, hi = q + 1 - 2 * w, q + 1 + 2 * w
    m = math.isqrt(hi - lo) + 1
    saved = (E.ctr.add, E.ctr.dbl, E.ctr.smul)
    # baby[±[j]P] = ±j, so that [lo + i m]P = [j]P  ⇒  k = lo + i m − j.
    baby: Dict[Point, int] = {}
    R: Point = None
    for j in range(m):
        baby.setdefault(R, j)
        baby.setdefault(E.neg(R), -j)
        R = E.add(R, P)
    step = E.mul(m, P)
    cur = E.mul(lo, P)
    found = set()
    for i in range(m + 1):
        if cur in baby:
            k = lo + i * m - baby[cur]
            if lo <= k <= hi:
                found.add(k)
        cur = E.add(cur, step)
    E.ctr.add, E.ctr.dbl, E.ctr.smul = saved
    return sorted(found)


def group_order_if_prime(E: Curve, rng: random.Random) -> Optional[int]:
    """Return #E(F_q) when it is prime, else None."""
    P = E.random_point(rng)
    ks = point_order_candidates(E, P)
    if len(ks) != 1:
        return None
    k = ks[0]
    if not is_probable_prime(k):
        return None
    # ord(P) | k prime ⇒ ord(P) = k > 4√q ⇒ #E = k.
    return k


@dataclass
class CurveInstance:
    q: int
    a: int
    b: int
    p: int  # prime group order
    G: Tuple[int, int]
    fac_p_minus_1: Dict[int, int]
    fac_p_plus_1: Dict[int, int]

    def curve(self, counter: Optional[Counter] = None) -> Curve:
        return Curve(self.q, self.a, self.b, counter)

    def to_json(self) -> dict:
        return {
            "q": self.q, "a": self.a, "b": self.b, "p": self.p, "G": list(self.G),
            "p_minus_1": {str(k): v for k, v in self.fac_p_minus_1.items()},
            "p_plus_1": {str(k): v for k, v in self.fac_p_plus_1.items()},
        }


def random_prime(bits: int, rng: random.Random) -> int:
    while True:
        c = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_probable_prime(c):
            return c


def best_divisor(fac: Dict[int, int], target: float) -> int:
    """Divisor closest to ``target`` in log scale."""
    return min(divisors(fac), key=lambda d: abs(math.log(d) - math.log(target)))


def find_curve(bits: int, seed: int, want: str, exponent: float,
               tolerance: float = 0.06, max_tries: int = 200000) -> CurveInstance:
    """Search random curves y² = x³ + ax + b over a ``bits``-bit prime until
    the group order p is prime and p∓1 (``want`` = 'p-1' or 'p+1') has a
    divisor d with |log_p d − exponent| ≤ tolerance."""
    rng = random.Random(seed)
    q = random_prime(bits, rng)
    for _ in range(max_tries):
        a, b = rng.randrange(q), rng.randrange(q)
        if (4 * a**3 + 27 * b**2) % q == 0:
            continue
        E = Curve(q, a, b)
        p = group_order_if_prime(E, rng)
        if p is None:
            continue
        n = p - 1 if want == "p-1" else p + 1
        fac = factorise(n)
        d = best_divisor(fac, p**exponent)
        if abs(math.log(d, p) - exponent) <= tolerance:
            G = E.random_point(rng)
            return CurveInstance(q, a, b, p, G, factorise(p - 1), factorise(p + 1))
    raise RuntimeError("no curve found")
