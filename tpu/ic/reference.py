"""Scalar oracle for the TPU index-calculus backend.

This module is deliberately plain Python with no JAX, NumPy, or hardware
dependency.  It is the *oracle*: the batched, MXU-shaped kernels in
``field.py`` / ``curve.py`` / ``pairtable.py`` / ``linalg.py`` are correct
iff they agree with this file element by element, exactly as the CUDA
backend in ``gpu/ecc2k`` is checked against a scalar reference before a
device ever runs it.

Every contract here is transcribed from the Rust library so the port
cannot drift from it, with the source cited at each definition:

* field ``F_2[t]/(irr)`` -- ``src/cryptanalysis/semaev_decomp.rs`` ``Gf2``
  (``mul`` = ``reduce(clmul(a,b))``, ``sqr`` = ``reduce(spread(a))``,
  ``inv`` = Fermat ``a^(2^n-2)``).
* binary Koblitz curve ``y^2 + xy = x^3 + a x^2 + b`` --
  ``src/cryptanalysis/koblitz_fast.rs`` ``FastCurve`` (``neg``/``double``/
  ``add``/``frobenius``).
* ``pack`` -- ``FastPoint::pack`` / ``pack_point`` in
  ``src/cryptanalysis/koblitz_index_calculus.rs``.
* ``pair_filter_hash`` -- same file.

None of this is a *performance* path; it exists to be obviously right.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

U64 = (1 << 64) - 1

# --------------------------------------------------------------------------
# field F_2[t] / (irr)
# --------------------------------------------------------------------------


def clmul(a: int, b: int) -> int:
    """Carry-less (GF(2)) product of two bit-polynomials, as ``Gf2::clmul``'s
    software path: XOR of ``a`` shifted by each set bit of ``b``."""
    w = 0
    aa = a
    while b:
        if b & 1:
            w ^= aa
        aa <<= 1
        b >>= 1
    return w


def gf2_reduce(w: int, irr: int, n: int) -> int:
    """Fold a product ``< 2^{2n-1}`` back mod ``irr`` (``Gf2::reduce``)."""
    for i in range(2 * n - 2, n - 1, -1):
        if (w >> i) & 1:
            w ^= irr << (i - n)
    return w & ((1 << n) - 1)


@dataclass(frozen=True)
class Field:
    """A binary field ``F_2[t]/(irr)`` of degree ``n`` (``n <= 62``)."""

    n: int
    irr: int  # includes the leading t^n bit

    @staticmethod
    def from_low_terms(n: int, low_terms: List[int]) -> "Field":
        irr = 1 << n
        for t in low_terms:
            irr |= 1 << t
        irr |= 1  # the constant term is always present for an irreducible
        return Field(n=n, irr=irr)

    @property
    def mask(self) -> int:
        return (1 << self.n) - 1

    def add(self, a: int, b: int) -> int:  # noqa: D401 - trivial
        return a ^ b

    def mul(self, a: int, b: int) -> int:
        return gf2_reduce(clmul(a, b), self.irr, self.n)

    def sqr(self, a: int) -> int:
        # bit spread: coefficient i -> 2i, then reduce.
        w = 0
        i = 0
        while a:
            if a & 1:
                w |= 1 << (2 * i)
            a >>= 1
            i += 1
        return gf2_reduce(w, self.irr, self.n)

    def sqr_k(self, a: int, k: int) -> int:
        for _ in range(k):
            a = self.sqr(a)
        return a

    def pow(self, a: int, e: int) -> int:
        r = 1
        base = a
        while e:
            if e & 1:
                r = self.mul(r, base)
            base = self.sqr(base)
            e >>= 1
        return r

    def inv(self, a: int) -> int:
        if a == 0:
            return 0
        return self.pow(a, (1 << self.n) - 2)

    # --- GF(2)-linear maps as bit matrices (the MXU reformulation) --------
    def squaring_matrix(self) -> List[int]:
        """Column ``j`` is ``sqr(t^j)`` as an n-bit integer.  Squaring is
        ``F_2``-linear, so ``sqr(a) = M @ a (mod 2)`` where ``M`` is this
        matrix.  This is what makes squaring (and hence the Fermat inverse's
        square-chain) a bit-matmul on the TPU's systolic array."""
        return [self.sqr(1 << j) for j in range(self.n)]

    def mul_by_matrix(self, c: int) -> List[int]:
        """Multiplication by the fixed element ``c`` is ``F_2``-linear.
        Column ``j`` is ``mul(c, t^j)``."""
        return [self.mul(c, 1 << j) for j in range(self.n)]

    def reduction_matrix(self) -> List[int]:
        """Maps a full carry-less product (``2n-1`` bits) to its reduced
        ``n``-bit residue.  Column ``k`` is ``t^k mod irr``, ``0<=k<2n-1``."""
        out = []
        for k in range(2 * self.n - 1):
            out.append(gf2_reduce(1 << k, self.irr, self.n))
        return out


# --------------------------------------------------------------------------
# binary Koblitz curve  y^2 + xy = x^3 + a x^2 + b
# --------------------------------------------------------------------------

# A point is (x, y, infinity).  Mirrors koblitz_fast::FastPoint.
INF = (0, 0, True)


@dataclass(frozen=True)
class Curve:
    field: Field
    a: int
    b: int

    def is_on_curve(self, p) -> bool:
        x, y, inf = p
        if inf:
            return True
        f = self.field
        lhs = f.sqr(y) ^ f.mul(x, y)
        x2 = f.sqr(x)
        rhs = f.mul(x2, x) ^ f.mul(self.a, x2) ^ self.b
        return lhs == rhs

    def neg(self, p):
        x, y, inf = p
        if inf:
            return p
        return (x, x ^ y, False)

    def double(self, p):
        x, y, inf = p
        if inf or x == 0:
            return INF
        f = self.field
        lam = x ^ f.mul(y, f.inv(x))
        x3 = f.sqr(lam) ^ lam ^ self.a
        y3 = f.sqr(x) ^ f.mul(lam ^ 1, x3)
        return (x3, y3, False)

    def add(self, p, q):
        px, py, pinf = p
        qx, qy, qinf = q
        if pinf:
            return q
        if qinf:
            return p
        if px == qx:
            if (py ^ qy) == px:
                return INF
            return self.double(p)
        f = self.field
        lam = f.mul(py ^ qy, f.inv(px ^ qx))
        x3 = f.sqr(lam) ^ lam ^ px ^ qx ^ self.a
        y3 = f.mul(lam, px ^ x3) ^ x3 ^ py
        return (x3, y3, False)

    def frobenius(self, p):
        x, y, inf = p
        if inf:
            return p
        f = self.field
        return (f.sqr(x), f.sqr(y), False)

    def mul_scalar(self, p, k: int):
        r = INF
        base = p
        while k:
            if k & 1:
                r = self.add(r, base)
            base = self.double(base)
            k >>= 1
        return r


def pack(p) -> int:
    """``FastPoint::pack``: ``0`` for O, else ``2(x+1) + s`` with sign bit
    ``s = (y > (x ^ y))``.  Valid for ``n <= 62``."""
    x, y, inf = p
    if inf:
        return 0
    sign = 1 if y > (x ^ y) else 0
    return ((x + 1) << 1) | sign


def pair_filter_hash(key: int) -> int:
    """``pair_filter_hash`` in ``koblitz_index_calculus.rs``: a Murmur3-style
    finalizer with the two 64-bit constants."""
    h = (key * 0xFF51AFD7ED558CCD) & U64
    h ^= h >> 33
    return (h * 0xC4CEB9FE1A85EC53) & U64
