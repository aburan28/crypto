"""Full-point odd-degree binary Vélu evaluation from an irreducible cubic kernel.

This independent research implementation works in A = Fq[U]/h(U), where
h(U) = U^3 + c2 U^2 + c1 U + c0 is the kernel-abscissa polynomial for a
degree-7 map. Summing the three conjugate kernel terms is the trace A/Fq.
It uses only the existing Fq field interface (mul, sqr, inv) and the
existing Koblitz group model; Sage is deliberately not imported here.
"""
from __future__ import annotations

from pathlib import Path
import sys

RESEARCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RESEARCH / "ecc2k130_relations"))
from relations import Koblitz  # noqa: E402


class CubicVeluMap:
    """Normalized map for y²+xy=x³+a*x²+b with a degree-7 cubic kernel."""

    def __init__(self, source: Koblitz, coeffs: tuple[int, int, int]):
        if len(coeffs) != 3 or any(type(z) is not int for z in coeffs):
            raise ValueError("expected three canonical field coefficients")
        self.source = source
        self.F = source.F
        self.c0, self.c1, self.c2 = coeffs
        if any(z < 0 or z >= (1 << self.F.deg) for z in coeffs):
            raise ValueError("noncanonical cubic coefficient")
        if not self.c0:
            raise ValueError("kernel polynomial has zero root")
        self.c2_sq = self.F.sqr(self.c2)
        self.red4 = (
            self.F.mul(self.c2, self.c0),
            self.F.mul(self.c2, self.c1) ^ self.c0,
            self.c2_sq ^ self.c1,
        )
        b_prime = source.b ^ self.c2 ^ self.c2_sq
        if not b_prime:
            raise ValueError("singular quotient")
        self.codomain = Koblitz(self.F, a=source.a, b=b_prime)

    @staticmethod
    def add(a: tuple[int, int, int], b: tuple[int, int, int]) -> tuple[int, int, int]:
        return (a[0] ^ b[0], a[1] ^ b[1], a[2] ^ b[2])

    def mul(self, a: tuple[int, int, int], b: tuple[int, int, int]) -> tuple[int, int, int]:
        f = self.F
        d0 = f.mul(a[0], b[0])
        d1 = f.mul(a[0], b[1]) ^ f.mul(a[1], b[0])
        d2 = f.mul(a[0], b[2]) ^ f.mul(a[1], b[1]) ^ f.mul(a[2], b[0])
        d3 = f.mul(a[1], b[2]) ^ f.mul(a[2], b[1])
        d4 = f.mul(a[2], b[2])
        return (
            d0 ^ f.mul(d3, self.c0) ^ f.mul(d4, self.red4[0]),
            d1 ^ f.mul(d3, self.c1) ^ f.mul(d4, self.red4[1]),
            d2 ^ f.mul(d3, self.c2) ^ f.mul(d4, self.red4[2]),
        )

    def trace(self, a: tuple[int, int, int]) -> int:
        # Tr(1)=1, Tr(U)=c2, Tr(U²)=c2² in characteristic two.
        return a[0] ^ self.F.mul(a[1], self.c2) ^ self.F.mul(a[2], self.c2_sq)

    def h(self, x: int) -> int:
        f = self.F
        x2 = f.sqr(x)
        return f.mul(x2, x) ^ f.mul(self.c2, x2) ^ f.mul(self.c1, x) ^ self.c0

    def __call__(self, point: tuple[int, int] | None) -> tuple[int, int] | None:
        if point is None:
            return None
        if (not isinstance(point, tuple) or len(point) != 2
                or any(type(z) is not int or z < 0 or z >= (1 << self.F.deg) for z in point)
                or not self.source.on_curve(point)):
            raise ValueError("input must be a canonical source-curve point")
        x, y = point
        f = self.F
        x2 = f.sqr(x)
        hx = self.h(x)
        if not hx:
            return None
        inv_hx = f.inv(hx)
        # (U+x)^-1 = [U²+(x+c2)U+(x²+c2*x+c1)]/h(x) modulo h(U).
        inv = (
            f.mul(x2 ^ f.mul(self.c2, x) ^ self.c1, inv_hx),
            f.mul(x ^ self.c2, inv_hx),
            inv_hx,
        )
        u = (0, 1, 0)
        c = self.mul(u, inv)
        c2 = self.mul(c, c)
        inv2 = self.mul(inv, inv)
        c3 = self.mul(c, c2)
        x_image = x ^ self.trace(c) ^ self.trace(c2)
        y_image = (y ^ f.mul(x2 ^ y, self.trace(self.mul(u, inv2)))
                   ^ self.trace(c3) ^ self.trace(c))
        image = (x_image, y_image)
        if not self.codomain.on_curve(image):
            raise ArithmeticError("cubic Vélu image failed codomain equation")
        return image
