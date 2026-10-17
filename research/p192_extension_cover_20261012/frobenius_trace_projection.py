"""Exact P-192 E(Fp6) to E(Fp)[r] cofactor projection through trace.

For #E(Fp6) = r*h, #E(Fp) = r, gcd(r,h) = 1, and gcd(r,6) = 1,
the r-primary component of every extension point is Fp-rational.  If
Tr(P) = sum_{i=0}^5 Frob_p^i(P), then

    [h]P = [h * 6^{-1} mod r] Tr(P).

The caller is responsible for checking the declared group orders.  The
method computes the trace in E(Fp6), coerces its coordinates to E(Fp), and
performs the remaining scalar multiplication in the prime field.
"""

from sage.all import EllipticCurve, Integer


class FrobeniusTraceProjector:
    def __init__(self, base_field, extension_field, a, b, r, h, extension_degree=6):
        self.p = Integer(base_field.cardinality())
        self.r = Integer(r)
        self.h = Integer(h)
        self.degree = Integer(extension_degree)
        assert self.degree == 6
        assert self.r.is_prime() and self.r % self.degree != 0
        assert self.h % self.r != 0
        assert self.r * self.h > 0
        self.E0 = EllipticCurve(base_field, [base_field(a), base_field(b)])
        self.E6 = EllipticCurve(extension_field, [extension_field(a), extension_field(b)])
        assert self.E0.cardinality() == self.r
        self.scalar = Integer(self.h * self.degree.inverse_mod(self.r) % self.r)

    def trace(self, point):
        """Six-step additive Frobenius trace, returned on E(Fp)."""
        P = self.E6(point)
        total = self.E6(0)
        x = P[0] if P else None
        y = P[1] if P else None
        for _ in range(6):
            total += self.E6(x, y) if x is not None else self.E6(0)
            if x is not None:
                x = x ** self.p
                y = y ** self.p
        if not total:
            return self.E0(0)
        assert total[0] ** self.p == total[0]
        assert total[1] ** self.p == total[1]
        return self.E0(self.E0.base_field()(total[0]),
                       self.E0.base_field()(total[1]))

    def project_delta(self, point, anchor_trace):
        T = self.trace(point) - anchor_trace
        Q = self.scalar * T
        assert self.r * Q == self.E0(0)
        return Q

    def direct_delta(self, point, anchor):
        Q = self.h * (self.E6(point) - self.E6(anchor))
        assert self.r * Q == self.E6(0)
        return Q

    @staticmethod
    def key(Q):
        if not Q:
            return None
        return int(Q[0]), int(Q[1])
