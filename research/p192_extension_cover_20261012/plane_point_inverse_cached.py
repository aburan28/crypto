"""Cached Type I plane-point inverse for repeated frozen support projections.

This keeps the same equations as ``PlanePointInverse`` while computing
Frobenius conjugates of model constants once per batch. It is a candidate
implementation until independently paired against the established inverse.
"""

from sage.all import PolynomialRing, prod, vector


class CachedPlanePointInverse:
    def __init__(self, model):
        self.model = model
        self.p = model.p
        self.L = model.L
        self.R = PolynomialRing(self.L, "t")
        t = self.R.gen()
        self.conjugates = [model.zeta ** (self.p ** i) for i in range(3)]
        self.quartic_roots = [
            [root ** (self.p ** i) for root in model.quartic_roots]
            for i in range(3)
        ]
        self.norm_roots = [
            root ** (self.p ** i)
            for root in (model.zeta, model.beta) for i in range(3)
        ]
        self.quartics = [prod(t - root for root in roots)
                         for roots in self.quartic_roots]
        self.quotient_constants = [
            tuple(value ** (self.p ** i)
                  for value in (model.inf, model.c0, model.d0, model.kappa))
            for i in range(3)
        ]

    def invert_affine(self, u, v):
        """Return checked cover coordinates and their three quartic values."""
        L = self.L
        A = [L(u) + self.conjugates[i] * L(v)
             + self.conjugates[i] ** 2 for i in range(3)]
        assert all(a != 0 for a in A)
        f = self.quartics
        common = (A[0] ** 2 * f[0] - A[1] ** 2 * f[1]).gcd(
            A[0] ** 2 * f[0] - A[2] ** 2 * f[2]).monic()
        assert common.degree() == 1
        X = -common[0]
        norm = prod(X - root for root in self.norm_roots)
        assert norm != 0
        lam = A[0] * f[0](X) / (A[1] * A[2] * norm)
        assert lam != 0
        y = [1 / (lam * a) for a in A]
        assert all(y[i] ** 2 == f[i](X) for i in range(3))
        assert prod(y) == norm
        Y = self.model.Vandermonde.solve_right(vector(L, y))
        point = (X, *Y)
        assert all(Y[0] + self.conjugates[i] * Y[1]
                   + self.conjugates[i] ** 2 * Y[2] == y[i]
                   for i in range(3))
        assert all(c ** self.p == c for c in point)
        roundtrip = self.model.Vandermonde.solve_right(
            vector(L, [1 / value for value in y]))
        assert roundtrip[0] / roundtrip[2] == L(u)
        assert roundtrip[1] / roundtrip[2] == L(v)
        return point, y

    def quotient_image(self, E, point, y, index):
        """Map a checked cover point through a cached quotient component."""
        assert index in (0, 1, 2)
        X = point[0]
        rho, c0, d0, kappa = self.quotient_constants[index]
        assert X != rho
        return E(c0 + d0 / (X - rho), kappa * y[index] / (X - rho) ** 2)
