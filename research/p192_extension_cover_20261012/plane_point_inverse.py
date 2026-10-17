"""Invert generic affine points of the descended Type I plane quartic.

For the three conjugate quartic equations y_i^2=f_i(X), the canonical plane
coordinate has A_i=1/(lambda*y_i). Two squared ratios recover X as the common
root of A_0^2*f_0-A_i^2*f_i. The norm equation fixes lambda and removes the
simultaneous-sign ambiguity. Exceptional points are reported by assertions;
the frozen line-support run checks that its inputs avoid them.
"""

from sage.all import prod, vector, PolynomialRing


class PlanePointInverse:
    def __init__(self, model):
        self.model = model
        self.p = model.p
        self.L = model.L
        self.R = PolynomialRing(self.L, "t")
        t = self.R.gen()
        self.quartics = [prod(t - root ** (self.p ** i)
                              for root in model.quartic_roots)
                         for i in range(3)]
        self.conjugates = [model.zeta ** (self.p ** i) for i in range(3)]

    def invert_affine(self, u, v):
        """Return a checked cover point above (u:v:1) over the prime field."""
        dt, p, L = self.model, self.p, self.L
        A = [L(u) + self.conjugates[i] * L(v) +
             self.conjugates[i] ** 2 for i in range(3)]
        assert all(a != 0 for a in A)
        f = self.quartics
        common = (A[0] ** 2 * f[0] - A[1] ** 2 * f[1]).gcd(
            A[0] ** 2 * f[0] - A[2] ** 2 * f[2]).monic()
        assert common.degree() == 1
        X = -common[0]
        norm = prod(X - root ** (p ** i)
                    for root in (dt.zeta, dt.beta) for i in range(3))
        assert norm != 0
        lam = A[0] * f[0](X) / (A[1] * A[2] * norm)
        assert lam != 0
        y = [1 / (lam * a) for a in A]
        assert all(y[i] ** 2 == f[i](X) for i in range(3))
        assert prod(y) == norm
        Y = dt.Vandermonde.solve_right(vector(L, y))
        point = (X, *Y)
        dt.check_model(point)
        assert all(c ** p == c for c in point)
        roundtrip = dt.Vandermonde.solve_right(vector(L, [1 / a for a in y]))
        assert roundtrip[0] / roundtrip[2] == L(u)
        assert roundtrip[1] / roundtrip[2] == L(v)
        return point


def elliptic_quotient_image(model, E, point, index):
    """Map one checked cover point through its indexed elliptic quotient."""
    assert index in (0, 1, 2)
    X = point[0]
    y = model.check_model(point)[index]
    power = model.p ** index
    rho, c0, d0, kappa = (value ** power for value in
                           (model.inf, model.c0, model.d0, model.kappa))
    assert X != rho
    return E(c0 + d0 / (X - rho), kappa * y / (X - rho) ** 2)
