"""Binary subfield curves and exact linear Artin--Schreier solving.

Field elements use the existing type-II ONB representation. The curve is
y^2+xy=x^3+A*x^2+B, B!=0; coefficients need not belong to F2.
"""
import curves


class LinearImage:
    """An F2-linear image with preimages; every vector update is charged."""
    def __init__(self, f, columns):
        self.f = f
        pivots = {}
        for value, pre in columns:
            for k, (basis, lift) in sorted(pivots.items(), reverse=True):
                if value >> k & 1:
                    value = f.add(value, basis)
                    pre = f.add(pre, lift)
            if value:
                pivots[value.bit_length() - 1] = (value, pre)
        self.rows = sorted(pivots.items(), reverse=True)

    def preimage(self, value):
        pre = 0
        for k, (basis, lift) in self.rows:
            if value >> k & 1:
                value = self.f.add(value, basis)
                pre = self.f.add(pre, lift)
        return None if value else pre


class Curve(curves.Curve):
    def __init__(self, f, a2, a6):
        super().__init__(f)
        if not a6:
            raise ValueError('singular curve: a6 must be nonzero')
        self.a2, self.a6 = a2, a6
        cols = []
        for i in range(f.m):
            pre = f.fromCoords(1 << i)
            cols.append((f.add(f.sqr(pre), pre), pre))
        self.artinSchreier = LinearImage(f, cols)
        if len(self.artinSchreier.rows) != f.m - 1:
            raise ArithmeticError('Artin--Schreier rank mismatch')

    def asRoots(self, value):
        w = self.artinSchreier.preimage(value)
        return [] if w is None else [w, self.f.add(w, self.one)]

    def onCurve(self, p):
        if p is None:
            return True
        f = self.f
        x, y = p
        x2 = f.sqr(x)
        return f.add(f.sqr(y), f.mul(x, y)) == f.add(
            f.add(f.mul(x2, x), f.mul(self.a2, x2)), self.a6)

    def dbl(self, p):
        if p is None or not p[0]:
            return None
        f = self.f
        x, y = p
        lam = f.add(x, f.mul(y, f.inv(x)))
        x3 = f.add(f.add(f.sqr(lam), lam), self.a2)
        return x3, f.add(f.sqr(x), f.mul(f.add(lam, self.one), x3))

    def add(self, p, q):
        if p is None:
            return q
        if q is None:
            return p
        if p[0] == q[0]:
            return self.dbl(p) if p == q else None
        f = self.f
        x1, y1 = p
        x2, y2 = q
        dx = f.add(x1, x2)
        lam = f.mul(f.add(y1, y2), f.inv(dx))
        x3 = f.add(f.add(f.add(f.sqr(lam), lam), dx), self.a2)
        return x3, f.add(f.add(f.mul(lam, f.add(x1, x3)), x3), y1)

    def pointFromX(self, x):
        f = self.f
        if not x:
            return 0, f.frob(self.a6, f.m - 1)
        rhs = f.add(f.add(x, self.a2), f.mul(self.a6, f.inv(f.sqr(x))))
        roots = self.asRoots(rhs)
        return None if not roots else (x, f.mul(x, roots[0]))

    def halfTrace(self, a):
        # Prevent accidental reuse of the old odd-degree-only implementation.
        raise ValueError('use asRoots: this curve supports even field degrees')


def subfieldGenerator(f, k):
    """Deterministic non-F2 relative trace of an ONB vector (k=2 or 3)."""
    if k not in (2, 3) or f.m % k or f.m <= k:
        raise ValueError('a proper F4 or F8 coefficient subfield is required')
    candidates = []
    for i in range(f.m):
        value = 0
        for j in range(f.m // k):
            value = f.add(value, f.frob(f.fromCoords(1 << i), j * k))
        if value not in (0, f.one()):
            if f.frob(value, k) != value:
                raise ArithmeticError('coefficient is not in the declared subfield')
            candidates.append(f.toCoords(value))
    if not candidates:
        raise ArithmeticError('relative trace did not span the subfield')
    return f.fromCoords(min(candidates))


def residualNorm(f, curve, target, a, b):
    if target is None or not b:
        raise ValueError('finite target and nonzero function y coefficient required')
    r, s = target
    c = f.add(f.add(f.sqr(r), f.mul(a, r)), f.mul(b, f.add(r, s)))
    b2 = f.sqr(b)
    norm = [f.add(f.sqr(c), f.mul(b2, curve.a6)), f.mul(b, c),
            f.add(f.add(f.sqr(a), f.mul(a, b)), f.mul(b2, curve.a2)),
            f.add(b, b2), f.one()]
    h = [0, 0, 0, f.one()]
    for i in (2, 1, 0):
        h[i] = f.add(norm[i + 1], f.mul(r, h[i + 1]))
    if f.add(norm[0], f.mul(r, h[0])):
        raise ArithmeticError('norm does not vanish at the pinned target')
    return h, c
