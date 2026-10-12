"""Early image rejection and a charged S3 pair-invariant table."""
import time
import algebra
import hybrid
import nagaoannihilator as scalar
from space import Space


class Field(scalar.CountedField):
    def __init__(self, n):
        super().__init__(n)
        self.binary = {}

    def word(self, key, count=1):
        row = self.binary.setdefault(self.phase, {})
        row[key] = row.get(key, 0) + count

    def parity(self, value, mask):
        self.word('word_and_popcount_parity')
        return (value & mask).bit_count() & 1

    def fullReport(self):
        return {'field_api': self.report(), 'binary_word_operations': self.binary}


class Curve(algebra.Curve):
    def add(self, p, q):
        self.f.bump('curveAdditionCalls')
        return super().add(p, q)

    def dbl(self, p):
        self.f.bump('curveDoublingCalls')
        return super().dbl(p)


class SupportImage(algebra.LinearImage):
    def __init__(self, f, columns, filtered):
        super().__init__(f, columns)
        self.checks = []
        if not filtered:
            return
        masks = {}
        # Compute the linear remainder of each ambient basis vector; each
        # independent residual bit supplies a parity check on the input.
        for j in range(f.m):
            value = f.fromCoords(1 << j)
            for k, (column, _) in self.rows:
                if value >> k & 1:
                    value = f.add(value, column)
            for i in range(1, f.m + 1):
                if value >> i & 1:
                    masks[i] = masks.get(i, 0) ^ (1 << (j + 1))
                    f.word('mask_construction_xor')
        pivots = {}
        for mask in masks.values():
            value = mask
            for k, column in sorted(pivots.items(), reverse=True):
                if value >> k & 1:
                    value ^= column
                    f.word('mask_construction_xor')
            if value:
                pivots[value.bit_length() - 1] = value
                self.checks.append(value)
        if len(self.checks) != f.m - len(self.rows):
            raise ArithmeticError('support parity-check rank mismatch')

    def contains(self, value):
        return all(not self.f.parity(value, mask) for mask in self.checks)


class ImageCache:
    def __init__(self, f, space, filtered, deadline):
        self.images = {}
        for u in space.values[1:]:
            if time.perf_counter() >= deadline:
                raise TimeoutError('support setup deadline')
            columns = [(f.add(f.sqr(w), f.mul(u, w)), w) for w in space.basis]
            image = SupportImage(f, columns, filtered)
            if len(image.rows) != space.d - 1:
                raise ArithmeticError('support image rank mismatch')
            self.images[u] = image


class Search:
    def __init__(self, f, curve, space, target, deadline, filtered=False, cache=None):
        self.f, self.curve, self.space, self.target = f, curve, space, target
        self.deadline, self.filtered = deadline, filtered
        self.cache = cache or ImageCache(f, space, filtered, deadline)
        r, _ = target
        zs = [z for z in space.values[1:] if z != r]
        invs = hybrid.batchInverse(f, zs + [f.add(r, z) for z in zs])
        self.geometry = [(z, f.mul(z, invs[len(zs) + i]), f.mul(f.add(r, z), f.sqr(invs[i])))
                         for i, z in enumerate(zs)]
        self.current = None
        self.generated = self.earlyRejected = 0

    def check(self):
        if time.perf_counter() >= self.deadline:
            raise TimeoutError('search deadline')

    def candidates(self):
        f = self.f
        r, _ = self.target
        for h2 in self.space.values:
            self.check()
            for b in self.curve.asRoots(f.add(h2, r)):
                if not b:
                    continue
                h, _ = algebra.residualNorm(f, self.curve, self.target, 0, b)
                self.current = b, h2, h[1]
                invB = f.inv(b)
                invB2 = f.sqr(invB)
                for z, alpha, beta in self.geometry:
                    self.check()
                    if z == h2:
                        continue
                    k = scalar.polyEval(f, h, z)
                    t = f.mul(b, alpha)
                    for w in self.curve.asRoots(f.mul(f.mul(k, beta), invB2)):
                        self.generated += 1
                        yield f.mul(t, w), b, invB, z

    def recover(self, a, b, invB, z):
        f = self.f
        r, _ = self.target
        if self.filtered:
            currentB, h2, h10 = self.current
            if b != currentB:
                raise ArithmeticError('stale coefficient cache')
            u = f.add(h2, z)
            h1 = f.add(f.add(f.sqr(a), f.mul(a, b)), h10)
            v = f.add(h1, f.mul(z, u))
            if not self.cache.images[u].contains(v):
                self.earlyRejected += 1
                return None
        h, c = algebra.residualNorm(f, self.curve, self.target, a, b)
        if scalar.polyEval(f, h, z):
            raise ArithmeticError('conditioned root missing')
        if not h[0] or not scalar.polyEval(f, h, r):
            return None
        u = f.add(h[2], z)
        v = f.add(h[1], f.mul(z, u))
        w = self.cache.images[u].preimage(v)
        if w is None:
            return None
        xs = [z, w, f.add(w, u)]
        if len(set(xs)) != 3 or any(x == 0 or x == r for x in xs):
            return None
        if any(x not in self.space.indices for x in xs):
            raise ArithmeticError('support violation')
        points = [(x, f.mul(f.add(f.add(f.sqr(x), f.mul(a, x)), c), invB)) for x in xs]
        return checked(f, self.curve, xs, points, self.target)


def checked(f, curve, xs, points, target):
    total = None
    for point in points:
        if not curve.onCurve(point):
            raise ArithmeticError('off-curve extraction')
        total = curve.add(total, point)
    if total != target:
        raise ArithmeticError('signed relation failed')
    return tuple(sorted(f.toCoords(x) for x in xs)), [[f.toCoords(v) for v in p] for p in points]


class PairTable:
    """S3 pair abscissas, built without group-law pair enumeration.

    u=x+y, v=xy: S3(x,y,q)=u²q²+vq+v²+B. Quadratic solving
    gives the two possible intermediate abscissas. All setup is charged.
    """
    def __init__(self, f, curve, space, deadline):
        self.f, self.curve, self.space = f, curve, space
        self.deadline = deadline
        self.table = {}
        self.base = {}
        invs = dict(zip(space.values[1:], hybrid.batchInverse(f, space.values[1:])))
        for x in space.values[1:]:
            self.check()
            rhs = f.add(f.add(x, curve.a2), f.mul(curve.a6, f.sqr(invs[x])))
            roots = curve.asRoots(rhs)
            if roots:
                point = x, f.mul(x, roots[0])
                if not curve.onCurve(point):
                    raise ArithmeticError('invalid factor point')
                self.base[x] = point
        xs = list(self.base)
        self.pairs = 0
        for i, x in enumerate(xs):
            for y in xs[i + 1:]:
                self.check()
                u, v = f.add(x, y), f.mul(x, y)
                invV = f.mul(invs[x], invs[y])
                rhs = f.mul(f.sqr(u), f.add(f.one(), f.mul(curve.a6, f.sqr(invV))))
                scale = f.mul(v, f.sqr(invs[u]))
                roots = curve.asRoots(rhs)
                if len(roots) != 2:
                    raise ArithmeticError('liftable pair must have two S3 roots')
                for w in roots:
                    q = f.mul(scale, w)
                    self.table.setdefault(q, []).append((x, y))
                self.pairs += 1

    def check(self):
        if time.perf_counter() >= self.deadline:
            raise TimeoutError('pair table deadline')

    def solve(self, target, mode, deadline):
        self.deadline = deadline
        f, c = self.f, self.curve
        solutions = {}
        r = target[0]
        for z, p3 in self.base.items():
            self.check()
            if z == r:
                continue
            for p in (p3, c.neg(p3)):
                q = c.add(target, c.neg(p))
                if q is None:
                    continue
                for x, y in self.table.get(q[0], []):
                    if z in (x, y) or r in (x, y):
                        continue
                    key = tuple(sorted(f.toCoords(v) for v in (x, y, z)))
                    if key in solutions:
                        continue
                    found = None
                    for px in (self.base[x], c.neg(self.base[x])):
                        for py in (self.base[y], c.neg(self.base[y])):
                            if c.add(px, py) == q:
                                found = checked(f, c, (x, y, z), (px, py, p), target)
                                break
                        if found:
                            break
                    if found is None:
                        raise ArithmeticError('S3 pair did not lift to the residual point')
                    solutions[found[0]] = found[1]
                    yield found
                    if mode == 'first':
                        return

    def metadata(self):
        return {'signed_base_size': 2 * len(self.base), 'unordered_pairs': self.pairs,
                'table_keys': len(self.table), 'table_entries': sum(len(v) for v in self.table.values())}
