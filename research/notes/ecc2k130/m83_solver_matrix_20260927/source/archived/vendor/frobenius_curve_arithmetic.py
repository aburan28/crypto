"""Pure-Python arithmetic copied unchanged from the earlier verified harness.
No native solver dependency. Field irreducibility check assumes prime degree.
The lifting benchmark constructs Curve without its factorization initializer.
"""

class Field:
    def __init__(self, n, modulus):
        self.n, self.modulus, self.mask = n, modulus, (1 << n)-1
        basis = [self.mul(1 << i, 1 << i) for i in range(n)]
        self.square_tables = []
        for offset in range(0, n, 8):
            table = []
            for b in range(256):
                v = 0
                for j in range(8):
                    if offset+j < n and (b >> j) & 1:
                        v ^= basis[offset+j]
                table.append(v)
            self.square_tables.append(table)
        # Rabin irreducibility criterion; n is prime in these experiments.
        x = 2
        for _ in range(n):
            x = self.sq(x)
        assert x == 2
        a, b = modulus, self.sq(2) ^ 2
        while b:
            a, b = b, self.remainder(a, b)
        assert a == 1

    @staticmethod
    def remainder(a, b):
        while a and a.bit_length() >= b.bit_length():
            a ^= b << (a.bit_length()-b.bit_length())
        return a

    def mul(self, a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> self.n:
                a ^= self.modulus
        return r

    def sq(self, a):
        r = 0
        for table in self.square_tables:
            r ^= table[a & 255]
            a >>= 8
        return r

    def inv(self, a):
        if not a:
            raise ZeroDivisionError
        u, v, g, h = a, self.modulus, 1, 0
        while u != 1:
            j = u.bit_length()-v.bit_length()
            if j < 0:
                u, v, g, h, j = v, u, h, g, -j
            u ^= v << j
            g ^= h << j
        return self.remainder(g, self.modulus)

    def trace(self, x):
        t, y = 0, x
        for _ in range(self.n):
            t ^= y
            y = self.sq(y)
        assert t in (0, 1)
        return t

class Curve:
    """E_0: y^2 + xy = x^3 + 1."""
    def __init__(self, field):
        self.f = field
        t0, t1 = 2, -1
        for _ in range(2, field.n+1):
            t0, t1 = t1, -t1-2*t0
        self.order = (1 << field.n)+1-t1
        v, p, factors = self.order, 2, []
        while p*p <= v:
            while v % p == 0:
                factors.append(p)
                v //= p
            p += 1
        if v > 1:
            factors.append(v)
        self.r = max(factors)
        self.cofactor = self.order//self.r
        for x in range(1, 1 << field.n):
            lifts = self.lift(x)
            if lifts:
                g = self.scale(lifts[0], self.cofactor)
                if g:
                    self.g = g
                    break
        assert self.scale(self.g, self.r) is None
        roots = [v for v in range(self.r) if (v*v+v+2) % self.r == 0]
        matches = [v for v in roots if self.scale(self.g, v) == self.frob(self.g)]
        assert len(matches) == 1
        self.lam = matches[0]
        assert pow(self.lam, field.n, self.r) == 1 and self.lam != 1

    def neg(self, p):
        return (p[0], p[0] ^ p[1]) if p else None

    def valid(self, p):
        if p is None:
            return True
        x, y = p
        f = self.f
        return f.sq(y) ^ f.mul(x, y) == f.mul(f.sq(x), x) ^ 1

    def add(self, p, q):
        if p is None:
            return q
        if q is None:
            return p
        x, y = p
        xx, yy = q
        f = self.f
        if x == xx:
            if y != yy or x == 0:
                return None
            slope = x ^ f.mul(y, f.inv(x))
            outx = f.sq(slope) ^ slope
            return outx, f.sq(x) ^ f.mul(slope ^ 1, outx)
        slope = f.mul(y ^ yy, f.inv(x ^ xx))
        outx = f.sq(slope) ^ slope ^ x ^ xx
        return outx, f.mul(slope, x ^ outx) ^ outx ^ y

    def scale(self, p, k):
        out = None
        while k:
            if k & 1:
                out = self.add(out, p)
            p = self.add(p, p)
            k >>= 1
        return out

    def frob(self, p, power=1):
        for _ in range(power % self.f.n):
            if p:
                p = self.f.sq(p[0]), self.f.sq(p[1])
        return p

    def lift(self, x):
        if x == 0:
            return [(0, 1)]
        f = self.f
        rhs = x ^ f.sq(f.inv(x))
        if f.trace(rhs):
            return []
        z, y = 0, rhs
        for _ in range((f.n+1)//2):
            z ^= y
            y = f.sq(f.sq(y))
        assert f.sq(z) ^ z == rhs
        p = x, f.mul(x, z)
        assert self.valid(p)
        return [p, self.neg(p)]

def binary_rank(values):
    pivots = {}
    for v in values:
        while v:
            j = v.bit_length()-1
            if j in pivots:
                v ^= pivots[j]
            else:
                pivots[j] = v
                break
    return len(pivots)
