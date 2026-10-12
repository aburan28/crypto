"""
Host-side (Python bigint) elliptic curve arithmetic over GF(p), short Weierstrass
y^2 = x^3 + a x + b. Used for seeding walks, verifying device results, and
generating small prime-order demo curves.
"""

import math
import random

INFINITY = None


class Curve:
    def __init__(self, p, a, b, order=None, generator=None):
        self.p = p
        self.a = a
        self.b = b
        self.order = order
        self.generator = generator

    def isOnCurve(self, point):
        if point is INFINITY:
            return True
        x, y = point
        return (y * y - (x * x * x + self.a * x + self.b)) % self.p == 0

    def negate(self, point):
        if point is INFINITY:
            return INFINITY
        x, y = point
        return (x, (-y) % self.p)

    def add(self, pointA, pointB):
        if pointA is INFINITY:
            return pointB
        if pointB is INFINITY:
            return pointA
        p = self.p
        x1, y1 = pointA
        x2, y2 = pointB
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return INFINITY
            lam = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def mul(self, scalar, point):
        result = INFINITY
        addend = point
        while scalar > 0:
            if scalar & 1:
                result = self.add(result, addend)
            addend = self.add(addend, addend)
            scalar >>= 1
        return result

    def randomPoint(self, rng):
        p = self.p
        assert p % 4 == 3, "sqrt shortcut needs p = 3 mod 4"
        while True:
            x = rng.randrange(p)
            rhs = (x * x * x + self.a * x + self.b) % p
            y = pow(rhs, (p + 1) // 4, p)
            if y * y % p == rhs:
                return (x, y)


# ----------------------------------------------------------------------------
# primality and curve generation
# ----------------------------------------------------------------------------

def isProbablePrime(n, rounds=32, rng=None):
    if n < 2:
        return False
    for small in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % small == 0:
            return n == small
    rng = rng or random.Random(0)
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for _ in range(rounds):
        a = rng.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def randomPrime(bits, rng, residueMod4=3):
    while True:
        candidate = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if candidate % 4 != residueMod4:
            continue
        if isProbablePrime(candidate, rng=rng):
            return candidate


def pointOrderInHasseInterval(curve, point):
    """Smallest m in the Hasse interval with m*P = O, via baby-step giant-step.

    O(p^(1/4)) curve operations. Returns None if no multiple is found (cannot
    happen for a point on the curve, but kept for safety).
    """
    p = curve.p
    sqrtP = math.isqrt(p)
    lo = p + 1 - 2 * sqrtP - 1
    hi = p + 1 + 2 * sqrtP + 1
    width = hi - lo + 1
    stride = math.isqrt(width) + 1

    babyTable = {}
    current = INFINITY
    for j in range(stride):
        key = current if current is INFINITY else current
        babyTable.setdefault(key, j)
        current = curve.add(current, point)

    giant = curve.mul(stride, point)
    base = curve.mul(lo, point)
    best = None
    for i in range(stride + 1):
        # base = (lo + i*stride) P ; want base + jP = O  <=>  -base = jP
        target = curve.negate(base)
        if target in babyTable:
            m = lo + i * stride + babyTable[target]
            if m > 0 and (best is None or m < best):
                best = m
        base = curve.add(base, giant)
    return best


def generatePrimeOrderCurve(bits, rng):
    """Random curve over a `bits`-bit prime whose group order is prime.

    For a random point P we find the smallest m in the Hasse interval with mP = O.
    If m is prime and m > hi/2 then ord(P) = m and #E = m (ord(P) | m prime, and
    #E is a multiple of ord(P) inside [lo, hi] with hi < 2m).
    """
    p = randomPrime(bits, rng)
    hi = p + 1 + 2 * math.isqrt(p) + 1
    while True:
        a = rng.randrange(p)
        b = rng.randrange(p)
        if (4 * a * a * a + 27 * b * b) % p == 0:
            continue
        curve = Curve(p, a, b)
        point = curve.randomPoint(rng)
        m = pointOrderInHasseInterval(curve, point)
        if m is None or m <= hi // 2:
            continue
        if not isProbablePrime(m, rng=rng):
            continue
        assert curve.mul(m, point) is INFINITY
        curve.order = m
        curve.generator = point
        return curve


SECP256K1 = Curve(
    p=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F,
    a=0,
    b=7,
    order=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141,
    generator=(
        0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
        0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
    ),
)
