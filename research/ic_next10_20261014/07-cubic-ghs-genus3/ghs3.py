#!/usr/bin/env python3
"""Item 7: cubic GHS (m=3, r=5, split arm) transfer E0(F_{q^3}) -> Jac(C)(F_q),
genus-3 double-large-prime index calculus, and matched rho on E0, all with
explicit F_q operation counters.

Counters: CNT[0] = F_q multiplications/squarings, CNT[1] = F_q inversions,
CNT[2] = multiplications mod ell (linear algebra), CNT[3] = inversions mod ell,
CNT[4] = group operations (Jacobian or elliptic additions/doublings).
"""
import json, math, random, subprocess, sys, time

CNT = [0, 0, 0, 0, 0]
RHO_RANDOM_MULT = bool(int(__import__('os').environ.get('RHO_RANDOM_MULT', '0')))


def reset():
    for i in range(5):
        CNT[i] = 0


def snap():
    return list(CNT)


def is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
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


# ---------------------------------------------------------------- F_p polys
# lists low -> high, trimmed (no trailing zeros); zero poly = []

def tr(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def padd(a, b, p):
    if len(a) < len(b):
        a, b = b, a
    r = a[:]
    for i, c in enumerate(b):
        r[i] = (r[i] + c) % p
    return tr(r)


def psub(a, b, p):
    n = max(len(a), len(b))
    r = [0] * n
    for i, c in enumerate(a):
        r[i] = c
    for i, c in enumerate(b):
        r[i] = (r[i] - c) % p
    return tr(r)


def pneg(a, p):
    return [(-c) % p for c in a]


def pmul(a, b, p):
    if not a or not b:
        return []
    CNT[0] += len(a) * len(b)
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return tr([c % p for c in r])


def pscal(a, c, p):
    if c == 1:
        return a[:]
    CNT[0] += len(a)
    return tr([x * c % p for x in a])


def pdivmod(a, b, p):
    lb = len(b)
    if len(a) < lb:
        return [], a[:]
    r = a[:]
    lead = b[-1]
    inv = 1
    if lead != 1:
        CNT[1] += 1
        inv = pow(lead, p - 2, p)
    q = [0] * (len(a) - lb + 1)
    for i in range(len(a) - lb, -1, -1):
        c = r[i + lb - 1]
        if c:
            if inv != 1:
                c = c * inv % p
                CNT[0] += 1
            q[i] = c
            CNT[0] += lb - 1
            for j in range(lb - 1):
                r[i + j] = (r[i + j] - c * b[j]) % p
            r[i + lb - 1] = 0
    return tr(q), tr(r[:lb - 1])


def pmod(a, b, p):
    return pdivmod(a, b, p)[1]


def pmonic(a, p):
    if not a or a[-1] == 1:
        return a[:]
    CNT[1] += 1
    inv = pow(a[-1], p - 2, p)
    CNT[0] += len(a) - 1
    return [c * inv % p for c in a[:-1]] + [1]


def pxgcd(a, b, p):
    """returns g, s, t with s a + t b = g, g monic."""
    r0, r1 = a[:], b[:]
    s0, s1 = [1], []
    t0, t1 = [], [1]
    while r1:
        q, r = pdivmod(r0, r1, p)
        r0, r1 = r1, r
        s0, s1 = s1, psub(s0, pmul(q, s1, p), p)
        t0, t1 = t1, psub(t0, pmul(q, t1, p), p)
    if not r0:
        return [], [], []
    lead = r0[-1]
    if lead != 1:
        CNT[1] += 1
        inv = pow(lead, p - 2, p)
        r0 = pscal(r0, inv, p); s0 = pscal(s0, inv, p); t0 = pscal(t0, inv, p)
    return r0, s0, t0


def pev(a, x, p):
    r = 0
    for c in reversed(a):
        r = (r * x + c) % p
    CNT[0] += max(len(a) - 1, 0)
    return r


# ---------------------------------------------------------- Cantor over F_p
G3 = 3


def creduce(u, v, h, p):
    while len(u) - 1 > G3:
        u2, rem = pdivmod(psub(h, pmul(v, v, p), p), u, p)
        assert not rem
        u = pmonic(u2, p)
        v = pmod(pneg(v, p), u, p)
    u = pmonic(u, p)
    return u, pmod(v, u, p)


def cadd(D1, D2, h, p):
    CNT[4] += 1
    u1, v1 = D1
    u2, v2 = D2
    if u1 == [1]:
        return (u2[:], v2[:])
    if u2 == [1]:
        return (u1[:], v1[:])
    d1, e1, e2 = pxgcd(u1, u2, p)
    if d1 == [1]:
        # CRT: v = v1 + u1 * (e1 (v2 - v1) mod u2)
        t = pmod(pmul(e1, psub(v2, v1, p), p), u2, p)
        u = pmul(u1, u2, p)
        v = padd(v1, pmul(u1, t, p), p)
        return creduce(u, v, h, p)
    d, c1, c2 = pxgcd(d1, padd(v1, v2, p), p)
    s1 = pmul(c1, e1, p)
    s2 = pmul(c1, e2, p)
    s3 = c2
    u = pdivmod(pmul(u1, u2, p), pmul(d, d, p), p)[0]
    num = padd(padd(pmul(pmul(s1, u1, p), v2, p), pmul(pmul(s2, u2, p), v1, p), p),
               pmul(s3, padd(pmul(v1, v2, p), h, p), p), p)
    v, rem = pdivmod(num, d, p)
    assert not rem
    v = pmod(v, u, p)
    return creduce(u, v, h, p)


def cneg(D, p):
    return (D[0], pneg(D[1], p))


def cmul(n, D, h, p):
    R = ([1], [])
    B = D
    while n > 0:
        if n & 1:
            R = cadd(R, B, h, p)
        n >>= 1
        if n:
            B = cadd(B, B, h, p)
    return R


# ------------------------------------------------------------ K = F_{p^3}
class KField:
    """K = F_p[t]/(t^3 - c1 t - c0)."""

    def __init__(self, p, rng):
        self.p = p
        for c1 in list(range(1, 20)) + [rng.randrange(1, p) for _ in range(200)]:
            for _ in range(200):
                c0 = rng.randrange(1, p)
                if self._irreducible(c1, c0):
                    self.c1, self.c0 = c1, c0
                    self._setup()
                    return
        raise RuntimeError("no modulus")

    def _irreducible(self, c1, c0):
        p = self.p
        f = [(-c0) % p, (-c1) % p, 0, 1]
        # cubic irreducible iff no root
        old = CNT[:]
        r = [0, 1]
        res = [1]
        e = p
        base = r
        while e:
            if e & 1:
                res = pmod(pmul(res, base, p), f, p)
            base = pmod(pmul(base, base, p), f, p)
            e >>= 1
        g = pxgcd(f, psub(res, [0, 1], p), p)[0]
        CNT[:] = old
        return g == [1]

    def _setup(self):
        p = self.p
        old = CNT[:]
        self.T1 = self.pw((0, 1, 0), p)
        self.T2 = self.mul(self.T1, self.T1)
        CNT[:] = old

    def mul(self, a, b):
        p = self.p
        a0, a1, a2 = a
        b0, b1, b2 = b
        d0 = a0 * b0
        d1 = a0 * b1 + a1 * b0
        d2 = a0 * b2 + a1 * b1 + a2 * b0
        d3 = a1 * b2 + a2 * b1
        d4 = a2 * b2
        c0, c1 = self.c0, self.c1
        CNT[0] += 9 + 2 + (0 if c1 == 1 else 2)
        return ((d0 + c0 * d3) % p, (d1 + c1 * d3 + c0 * d4) % p, (d2 + c1 * d4) % p)

    def sqr(self, a):
        p = self.p
        a0, a1, a2 = a
        d0 = a0 * a0
        d1 = 2 * a0 * a1
        d2 = 2 * a0 * a2 + a1 * a1
        d3 = 2 * a1 * a2
        d4 = a2 * a2
        c0, c1 = self.c0, self.c1
        CNT[0] += 6 + 2 + (0 if c1 == 1 else 2)
        return ((d0 + c0 * d3) % p, (d1 + c1 * d3 + c0 * d4) % p, (d2 + c1 * d4) % p)

    def add(self, a, b):
        p = self.p
        return ((a[0] + b[0]) % p, (a[1] + b[1]) % p, (a[2] + b[2]) % p)

    def sub(self, a, b):
        p = self.p
        return ((a[0] - b[0]) % p, (a[1] - b[1]) % p, (a[2] - b[2]) % p)

    def neg(self, a):
        p = self.p
        return ((-a[0]) % p, (-a[1]) % p, (-a[2]) % p)

    def smul(self, a, c):
        p = self.p
        CNT[0] += 3
        return (a[0] * c % p, a[1] * c % p, a[2] * c % p)

    def frob(self, a):
        p = self.p
        T1, T2 = self.T1, self.T2
        CNT[0] += 6
        return ((a[0] + a[1] * T1[0] + a[2] * T2[0]) % p,
                (a[1] * T1[1] + a[2] * T2[1]) % p,
                (a[1] * T1[2] + a[2] * T2[2]) % p)

    def inv(self, a):
        p = self.p
        b1 = self.frob(a)
        b2 = self.frob(b1)
        b = self.mul(b1, b2)
        N = (a[0] * b[0] + self.c0 * (a[1] * b[2] + a[2] * b[1])) % p
        CNT[0] += 4
        CNT[1] += 1
        Ni = pow(N, p - 2, p)
        return self.smul(b, Ni)

    def pw(self, a, e):
        r = (1, 0, 0)
        b = a
        while e:
            if e & 1:
                r = self.mul(r, b)
            e >>= 1
            if e:
                b = self.sqr(b)
        return r

    def is_zero(self, a):
        return a == (0, 0, 0)

    def emb(self, c):
        return (c % self.p, 0, 0)

    def in_fp(self, a):
        return a[1] == 0 and a[2] == 0

    def rand(self, rng):
        p = self.p
        return (rng.randrange(p), rng.randrange(p), rng.randrange(p))

    def sqrt(self, a, rng):
        """Tonelli-Shanks in K (not counted: instance generation only)."""
        old = CNT[:]
        try:
            if self.is_zero(a):
                return a
            n = self.p ** 3
            if self.pw(a, (n - 1) // 2) != (1, 0, 0):
                return None
            Q, S = n - 1, 0
            while Q % 2 == 0:
                Q //= 2; S += 1
            while True:
                z = self.rand(rng)
                if not self.is_zero(z) and self.pw(z, (n - 1) // 2) != (1, 0, 0):
                    break
            M, c, t, R = S, self.pw(z, Q), self.pw(a, Q), self.pw(a, (Q + 1) // 2)
            while t != (1, 0, 0):
                i, tt = 0, t
                while tt != (1, 0, 0):
                    tt = self.sqr(tt); i += 1
                b = c
                for _ in range(M - i - 1):
                    b = self.sqr(b)
                M, c = i, self.sqr(b)
                t, R = self.mul(t, c), self.mul(R, b)
            return R
        finally:
            CNT[:] = old


# generic polys over K (tuples), lists low->high
KZ = (0, 0, 0)


def ktr(a):
    while a and a[-1] == KZ:
        a.pop()
    return a


def kpadd(K, a, b):
    n = max(len(a), len(b))
    return ktr([K.add(a[i] if i < len(a) else KZ, b[i] if i < len(b) else KZ) for i in range(n)])


def kpsub(K, a, b):
    n = max(len(a), len(b))
    return ktr([K.sub(a[i] if i < len(a) else KZ, b[i] if i < len(b) else KZ) for i in range(n)])


def kpneg(K, a):
    return [K.neg(c) for c in a]


def kpmul(K, a, b):
    if not a or not b:
        return []
    r = [KZ] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x != KZ:
            for j, y in enumerate(b):
                if y != KZ:
                    r[i + j] = K.add(r[i + j], K.mul(x, y))
    return ktr(r)


def kpscal(K, a, c):
    return ktr([K.mul(x, c) for x in a])


def kpdivmod(K, a, b):
    lb = len(b)
    if len(a) < lb:
        return [], a[:]
    r = a[:]
    inv = K.inv(b[-1])
    q = [KZ] * (len(a) - lb + 1)
    for i in range(len(a) - lb, -1, -1):
        c = r[i + lb - 1]
        if c != KZ:
            c = K.mul(c, inv)
            q[i] = c
            for j in range(lb - 1):
                r[i + j] = K.sub(r[i + j], K.mul(c, b[j]))
            r[i + lb - 1] = KZ
    return ktr(q), ktr(r[:lb - 1])


def kpmod(K, a, b):
    return kpdivmod(K, a, b)[1]


def kpmonic(K, a):
    if not a:
        return a
    return kpscal(K, a, K.inv(a[-1]))


def kpxgcd(K, a, b):
    r0, r1 = a[:], b[:]
    s0, s1 = [(1, 0, 0)], []
    t0, t1 = [], [(1, 0, 0)]
    while r1:
        q, r = kpdivmod(K, r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, kpsub(K, s0, kpmul(K, q, s1))
        t0, t1 = t1, kpsub(K, t0, kpmul(K, q, t1))
    inv = K.inv(r0[-1])
    return kpscal(K, r0, inv), kpscal(K, s0, inv), kpscal(K, t0, inv)


def kcreduce(K, u, v, h):
    while len(u) - 1 > G3:
        u2, rem = kpdivmod(K, kpsub(K, h, kpmul(K, v, v)), u)
        assert not rem
        u = kpmonic(K, u2)
        v = kpmod(K, kpneg(K, v), u)
    u = kpmonic(K, u)
    return u, kpmod(K, v, u)


def kcadd(K, D1, D2, h):
    CNT[4] += 1
    u1, v1 = D1
    u2, v2 = D2
    one = [(1, 0, 0)]
    if u1 == one:
        return D2
    if u2 == one:
        return D1
    d1, e1, e2 = kpxgcd(K, u1, u2)
    if d1 == one:
        t = kpmod(K, kpmul(K, e1, kpsub(K, v2, v1)), u2)
        return kcreduce(K, kpmul(K, u1, u2), kpadd(K, v1, kpmul(K, u1, t)), h)
    d, c1, c2 = kpxgcd(K, d1, kpadd(K, v1, v2))
    s1 = kpmul(K, c1, e1); s2 = kpmul(K, c1, e2); s3 = c2
    u = kpdivmod(K, kpmul(K, u1, u2), kpmul(K, d, d))[0]
    num = kpadd(K, kpadd(K, kpmul(K, kpmul(K, s1, u1), v2), kpmul(K, kpmul(K, s2, u2), v1)),
                kpmul(K, s3, kpadd(K, kpmul(K, v1, v2), h)))
    v, rem = kpdivmod(K, num, d)
    assert not rem
    return kcreduce(K, u, kpmod(K, v, u), h)


# -------------------------------------------- elliptic curve over K (short W.)
class ECW:
    def __init__(self, K, A, B):
        self.K, self.A, self.B = K, A, B

    def add(self, P, Q):
        CNT[4] += 1
        K = self.K
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if K.add(y1, y2) == KZ:
                return None
            lam = K.mul(K.add(K.smul(K.sqr(x1), 3), self.A), K.inv(K.add(y1, y1)))
        else:
            lam = K.mul(K.sub(y2, y1), K.inv(K.sub(x2, x1)))
        x3 = K.sub(K.sub(K.sqr(lam), x1), x2)
        y3 = K.sub(K.mul(lam, K.sub(x1, x3)), y1)
        return (x3, y3)

    def mul(self, n, P):
        R = None
        B = P
        while n > 0:
            if n & 1:
                R = self.add(R, B)
            n >>= 1
            if n:
                B = self.add(B, B)
        return R

    def on(self, P):
        if P is None:
            return True
        K = self.K
        x, y = P
        return K.sqr(y) == K.add(K.add(K.mul(K.sqr(x), x), K.mul(self.A, x)), self.B)


# ---------------------------------------------------------- instance build
def gp_order(h, p):
    coeffs = ",".join(str(c) for c in reversed(h))
    cmd = f"h=Pol([{coeffs}])*Mod(1,{p});P=hyperellcharpoly(h);print(subst(P,'x,1))"
    out = subprocess.run(["gp", "-q", "-f", "-s", "400000000"], input=cmd, capture_output=True, text=True, timeout=600)
    return int(out.stdout.strip())


def gp_order_E(K, A, B):
    p = K.p
    def kstr(a):
        return f"({a[0]}+{a[1]}*t+{a[2]}*t^2)"
    cmd = (f"t=ffgen(Mod(1,{p})*(x^3-{K.c1}*x-{K.c0}),'t);"
           f"E=ellinit([{kstr(A)},{kstr(B)}]);print(ellcard(E))")
    out = subprocess.run(["gp", "-q", "-f", "-s", "400000000"], input=cmd, capture_output=True, text=True, timeout=600)
    return int(out.stdout.strip())


class Instance:
    pass


def build_instance(p, rng, cmax=64):
    """Return an Instance or None (rejected)."""
    old = CNT[:]
    K = KField(p, rng)
    while True:
        pr = K.rand(rng)
        if not K.in_fp(pr):
            break
    r0 = pr
    r1 = K.frob(r0)
    r2 = K.frob(r1)
    a, b = rng.sample(range(p), 2)
    inst = Instance()
    inst.p, inst.K, inst.a, inst.b, inst.r = p, K, a, b, (r0, r1, r2)
    # conic Tr(c0 s^2) = 0, c0 = r2 - r0
    c0 = K.sub(r2, r0)
    trb = []
    for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        x = e
        s = KZ
        for _ in range(3):
            s = K.add(s, x)
            x = K.frob(x)
        assert K.in_fp(s)
        trb.append(s[0])

    def Tr(x):
        return (x[0] * trb[0] + x[1] * trb[1] + x[2] * trb[2]) % p

    def Qf(s):
        return Tr(K.mul(c0, K.sqr(s)))

    def Bf(s1, s2):
        return Tr(K.mul(c0, K.mul(s1, s2)))

    # find a conic point s*
    while True:
        s1, s2 = K.rand(rng), K.rand(rng)
        A2, A1, A0 = Qf(s2), 2 * Bf(s1, s2) % p, Qf(s1)
        if A2 == 0:
            continue
        disc = (A1 * A1 - 4 * A2 * A0) % p
        if disc == 0 or pow(disc, (p - 1) // 2, p) != 1:
            continue
        sq = sqrt_fp(disc, p)
        lam = (-A1 + sq) * pow(2 * A2, p - 2, p) % p
        sstar = K.add(s1, K.smul(s2, lam))
        if sstar != KZ:
            break
    assert Qf(sstar) == 0

    def param(e1, e2):
        # s(z) = -Q(e(z)) s* + 2 B(s*, e(z)) e(z), e(z) = e1 + z e2
        q0, q1, q2 = Qf(e1), 2 * Bf(e1, e2) % p, Qf(e2)
        b0, b1 = Bf(sstar, e1), Bf(sstar, e2)
        # -Q s*: coefficients z^0..z^2 ; 2B e: (b0 + b1 z)(e1 + z e2)*2
        c = [K.sub(K.smul(e1, 2 * b0 % p), K.smul(sstar, q0)),
             K.sub(K.add(K.smul(e1, 2 * b1 % p), K.smul(e2, 2 * b0 % p)), K.smul(sstar, q1)),
             K.sub(K.smul(e2, 2 * b1 % p), K.smul(sstar, q2))]
        return ktr(c)

    def build_from(sz):
        X0 = sz
        X1 = [K.frob(c) for c in X0]
        X2 = [K.frob(c) for c in X1]
        X0s, X2s = kpmul(K, X0, X0), kpmul(K, X2, X2)
        Delta = kpsub(K, X0s, X2s)
        Num = kpsub(K, kpscal(K, X0s, r0), kpscal(K, X2s, r1))
        kap = K.inv(Delta[-1])
        xd = kpscal(K, Delta, kap)
        xn = kpscal(K, Num, kap)
        if not all(K.in_fp(c) for c in xd + xn):
            return None
        xdp = tr([c[0] for c in xd]); xnp = tr([c[0] for c in xn])
        h = pmul(psub(xnp, pscal(xdp, a, p), p), psub(xnp, pscal(xdp, b, p), p), p)
        return dict(X0=X0, X2=X2, Delta=Delta, kap=kap, xd=xdp, xn=xnp, h=h)

    while True:
        e1, e2 = K.rand(rng), K.rand(rng)
        sz = param(e1, e2)
        bf = build_from(sz)
        if bf is None:
            CNT[:] = old
            raise RuntimeError("x(z) not F_q-rational: Frobenius ordering wrong")
        h = bf["h"]
        roots = [z for z in range(p) if pev(h, z, p) == 0]
        if roots:
            break
        CNT[:] = old
        return None  # no F_q Weierstrass point on this conic (should not happen in split arm)
    # s(z0) for a root z0 whose point is not s* itself
    sW = None
    for z0 in roots:
        acc = KZ
        for c in reversed(sz):
            acc = K.add(K.smul(acc, z0) if acc != KZ else KZ, c)
        if Bf(sstar, acc) != 0:
            sW = acc
            break
    if sW is None:
        CNT[:] = old
        return None
    assert Qf(sW) == 0
    for _ in range(1000):
        e1 = K.rand(rng)
        sz2 = param(e1, sW)
        if len(sz2) == 3:
            bf = build_from(sz2)
            if bf is not None and len(bf["h"]) == 8:
                break
    else:
        CNT[:] = old
        return None
    h = bf["h"]
    # squarefree check
    if pxgcd(h, pder(h, p), p)[0] != [1]:
        CNT[:] = old
        return None
    inst.h, inst.xn, inst.xd = h, bf["xn"], bf["xd"]
    X0, X2, Delta, kap = bf["X0"], bf["X2"], bf["Delta"], bf["kap"]
    inst.numR = kpscal(K, kpmul(K, Delta, Delta), kap)
    inst.denR = kpscal(K, kpmul(K, X0, X2), K.sub(r0, r1))
    # Weierstrass model of E0: y^2 = (x-r0)(x-r1)(x-a)(x-b)
    rho = [r1, K.emb(a), K.emb(b)]
    # Y^2 = prod (1 + (r0 - rho) X) = l3 X^3 + l2 X^2 + l1 X + 1
    poly = [(1, 0, 0)]
    for rr in rho:
        poly = kpmul(K, poly, [(1, 0, 0), K.sub(r0, rr)])
    l1, l2, l3 = poly[1], poly[2], poly[3]
    # Y'^2 = X'^3 + l2 X'^2 + l1 l3 X' + l3^2, X' = l3 X, Y' = l3 Y
    # depress: X' = x - l2/3
    i3 = pow(3, p - 2, p)
    s = K.smul(l2, i3)
    A = K.sub(K.mul(l1, l3), K.mul(K.smul(l2, i3), l2))
    # B = l3^2 - s*(l1 l3) + 2 s^3 ... compute directly
    # f(x - s) with f(X) = X^3 + l2 X^2 + c1 X + c0
    c1, cc0 = K.mul(l1, l3), K.sqr(l3)
    ms = K.neg(s)
    B = K.add(K.add(K.add(K.mul(K.sqr(ms), ms), K.mul(l2, K.sqr(ms))), K.mul(c1, ms)), cc0)
    inst.l3, inst.shift = l3, s
    inst.E = ECW(K, A, B)
    CNT[:] = old
    # order: #Jac(C)(F_q) = #E0(K), counted on E0 by SEA in PARI
    N = gp_order_E(K, A, B)
    inst.N = N
    ell, c = None, None
    for cc in range(4, cmax + 1, 4):
        if N % cc == 0 and is_prime(N // cc) and (N // cc) % 2 == 1 and N // cc > cmax:
            ell, c = N // cc, cc
            break
    if ell is None:
        return None
    inst.ell, inst.c = ell, c
    return inst


def pder(a, p):
    return tr([(i * a[i]) % p for i in range(1, len(a))])


def sqrt_fp(a, p):
    a %= p
    if a == 0:
        return 0
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    Q, S = p - 1, 0
    while Q % 2 == 0:
        Q //= 2; S += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    M, c, t, R = S, pow(z, Q, p), pow(a, Q, p), pow(a, (Q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p; i += 1
        b = pow(c, 1 << (M - i - 1), p)
        M, c, t, R = i, b * b % p, t * b * b % p, R * b % p
    return R


def w_to_quartic(inst, P):
    """short-W point -> quartic point (x, y) of E0 : y^2 = F0(x)."""
    K = inst.K
    x, y = P
    Xp = K.sub(x, inst.shift)          # X'
    if Xp == KZ:
        return None
    X = K.mul(Xp, K.inv(inst.l3))       # X = X'/l3
    Y = K.mul(y, K.inv(inst.l3))
    Xi = K.inv(X)
    return (K.add(inst.r[0], Xi), K.mul(Y, K.sqr(Xi)))


def random_point(inst, rng):
    K, E = inst.K, inst.E
    old = CNT[:]
    while True:
        x = K.rand(rng)
        rhs = K.add(K.add(K.mul(K.sqr(x), x), K.mul(E.A, x)), E.B)
        y = K.sqrt(rhs, rng)
        if y is not None:
            CNT[:] = old
            return (x, y)


def transfer(inst, P):
    """phi(P) in Jac(C)(F_q) as an F_p Mumford pair; P on short-W model."""
    K, p = inst.K, inst.p
    Q = w_to_quartic(inst, P)
    if Q is None:
        raise ValueError("exceptional")
    xP, yP = Q
    xn = [(c, 0, 0) for c in inst.xn]
    xd = [(c, 0, 0) for c in inst.xd]
    U = kpsub(K, xn, kpscal(K, xd, xP))
    if len(U) != 5:
        raise ValueError("exceptional")
    U = kpmonic(K, U)
    num = kpmod(K, inst.numR, U)
    den = kpmod(K, inst.denR, U)
    g, s, _ = kpxgcd(K, den, U)
    if g != [(1, 0, 0)]:
        raise ValueError("exceptional")
    v = kpmod(K, kpscal(K, kpmul(K, num, s), yP), U)
    hK = [(c, 0, 0) for c in inst.h]
    assert kpmod(K, kpsub(K, kpmul(K, v, v), hK), U) == [], "v^2 != h mod U"
    D = kcreduce(K, U, v, hK)
    D1 = ([K.frob(c) for c in D[0]], [K.frob(c) for c in D[1]])
    D2 = ([K.frob(c) for c in D1[0]], [K.frob(c) for c in D1[1]])
    S = kcadd(K, kcadd(K, D, D1, hK), D2, hK)
    assert all(K.in_fp(c) for c in S[0] + S[1]), "phi(P) not F_q-rational"
    return (tr([c[0] for c in S[0]]), tr([c[0] for c in S[1]]))


# ------------------------------------------------------------------- rho
def rho(inst, G, T, rng, r=20, max_steps=None):
    """single-walk r-adding rho with distinguished points; returns k or None."""
    E, K, ell = inst.E, inst.K, inst.ell
    # multipliers by a Fibonacci chain from two random combinations (4 scalar mults)
    nrand = r if RHO_RANDOM_MULT else 2
    coeffs = [(rng.randrange(ell), rng.randrange(1, ell)) for _ in range(nrand)]
    Ms = [E.add(E.mul(a, G), E.mul(b, T)) for a, b in coeffs]
    for j in range(nrand, r):
        Ms.append(E.add(Ms[-1], Ms[-2]))
        coeffs.append(((coeffs[-1][0] + coeffs[-2][0]) % ell, (coeffs[-1][1] + coeffs[-2][1]) % ell))
    R = E.add(Ms[0], Ms[r - 1])
    a, b = (coeffs[0][0] + coeffs[r - 1][0]) % ell, (coeffs[0][1] + coeffs[r - 1][1]) % ell
    dpbits = max(2, int(math.log2(math.isqrt(ell))) - 6)
    mask = (1 << dpbits) - 1
    seen = {}
    steps = 0
    while True:
        if R is None:
            return None, steps
        hx = R[0][0] * 1000003 + R[0][1] * 7919 + R[0][2]
        if (hx >> 3) & mask == 0:
            key = R[0] + R[1]
            if key in seen:
                a2, b2 = seen[key]
                if (b2 - b) % ell:
                    k = (a - a2) * pow((b2 - b) % ell, ell - 2, ell) % ell
                    return k, steps
                return None, steps
            seen[key] = (a, b)
        j = (hx >> 1) % r
        R = E.add(R, Ms[j])
        a = (a + coeffs[j][0]) % ell
        b = (b + coeffs[j][1]) % ell
        steps += 1
        if max_steps and steps > max_steps:
            return None, steps


# --------------------------------------------------------------- IC (2LP)
def mulmod3(a, b, u, p):
    """a, b residues mod monic cubic u (len-3 lists); counted 9M + 6M."""
    a0, a1, a2 = a
    b0, b1, b2 = b
    d0 = a0 * b0
    d1 = a0 * b1 + a1 * b0
    d2 = a0 * b2 + a1 * b1 + a2 * b0
    d3 = a1 * b2 + a2 * b1
    d4 = a2 * b2
    u0, u1, u2 = u[0], u[1], u[2]
    # x^4 and x^3 reduction: x^3 = -(u2 x^2 + u1 x + u0)
    d3 -= d4 * u2; d2 -= d4 * u1; d1 -= d4 * u0
    d2 -= d3 * u2; d1 -= d3 * u1; d0 -= d3 * u0
    CNT[0] += 15
    return [d0 % p, d1 % p, d2 % p]


def sqrmod3(a, u, p):
    a0, a1, a2 = a
    d0 = a0 * a0
    d1 = 2 * a0 * a1
    d2 = 2 * a0 * a2 + a1 * a1
    d3 = 2 * a1 * a2
    d4 = a2 * a2
    u0, u1, u2 = u[0], u[1], u[2]
    d3 -= d4 * u2; d2 -= d4 * u1; d1 -= d4 * u0
    d2 -= d3 * u2; d1 -= d3 * u1; d0 -= d3 * u0
    CNT[0] += 12
    return [d0 % p, d1 % p, d2 % p]


def xmul3(a, u, p):
    """a * x mod u."""
    a0, a1, a2 = a
    CNT[0] += 3
    return [(-a2 * u[0]) % p, (a0 - a2 * u[1]) % p, (a1 - a2 * u[2]) % p]


def splits_cubic(u, p):
    """True iff monic cubic u divides x^p - x."""
    r = [0, 1, 0]
    bits = bin(p)[3:]
    for bit in bits:
        r = sqrmod3(r, u, p)
        if bit == '1':
            r = xmul3(r, u, p)
    return r == [0, 1, 0]


def splits(u, p):
    d = len(u) - 1
    if d <= 1:
        return True
    if d == 3:
        return splits_cubic(u, p)
    # deg 2: discriminant
    disc = (u[1] * u[1] - 4 * u[0]) % p
    CNT[0] += 2
    CNT[0] += int(1.5 * p.bit_length())
    return disc != 0 and pow(disc, (p - 1) // 2, p) == 1


def ppowmod(base, e, m, p):
    r = [1]
    b = pmod(base, m, p)
    while e:
        if e & 1:
            r = pmod(pmul(r, b, p), m, p)
        e >>= 1
        if e:
            b = pmod(pmul(b, b, p), m, p)
    return r


def roots_split(u, p, rng):
    """roots of a split squarefree monic u (Cantor-Zassenhaus), counted."""
    d = len(u) - 1
    if d == 0:
        return []
    if d == 1:
        return [(-u[0]) % p]
    if d == 2:
        disc = (u[1] * u[1] - 4 * u[0]) % p
        sq = sqrt_fp(disc, p)
        CNT[0] += 2 + int(1.5 * p.bit_length()) + 2
        CNT[1] += 1
        i2 = (p + 1) // 2
        return [((-u[1] + sq) * i2) % p, ((-u[1] - sq) * i2) % p]
    while True:
        delta = rng.randrange(p)
        w = ppowmod([delta, 1], (p - 1) // 2, u, p)
        g = pxgcd(u, psub(w, [1], p), p)[0]
        if 0 < len(g) - 1 < d:
            g2 = pdivmod(u, g, p)[0]
            return roots_split(g, p, rng) + roots_split(g2, p, rng)


def lanczos(rows, rhs, ncols, ell):
    """solve S x = b in least squares via Lanczos on A = S^T S. counted mod ell."""
    def Sv(v):
        out = []
        for r in rows:
            s = 0
            for j, c in r:
                s += c * v[j]
            out.append(s % ell)
        CNT[2] += nnz
        return out

    def STv(w):
        out = [0] * ncols
        for r, wi in zip(rows, w):
            if wi:
                for j, c in r:
                    out[j] += c * wi
        CNT[2] += nnz
        return [x % ell for x in out]

    def dot(a, b):
        CNT[2] += ncols
        return sum(x * y for x, y in zip(a, b)) % ell

    def inv(x):
        CNT[3] += 1
        return pow(x, ell - 2, ell)

    nnz = sum(len(r) for r in rows)
    b = STv(rhs)
    x = [0] * ncols
    w_prev, Aw_prev, wAw_prev = None, None, None
    w = b[:]
    for it in range(ncols + 5):
        Aw = STv(Sv(w))
        wAw = dot(w, Aw)
        if wAw == 0:
            break
        iw = inv(wAw)
        coef = dot(w, b) * iw % ell
        x = [(xi + coef * wi) % ell for xi, wi in zip(x, w)]
        CNT[2] += ncols
        AwAw = dot(Aw, Aw)
        alpha = AwAw * iw % ell
        nw = [(a - alpha * wi) % ell for a, wi in zip(Aw, w)]
        CNT[2] += ncols
        if w_prev is not None:
            beta = dot(Aw, Aw_prev) * inv(wAw_prev) % ell
            nw = [(a - beta * wi) % ell for a, wi in zip(nw, w_prev)]
            CNT[2] += ncols
        w_prev, Aw_prev, wAw_prev = w, Aw, wAw
        w = nw
        if not any(w):
            break
    return x


def ic_solve(inst, D1, D2, k_true, rng, kappa=2.0, variant="2lp", r=20, log=None, max_steps=None):
    """return dict with k, correctness, counters breakdown."""
    p, h, ell = inst.p, inst.h, inst.ell
    t0 = time.time()
    c_before = snap()
    # factor base: degree-1 primes with smallest x, h(x) nonzero square
    nfb_target = math.ceil(kappa * p ** (2.0 / 3.0)) if variant == "2lp" else None
    fb = {}
    ycan = {}
    x = 0
    while x < p:
        hx = pev(h, x, p)
        CNT[0] += int(1.5 * p.bit_length())
        if hx and pow(hx, (p - 1) // 2, p) == 1:
            if nfb_target is None or len(fb) < nfb_target:
                fb[x] = len(fb)
            elif variant == "2lp":
                break
        x += 1
    nfb = len(fb)
    kcol = nfb
    use_lp = variant == "2lp"
    # walk setup
    a_s, b_s = rng.randrange(1, ell), rng.randrange(1, ell)
    M0 = cmul(a_s, D1, h, p)
    M1 = cmul(b_s, D1, h, p)
    Ms, Ma = [M0, M1], [a_s, b_s]
    for i in range(2, r):
        Ms.append(cadd(Ms[-1], Ms[-2], h, p))
        Ma.append((Ma[-1] + Ma[-2]) % ell)
    al, be = rng.randrange(ell), rng.randrange(1, ell)
    R = cadd(cmul(al, D1, h, p), cmul(be, D2, h, p), h, p)
    c_setup = snap()
    seen = set()
    full_rows = []   # list of (dict col->coef, rhs)
    expr = {}        # LP x -> (dict col->coef, const)  meaning lambda_L = sum coef*x_col + const
    pend = {}        # LP x -> list of (other x, rel)
    stats = dict(steps=0, cand_split=0, full=0, lp1=0, lp2=0, lp3=0, restarts=0, cycles=0, weier=0)
    next_check = nfb + 10
    kfound = None
    last_new = 0
    checks = 0

    def emit(rowd, rhs):
        rowd = {j: c % ell for j, c in rowd.items() if c % ell}
        full_rows.append((rowd, rhs % ell))

    def sub_expr(rowd, rhs, s, e):
        # rowd: sum + s*lambda_L = rhs -> substitute lambda_L = e
        d = dict(rowd)
        for j, c in e[0].items():
            d[j] = (d.get(j, 0) + s * c) % ell
        return d, (rhs - s * e[1]) % ell

    def solve_for(rowd, rhs, s):
        # rowd + s*lambda = rhs -> lambda = (rhs - rowd)/s, s = +-1
        si = s % ell if s == 1 else ell - 1
        return ({j: (-c * si) % ell for j, c in rowd.items()}, rhs * si % ell)

    def propagate(start):
        stack = [start]
        while stack:
            v = stack.pop()
            for (w, rel) in pend.pop(v, []):
                # rel = (rowd, rhs, [(xa, sa), (xb, sb)])
                rowd, rhs, lps = rel
                if rel[0] is None:
                    continue
                d, rr = rowd, rhs
                for (xl, sl) in lps:
                    if xl == v:
                        d, rr = sub_expr(d, rr, sl, expr[v])
                sw = [sl for (xl, sl) in lps if xl == w][0]
                if w in expr:
                    d2, r2 = sub_expr(d, rr, sw, expr[w])
                    emit(d2, r2)
                    stats["cycles"] += 1
                else:
                    expr[w] = solve_for(d, rr, sw)
                    stack.append(w)
                # remove the reverse entry
                lst = pend.get(w)
                if lst is not None:
                    pend[w] = [t for t in lst if t[1] is not rel]

    while kfound is None:
        stats["steps"] += 1
        if max_steps and stats["steps"] > max_steps:
            break
        u, v = R
        key = None
        if splits(u, p):
            key = hash((tuple(u), tuple(v)))
            if key in seen:
                stats["restarts"] += 1
                R = cadd(cadd(R, Ms[0], h, p), D2, h, p)
                al, be = (al + Ma[0]) % ell, (be + 1) % ell
                last_new = stats["steps"]
                continue
            seen.add(key)
            last_new = stats["steps"]
            stats["cand_split"] += 1
            rts = roots_split(u, p, rng)
            rowd = {kcol: (-be) % ell}
            lps = []
            ok = True
            for x0 in rts:
                y0 = pev(v, x0, p)
                if y0 == 0:
                    stats["weier"] += 1
                    continue
                s = 1 if y0 <= p - y0 else -1
                if x0 in fb:
                    j = fb[x0]
                    rowd[j] = (rowd.get(j, 0) + s) % ell
                else:
                    if not use_lp:
                        ok = False
                        break
                    lps.append((x0, s))
            if ok:
                nl = len(lps)
                if nl == 0:
                    stats["full"] += 1
                    emit(rowd, al)
                elif nl == 1:
                    stats["lp1"] += 1
                    xl, sl = lps[0]
                    if xl in expr:
                        d, rr = sub_expr(rowd, al, sl, expr[xl])
                        emit(d, rr)
                        stats["cycles"] += 1
                    else:
                        expr[xl] = solve_for(rowd, al, sl)
                        propagate(xl)
                elif nl == 2:
                    stats["lp2"] += 1
                    (xa, sa), (xb, sb) = lps
                    ina, inb = xa in expr, xb in expr
                    if ina and inb:
                        d, rr = sub_expr(rowd, al, sa, expr[xa])
                        d, rr = sub_expr(d, rr, sb, expr[xb])
                        emit(d, rr)
                        stats["cycles"] += 1
                    elif ina or inb:
                        if inb:
                            (xa, sa), (xb, sb) = (xb, sb), (xa, sa)
                        d, rr = sub_expr(rowd, al, sa, expr[xa])
                        expr[xb] = solve_for(d, rr, sb)
                        propagate(xb)
                    else:
                        rel = (rowd, al, lps)
                        pend.setdefault(xa, []).append((xb, rel))
                        pend.setdefault(xb, []).append((xa, rel))
                else:
                    stats["lp3"] += 1
        if stats["steps"] - last_new > 60:
            stats["restarts"] += 1
            R = cadd(cadd(R, Ms[0], h, p), D2, h, p)
            al, be = (al + Ma[0]) % ell, (be + 1) % ell
            last_new = stats["steps"]
            continue
        if len(full_rows) >= next_check:
            # filter: drop singleton columns (never the k column)
            rows = [dict(rw) for rw, _ in full_rows]
            rhs = [rr for _, rr in full_rows]
            alive = [True] * len(rows)
            changed = True
            while changed:
                changed = False
                cnt = {}
                for i, rw in enumerate(rows):
                    if alive[i]:
                        for j in rw:
                            cnt[j] = cnt.get(j, 0) + 1
                for i, rw in enumerate(rows):
                    if alive[i] and any(cnt[j] == 1 and j != kcol for j in rw):
                        alive[i] = False
                        changed = True
            cols = {}
            for i, rw in enumerate(rows):
                if alive[i]:
                    for j in rw:
                        if j not in cols:
                            cols[j] = len(cols)
            nrows = sum(alive)
            if kcol in cols and nrows >= len(cols) + 5:
                checks += 1
                srows = [[(cols[j], c) for j, c in rows[i].items()] for i in range(len(rows)) if alive[i]]
                srhs = [rhs[i] for i in range(len(rows)) if alive[i]]
                c_la0 = snap()
                sol = lanczos(srows, srhs, len(cols), ell)
                kk = sol[cols[kcol]]
                # verification (not counted): k*D1 == D2
                old = CNT[:]
                okk = cmul(kk, D1, h, p) == D2
                CNT[:] = old
                if okk:
                    kfound = kk
                    break
            next_check = int(len(full_rows) * 1.1) + 5
        if kfound is not None:
            break
        # next walk step
        j = (u[0] * 31 + (u[1] if len(u) > 1 else 0) * 7 + (v[0] if v else 0)) % r
        R = cadd(R, Ms[j], h, p)
        al = (al + Ma[j]) % ell
    c_end = snap()
    dt = time.time() - t0
    return dict(k=kfound, correct=(kfound == k_true % ell), nfb=nfb, rows=len(full_rows),
                lp_vertices=len(expr), checks=checks, wall=dt,
                counters=[c_end[i] - c_before[i] for i in range(5)],
                setup_counters=[c_setup[i] - c_before[i] for i in range(5)], **stats)


# ------------------------------------------------------------------ driver
def next_prime(n):
    while not is_prime(n):
        n += 1
    return n


def cost_units(c, p, cl=1):
    """F_q-mult equivalents: M + ceil(log2 q) I + cl*M_ell + cl*log2(ell) I_ell (ell inv by Fermat)."""
    return c[0] + math.ceil(math.log2(p)) * c[1] + cl * (c[2]) + cl * c[3] * 0


def run(bits, curve, arms, rho_targets, kappa, out, basic_max_bits=99):
    p = next_prime((1 << bits) + 1)
    seed = 20261014 * 1000 + bits * 10 + curve
    rng = random.Random(seed)
    t0 = time.time()
    tries = 0
    while True:
        tries += 1
        inst = build_instance(p, rng)
        if inst is not None:
            break
    rec = dict(bits=bits, q=p, curve=curve, seed=seed, build_tries=tries, N=inst.N, ell=inst.ell, c=inst.c,
               h=inst.h, a=inst.a, b=inst.b, r0=list(inst.r[0]), Kmod=[inst.K.c1, inst.K.c0],
               build_wall=time.time() - t0)
    E = inst.E
    # generator of E0(K)[ell]
    while True:
        P = random_point(inst, rng)
        old = CNT[:]
        G = E.mul(inst.c, P)
        CNT[:] = old
        if G is not None:
            break
    old = CNT[:]
    assert E.on(G) and E.mul(inst.ell, G) is None
    # checks (not counted)
    phiG = transfer(inst, G)
    chk = dict(phiG_nonzero=phiG != ([1], []),
               ell_phiG_zero=cmul(inst.ell, phiG, inst.h, p) == ([1], []),
               N_phiG_zero=cmul(inst.N, phiG, inst.h, p) == ([1], []))
    hom = 0
    for _ in range(3):
        P1, P2 = random_point(inst, rng), random_point(inst, rng)
        P3 = E.add(P1, P2)
        try:
            hom += transfer(inst, P3) == cadd(transfer(inst, P1), transfer(inst, P2), inst.h, p)
        except ValueError:
            pass
    chk["hom_ok"] = hom
    if bits <= 11:
        chk["jac_order_hyperellcharpoly_equals_N"] = gp_order(inst.h, p) == inst.N
    deg_h = len(inst.h) - 1
    chk["deg_h"] = deg_h
    CNT[:] = old
    rec["checks"] = chk
    results = []
    if "ic2lp" in arms or "basic" in arms:
        for arm in [a for a in arms if a in ("ic2lp", "basic")]:
            if arm == "basic" and bits > basic_max_bits:
                continue
            k = rng.randrange(1, inst.ell)
            old = CNT[:]
            T = E.mul(k, G)
            CNT[:] = old
            reset()
            tt = time.time()
            # transfer cost counted (both points)
            D1 = transfer(inst, G)
            D2 = transfer(inst, T)
            ctr = snap()
            res = ic_solve(inst, D1, D2, k, rng, kappa=kappa, variant="2lp" if arm == "ic2lp" else "basic")
            res["transfer_counters"] = ctr
            tot = [ctr[i] + res["counters"][i] for i in range(5)]
            res["total_counters"] = tot
            res["cost"] = cost_units(tot, p)
            res["cost_cl9"] = cost_units(tot, p, 9)
            res["arm"] = arm
            res["wall_total"] = time.time() - tt
            results.append(res)
            print(json.dumps(dict(bits=bits, curve=curve, arm=arm, correct=res["correct"], cost=res["cost"],
                                  steps=res["steps"], wall=round(res["wall_total"], 1))), file=sys.stderr, flush=True)
    if "rho" in arms:
        for t in range(rho_targets):
            k = rng.randrange(1, inst.ell)
            old = CNT[:]
            T = E.mul(k, G)
            CNT[:] = old
            reset()
            tt = time.time()
            kk, steps = rho(inst, G, T, rng, max_steps=40 * math.isqrt(inst.ell) + 1000)
            c = snap()
            res = dict(arm="rho", target=t, k=kk, correct=(kk == k), steps=steps, counters=c,
                       cost=cost_units(c, p), cost_cl9=cost_units(c, p, 9), wall_total=time.time() - tt)
            results.append(res)
            print(json.dumps(dict(bits=bits, curve=curve, arm="rho", correct=res["correct"], cost=res["cost"],
                                  steps=steps, wall=round(res["wall_total"], 1))), file=sys.stderr, flush=True)
    rec["results"] = results
    with open(out, "a") as f:
        f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, required=True)
    ap.add_argument("--curve", type=int, required=True)
    ap.add_argument("--arms", default="ic2lp,rho")
    ap.add_argument("--rho-targets", type=int, default=3)
    ap.add_argument("--kappa", type=float, default=2.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    run(a.bits, a.curve, a.arms.split(","), a.rho_targets, a.kappa, a.out)
