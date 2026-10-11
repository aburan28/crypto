"""Image-contained versus fiber-saturated factor bases: an exact census.

Claim under test (RESEARCH_ISOGENY_FIBERS_ECDLP.md section 14.1).  For an
m-term decomposition over a base F, the number D_m of distinct points of the
prime-order subgroup G that are sums of m signed base points scales as
    D_m ~ c * |F|^m / m!,
and the class factor c depends on how F sits relative to the isogeny
structure of E(F_q):
  * a base inside the image of an isogeny of index i (a subgroup of index
    i in the h-part) has c multiplied by i (the cofactor classes cancel
    more often);
  * a base saturated under the kernel M of an isogeny (closed under
    P -> P + M, |M| = mu) has c divided by mu^(m-1) (each sum is
    represented mu^(m-1) times) while its projected column count is
    divided by mu.
With the mitm pipeline cost  C_m = U * (r / D_m) * L_m  (U columns,
L_2 = 1 lookup per target, L_3 = |F| lookups) this predicts: image
containment gains i at every m; kernel saturation is neutral at m = 2 and
loses mu^(m-2) at m >= 3.

Curve: Koblitz K_0: y^2 + x y = x^3 + 1 over GF(2^n), small n, exact
enumeration of every relevant sum.  Pure Python, no dependencies.
"""
import json, random, sys
from itertools import combinations_with_replacement

# ------------------------------------------------------------------ GF(2^n)
def poly_mod(a, f, n):
    while a.bit_length() > n:
        a ^= f << (a.bit_length() - 1 - n)
    return a

def gmul(a, b, f, n):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> n:
            a ^= f
    return r

def ginv(a, f, n):
    # extended Euclid over GF(2)[x]
    u, v, g1, g2 = a, f, 1, 0
    while u != 1:
        j = u.bit_length() - v.bit_length()
        if j < 0:
            u, v, g1, g2, j = v, u, g2, g1, -j
        u ^= v << j
        g1 ^= g2 << j
    return poly_mod(g1, f, n)

def gsq(a, f, n):
    return gmul(a, a, f, n)

def trace(a, f, n):
    t, c = 0, a
    for _ in range(n):
        t ^= c
        c = gsq(c, f, n)
    return t & 1 if t in (0, 1) else None   # trace lies in GF(2)

def half_trace(c, f, n):
    h, s = c, c
    for _ in range((n - 1) // 2):
        s = gsq(gsq(s, f, n), f, n)
        h ^= s
    return h

def irreducible(f, n):
    # n prime: irreducible iff x^(2^n) = x mod f and f has no root in GF(2)
    if not (f & 1) or bin(f).count("1") % 2 == 0:
        return False
    x = 2
    for _ in range(n):
        x = gsq(x, f, n)
    return x == 2

def find_modulus(n):
    for k in range(1, n):
        f = (1 << n) | (1 << k) | 1
        if irreducible(f, n):
            return f
    for a in range(1, n):
        for b in range(1, a):
            for c in range(1, b):
                f = (1 << n) | (1 << a) | (1 << b) | (1 << c) | 1
                if irreducible(f, n):
                    return f
    raise ValueError(n)

# ------------------------------------------------------------------ curve
class K0:
    """y^2 + x y = x^3 + a x^2 + 1 over GF(2^n)."""
    def __init__(self, n, a=0):
        self.n, self.a = n, a
        self.f = find_modulus(n)
        t = -1 if a == 0 else 1
        v0, v1 = 2, t
        for _ in range(n - 1):
            v0, v1 = v1, t * v1 - 2 * v0
        self.N = (1 << n) + 1 - v1

    def mul(self, a, b): return gmul(a, b, self.f, self.n)
    def inv(self, a): return ginv(a, self.f, self.n)

    def lift(self, x):
        """Both points over x (x != 0), or []."""
        if x == 0:
            return [(0, 1)]
        xi = self.inv(x)
        c = x ^ (self.a) ^ self.mul(xi, xi)        # z^2 + z = x + a + 1/x^2
        z = half_trace(c, self.f, self.n)
        if self.mul(z, z) ^ z != c:
            return []
        y = self.mul(x, z)
        return [(x, y), (x, y ^ x)]

    def neg(self, P):
        return None if P is None else (P[0], P[0] ^ P[1])

    def add(self, P, Q):
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2:
            if y1 ^ y2 == x1 or x1 == 0 and y1 != y2:   # Q = -P
                return None
            if x1 == 0:
                return None
            lam = x1 ^ self.mul(y1, self.inv(x1))
            x3 = self.mul(lam, lam) ^ lam ^ self.a
            y3 = self.mul(x1, x1) ^ self.mul(lam ^ 1, x3)
            return (x3, y3)
        lam = self.mul(y1 ^ y2, self.inv(x1 ^ x2))
        x3 = self.mul(lam, lam) ^ lam ^ x1 ^ x2 ^ self.a
        y3 = self.mul(lam, x1 ^ x3) ^ x3 ^ y1
        return (x3, y3)

    def smul(self, k, P):
        R = None
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R

def factor(N):
    out, d = {}, 2
    while d * d <= N:
        while N % d == 0:
            out[d] = out.get(d, 0) + 1
            N //= d
        d += 1
    if N > 1:
        out[N] = out.get(N, 0) + 1
    return out

# ------------------------------------------------------------------ census
def census(E, base, r, h, m_list=(2, 3)):
    """Exact D_m (distinct points of G that are sums of m signed base points),
    projected column count U, and the mitm cost model, for one base."""
    # class of each point modulo G, as the point [r]P in E[h]
    cls = [E.smul(r, P) for P in base]
    keys = {}
    cid = [keys.setdefault(c, len(keys)) for c in cls]
    # class group arithmetic via the points themselves
    reps = list(keys)
    add_cls = {}
    for i, a in enumerate(reps):
        for j, b in enumerate(reps):
            add_cls[(i, j)] = keys.setdefault(E.add(a, b), len(keys)) if E.add(a, b) in keys else None
    zero = keys.get(None)
    # projected columns: [h]P up to sign
    cols = set()
    for P in base:
        Q = E.smul(h, P)
        if Q is None:
            continue
        cols.add(min(Q, E.neg(Q)))
    U = len(cols)
    out = dict(points=len(base), columns=U, points_per_column=len(base) / max(U, 1),
               classes=len(set(cid)))
    # pair sums, bucketed by class
    pairs = {}
    for i, j in combinations_with_replacement(range(len(base)), 2):
        s = E.add(base[i], base[j])
        k = add_cls.get((cid[i], cid[j]))
        if k is None:
            k = keys.setdefault(E.smul(r, s) if s is not None else None, len(keys))
        pairs.setdefault(k, set()).add(s)
    zero = keys.get(None)
    for m in m_list:
        if m == 2:
            hits = {s for s in pairs.get(zero, set()) if s is not None}
        else:
            hits = set()
            for i, P in enumerate(base):
                # need class(pair) = -class(P): find pair classes k with k + cid[i] = zero
                for k, ss in pairs.items():
                    rep_pair = next(iter(ss))
                    c_pair = E.smul(r, rep_pair) if rep_pair is not None else None
                    if E.add(c_pair, cls[i]) is not None:
                        continue
                    for s in ss:
                        if s is None:          # P + Q - Q: degenerate, not a 3-term sum
                            continue
                        t = E.add(s, P)
                        if t is not None:
                            hits.add(t)
        D = len(hits)
        Fm = len(base) ** m
        fact = 2 if m == 2 else 6
        c = D * fact / Fm
        L = 1 if m == 2 else len(base)
        cost = U * (r / D) * L if D else float("inf")
        out[f"m{m}"] = dict(D=D, hit_rate=D / r, c=c, cost_model=cost)
    return out

def main(n=19, d=6, seed=1):
    E = K0(n, 0)
    fac = factor(E.N)
    r = max(fac)
    h = E.N // r
    rng = random.Random(seed)
    two_part = 1
    while E.N % (two_part * 2) == 0:
        two_part *= 2
    odd = E.N // two_part

    def in_2jE(P, j):
        P2 = E.smul(odd, P)                 # 2-primary component
        Q = P2
        for _ in range((two_part.bit_length() - 1) - j):
            Q = E.add(Q, Q)
        return Q is None

    def xbase(xs):
        pts = []
        for x in xs:
            if x == 0:
                continue
            pts += E.lift(x)
        return pts

    def span(gens):
        S = {0}
        for g in gens:
            S |= {s ^ g for s in S}
        return S

    V = span([1 << i for i in range(d)])                      # contains 1
    tr = lambda x: trace(x, E.f, E.n)
    Vtz = {x for x in span([1 << i for i in range(d + 1)]) if tr(x) == 0}      # dim d
    Vtz2 = {x for x in span([1 << i for i in range(d + 2)]) if tr(x) == 0}     # dim d+1

    bases = {}
    bases["x (plain subspace, i = 1)"] = xbase(V)
    bases["tz (trace zero = inside 2E, i = 2)"] = xbase(Vtz)
    four = [P for P in xbase(Vtz2) if in_2jE(P, 2)]
    bases["4E (inside 4E, i = 4)"] = four
    # u-frame: closed under +T, T = (0, 1)
    u_x = []
    for u in V:
        if u in (0, 1):
            continue
        u_x.append(E.inv(u) ^ 1)
    bases["u (+T saturated, mu = 2)"] = xbase(u_x)
    # subgroup base: random points of G, same size as the x base
    target = len(bases["x (plain subspace, i = 1)"])
    G = set()
    while len(G) < target:
        x = rng.randrange(1, 1 << n)
        L = E.lift(x)
        if not L:
            continue
        Q = E.smul(h, L[0])
        if Q is None:
            continue
        G.add(Q); G.add(E.neg(Q))
    bases["G (inside the prime-order subgroup, i = h)"] = sorted(G)

    info = dict(n=n, N=E.N, factorisation={str(k): v for k, v in fac.items()}, r=r, h=h,
                two_part=two_part, modulus=hex(E.f), subspace_dim=d)
    # sanity: group order
    for P in bases["x (plain subspace, i = 1)"][:4]:
        assert E.smul(E.N, P) is None
    res = {}
    for name, B in bases.items():
        res[name] = census(E, B, r, h)
    return info, res

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 19
    d = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    info, res = main(n, d)
    print(json.dumps(dict(curve=info, bases=res), indent=1))
