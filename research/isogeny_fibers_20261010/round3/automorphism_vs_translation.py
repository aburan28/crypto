"""Automorphism orbits versus isogeny-kernel translations, on one curve.

A j = 0 curve E: y^2 = x^3 + b over F_p, p = 1 mod 3, with N = 3 r
(r prime), carries two order-3 actions on its points:
  * the automorphism zeta: (x, y) -> (w x, y), w^3 = 1   (fixes O, an
    isogeny of degree 1, acts on G as a scalar);
  * translation by a rational 3-torsion point T          (the kernel of the
    3-isogeny E -> E / <T>; a fiber of that isogeny is P + <T>).
Same group order, one acts as an automorphism, one as a kernel translation.

Bases of equal size (about 72 signed points):
  plain       small-x points
  zeta-closed S u zeta S u zeta^2 S       (S small-x, |S| = |F| / 3)
  T-closed    S u (S + T) u (S + 2T)      (fiber-saturated, mu = 3)
  G           random points of the prime-order subgroup (image of [3])
For each: exact D_m = number of distinct points of G that are sums of m
signed base points (m = 2, 3); columns U = orbits of the projections [h]P
under the automorphism group <-1, zeta> (the folding any index calculus
applies); class factor c = D_m * m! / |F|^m; mitm cost model
U * (r / D_m) * L_m with L_2 = 1, L_3 = |F|.

Prediction (RESEARCH_ISOGENY_FIBERS_ECDLP.md section 15.2): relative to
plain, zeta-closed keeps c and divides U by 3; T-closed divides U by 3
and c by 3^(m-1); G multiplies c by h = 3.
"""
import json, random, sys
from itertools import combinations_with_replacement

def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

def sqrt_mod(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1: tt = tt * tt % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r

class Curve:
    def __init__(self, p, b):
        self.p, self.b = p, b
    def on(self, P):
        return P is None or (P[1] * P[1] - P[0] ** 3 - self.b) % self.p == 0
    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.p)
    def add(self, P, Q):
        p = self.p
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % p == 0: return None
            lam = 3 * P[0] * P[0] * pow(2 * P[1], -1, p) % p
        else:
            lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
        x3 = (lam * lam - P[0] - Q[0]) % p
        return (x3, (lam * (P[0] - x3) - P[1]) % p)
    def smul(self, k, P):
        R = None
        while k:
            if k & 1: R = self.add(R, P)
            P = self.add(P, P); k >>= 1
        return R
    def lift(self, x):
        s = sqrt_mod(x ** 3 + self.b, self.p)
        if s is None: return []
        return [(x, s)] if s == 0 else [(x, s), (x, self.p - s)]
    def order(self):
        p = self.p
        return p + 1 + sum(1 if pow((x ** 3 + self.b) % p, (p - 1) // 2, p) == 1 else
                           (-1 if (x ** 3 + self.b) % p else 0) for x in range(p))

def find_curve(bits, seed):
    rng = random.Random(seed)
    while True:
        p = rng.randrange(1 << (bits - 1), 1 << bits)
        if p % 3 != 1 or not is_prime(p): continue
        for b in range(1, 12):
            E = Curve(p, b)
            N = E.order()
            if N % 3 == 0 and N % 9 != 0 and is_prime(N // 3):
                return E, N

def census(E, base, r, h, w):
    p = E.p
    zeta = lambda P: None if P is None else (w * P[0] % p, P[1])
    cls = [E.smul(r, P) for P in base]
    # columns: orbits of [h]P under <-1, zeta>
    seen, U = set(), 0
    for P in base:
        Q = E.smul(h, P)
        if Q is None or Q in seen: continue
        U += 1
        R = Q
        for _ in range(3):
            seen.add(R); seen.add(E.neg(R)); R = zeta(R)
    out = dict(points=len(base), columns=U, classes=len(set(cls)))
    pairs = {}
    for i, j in combinations_with_replacement(range(len(base)), 2):
        s = E.add(base[i], base[j])
        k = E.add(cls[i], cls[j])
        pairs.setdefault(k, set()).add(s)
    for m in (2, 3):
        if m == 2:
            hits = {s for s in pairs.get(None, set()) if s is not None}
        else:
            hits = set()
            for i, P in enumerate(base):
                need = E.neg(cls[i])
                for s in pairs.get(need, ()):
                    if s is None: continue
                    t = E.add(s, P)
                    if t is not None: hits.add(t)
        D = len(hits)
        fact = 2 if m == 2 else 6
        L = 1 if m == 2 else len(base)
        out[f"m{m}"] = dict(D=D, hit_rate=D / r, c=D * fact / len(base) ** m,
                            cost_model=U * (r / D) * L if D else None)
    return out

def main(bits=18, seed=5, size=72):
    E, N = find_curve(bits, seed)
    p, r, h = E.p, N // 3, 3
    g = 2
    while pow(g, (p - 1) // 3, p) == 1: g += 1
    w = pow(g, (p - 1) // 3, p)
    rng = random.Random(seed)
    # rational 3-torsion point
    T = None
    x = 1
    while T is None:
        for P in E.lift(x):
            Q = E.smul(r, P)
            if Q is not None: T = Q; break
        x += 1
    assert E.smul(3, T) is None
    zeta = lambda P: (w * P[0] % p, P[1])
    def small(k):
        pts, x = [], 1
        while len(pts) < k:
            pts += E.lift(x); x += 1
        return pts[:k] if len(pts[:k]) % 2 == 0 else pts[:k + 1]
    bases = {}
    bases["plain (small x)"] = small(size)
    S = small(size // 3)
    bases["zeta-closed (automorphism orbits, a = 3)"] = S + [zeta(P) for P in S] + [zeta(zeta(P)) for P in S]
    bases["T-closed (3-isogeny fibers, mu = 3)"] = S + [E.add(P, T) for P in S] + [E.add(P, E.add(T, T)) for P in S]
    G = set()
    while len(G) < size:
        L = E.lift(rng.randrange(1, p))
        if not L: continue
        Q = E.smul(h, L[0])
        if Q is None: continue
        G.add(Q); G.add(E.neg(Q))
    bases["G (inside the image of [3], i = 3)"] = sorted(G)
    for name, B in bases.items():
        assert len(set(B)) == len(B), name
        assert all(E.on(P) for P in B)
    res = {name: census(E, B, r, h, w) for name, B in bases.items()}
    return dict(p=p, b=E.b, N=N, r=r, h=h, omega=w, T=T), res

if __name__ == "__main__":
    bits = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    info, res = main(bits, seed)
    print(json.dumps(dict(curve=info, bases=res), indent=1))
