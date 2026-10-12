"""Linear (A, P) oracle for the S3 core over a geometric factor base.

Curve y^2 + xy = x^3 + 1 over F_{2^n}: S3(x1, x2, x3) = e2^2 + e3 + 1.  Fix
x1 = t and put A = X2 + X3, P = X2 X3.  Then e2 = tA + P, e3 = tP and

    S3(t, X2, X3) = P^2 + t^2 A^2 + t P + 1,

which is F_2-linear in the bits of (A, P) because squaring is F_2-linear.
X2 and X3 are the two roots of Z^2 + A Z + P.

For X2, X3 in V, A lies in V (l bits) and P lies in span(V * V).
  geometric V = theta <1, g, ..., g^(l-1)>: span(V*V) <= theta^2 <1..g^(2l-2)>, 2l-1 bits
  random V:                                 span(V*V) = span{v_i v_j}, l(l+1)/2 bits
One trial = one n x (l + dim span(V*V)) linear solve for a target abscissa t.
Accept a solution only if both roots of Z^2 + A Z + P lie in V, lift to curve
points, and some +-P2 +-P3 equals the target point T.
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "linearization_reach_20260930"))
import lr  # noqa: E402  (frozen field, curve and quadratic code, reused unchanged)

MAX_FAMILY_BITS = 6   # enumerate at most 2^6 points of an underdetermined solution family


def rank_of(vectors):
    piv = {}
    for v in vectors:
        while v:
            low = (v & -v).bit_length() - 1
            if low in piv:
                v ^= piv[low]
            else:
                piv[low] = v
                break
    return len(piv)


class Base:
    """Factor base V and a basis of span(V*V); membership test for V."""

    def __init__(self, F, l, kind, rng):
        n = F.n
        while True:
            if kind == "geometric":
                theta, g = rng.getrandbits(n) or 1, rng.getrandbits(n) or 2
                V = [F.mul(theta, self._pw(F, g, i)) for i in range(l)]
                prods = [F.mul(F.sq(theta), self._pw(F, g, k)) for k in range(2 * l - 1)]
            else:
                V = [rng.getrandbits(n) for _ in range(l)]
                prods = [F.mul(V[i], V[j]) for i in range(l) for j in range(i, l)]
            if rank_of(V) == l:
                break
        self.F, self.l, self.V = F, l, V
        # reduce the product spanning set to a basis
        basis, piv = [], {}
        for p in prods:
            v = p
            while v:
                low = (v & -v).bit_length() - 1
                if low in piv:
                    v ^= piv[low]
                else:
                    piv[low] = v
                    basis.append(p)
                    break
        self.P = basis
        # membership in V: echelon form of V
        self.vpiv = {}
        for v in V:
            w = v
            while w:
                low = (w & -w).bit_length() - 1
                if low in self.vpiv:
                    w ^= self.vpiv[low]
                else:
                    self.vpiv[low] = w
                    break

    @staticmethod
    def _pw(F, g, e):
        r = 1
        for _ in range(e):
            r = F.mul(r, g)
        return r

    def in_V(self, x):
        while x:
            low = (x & -x).bit_length() - 1
            if low not in self.vpiv:
                return False
            x ^= self.vpiv[low]
        return True

    def element(self, bits):
        x = 0
        for v, b in zip(self.V, bits):
            if b:
                x ^= v
        return x

    @property
    def unknowns(self):
        return self.l + len(self.P)


def solve(base, t):
    """All (A, P) of the linear system for target abscissa t; None if too many."""
    F, n = base.F, base.F.n
    t2 = F.sq(t)
    cols = [F.mul(t2, F.sq(v)) for v in base.V]                      # A coordinates
    cols += [F.sq(p) ^ F.mul(t, p) for p in base.P]                   # P coordinates
    N = len(cols)
    rows = []
    for bit in range(n):
        r = 0
        for k, cv in enumerate(cols):
            if (cv >> bit) & 1:
                r |= 1 << k
        if bit == 0:          # right-hand side: the constant 1
            r |= 1 << N
        rows.append(r)
    rk, cons, piv = lr.eliminate(rows, N)
    if not cons:
        return rk, []
    free = [c for c in range(N) if c not in piv]
    if len(free) > MAX_FAMILY_BITS:
        return rk, None
    out = []
    for mask in range(1 << len(free)):
        sol = [0] * N
        for k, c in enumerate(free):
            sol[c] = (mask >> k) & 1
        for c in sorted(piv, reverse=True):
            r, v = piv[c], (piv[c] >> N) & 1
            for c2 in range(c + 1, N):
                if (r >> c2) & 1:
                    v ^= sol[c2]
            sol[c] = v
        A = 0
        for v, b in zip(base.V, sol[:base.l]):
            if b:
                A ^= v
        P = 0
        for p, b in zip(base.P, sol[base.l:]):
            if b:
                P ^= p
        out.append((A, P))
    return rk, out


def decompositions(base, T, cand):
    """Verified pairs (X2, X3) from candidate (A, P); counts fakes."""
    F = base.F
    good, fake = [], 0
    for A, P in cand:
        if A == 0:
            fake += 1           # X2 = X3: not a two-point decomposition
            continue
        roots = lr.quad_roots(F, 1, A, P)
        if len(roots) != 2 or not all(base.in_V(z) for z in roots):
            fake += 1
            continue
        pts = [lr.lift(F, z) for z in roots]
        if any(p is None for p in pts):
            fake += 1
            continue
        hit = False
        for s in range(4):
            q2 = pts[0] if s & 1 else lr.neg(pts[0])
            q3 = pts[1] if s & 2 else lr.neg(pts[1])
            S = lr.add(F, q2, q3)
            if S is not None and (S == T or S == lr.neg(T)):
                hit = True
        if hit:
            good.append(tuple(roots))
        else:
            fake += 1
    return good, fake


def random_target(F, rng):
    while True:
        T = lr.lift(F, rng.getrandbits(F.n))
        if T:
            return T


def planted_target(F, base, rng):
    while True:
        x2 = base.element([rng.getrandbits(1) for _ in range(base.l)])
        x3 = base.element([rng.getrandbits(1) for _ in range(base.l)])
        if x2 == x3:
            continue
        p2, p3 = lr.lift(F, x2), lr.lift(F, x3)
        if p2 is None or p3 is None:
            continue
        T = lr.add(F, p2, p3)
        if T is not None:
            return T, {x2, x3}


def cell(n, l, kind, trials, planted, rng):
    F = lr.Field(n)
    res = dict(n=n, l=l, kind=kind, trials=trials, success=0, fake=0, too_big=0,
               rank_hist={}, unknowns=None, planted=planted, planted_recovered=0,
               planted_too_big=0)
    for _ in range(trials):
        base = Base(F, l, kind, rng)
        res["unknowns"] = base.unknowns
        T = random_target(F, rng)
        rk, cand = solve(base, T[0])
        res["rank_hist"][str(rk)] = res["rank_hist"].get(str(rk), 0) + 1
        if cand is None:
            res["too_big"] += 1
            continue
        good, fake = decompositions(base, T, cand)
        res["success"] += bool(good)
        res["fake"] += fake
    for _ in range(planted):
        base = Base(F, l, kind, rng)
        T, pair = planted_target(F, base, rng)
        rk, cand = solve(base, T[0])
        if cand is None:
            res["planted_too_big"] += 1
            continue
        good, _ = decompositions(base, T, cand)
        res["planted_recovered"] += any(set(g) == pair for g in good)
    res["log2p"] = math.log2(res["success"] / trials) if res["success"] else None
    res["law"] = 2 * l - n - 1
    return res


if __name__ == "__main__":
    n, l, T, Pl, seed = map(int, sys.argv[1:6])
    kind = sys.argv[6]
    print(json.dumps(cell(n, l, kind, T, Pl, random.Random(seed))))
