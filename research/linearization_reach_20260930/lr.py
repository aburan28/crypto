"""Linearization reach of the bilinear S3 core in binary point decomposition.

Curve E: y^2 + xy = x^3 + 1 over F_{2^n} (a2 = 0, a6 = 1, the ECC2K-130 family).
S3(x1, x2, x3) = (x1 x2 + x1 x3 + x2 x3)^2 + x1 x2 x3 + 1.

With x3 = t fixed, S3(X2, X3, t) = 0 is F_2-bilinear in the coordinates of X2 and
X3, because squaring is F_2-linear.  A trial of the guess-and-linearize oracle:

  1. guess X1 in V, solve S3(X1, t, xR) = 0 for t (a quadratic, half-trace);
  2. guess a coset w0 + W of a b-dimensional W < V for X3;
  3. linearize S3(X2, X3, t) = 0 over the monomials {c_i d_j, c_i, d_j}, with
     X2 = sum c_i v_i in V and X3 = w0 + sum d_j w_j;
  4. enumerate the affine solution family, keep product-consistent points,
     rebuild the points and check +-P1 +-P2 +-P3 = R on the curve.

Pure Python 3, no dependencies.  See PROTOCOL.md for the frozen predictions.
"""
import json
import math
import random
import sys

# Irreducible polynomials; checked at start-up (x^(2^n) = x, n prime).
IRR = {13: (1 << 13) | 0b11011, 17: (1 << 17) | 0b1001, 19: (1 << 19) | 0b100111,
       23: (1 << 23) | 0b100001, 29: (1 << 29) | 0b101, 31: (1 << 31) | 0b1001}


class Field:
    def __init__(self, n):
        self.n, self.P = n, IRR[n]
        x = 2
        for _ in range(n):
            x = self.mul(x, x)
        assert x == 2, f"modulus for n={n} is not irreducible"

    def mul(self, a, b):
        n, P, r = self.n, self.P, 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if (a >> n) & 1:
                a ^= P
        return r

    def sq(self, a):
        return self.mul(a, a)

    def inv(self, a):
        r, e = 1, (1 << self.n) - 2
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            e >>= 1
        return r

    def tr(self, a):
        s = 0
        for _ in range(self.n):
            s ^= a
            a = self.sq(a)
        return s

    def htr(self, a):  # n odd: z^2 + z = a whenever Tr(a) = 0
        s = 0
        for _ in range((self.n + 1) // 2):
            s ^= a
            a = self.sq(self.sq(a))
        return s


def S3(F, x1, x2, x3):
    s = F.mul(x1, x2) ^ F.mul(x1, x3) ^ F.mul(x2, x3)
    return F.sq(s) ^ F.mul(F.mul(x1, x2), x3) ^ 1


def quad_roots(F, A, B, C):
    """Roots of A t^2 + B t + C = 0 (B != 0 path; A = 0 linear)."""
    if A == 0:
        return [F.mul(C, F.inv(B))] if B else []
    if B == 0:
        return []  # measure-zero here; skipped, counted as a failed trial
    k = F.mul(F.mul(A, C), F.inv(F.sq(B)))
    if F.tr(k):
        return []
    z, f = F.htr(k), F.mul(B, F.inv(A))
    return [F.mul(f, z), F.mul(f, z ^ 1)]


def lift(F, x):
    if x == 0:
        return None
    k = x ^ F.inv(F.sq(x))  # (y/x)^2 + (y/x) = x + 1/x^2
    if F.tr(k):
        return None
    return (x, F.mul(x, F.htr(k)))


def add(F, P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2:
        if y1 ^ y2 == x1 or x1 == 0:
            return None
        lam = x1 ^ F.mul(y1, F.inv(x1))
        x3 = F.sq(lam) ^ lam
        return (x3, F.sq(x1) ^ F.mul(lam ^ 1, x3))
    lam = F.mul(y1 ^ y2, F.inv(x1 ^ x2))
    x3 = F.sq(lam) ^ lam ^ x1 ^ x2
    return (x3, F.mul(lam, x1 ^ x3) ^ x3 ^ y1)


def neg(P):
    return None if P is None else (P[0], P[0] ^ P[1])


def eliminate(rows, ncols):
    """rows: bitmasks, bit ncols = right-hand side. Returns (rank, consistent, pivots)."""
    piv, cons = {}, True
    for r in rows:
        for c in range(ncols):
            if not (r >> c) & 1:
                continue
            if c in piv:
                r ^= piv[c]
            else:
                piv[c] = r
                break
        else:
            if (r >> ncols) & 1:
                cons = False
    return len(piv), cons, piv


def columns(F, V, W, w0, t):
    """Field-element coefficient of each linearized monomial, and the constant."""
    t2, w02 = F.sq(t), F.sq(w0)
    cols = []
    for vi in V:
        for wj in W:
            cols.append(F.mul(F.sq(vi), F.sq(wj)) ^ F.mul(F.mul(vi, wj), t))
    for vi in V:
        cols.append(F.mul(F.sq(vi), w02) ^ F.mul(F.mul(vi, w0), t) ^ F.mul(F.sq(vi), t2))
    for wj in W:
        cols.append(F.mul(F.sq(wj), t2))
    return cols, F.mul(w02, t2) ^ 1


def to_rows(n, cols, rhs):
    N, rows = len(cols), []
    for bit in range(n):
        r = 0
        for k, cv in enumerate(cols):
            if (cv >> bit) & 1:
                r |= 1 << k
        if (rhs >> bit) & 1:
            r |= 1 << N
        rows.append(r)
    return rows


def rand_in(V, rng):
    x = 0
    for v in V:
        if rng.getrandbits(1):
            x ^= v
    return x


def rank_defect_cell(n, l, b, trials, rng):
    """Rank defect N - rank of the linearized system with W < V (the real setting)."""
    F, hist = Field(n), {}
    for _ in range(trials):
        V = [rng.getrandbits(n) for _ in range(l)]
        W = [V[i] for i in rng.sample(range(l), b)]
        cols, rhs = columns(F, V, W, rand_in(V, rng), rng.getrandbits(n))
        rk, _, _ = eliminate(to_rows(n, cols, rhs), len(cols))
        d = len(cols) - rk
        hist[d] = hist.get(d, 0) + 1
    return hist


MAX_FAMILY = 8  # enumerate at most 2^8 points of the affine solution family


def oracle_trial(F, l, b, rng, stats):
    n = F.n
    V = [rng.getrandbits(n) for _ in range(l)]
    while True:
        R = lift(F, rng.getrandbits(n))
        if R:
            break
    xR = R[0]
    X1 = rand_in(V, rng)
    A = F.sq(X1) ^ F.sq(xR)
    B = F.mul(X1, xR)
    C = F.mul(F.sq(X1), F.sq(xR)) ^ 1
    W = [V[i] for i in rng.sample(range(l), b)]
    w0 = rand_in(V, rng)
    for t in quad_roots(F, A, B, C):
        cols, rhs = columns(F, V, W, w0, t)
        N = len(cols)
        rk, cons, piv = eliminate(to_rows(n, cols, rhs), N)
        stats["solves"] += 1
        if not cons or N - rk > MAX_FAMILY:
            if cons:
                stats["family_too_big"] += 1
            continue
        free = [c for c in range(N) if c not in piv]
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
            cv, dv = sol[l * b:l * b + l], sol[l * b + l:]
            if any(sol[i * b + j] != (cv[i] & dv[j]) for i in range(l) for j in range(b)):
                continue
            X2 = 0
            for i, v in enumerate(V):
                if cv[i]:
                    X2 ^= v
            X3 = w0
            for j, w in enumerate(W):
                if dv[j]:
                    X3 ^= w
            if S3(F, X2, X3, t) or S3(F, X1, t, xR):
                stats["algebra_mismatch"] += 1
                continue
            pts = [lift(F, x) for x in (X1, X2, X3)]
            if any(p is None for p in pts):
                continue
            for s in range(8):
                q = [pts[k] if (s >> k) & 1 else neg(pts[k]) for k in range(3)]
                if add(F, add(F, q[0], q[1]), q[2]) == R:
                    return True
            stats["curve_reject"] += 1
    return False


def main(out_dir, seed):
    rng = random.Random(seed)
    res = {"seed": seed, "rank_defect": [], "oracle": []}
    for n in (23, 29, 31):
        for l in range(3, n):
            for b in range(1, min(l, 4) + 1):
                if not (n - 8 <= l * b + l <= n + 4):
                    continue
                hist = rank_defect_cell(n, l, b, 200, rng)
                res["rank_defect"].append(dict(n=n, l=l, b=b, lb_plus_l=l * b + l,
                                               hist={str(k): v for k, v in sorted(hist.items())}))
    for n, l, trials in [(13, 3, 40000), (13, 4, 40000), (17, 4, 60000),
                         (17, 5, 60000), (19, 5, 60000)]:
        F = Field(n)
        bmax = (n - l) // l
        for b in range(0, bmax + 1):
            stats = dict(solves=0, family_too_big=0, algebra_mismatch=0, curve_reject=0)
            ok = sum(oracle_trial(F, l, b, rng, stats) for _ in range(trials))
            row = dict(n=n, l=l, b=b, bmax=bmax, success=ok, trials=trials,
                       log2p=(math.log2(ok / trials) if ok else None),
                       law=l + b - n, **stats)
            res["oracle"].append(row)
            print(json.dumps(row), flush=True)
    with open(f"{out_dir}/results.json", "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
