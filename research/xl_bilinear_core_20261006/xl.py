"""Bidegree XL on the bilinear S3 core.

System (from ../linearization_reach_20260930): with t fixed, S3(X2, X3, t) = 0 is
n Boolean equations, bilinear in c (X2 in V, dim l) and d (X3 in w0 + W, W < V,
dim b).  XL at bidegree (D1, D2) multiplies every equation by every Boolean
monomial c^alpha d^beta with |alpha| <= D1 - 1, |beta| <= D2 - 1, reduces
multilinearly (x^2 = x), and row-reduces the Macaulay matrix with the linear
and constant columns ordered last.  The rows whose pivot falls in that block
span every affine-linear polynomial in the degree-(D1, D2) row space.

Instances are PLANTED: X2 in V and X3 in w0 + W are drawn first and t is a
root of S3(X2, X3, t) = 0, so a solution exists.  XL "resolves" the instance
when the affine solution set of the extracted linear polynomials has at most
4 points and contains the planted solution.  (The X2 <-> X3 swap on W gives a
second genuine solution, so 2 points is the expected floor.)
"""
import itertools
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "linearization_reach_20260930"))
import lr  # noqa: E402  (frozen field, curve and column code, reused unchanged)


def boolean_monomials(k, dmax):
    """All subsets of range(k) of size <= dmax, as bitmasks, by degree."""
    out = []
    for d in range(dmax + 1):
        for comb in itertools.combinations(range(k), d):
            m = 0
            for i in comb:
                m |= 1 << i
            out.append(m)
    return out


def equation_polys(F, V, W, w0, t, l, b):
    """n bit-equations as dicts {(cmask, dmask): 1} over F_2."""
    cols, rhs = lr.columns(F, V, W, w0, t)
    monos = [(1 << i, 1 << j) for i in range(l) for j in range(b)]
    monos += [(1 << i, 0) for i in range(l)] + [(0, 1 << j) for j in range(b)]
    polys = []
    for bit in range(F.n):
        p = set()
        for m, cv in zip(monos, cols):
            if (cv >> bit) & 1:
                p ^= {m}
        if (rhs >> bit) & 1:
            p ^= {(0, 0)}
        polys.append(p)
    return polys


def popcount(x):
    return bin(x).count("1")


def macaulay(polys, l, b, D1, D2):
    """Rows as Python ints over a column index with linear/constant columns last."""
    cm = boolean_monomials(l, D1)
    dm = boolean_monomials(b, D2)
    monos = [(x, y) for x in cm for y in dm]
    lin = [(1 << i, 0) for i in range(l)] + [(0, 1 << j) for j in range(b)] + [(0, 0)]
    linset = set(lin)
    nonlin = [m for m in monos if m not in linset]
    # higher total degree first so elimination removes it first
    nonlin.sort(key=lambda m: -(popcount(m[0]) + popcount(m[1])))
    order = nonlin + lin
    col = {m: k for k, m in enumerate(order)}
    mult = [(x, y) for x in boolean_monomials(l, D1 - 1) for y in boolean_monomials(b, D2 - 1)]
    rows = []
    for p in polys:
        for (x, y) in mult:
            r = 0
            for (pc, pd) in p:
                r ^= 1 << col[(pc | x, pd | y)]
            if r:
                rows.append(r)
    return rows, order, len(nonlin)


def reduce_rows(rows):
    piv = {}
    for r in rows:
        while r:
            low = (r & -r).bit_length() - 1
            if low in piv:
                r ^= piv[low]
            else:
                piv[low] = r
                break
    return piv


def extract_linear(piv, n_nonlin, l, b):
    """Linear polynomials in the row space: pivots in the linear block."""
    out = []
    for c, r in piv.items():
        if c >= n_nonlin:
            out.append(r >> n_nonlin)  # bits 0..l+b-1 variables, bit l+b constant
    return out


def affine_solutions(lin_rows, nvar, cap=4):
    """Solve the linear system; None if inconsistent, else (dim, points<=cap)."""
    piv = {}
    for r in lin_rows:
        while r:
            low = (r & -r).bit_length() - 1
            if low == nvar:          # 1 = 0
                return None
            if low in piv:
                r ^= piv[low]
            else:
                piv[low] = r
                break
    free = [v for v in range(nvar) if v not in piv]
    dim = len(free)
    if dim > 2:
        return dim, None
    pts = []
    for mask in range(1 << dim):
        x = [0] * nvar
        for k, v in enumerate(free):
            x[v] = (mask >> k) & 1
        for v in sorted(piv, reverse=True):
            r = piv[v]
            s = (r >> nvar) & 1
            for u in range(v + 1, nvar):
                if (r >> u) & 1:
                    s ^= x[u]
            x[v] = s
        pts.append(x)
    return dim, pts


def planted_instance(F, l, b, rng):
    while True:
        V = [rng.getrandbits(F.n) for _ in range(l)]
        idx = rng.sample(range(l), b)
        W = [V[i] for i in idx]
        w0vec = [rng.getrandbits(1) for _ in range(l)]
        w0 = 0
        for i, bit in enumerate(w0vec):
            if bit:
                w0 ^= V[i]
        c = [rng.getrandbits(1) for _ in range(l)]
        d = [rng.getrandbits(1) for _ in range(b)]
        X2 = 0
        for i in range(l):
            if c[i]:
                X2 ^= V[i]
        X3 = w0
        for j in range(b):
            if d[j]:
                X3 ^= W[j]
        # S3(X2, X3, t) = (X2^2 + X3^2) t^2 + X2 X3 t + X2^2 X3^2 + 1 in t
        A = F.sq(X2) ^ F.sq(X3)
        B = F.mul(X2, X3)
        C = F.mul(F.sq(X2), F.sq(X3)) ^ 1
        roots = lr.quad_roots(F, A, B, C)
        if roots:
            t = roots[rng.getrandbits(1) % len(roots)]
            assert lr.S3(F, X2, X3, t) == 0
            return V, W, w0, t, c + d


def xl_cell(n, l, b, D1, D2, trials, rng):
    F = lr.Field(n)
    res = dict(n=n, l=l, b=b, D1=D1, D2=D2, trials=trials, resolved=0, inconsistent=0,
               dims={}, rows=0, cols=0, nonlin=0, rank=0, lin_found=0)
    for _ in range(trials):
        V, W, w0, t, planted = planted_instance(F, l, b, rng)
        polys = equation_polys(F, V, W, w0, t, l, b)
        rows, order, nn = macaulay(polys, l, b, D1, D2)
        piv = reduce_rows(rows)
        lins = extract_linear(piv, nn, l, b)
        sol = affine_solutions(lins, l + b)
        res["rows"], res["cols"], res["nonlin"] = len(rows), len(order), nn
        res["rank"] += len(piv)
        res["lin_found"] += len(lins)
        if sol is None:
            res["inconsistent"] += 1     # impossible on a planted instance: a bug
            continue
        dim, pts = sol
        res["dims"][str(dim)] = res["dims"].get(str(dim), 0) + 1
        if pts is not None and planted in pts:
            res["resolved"] += 1
    res["rank"] /= trials
    res["lin_found"] /= trials
    return res


if __name__ == "__main__":
    n, l, b, D1, D2, T, seed = map(int, sys.argv[1:8])
    print(json.dumps(xl_cell(n, l, b, D1, D2, T, random.Random(seed))))
