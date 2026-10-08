"""Multihomogeneous XL on the trilinear S3 chain with t unknown.

Variables (Boolean, four blocks):
  A = c1 (a1 bits): X1 = w1 + sum c1_i u_i, W1 < V spanned by a1 basis vectors of V
  B = c2 (l bits):  X2 = sum c2_i v_i, V the factor base
  C = d  (b bits):  X3 = w0 + sum d_j w_j, W < V
  T = tau (n bits): t = sum tau_k z^k, polynomial basis of F_{2^n}
Equations (2n bits):
  E1 = S3(X1, t, xR)   multidegree <= (1, 0, 0, 1)  bilinear
  E2 = S3(X2, X3, t)   multidegree <= (0, 1, 1, 1)  trilinear
XL at block caps (dA, dB, dC, dT): multiply each E1 bit by every Boolean monomial
of block degree <= (dA-1, dB, dC, dT-1), each E2 bit by <= (dA, dB-1, dC-1, dT-1),
keep only products inside the caps' monomial set (all are, by construction),
row-reduce with the affine-linear columns last, extract the linear part.

With a1 = 0 (X1 fixed) E1 is linear in tau and the system is the bilinear core
of ../xl_bilinear_core_20261006.  The question is what a1 > 0 costs.

Instances are planted: draw X1, X2, X3; t is a root of S3(X2, X3, t) = 0 and
xR a root of S3(X1, t, xR) = 0.  Resolved = the extracted linear system leaves
at most 8 affine points (dimension <= 3), one of them the planted point.  Up
to four genuine solutions are expected: two roots t times the X2 <-> X3 swap.
"""
import itertools
import json
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "linearization_reach_20260930"))
import lr  # noqa: E402  (frozen field and quadratic solver, reused unchanged)

# x^11 + x^2 + 1 is irreducible over F_2; lr.Field checks it at construction.
lr.IRR.setdefault(11, (1 << 11) | 0b101)


# ---- polynomials over F_{2^n} in Boolean variables: {mask: field element} ----

def padd(p, q):
    r = dict(p)
    for m, a in q.items():
        v = r.get(m, 0) ^ a
        if v:
            r[m] = v
        else:
            r.pop(m, None)
    return r


def pmul(F, p, q):
    r = {}
    for m1, a in p.items():
        for m2, b in q.items():
            m = m1 | m2
            v = r.get(m, 0) ^ F.mul(a, b)
            if v:
                r[m] = v
            else:
                r.pop(m, None)
    return r


def psq(F, p):
    # char 2 and x^2 = x for Boolean x: (sum a_m m)^2 = sum a_m^2 m
    r = {}
    for m, a in p.items():
        v = r.get(m, 0) ^ F.sq(a)
        if v:
            r[m] = v
        else:
            r.pop(m, None)
    return r


def S3poly(F, x1, x2, x3):
    s = padd(padd(pmul(F, x1, x2), pmul(F, x1, x3)), pmul(F, x2, x3))
    return padd(padd(psq(F, s), pmul(F, pmul(F, x1, x2), x3)), {0: 1})


def affine(const, basis, offset):
    """const + sum var_{offset+i} * basis[i] as a polynomial."""
    p = {0: const} if const else {}
    for i, v in enumerate(basis):
        if v:
            p = padd(p, {1 << (offset + i): v})
    return p


def bit_equations(F, poly):
    """One F_2 equation per coordinate bit: the set of monomials present in it."""
    return [{m for m, a in poly.items() if (a >> bit) & 1} for bit in range(F.n)]


# ---- blocks and monomials ----

class Blocks:
    def __init__(self, a1, l, b, n):
        self.sizes = (a1, l, b, n)
        self.off = (0, a1, a1 + l, a1 + l + b)
        self.nvar = a1 + l + b + n
        self.masks = tuple(((1 << s) - 1) << o for s, o in zip(self.sizes, self.off))

    def block_subsets(self, k, cap):
        size, off = self.sizes[k], self.off[k]
        out = []
        for d in range(min(cap, size) + 1):
            for comb in itertools.combinations(range(size), d):
                m = 0
                for i in comb:
                    m |= 1 << (off + i)
                out.append(m)
        return out

    def monomials(self, caps):
        if any(c < 0 for c in caps):
            return []
        parts = [self.block_subsets(k, caps[k]) for k in range(4)]
        return [a | b | c | d for a in parts[0] for b in parts[1] for c in parts[2] for d in parts[3]]


def popcount(x):
    return bin(x).count("1")


def macaulay(B, eqs1, eqs2, caps):
    dA, dB, dC, dT = caps
    cols = B.monomials(caps)
    lin = [1 << v for v in range(B.nvar)] + [0]
    linset = set(lin)
    nonlin = [m for m in cols if m not in linset]
    nonlin.sort(key=lambda m: -popcount(m))
    order = nonlin + lin
    col = {m: k for k, m in enumerate(order)}
    rows = []
    for eqs in (eqs1, eqs2):
        # multiplier caps = caps minus the equations' actual block degrees
        deg = [max((popcount(m & B.masks[k]) for e in eqs for m in e), default=0) for k in range(4)]
        mcaps = tuple(caps[k] - deg[k] for k in range(4))
        mults = B.monomials(mcaps)
        for e in eqs:
            for mu in mults:
                r = 0
                for m in e:
                    r ^= 1 << col[m | mu]
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


def solve_linear(piv, n_nonlin, nvar, cap_dim=3):
    lins = [r >> n_nonlin for c, r in piv.items() if c >= n_nonlin]
    p2 = {}
    for r in lins:
        while r:
            low = (r & -r).bit_length() - 1
            if low == nvar:
                return None, None
            if low in p2:
                r ^= p2[low]
            else:
                p2[low] = r
                break
    free = [v for v in range(nvar) if v not in p2]
    if len(free) > cap_dim:
        return len(free), None
    pts = []
    for mask in range(1 << len(free)):
        x = [0] * nvar
        for k, v in enumerate(free):
            x[v] = (mask >> k) & 1
        for v in sorted(p2, reverse=True):
            r = p2[v]
            s = (r >> nvar) & 1
            for u in range(v + 1, nvar):
                if (r >> u) & 1:
                    s ^= x[u]
            x[v] = s
        pts.append(x)
    return len(free), pts


# ---- planted instances ----

def rand_bits(rng, k):
    return [rng.getrandbits(1) for _ in range(k)]


def combo(const, basis, bits):
    x = const
    for v, b in zip(basis, bits):
        if b:
            x ^= v
    return x


def planted(F, a1, l, b, rng):
    n = F.n
    while True:
        V = [rng.getrandbits(n) for _ in range(l)]
        W = [V[i] for i in rng.sample(range(l), b)]
        w0 = combo(0, V, rand_bits(rng, l))
        # X1 is a factor-base point: its coset w1 + W1 lies in V, so a1 <= l
        U = [V[i] for i in rng.sample(range(l), a1)]
        w1 = combo(0, V, rand_bits(rng, l))
        c1, c2, d = rand_bits(rng, a1), rand_bits(rng, l), rand_bits(rng, b)
        X1, X2, X3 = combo(w1, U, c1), combo(0, V, c2), combo(w0, W, d)
        ts = lr.quad_roots(F, F.sq(X2) ^ F.sq(X3), F.mul(X2, X3), F.mul(F.sq(X2), F.sq(X3)) ^ 1)
        if not ts:
            continue
        t = ts[rng.getrandbits(1) % len(ts)]
        # S3(X1, t, xR) as a quadratic in xR
        xs = lr.quad_roots(F, F.sq(X1) ^ F.sq(t), F.mul(X1, t), F.mul(F.sq(X1), F.sq(t)) ^ 1)
        if not xs:
            continue
        xR = xs[rng.getrandbits(1) % len(xs)]
        assert lr.S3(F, X1, t, xR) == 0 and lr.S3(F, X2, X3, t) == 0
        tau = [(t >> k) & 1 for k in range(n)]
        return dict(U=U, w1=w1, V=V, W=W, w0=w0, xR=xR, sol=c1 + c2 + d + tau)


def build(F, inst, a1, l, b):
    n = F.n
    Bk = Blocks(a1, l, b, n)
    X1 = affine(inst["w1"], inst["U"], Bk.off[0])
    X2 = affine(0, inst["V"], Bk.off[1])
    X3 = affine(inst["w0"], inst["W"], Bk.off[2])
    t = affine(0, [1 << k for k in range(n)], Bk.off[3])
    xR = {0: inst["xR"]} if inst["xR"] else {}
    E1 = S3poly(F, X1, t, xR)
    E2 = S3poly(F, X2, X3, t)
    return Bk, bit_equations(F, E1), bit_equations(F, E2)


def cell(n, a1, l, b, caps, trials, rng):
    F = lr.Field(n)
    out = dict(n=n, a1=a1, l=l, b=b, caps=list(caps), trials=trials, resolved=0,
               inconsistent=0, dims={}, rows=0, cols=0)
    for _ in range(trials):
        inst = planted(F, a1, l, b, rng)
        Bk, e1, e2 = build(F, inst, a1, l, b)
        rows, order, nn = macaulay(Bk, e1, e2, caps)
        out["rows"], out["cols"] = len(rows), len(order)
        piv = reduce_rows(rows)
        dim, pts = solve_linear(piv, nn, Bk.nvar)
        if dim is None:
            out["inconsistent"] += 1
            continue
        out["dims"][str(dim)] = out["dims"].get(str(dim), 0) + 1
        if pts is not None and inst["sol"] in pts:
            out["resolved"] += 1
    return out


if __name__ == "__main__":
    n, a1, l, b, dA, dB, dC, dT, T, seed = map(int, sys.argv[1:11])
    print(json.dumps(cell(n, a1, l, b, (dA, dB, dC, dT), T, random.Random(seed))))
