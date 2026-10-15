#!/usr/bin/env python3
"""Exact degree-2 refutation criterion for the m = 2 system S3(x1, x2, X) = 0, x_i in V.

The Weil-descended equations are the n coordinates of S3, so their F_2-span is
{Tr(theta * S3) : theta in F_{2^n}}.  On y^2 + xy = x^3 + 1,
    S3 = (x1 x2 + x1 X + x2 X)^2 + x1 x2 X + 1
and, using Tr(a^2) = Tr(a),
    Tr(theta S3) = Tr(c x1 x2) + Tr(s X x1) + Tr(s X x2) + Tr(theta),
    s = sqrt(theta),  c = s + theta X.
x1 x2 is bilinear in the Boolean coordinates of x1 and x2, so no Boolean
reduction mixes the parts.  The degree-2 Macaulay matrix (just the equations)
therefore contains 1 iff some theta satisfies
    Tr(c v_i v_j) = 0 for all basis pairs,  Tr(s X v_i) = 0 for all i,  Tr(theta) = 1.
That is linear algebra over F_2 in theta.  This script decides it exactly for
a given V and target X.  It needs no Groebner run, so it is an independent check
of what arm D measured.  It also reports dim span(V*V), the quantity that makes
the quadratic conditions collapse for Frobenius-stable V.

Usage: python3 -I mechanism.py <outdir>
"""
import json
import random
import subprocess
import sys


def irreducible(n):
    out = subprocess.run(["gp", "-q", "-f"], input=f"print(Vec(lift(ffinit(2,{n}))))",
                         capture_output=True, text=True, check=True).stdout
    coeffs = [int(c) & 1 for c in out.strip().strip("[]").split(",")]
    m = 0
    for c in coeffs:
        m = (m << 1) | c
    return m


class F:
    def __init__(self, n):
        self.n, self.mod = n, irreducible(n)
        self.tmask = 0
        for t in range(n):
            self.tmask |= self._tr_slow(1 << t) << t

    def mul(self, a, b):
        r = 0
        n, mod = self.n, self.mod
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if (a >> n) & 1:
                a ^= mod
        return r

    def sqr(self, a):
        return self.mul(a, a)

    def pw2(self, a, k):
        for _ in range(k % self.n):
            a = self.sqr(a)
        return a

    def tr(self, a):
        return bin(a & self.tmask).count("1") & 1  # trace is F_2-linear

    def _tr_slow(self, a):
        t, acc = a, 0
        for _ in range(self.n):
            acc ^= t
            t = self.sqr(t)
        return acc & 1  # trace lies in F_2 = {0, 1}

    def sqrt(self, a):
        return self.pw2(a, self.n - 1)

    def inv(self, a):
        # a^(2^n - 2)
        r, e = 1, (1 << self.n) - 2
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.sqr(a)
            e >>= 1
        return r


def rank(vecs):
    piv = {}
    r = 0
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v
                r += 1
                break
    return r


def kernel_of_rows(rows, nbits):
    """Solutions x in F_2^nbits of <row, x> = 0 for all rows (rows as int masks).
    Returns a basis of the solution space."""
    piv = {}  # pivot bit -> row (reduced)
    for r in rows:
        for p, pr in piv.items():
            if (r >> p) & 1:
                r ^= pr
        if r:
            p = r.bit_length() - 1
            for q in list(piv):
                if (piv[q] >> p) & 1:
                    piv[q] ^= r
            piv[p] = r
    free = [b for b in range(nbits) if b not in piv]
    basis = []
    for f in free:
        x = 1 << f
        for p, pr in piv.items():
            if (pr >> f) & 1:
                x |= 1 << p
        basis.append(x)
    return basis


def stable_basis(fld, g):
    n = fld.n
    cols = []
    for k in range(n):
        x, acc = 1 << k, 0
        for i in range(g.bit_length()):
            if (g >> i) & 1:
                acc ^= fld.pw2(x, i)
        cols.append(acc)
    # kernel of the map e_k -> cols[k]: rows of the transposed system
    rows = []
    for bit in range(n):
        r = 0
        for k in range(n):
            if (cols[k] >> bit) & 1:
                r |= 1 << k
        rows.append(r)
    return kernel_of_rows(rows, n)


def random_basis(n, k, rng, fld=None):
    """Random k-dim subspace; with fld given, drawn inside the trace-zero hyperplane."""
    b = []
    while rank(b) < k:
        v = rng.getrandbits(n)
        if fld is not None and fld.tr(v):
            continue
        if v and rank(b + [v]) > len(b):
            b.append(v)
    return b


def degree2_refutes(fld, basis, X):
    """Exact test for some theta with the three conditions.  theta is written in
    the monomial basis e_t = z^t, t < n.  Each condition is linear in theta."""
    n = fld.n
    E = [1 << t for t in range(n)]
    S = [fld.sqrt(e) for e in E]  # sqrt is F_2-linear
    C = [S[t] ^ fld.mul(E[t], X) for t in range(n)]
    rows = []
    prods = [fld.mul(basis[i], basis[j]) for i in range(len(basis)) for j in range(len(basis))]
    W = []  # basis of span(V*V); quadratic conditions only depend on it
    for p in prods:
        if rank(W + [p]) > len(W):
            W.append(p)
    for w in W:
        rows.append(sum(fld.tr(fld.mul(C[t], w)) << t for t in range(n)))
    for v in basis:
        xv = fld.mul(X, v)
        rows.append(sum(fld.tr(fld.mul(S[t], xv)) << t for t in range(n)))
    sol = kernel_of_rows(rows, n)
    tr_e = [fld.tr(e) for e in E]
    # some solution with Tr(theta) = 1 ?
    return any(sum(tr_e[t] for t in range(n) if (x >> t) & 1) & 1 for x in sol), len(W), len(sol)


def on_curve(fld, x):
    return x == 0 or fld.tr(x ^ fld.inv(fld.sqr(x))) == 0


def main():
    out = sys.argv[1]
    rng = random.Random(20261014)
    cases = [(7, 4, "g=0x17"), (18, 9, "g=0x201"), (17, 9, "g=0x279"), (31, 15, "g=0x808b"),
             (31, 16, "g=0x10117"), (31, 15, "random"), (31, 16, "random"),
             (17, 9, "random"), (18, 9, "random"), (31, 15, "random-tr0"), (31, 16, "random-tr0"),
             (17, 9, "random-tr0")]
    res = []
    for n, k, base in cases:
        fld = F(n)
        for rep in range(1 if base.startswith("g=") else 3):
            if base.startswith("g="):
                basis = stable_basis(fld, int(base[2:], 16))
            else:
                basis = random_basis(n, k, rng, fld if base == "random-tr0" else None)
            v_in_tr0 = all(fld.tr(v) == 0 for v in basis)
            pred = 0
            assert len(basis) == k, (n, k, base, len(basis))
            hits = tot = 0
            wdim = None
            for _ in range(200):
                X = rng.getrandbits(n)
                if X == 0 or not on_curve(fld, X):
                    continue
                r, wdim, _ = degree2_refutes(fld, basis, X)
                hits += r
                pred += int(v_in_tr0 and fld.tr(fld.inv(X)) == 1)
                tot += 1
            row = dict(n=n, k=k, base=base, rep=rep, span_VV_dim=wdim, targets=tot,
                       degree2_refutable=hits, fraction=hits / tot, V_in_trace0=v_in_tr0,
                       predicted_by_trace_rule=pred, rule_matches=(pred == hits))
            print(row, flush=True)
            res.append(row)
    json.dump(res, open(f"{out}/mechanism.json", "w"), indent=1)


if __name__ == "__main__":
    main()
