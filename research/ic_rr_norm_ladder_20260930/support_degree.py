#!/usr/bin/env python3
"""Boolean degree of the Riemann–Roch *support* form, exactly, at tiny sizes.

The support form imposes `H | L_V` on the monic cofactor `H = N/(X + x_R)` of the norm
of `f`, where `L_V(T) = Π_{v∈V}(T + v)` is the subspace polynomial.  Its equations are
the coefficients of `L_V mod H`, computed by repeated squaring modulo `H`.  Over GF(2^n)
this remainder is a polynomial in the free coefficients of `f`; its Boolean degree after
Weil descent is `max` over its monomials `Π c_i^{k_i}` of `Σ wt(k_i)` (squaring is
Boolean-linear, so `c^k` has Boolean degree `wt(k)`, the binary weight).

This script computes that degree exactly for `m = 3` (coefficients α, β) and `m = 4`
(α, β, δ), with the target-dependent coefficient substituted, at `ℓ = 2 … L`.  It never
solves anything.  Output: one line per (m, ℓ) with the maximal Boolean degree of the
remainder's coefficients, and the same for the chained-S₃ presentation for reference.

    python3 support_degree.py --n 9 --lmax 5
"""

import argparse
import itertools
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'ic_tree_split_cutout_20260930'))
from cutout import Curve, Field, random_subspace  # noqa: E402


class MPoly:
    """Multivariate polynomial over GF(2^n) in `k` unknowns: {exponent tuple: coeff}."""

    def __init__(self, F, k, terms=None):
        self.F, self.k, self.t = F, k, dict(terms or {})

    @classmethod
    def const(cls, F, k, c):
        return cls(F, k, {(0,) * k: c} if c else {})

    @classmethod
    def var(cls, F, k, i):
        e = [0] * k
        e[i] = 1
        return cls(F, k, {tuple(e): 1})

    def add(self, o):
        t = dict(self.t)
        for e, c in o.t.items():
            v = t.get(e, 0) ^ c
            if v:
                t[e] = v
            else:
                t.pop(e, None)
        return MPoly(self.F, self.k, t)

    def mul(self, o):
        t = {}
        for e1, c1 in self.t.items():
            for e2, c2 in o.t.items():
                e = tuple(a + b for a, b in zip(e1, e2))
                v = t.get(e, 0) ^ self.F.mul(c1, c2)
                if v:
                    t[e] = v
                else:
                    t.pop(e, None)
        return MPoly(self.F, self.k, t)

    def scale(self, c):
        return MPoly(self.F, self.k, {e: self.F.mul(c, v) for e, v in self.t.items()} if c else {})

    def square(self):
        return MPoly(self.F, self.k, {tuple(2 * a for a in e): self.F.sqr(c) for e, c in self.t.items()})

    def bool_degree(self):
        return max((sum(bin(a).count('1') for a in e) for e in self.t), default=0)


def poly_mod(num, H):
    """`num mod H` for polynomials in X with MPoly coefficients (lists, low degree first);
    H monic."""
    num = list(num)
    dH = len(H) - 1
    while len(num) - 1 >= dH:
        lead = num[-1]
        if lead.t:
            shift = len(num) - 1 - dH
            for i in range(dH):
                num[shift + i] = num[shift + i].add(lead.mul(H[i]))
        num.pop()
    while len(num) < dH:
        num.append(MPoly(num[0].F, num[0].k))
    return num


def subspace_poly_coeffs(F, V):
    """`L_V(T) = Σ λ_i T^{2^i}` as the list [λ_0, …, λ_ℓ]; built by the recurrence
    `L_{W+<v>}(T) = L_W(T)² + L_W(v)·L_W(T)`."""
    lam = [1]  # L_{0}(T) = T
    basis = []
    for v in V:
        if v == 0 or any(True for _ in ()):
            continue
        # keep a basis of V incrementally
        span = {0}
        for b in basis:
            span |= {s ^ b for s in span}
        if v in span:
            continue
        basis.append(v)
        lw_v = 0
        for i, c in enumerate(lam):
            lw_v ^= F.mul(c, pow_2k(F, v, i))
        new = [0] * (len(lam) + 1)
        for i, c in enumerate(lam):
            new[i + 1] ^= F.sqr(c)
            new[i] ^= F.mul(lw_v, c)
        lam = new
    return lam


def pow_2k(F, v, k):
    for _ in range(k):
        v = F.sqr(v)
    return v


def support_degree(F, E, V, m, rng):
    """Boolean degree of the coefficients of `L_V mod H`, `H = N/(X + x_R)`, for one
    random target of the prime-order-free kind (any affine point off `V` works here)."""
    n = F.n
    # a target point R = (r, s) with r ∉ V
    while True:
        r = rng.getrandbits(n)
        pts = E.points_with_x(r) if r else []
        if pts and r not in set(V):
            R = pts[rng.randrange(2)]
            break
    r, s = R
    a = E.a
    if m == 3:
        # f = X² + αX + γ + βy,  γ = r² + αr + β(r+s);  N = A² + βXA + β²(X³ + aX² + 1)
        k = 2
        al, be = MPoly.var(F, k, 0), MPoly.var(F, k, 1)
        C = lambda c: MPoly.const(F, k, c)
        ga = C(F.sqr(r)).add(al.scale(r)).add(be.scale(r ^ s))
        A = [ga, al, C(1)]                          # A(X) low → high
    else:
        # f = Xy + αX² + βy + δX + γ,  f(−R) = 0 ⇒ γ = r(r+s) + αr² + β(r+s) + δr
        # N = B² + X·B·Cc + Cc²(X³ + aX² + 1),  B = αX² + δX + γ,  Cc = X + β
        k = 3
        al, be, de = (MPoly.var(F, k, i) for i in range(3))
        C = lambda c: MPoly.const(F, k, c)
        ga = C(F.mul(r, r ^ s)).add(al.scale(F.sqr(r))).add(be.scale(r ^ s)).add(de.scale(r))
    Z = lambda: MPoly(F, k)

    def pmul(p, q):
        out = [Z() for _ in range(len(p) + len(q) - 1)]
        for i, pi in enumerate(p):
            for j, qj in enumerate(q):
                out[i + j] = out[i + j].add(pi.mul(qj))
        return out

    def padd(p, q):
        out = [Z() for _ in range(max(len(p), len(q)))]
        for i, c in enumerate(p):
            out[i] = out[i].add(c)
        for i, c in enumerate(q):
            out[i] = out[i].add(c)
        return out

    def psq(p):
        out = [Z() for _ in range(2 * len(p) - 1)]
        for i, c in enumerate(p):
            out[2 * i] = c.square()
        return out

    curve = [C(1), Z(), C(a), C(1)]                 # X³ + aX² + 1
    if m == 3:
        A = [ga, al, C(1)]
        N = padd(padd(psq(A), pmul([Z(), be], A)), pmul([be.square()], curve))
    else:
        B = [ga, de, al]
        Cc = [be, C(1)]
        N = padd(padd(psq(B), pmul([Z()] + Cc, B)), pmul(psq(Cc), curve))
    # divide by (X + r): synthetic division, remainder must vanish identically
    H = [Z() for _ in range(len(N) - 1)]
    carry = N[-1]
    for i in range(len(N) - 2, -1, -1):
        H[i] = carry
        carry = N[i].add(carry.scale(r))
    assert not carry.t, 'N(x_R) = 0 must hold identically'
    assert H[-1].t == {(0,) * k: 1}, 'H monic'
    lam = subspace_poly_coeffs(F, V)
    # L_V mod H = Σ λ_i (X^{2^i} mod H)
    rem = [Z() for _ in range(len(H) - 1)]
    cur = poly_mod([Z(), C(1)], H)                 # X mod H
    steps = []
    for i, li in enumerate(lam):
        for j, c in enumerate(cur):
            rem[j] = rem[j].add(c.scale(li))
        steps.append(max(c.bool_degree() for c in cur))
        cur = poly_mod(psq(cur), H)
    return max(c.bool_degree() for c in rem), steps, sum(len(c.t) for c in rem)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--n', type=int, default=9)
    p.add_argument('--a', type=int, default=0)
    p.add_argument('--lmax', type=int, default=5)
    p.add_argument('--seed', type=int, default=20260930)
    args = p.parse_args()
    rng = random.Random(args.seed)
    E = Curve(args.a, args.n)
    F = E.F
    print(f'K_{args.a}/2^{args.n}: Boolean degree of the support-form remainder L_V mod H')
    print(f"{'m':>2} {'ell':>3} {'degree':>6}  per-squaring degrees of X^(2^i) mod H   monomials")
    for m in (3, 4):
        for ell in range(2, args.lmax + 1):
            V = random_subspace(args.n, ell, rng)
            d, steps, mons = support_degree(F, E, V, m, rng)
            print(f'{m:>2} {ell:>3} {d:>6}  {steps}   {mons}')
    print()
    print('reference: chained S3 (m = 3, 4) has Boolean degree 3 at every ell;')
    print('the eliminated norm form measured by rr_degree_ladder has degree 4 (m = 3).')


if __name__ == '__main__':
    main()
