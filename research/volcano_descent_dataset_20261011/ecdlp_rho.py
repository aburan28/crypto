#!/usr/bin/env python3
"""Counted-operation ECDLP solver for toy prime-field curves y^2 = x^3 + a x + b.

Pohlig-Hellman over the factorisation of the base point's order, with
distinguished-point Pollard rho (r-adding walk) in every prime-order layer
larger than SMALL, brute force below.  The counted unit is one affine group
operation (addition or doubling, each one field inversion); scalar
multiplications inside Pohlig-Hellman and walk set-up are counted too and
reported separately from the rho walk itself.

Equivalence-class modes (all use the same walk and distinguished-point code):

  plain  walk on points
  neg    walk on {P, -P}                    (every curve)
  aut    walk on the full F_p-rational automorphism orbit
         (mu_4 for j = 1728, mu_6 for j = 0, otherwise same as neg)

This is a pure-Python reference implementation, not the repository's ecbench
harness, and its counted unit is not ecbench's.  Results from here must not be
mixed with ecbench rows.
"""
from __future__ import annotations

import math
import random

from cypari import pari

SMALL = 1 << 12   # prime-order subgroups up to this size are solved by brute force


class Curve:
    __slots__ = ("p", "a", "b", "ops")

    def __init__(self, p, a, b):
        self.p, self.a, self.b, self.ops = p, a, b, 0

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        p = self.p
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        self.ops += 1
        x3 = (lam * lam - x1 - x2) % p
        return x3, (lam * (x1 - x3) - y1) % p

    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.p)

    def mul(self, k, P):
        if k < 0:
            k, P = -k, self.neg(P)
        R = None
        while k:
            if k & 1:
                R = self.add(R, P)
            k >>= 1
            if k:
                P = self.add(P, P)
        return R

    def on_curve(self, P):
        x, y = P
        return (y * y - x * x * x - self.a * x - self.b) % self.p == 0


def _sqrt_minus_one(p):
    for g in range(2, 1000):
        z = pow(g, (p - 1) // 4, p)
        if z * z % p == p - 1:
            return z
    raise ValueError


def _cube_root_unity(p):
    for g in range(2, 1000):
        z = pow(g, (p - 1) // 3, p)
        if z != 1:
            return z
    raise ValueError


def automorphisms(C: Curve, G, r, mode):
    """List of (map, lam) with map(P) = [lam]P on the order-r subgroup <G>."""
    p = C.p
    out = [(lambda P: C.neg(P), r - 1)]
    if mode != "aut":
        return out
    if C.b % p == 0 and C.a % p != 0:      # j = 1728, mu_4
        i = _sqrt_minus_one(p)
        roots = [u for u in _roots(r, [1, 0, 1])]
        iota = lambda P: ((-P[0]) % p, (i * P[1]) % p)
        lam = next(u for u in roots if C.mul(u, G) == iota(G))
        out.append((iota, lam))
        out.append((lambda P: iota(iota(iota(P))), (-lam) % r))   # iota^3 = -iota
    elif C.a % p == 0 and C.b % p != 0:    # j = 0, mu_6
        z = _cube_root_unity(p)
        roots = [u for u in _roots(r, [1, 1, 1])]
        zeta = lambda P: ((z * P[0]) % p, P[1])
        lam = next(u for u in roots if C.mul(u, G) == zeta(G))
        zeta2 = lambda P: zeta(zeta(P))
        out.append((zeta, lam))
        out.append((zeta2, lam * lam % r))
        out.append((lambda P: C.neg(zeta(P)), (-lam) % r))
        out.append((lambda P: C.neg(zeta2(P)), (-lam * lam) % r))
    return out


def _roots(r, coeffs):
    """Roots mod prime r of c0 u^2 + c1 u + c2 (coeffs = [c0, c1, c2])."""
    c0, c1, c2 = coeffs
    pol = pari(f"Mod({c0},{r})*x^2+Mod({c1},{r})*x+Mod({c2},{r})")
    return [int(z) for z in pari.polrootsmod(pol, r).lift()] if True else []


def rho_prime(C: Curve, G, Q, r, mode, rng, S=8192):
    """DLP in the subgroup of prime order r: Q = [x]G.  Returns (x, counted ops).

    Distinguished-point rho with an S-entry additive step table.  Trails start
    from the sum of two table points (one group operation).  For the equivalence
    modes (neg, aut) fruitless short cycles are detected against the last HIST
    points and left by one deterministic jump from the smallest point of the
    cycle, so every trail that enters a cycle leaves it at the same place.
    """
    if r <= SMALL:
        R = None
        for x in range(r):
            if R == Q:
                return x, 0
            R = C.add(R, G)
        raise ArithmeticError("no solution")
    auts = automorphisms(C, G, r, mode) if mode != "plain" else []
    dbits = max(1, int(round(math.log2(r) / 2)) - 5)
    mask = (1 << dbits) - 1
    cap = 20 << dbits
    HIST = 12
    # Step table: 16 random combinations u G + v Q, then each further entry is the sum of the
    # previous entry and a random earlier one (one group operation per entry).
    table = []
    for _ in range(16):
        u, v = rng.randrange(r), rng.randrange(r)
        table.append((u, v, C.add(C.mul(u, G), C.mul(v, Q))))
    while len(table) < S:
        u1, v1, R1 = table[-1]
        u2, v2, R2 = table[rng.randrange(len(table) - 1)]
        R3 = C.add(R1, R2)
        if R3 is None:
            continue
        table.append(((u1 + u2) % r, (v1 + v2) % r, R3))
    start_ops = C.ops

    def canon(P, a, b):
        if not auts:
            return P, a, b
        best, ba, bb = P, a, b
        for f, lam in auts:
            Pi = f(P)
            if Pi < best:
                best, ba, bb = Pi, lam * a % r, lam * b % r
        return best, ba, bb

    seen = {}
    trail = 0
    while True:
        trail += 1
        i, j = rng.randrange(S), rng.randrange(S)
        while j == i:
            j = rng.randrange(S)
        ui, vi, Ri = table[i]
        uj, vj, Rj = table[j]
        P = C.add(Ri, Rj)
        a, b = (ui + uj) % r, (vi + vj) % r
        if P is None:
            continue
        P, a, b = canon(P, a, b)
        hist = []
        for n in range(cap):
            if P is None:
                break
            if (P[0] & mask) == 0:
                prev = seen.get(P)
                if prev is not None and prev[2] != trail:
                    a2, b2, _ = prev
                    if (b - b2) % r:
                        x = (a2 - a) * pow(b - b2, -1, r) % r
                        if C.mul(x, G) == Q:
                            return x, C.ops - start_ops
                elif prev is None:
                    seen[P] = (a, b, trail)
                break
            if auts:
                # fruitless-cycle escape
                k = next((t for t in range(len(hist)) if hist[t][0] == P), None)
                if k is not None:
                    cyc = hist[k:]
                    P, a, b = min(cyc, key=lambda e: e[0])
                    u, v, R = table[(P[0] + 1) % S]
                    P = C.add(P, R)
                    a, b = (a + u) % r, (b + v) % r
                    hist = []
                    if P is None:
                        break
                    P, a, b = canon(P, a, b)
                    continue
                hist.append((P, a, b))
                if len(hist) > HIST:
                    hist.pop(0)
            u, v, R = table[P[0] % S]
            P = C.add(P, R)
            a, b = (a + u) % r, (b + v) % r
            if P is None:
                break
            P, a, b = canon(P, a, b)


def solve_ph(C: Curve, G, Q, n_fac, mode, rng):
    """Pohlig-Hellman.  n_fac: list of (prime, exponent) of ord(G).
    Returns dict(x, ops_total, ops_rho, r_max, steps_max)."""
    n = 1
    for q, e in n_fac:
        n *= q ** e
    residues, mods = [], []
    rho_ops = 0
    r_max, steps_max = 0, 0
    C.ops = 0
    for q, e in n_fac:
        qe = q ** e
        Gq = C.mul(n // qe, G)
        Qq = C.mul(n // qe, Q)
        Gq_pow = [Gq]                       # Gq * q^i helper
        # digit-by-digit lift
        x_q = 0
        Gq_qe1 = C.mul(q ** (e - 1), Gq)     # element of order q
        for i in range(e):
            T = C.add(Qq, C.neg(C.mul(x_q, Gq)))
            T = C.mul(q ** (e - 1 - i), T)
            d, st = rho_prime(C, Gq_qe1, T, q, mode, rng)
            rho_ops += st
            if q > r_max:
                r_max, steps_max = q, st
            x_q += d * q ** i
        residues.append(x_q % qe)
        mods.append(qe)
    x, M = 0, 1
    for r_, m_ in zip(residues, mods):          # incremental CRT, moduli pairwise coprime
        x += M * ((r_ - x) * pow(M, -1, m_) % m_)
        M *= m_
    x %= M
    assert C.mul(x, G) == Q
    return dict(x=x, ops_total=C.ops, ops_rho=rho_ops, r_max=r_max, steps_max=steps_max)


# ------------------------------------------------------------------ set-up
def pick_base_point(p, a, b, N, Nfac, rng, tries=8):
    """Deterministic (seeded) base point of near-maximal order.

    Returns (point, order, factorisation [(q, e), ...]); the point with the
    largest order among `tries` random points is used.
    """
    from volcano_descent import rand_point
    E = pari.ellinit([a, b], p)
    best = None
    for _ in range(tries):
        P = rand_point(E, p, rng)
        o = int(pari.ellorder(E, list(P), Nfac))
        if best is None or o > best[1]:
            best = (P, o)
    P, o = best
    fac = pari.factor(o)
    return P, o, [(int(fac[i, 0]), int(fac[i, 1])) for i in range(int(fac.nrows()))]
