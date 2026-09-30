#!/usr/bin/env python3
"""Cut-out degree of the m = 4 tree-split set, and controls (PREREGISTRATION.md).

Objects, on K_a : y^2 + xy = x^3 + a x^2 + 1 over GF(2^n):
  F   factor base: points with x in V \\ {0}, V a random dim-l subspace (seeded)
  T2  {x(P+Q) : P, Q in F, P+Q != O}                                    subset of GF(2)^n
  Z   {(e+f, e*f) : e, f in T2} (unordered, e = f allowed)               subset of GF(2)^{2n}
Controls of the same size and ambient dimension:
  R   uniform random subset
  ZR  (for Z only) {(e+f, e*f) : e, f in a uniform random set of |T2| field elements}

Statistics per set S in GF(2)^N:
  h_S(D)     GF(2) rank of the evaluation matrix of multilinear monomials of degree <= D
  deficit(D) = min(|S|, M(D)) - h_S(D)      (0 for a generic set)
  D_s        sampled cut-out degree: smallest D <= D_MAX at which none of K points drawn
             uniformly from GF(2)^N \\ S has its monomial vector in the row space of S's
             (i.e. is a common zero of every degree-<=D polynomial vanishing on S)
"""

import argparse
import itertools
import json
import random
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- GF(2^n)


def is_irreducible(poly, n):
    # Ben-Or / Rabin on small n: x^(2^i) mod f for i <= n/2 has gcd 1 with x^(2^i)-x
    def mulmod(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> n & 1:
                a ^= poly
        return r

    def gcd(a, b):
        while b:
            while a and a.bit_length() >= b.bit_length():
                a ^= b << (a.bit_length() - b.bit_length())
            a, b = b, a
        return a

    x = 2
    t = x
    for _ in range(1, n // 2 + 1):
        t = mulmod(t, t)
        if gcd(poly, t ^ x) != 1:
            return False
    return True


def find_irreducible(n):
    low = [1 << k for k in range(1, n)]
    for k in range(1, n):
        p = (1 << n) | (1 << k) | 1
        if is_irreducible(p, n):
            return p
    for a, b, c in itertools.combinations(range(1, n), 3):
        p = (1 << n) | (1 << a) | (1 << b) | (1 << c) | 1
        if is_irreducible(p, n):
            return p
    raise ValueError('no irreducible found')


class Field:
    def __init__(self, n):
        self.n, self.poly = n, find_irreducible(n)

    def mul(self, a, b):
        r, n, poly = 0, self.n, self.poly
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> n & 1:
                a ^= poly
        return r

    def sqr(self, a):
        return self.mul(a, a)

    def inv(self, a):
        # a^(2^n - 2)
        r, e, base = 1, (1 << self.n) - 2, a
        while e:
            if e & 1:
                r = self.mul(r, base)
            base = self.mul(base, base)
            e >>= 1
        return r

    def trace(self, a):
        t, x = 0, a
        for _ in range(self.n):
            t ^= x
            x = self.sqr(x)
        return t & 1 if t in (0, 1) else None

    def half_trace(self, c):  # n odd: z with z^2 + z = c when Tr(c) = 0
        z, x = 0, c
        for i in range(0, self.n, 2):
            z ^= x
            x = self.sqr(self.sqr(x))
        return z


# ---------------------------------------------------------------- Koblitz curve


class Curve:
    def __init__(self, a, n):
        self.F, self.a = Field(n), a

    def points_with_x(self, x):
        F = self.F
        if x == 0:
            return []
        # y = x z:  z^2 + z = x + a + 1/x^2
        c = x ^ self.a ^ F.inv(F.sqr(x))
        if F.trace(c) != 0:
            return []
        z = F.half_trace(c)
        y = F.mul(x, z)
        return [(x, y), (x, y ^ x)]

    def add(self, P, Q):
        F = self.F
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if y1 ^ y2 == x1 or x1 == 0:  # Q = -P  (-(x,y) = (x, x+y))
                return None
            lam = x1 ^ F.mul(y1, F.inv(x1))
            x3 = F.sqr(lam) ^ lam ^ self.a
            y3 = F.sqr(x1) ^ F.mul(lam ^ 1, x3)
            return (x3, y3)
        lam = F.mul(y1 ^ y2, F.inv(x1 ^ x2))
        x3 = F.sqr(lam) ^ lam ^ x1 ^ x2 ^ self.a
        y3 = F.mul(lam, x1 ^ x3) ^ x3 ^ y1
        return (x3, y3)

    def on_curve(self, P):
        F = self.F
        x, y = P
        return F.sqr(y) ^ F.mul(x, y) == F.mul(F.sqr(x), x) ^ F.mul(self.a, F.sqr(x)) ^ 1


def random_subspace(n, ell, rng):
    basis, span = [], {0}
    while len(basis) < ell:
        v = rng.getrandbits(n)
        if v not in span:
            basis.append(v)
            span |= {s ^ v for s in span}
    return sorted(span)


def objects(a, n, ell, seed):
    rng = random.Random(seed ^ (a << 40) ^ (n << 32) ^ (ell << 24))
    E = Curve(a, n)
    V = random_subspace(n, ell, rng)
    base = [P for x in V for P in E.points_with_x(x)]
    assert all(E.on_curve(P) for P in base)
    T2 = set()
    for P, Q in itertools.combinations_with_replacement(base, 2):
        R = E.add(P, Q)
        if R is not None:
            T2.add(R[0])
    T2 = sorted(T2)
    return E, rng, base, T2


def sym2(F, T):
    n = F.n
    out = set()
    for e, f in itertools.combinations_with_replacement(T, 2):
        out.add(((e ^ f) << n) | F.mul(e, f))
    return sorted(out)


# ---------------------------------------------------------------- GF(2) linear algebra


def monomials(N, D):
    return [m for d in range(D + 1) for m in itertools.combinations(range(N), d)]


def eval_matrix(points, N, monos):
    """Packed bit matrix: row per point, bit j = monomial j evaluated at the point."""
    P = np.array(points, dtype=np.uint64)
    bits = [(P >> np.uint64(i)) & np.uint64(1) for i in range(N)]
    W = (len(monos) + 63) // 64
    out = np.zeros((len(points), W), dtype=np.uint64)
    for j, m in enumerate(monos):
        col = np.ones(len(points), dtype=np.uint64)
        for i in m:
            col &= bits[i]
        out[:, j >> 6] |= col << np.uint64(j & 63)
    return out


SPAN_BIN = Path(__file__).resolve().parents[2] / 'target/release/examples/gf2_span'


def rank_and_span(A, samples):
    """(GF(2) rank of A's rows, number of sample rows in their span), by examples/gf2_span.rs."""
    head = np.array([A.shape[0], samples.shape[0], A.shape[1]], dtype='<u8')
    blob = head.tobytes() + A.astype('<u8').tobytes() + samples.astype('<u8').tobytes()
    out = subprocess.run([str(SPAN_BIN)], input=blob, capture_output=True, check=True).stdout
    r = json.loads(out)
    return r['rank'], r['in_span']


def profile(points, N, d_max, k_samples, rng, max_cols):
    s = len(points)
    pset = set(points)
    samples = []
    while len(samples) < k_samples:
        p = rng.getrandbits(N)
        if p not in pset:
            samples.append(p)
    rows, cut = [], None
    for D in range(1, d_max + 1):
        monos = monomials(N, D)
        if len(monos) > max_cols:
            rows.append({'D': D, 'skipped': f'{len(monos)} monomials > max_cols'})
            break
        t0 = time.time()
        h, false_zeros = rank_and_span(eval_matrix(points, N, monos), eval_matrix(samples, N, monos))
        rows.append({'D': D, 'M': len(monos), 'rank': h, 'deficit': min(s, len(monos)) - h,
                     'false_zeros': false_zeros, 'samples': k_samples, 'secs': round(time.time() - t0, 2)})
        if false_zeros == 0 and cut is None:
            cut = D
            break
    return {'size': s, 'N': N, 'D_s': cut, 'rows': rows}


# ---------------------------------------------------------------- driver


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('--a', type=int, required=True)
    p.add_argument('--n', type=int, required=True)
    p.add_argument('--ell', type=int, required=True)
    p.add_argument('--object', choices=['T2', 'Z'], required=True)
    p.add_argument('--seed', type=int, default=20260930)
    p.add_argument('--d-max', type=int, default=6)
    p.add_argument('--samples', type=int, default=2000)
    p.add_argument('--max-cols', type=int, default=90000)
    args = p.parse_args(argv)
    E, rng, base, T2 = objects(args.a, args.n, args.ell, args.seed)
    F, n = E.F, args.n
    out = {'a': args.a, 'n': n, 'ell': args.ell, 'object': args.object, 'seed': args.seed,
           'field_poly': F.poly, 'factor_base_points': len(base), 'T2_size': len(T2)}
    srng = random.Random(args.seed + 1)
    if args.object == 'T2':
        arms = {'curve': (T2, n),
                'R': (sorted(random.Random(args.seed + 2).sample(range(1, 1 << n), len(T2))), n)}
    else:
        Z = sym2(F, T2)
        rand_T = random.Random(args.seed + 3).sample(range(1, 1 << n), len(T2))
        ZR = sym2(F, rand_T)
        R = sorted(random.Random(args.seed + 2).sample(range(1, 1 << (2 * n)), len(Z)))
        arms = {'curve': (Z, 2 * n), 'ZR': (ZR, 2 * n), 'R': (R, 2 * n)}
    for name, (pts, N) in arms.items():
        out[name] = profile(pts, N, args.d_max, args.samples, srng, args.max_cols)
    print(json.dumps(out, sort_keys=True))


if __name__ == '__main__':
    sys.exit(main())
