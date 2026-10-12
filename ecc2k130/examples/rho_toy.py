#!/usr/bin/env python3
"""A readable Pollard rho for a tiny Koblitz curve.

This is the same iteration the ECC2K-130 GPU client runs, stripped of
bitslicing, CUDA, checkpoints and S3. It recovers a planted discrete
logarithm on GF(2^23) in under a second so a reader can watch the whole
pipeline: walk, distinguished point, collision, k, and [k]P = Q.

    python3 ecc2k130/examples/rho_toy.py

The production walk is in include/ref.h (scalar) and the packed CUDA
kernels. This file is the version meant to be read.
"""

from __future__ import annotations

import random
import sys

# y^2 + x y = x^3 + 1 over GF(2^23), f = t^23 + t^5 + 1.
# Prime subgroup order r is 21 bits: small enough to finish, large enough
# that the walk is not a single step.
M = 23
MOD = (1 << 23) | (1 << 5) | 1
DP_WEIGHT = 8
J_MIN = 3
N_BRANCHES = 8


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & (1 << M):
            a ^= MOD
    return r


def sqr(a):
    return mul(a, a)


def pow_field(a, e):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a)
        a = sqr(a)
        e >>= 1
    return r


def inv(a):
    return pow_field(a, (1 << M) - 2)


def frob(a, n=1):
    for _ in range(n % M):
        a = sqr(a)
    return a


def trace(a):
    t, cur = 0, a
    for _ in range(M):
        t ^= cur
        cur = sqr(cur)
    return t


def half_trace(c):
    z, cur = 0, c
    for _ in range((M - 1) // 2 + 1):
        z ^= cur
        cur = sqr(sqr(cur))
    return z


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 ^ y2 == x1:
            return None
        lam = x1 ^ mul(y1, inv(x1))
        x3 = sqr(lam) ^ lam
        y3 = sqr(x1) ^ mul(lam ^ 1, x3)
        return (x3, y3)
    lam = mul(y1 ^ y2, inv(x1 ^ x2))
    x3 = sqr(lam) ^ lam ^ x1 ^ x2
    y3 = mul(lam, x1 ^ x3) ^ x3 ^ y1
    return (x3, y3)


def neg(P):
    return None if P is None else (P[0], P[0] ^ P[1])


def scalar_mul(k, P):
    if k < 0:
        return scalar_mul(-k, neg(P))
    R = None
    while k:
        if k & 1:
            R = add(R, P)
        P = add(P, P)
        k >>= 1
    return R


def frob_pt(P, n=1):
    if P is None:
        return None
    return (frob(P[0], n), frob(P[1], n))


def lift_x(x):
    if x == 0:
        return None
    # y^2 + x y = x^3 + 1  =>  z^2 + z = x + 1/x^2, y = x z
    c = x ^ inv(sqr(x))
    if trace(c):
        return None
    return (x, mul(x, half_trace(c)))


def on_curve(P):
    if P is None:
        return True
    x, y = P
    return (sqr(y) ^ mul(x, y)) == (mul(sqr(x), x) ^ 1)


def group_order():
    t0, t1 = 2, -1
    for _ in range(M - 1):
        t0, t1 = t1, -t1 - 2 * t0
    return (1 << M) + 1 - t1


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m_, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        b = pow(c, 1 << (m_ - i - 1), p)
        m_, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def gf2_rank(rows):
    rows = list(rows)
    rank = 0
    for col in range(M):
        piv = None
        for i in range(rank, len(rows)):
            if rows[i] >> col & 1:
                piv = i
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(len(rows)):
            if i != rank and (rows[i] >> col & 1):
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def normal_masks():
    gamma = None
    for cand in range(2, 1 << 12):
        rows, g = [], cand
        for _ in range(M):
            rows.append(g)
            g = sqr(g)
        if gf2_rank(rows) == M:
            gamma = cand
            break
    if gamma is None:
        raise RuntimeError("no normal element")
    masks, g = [], gamma
    for _ in range(M):
        mask = 0
        for k in range(M):
            if trace(mul(1 << k, g)):
                mask |= 1 << k
        masks.append(mask)
        g = sqr(g)
    return masks


def hw(x, masks):
    return sum(1 for mask in masks if bin(x & mask).count("1") & 1)


def j_of(P, masks):
    return ((hw(P[0], masks) // 2) % N_BRANCHES) + J_MIN


def step(P, masks):
    j = j_of(P, masks)
    return add(P, frob_pt(P, j)), j


def canonical_x(P):
    best, cur = P[0], P[0]
    for _ in range(1, M):
        cur = sqr(cur)
        if cur < best:
            best = cur
    return best


def setup(rng):
    n = group_order()
    cofactor = 4
    r = n // cofactor
    if n % cofactor or not is_prime(r):
        raise RuntimeError("unexpected group order %d" % n)
    while True:
        P = lift_x(rng.getrandbits(M))
        if P is None or not on_curve(P):
            continue
        G = scalar_mul(cofactor, P)
        if G is not None and on_curve(G) and scalar_mul(r, G) is None:
            break
    disc = sqrt_mod((1 - 8) % r, r)
    inv2 = pow(2, -1, r)
    s = None
    for cand in ((-1 + disc) * inv2 % r, (-1 - disc) * inv2 % r):
        if frob_pt(G) == scalar_mul(cand, G):
            s = cand
            break
    if s is None:
        raise RuntimeError("Frobenius eigenvalue mismatch")
    return r, G, s


def walk_to_dp(start, masks, max_iters=1 << 16):
    P = start
    counts = [0] * N_BRANCHES
    for it in range(max_iters):
        if hw(P[0], masks) <= DP_WEIGHT:
            return P, it, counts
        j = j_of(P, masks)
        counts[j - J_MIN] += 1
        P, _ = step(P, masks)
        if P is None:
            return None, it, counts
    return None, max_iters, counts


def multiplier(counts, s, r):
    mu = 1
    for j, n_j in enumerate(counts):
        if not n_j:
            continue
        factor = (1 + pow(s, j + J_MIN, r)) % r
        mu = mu * pow(factor, n_j, r) % r
    return mu


def recover_k(end_a, a_coeff, b_coeff, end_b, a_b, b_b, s, r, G, Q):
    for c in range(M):
        rot = frob_pt(end_b, c)
        if rot == end_a:
            sc = pow(s, c, r)
        elif rot == neg(end_a):
            sc = (-pow(s, c, r)) % r
        else:
            continue
        num = (a_coeff - sc * a_b) % r
        den = (sc * b_b - b_coeff) % r
        if den == 0:
            continue
        k = num * pow(den, -1, r) % r
        if scalar_mul(k, G) == Q:
            return k
    return None


def search(r, G, Q, s, masks, rng, walks=64):
    store = {}
    stats = {"steps": 0, "dps": 0, "walks": 0}
    for _ in range(walks):
        seed = rng.randrange(1, r)
        start = add(scalar_mul(seed, G), Q)  # a=seed, b=1
        stats["walks"] += 1
        end, steps, counts = walk_to_dp(start, masks)
        stats["steps"] += steps
        if end is None:
            continue
        stats["dps"] += 1
        mu = multiplier(counts, s, r)
        a, b = mu * seed % r, mu % r
        key = canonical_x(end)
        prev = store.get(key)
        if prev is None:
            store[key] = (end, a, b, seed, steps)
            continue
        if prev[3] == seed:
            continue
        k = recover_k(prev[0], prev[1], prev[2], end, a, b, s, r, G, Q)
        if k is not None:
            return k, stats, (prev[3], seed, key, prev[4], steps)
    return None, stats, None


def main():
    rng = random.Random(1)
    r, G, s = setup(rng)
    planted = rng.randrange(1, r)
    Q = scalar_mul(planted, G)
    masks = normal_masks()
    print("curve  y^2 + xy = x^3 + 1  over GF(2^%d)" % M)
    print("prime subgroup order r = %d (~2^%d)" % (r, r.bit_length()))
    print("walk    R <- R + tau^j(R),  j = 3 + ((HW(x)/2) mod 8)")
    print("distinguish HW(x) <= %d, same rule as ECC2K-130" % DP_WEIGHT)
    print("planted k = %d" % planted)
    k, stats, hit = search(r, G, Q, s, masks, rng)
    print("walked %d steps across %d seeds, %d distinguished points"
          % (stats["steps"], stats["walks"], stats["dps"]))
    if k is None:
        print("no collision this trial", file=sys.stderr)
        return 1
    seed_a, seed_b, key, steps_a, steps_b = hit
    print("collision on orbit %d" % key)
    print("  walk A seed %d reached it in %d steps" % (seed_a, steps_a))
    print("  walk B seed %d reached it in %d steps" % (seed_b, steps_b))
    print("recovered k = %d" % k)
    print("verified [k]P = Q" if k == planted else "MISMATCH", )
    return 0 if k == planted else 1


if __name__ == "__main__":
    raise SystemExit(main())
