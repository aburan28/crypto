#!/usr/bin/env python3
"""Can a point of ECC2K-130 be written as a sum of factor-base points?

The index-calculus question, at the challenge parameters: for a target
`R` in the order-`r` subgroup of `K_0 : y^2 + xy = x^3 + 1` over `F_2^131`,
and a factor base `F = {P : x(P) in V}` over an `F_2`-subspace `V` of
dimension `l`, does `R = P_1 + ... + P_m` with every `P_i` in `F`?

Three separate things get answered separately, because they have different
answers:

  1. do such decompositions EXIST      -- yes, and the yield law is measured
  2. are they ADMISSIBLE at n = 131    -- yes; the cofactor-4 class is not an
                                          obstruction, measured on the curve
  3. can they be FOUND                 -- this is the whole difficulty, and
                                          the cost surface is derived here

The boundary every cost is quoted against is Pollard rho with the
`<-1> x <pi>` speed-up on the same subgroup, `2^60.8090`, re-derived here and
matching `experiments/ecc2k130_extension_field_boundary.json` to the digit.

    python3 scripts/ecc2k130_point_decomposition.py [--targets N] [--sample N]

Writes `experiments/ecc2k130_point_decomposition.json`.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "experiments/ecc2k130_point_decomposition.json"

N131 = 131
# The challenge reduction polynomial, x^131 + x^13 + x^2 + x + 1 (checked below).
IRR131 = (1 << 131) | (1 << 13) | (1 << 2) | (1 << 1) | 1
SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


# ── finite field F_2^n ────────────────────────────────────────────────────

class GF2m:
    """`F_2[z]/(irr)` on Python ints, with optional log tables for small n."""

    def __init__(self, deg: int, irr: int, tables: bool = False):
        self.deg, self.irr, self.mask = deg, irr, (1 << deg) - 1
        # the reduction fold: irr = z^deg + tail, so z^deg = tail
        self.tail = irr ^ (1 << deg)
        self.taps = [i for i in range(deg) if (self.tail >> i) & 1]
        trace_bits = 0
        for i in range(deg):
            v, s = 1 << i, 0
            for _ in range(deg):
                s ^= v
                v = self.sqr(v)
            assert s in (0, 1), "trace must land in F_2"
            trace_bits |= (s & 1) << i
        self.trace_bits = trace_bits
        self.exp = self.log = None
        if tables:
            self._build_tables()

    # -- reduction: fold the high part down through the tail, to fixpoint
    def _reduce(self, a: int) -> int:
        while a > self.mask:
            hi = a >> self.deg
            a = (a & self.mask) ^ self._spread(hi)
        return a

    def _spread(self, h: int) -> int:
        r = 0
        for t in self.taps:
            r ^= h << t
        return r

    def mul(self, a: int, b: int) -> int:
        if self.log is not None:
            if a == 0 or b == 0:
                return 0
            return self.exp[self.log[a] + self.log[b]]
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if (a >> self.deg) & 1:
                a ^= self.irr
        return r

    def sqr(self, a: int) -> int:
        s = 0
        i = 0
        while a:
            if a & 1:
                s |= 1 << (2 * i)
            a >>= 1
            i += 1
        return self._reduce(s)

    def inv(self, a: int) -> int:
        assert a, "zero has no inverse"
        if self.log is not None:
            return self.exp[(self.mask - self.log[a]) % self.mask]
        r, base, e = 1, a, (1 << self.deg) - 2
        while e:
            if e & 1:
                r = self.mul(r, base)
            base = self.sqr(base)
            e >>= 1
        return r

    def sqrt(self, a: int) -> int:
        for _ in range(self.deg - 1):
            a = self.sqr(a)
        return a

    def trace(self, a: int) -> int:
        return bin(a & self.trace_bits).count("1") & 1

    def half_trace(self, c: int) -> int:
        """`t` with `t^2 + t = c`, for odd deg and `Tr(c) = 0`."""
        assert self.deg % 2 == 1
        t, acc = c, c
        for _ in range((self.deg - 1) // 2):
            acc = self.sqr(self.sqr(acc))
            t ^= acc
        return t

    def solve_artin_schreier(self, c: int):
        if self.trace(c):
            return None
        return self.half_trace(c)

    def _build_tables(self):
        g = self._find_generator()
        exp = [0] * (2 * self.mask)
        log = [0] * (self.mask + 1)
        x = 1
        for i in range(self.mask):
            exp[i] = x
            log[x] = i
            x = self.mul(x, g)
        assert x == 1, "generator did not have full order"
        for i in range(self.mask, 2 * self.mask):
            exp[i] = exp[i - self.mask]
        self.exp, self.log = exp, log

    def _find_generator(self) -> int:
        order = self.mask
        facs = sorted(set(prime_factors(order)))
        for cand in range(2, 1 << self.deg):
            if all(self.pow(cand, order // p) != 1 for p in facs):
                return cand
        raise AssertionError("no generator")

    def pow(self, a: int, e: int) -> int:
        r, base = 1, a
        while e:
            if e & 1:
                r = self.mul(r, base)
            base = self.mul(base, base)
            e >>= 1
        return r


def prime_factors(n: int):
    out, d = [], 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin over the first twelve primes."""
    if n < 2:
        return False
    for p in SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in SMALL_PRIMES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def poly_is_irreducible(f: int, deg: int) -> bool:
    """Rabin's test for `f` over `F_2`."""
    def mulmod(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a.bit_length() - 1 >= deg:
                a ^= f << (a.bit_length() - 1 - deg)
        return r

    def powx(k):                      # x^(2^k) mod f
        h = 2
        for _ in range(k):
            h = mulmod(h, h)
        return h

    def gcd(a, b):
        while b:
            while a and a.bit_length() >= b.bit_length():
                a ^= b << (a.bit_length() - b.bit_length())
            a, b = b, a
        return a

    if powx(deg) != 2:
        return False
    return all(gcd(powx(deg // p) ^ 2, f) == 1 for p in set(prime_factors(deg)))


def find_irreducible(deg: int) -> int:
    for extra in range(1, 1 << min(deg, 20)):
        f = (1 << deg) | (extra << 1) | 1
        if f.bit_length() - 1 == deg and poly_is_irreducible(f, deg):
            return f
    raise AssertionError(f"no irreducible of degree {deg}")


# ── the Koblitz curve y^2 + xy = x^3 + a x^2 + b, a = 0, b = 1 ────────────

class Koblitz:
    def __init__(self, field: GF2m, b: int = 1):
        self.F, self.b = field, b
        self.sqrt_b = field.sqrt(b)

    def points_over(self, x: int):
        """The zero, one or two affine points with abscissa `x`."""
        F = self.F
        if x == 0:
            return [(0, self.sqrt_b)]
        c = x ^ F.mul(self.b, F.sqr(F.inv(x)))       # x + b/x^2
        z = F.solve_artin_schreier(c)
        if z is None:
            return []
        return [(x, F.mul(x, z)), (x, F.mul(x, z ^ 1))]

    def neg(self, P):
        return None if P is None else (P[0], P[1] ^ P[0])

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        F = self.F
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if y2 == (y1 ^ x1):                       # Q = -P
                return None
            if x1 == 0:
                return None
            lam = x1 ^ F.mul(y1, F.inv(x1))
            x3 = F.sqr(lam) ^ lam
            return (x3, F.sqr(x1) ^ F.mul(lam ^ 1, x3))
        lam = F.mul(y1 ^ y2, F.inv(x1 ^ x2))
        x3 = F.sqr(lam) ^ lam ^ x1 ^ x2
        return (x3, F.mul(lam, x1 ^ x3) ^ x3 ^ y1)

    def mul(self, P, k: int):
        R, Q = None, P
        while k:
            if k & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

    def on_curve(self, P) -> bool:
        if P is None:
            return True
        F = self.F
        x, y = P
        lhs = F.sqr(y) ^ F.mul(x, y)
        rhs = F.mul(x, F.sqr(x)) ^ self.b
        return lhs == rhs

    # -- Semaev's fourth summation polynomial, b as above
    def s3_in_last(self, x1: int, x2: int):
        F = self.F
        return (F.sqr(x1 ^ x2), F.mul(x1, x2), F.sqr(F.mul(x1, x2)) ^ self.b)

    def s4(self, x1: int, x2: int, x3: int, x4: int) -> int:
        F = self.F
        a1, b1, c1 = self.s3_in_last(x1, x2)
        a2, b2, c2 = self.s3_in_last(x3, x4)
        t1 = F.sqr(F.mul(a1, c2) ^ F.mul(a2, c1))
        t2 = F.mul(F.mul(a1, b2) ^ F.mul(a2, b1), F.mul(b1, c2) ^ F.mul(b2, c1))
        return t1 ^ t2

    # -- lockstep scalar multiplication with Montgomery batch inversion
    def batch_mul(self, points, k: int):
        acc = [None] * len(points)
        cur = list(points)
        while k:
            if k & 1:
                acc = self._batch_add(acc, cur)
            k >>= 1
            if k:
                cur = self._batch_add(cur, cur)
        return acc

    def _batch_add(self, Ps, Qs):
        F = self.F
        dens, slots = [], []
        out = [None] * len(Ps)
        for i, (P, Q) in enumerate(zip(Ps, Qs)):
            if P is None or Q is None:
                out[i] = Q if P is None else P
                continue
            x1, y1 = P
            x2, y2 = Q
            if x1 == x2:
                if y2 == (y1 ^ x1) or x1 == 0:
                    out[i] = None
                    continue
                dens.append(x1)
            else:
                dens.append(x1 ^ x2)
            slots.append(i)
        for inv, i in zip(self._batch_inv(dens), slots):
            x1, y1 = Ps[i]
            x2, y2 = Qs[i]
            if x1 == x2:
                lam = x1 ^ F.mul(y1, inv)
                x3 = F.sqr(lam) ^ lam
                out[i] = (x3, F.sqr(x1) ^ F.mul(lam ^ 1, x3))
            else:
                lam = F.mul(y1 ^ y2, inv)
                x3 = F.sqr(lam) ^ lam ^ x1 ^ x2
                out[i] = (x3, F.mul(lam, x1 ^ x3) ^ x3 ^ y1)
        return out

    def _batch_inv(self, xs):
        F = self.F
        if not xs:
            return []
        pref = [1] * (len(xs) + 1)
        for i, v in enumerate(xs):
            pref[i + 1] = F.mul(pref[i], v)
        run = F.inv(pref[-1])
        out = [0] * len(xs)
        for i in range(len(xs) - 1, -1, -1):
            out[i] = F.mul(run, pref[i])
            run = F.mul(run, xs[i])
        return out


def curve_order(n: int, trace: int = -1, q: int = 2) -> int:
    """`#E(F_q^n)` from the Koblitz recurrence `s_i = t s_{i-1} - q s_{i-2}`."""
    s0, s1 = 2, trace
    for _ in range(n - 1):
        s0, s1 = s1, trace * s1 - q * s0
    return q ** n + 1 - s1


def subspace(field: GF2m, l: int):
    """The `F_2`-span of `1, z, ..., z^(l-1)`: every `l`-bit pattern."""
    return list(range(1 << l))


def factor_base(curve: Koblitz, l: int):
    pts = []
    for x in subspace(curve.F, l):
        pts.extend(curve.points_over(x))
    return pts


# ── 1. does a decomposition exist?  the yield law, measured on toys ───────

def yield_law(base_size: int, m: int, group_order: int) -> float:
    """Expected number of `m`-subsets of the base summing to a fixed target.

    The `m`-subsets number `C(|F|, m)` and their sums are spread over the whole
    group, so a fixed `R` is hit `C(|F|, m) / #E` times on average -- which for
    `|F| = 2^l` is `2^(m l) / (m! 2^n)` to leading order.
    """
    return math.comb(base_size, m) / group_order


_TOY_CACHE: dict = {}


def toy_instance(n: int):
    """`K_0 / F_2^n` with log tables, every affine point, and the cofactor.

    The full point list doubles as a check on the Koblitz recurrence: the count
    must match `#E(F_2^n)` exactly.
    """
    if n not in _TOY_CACHE:
        field = GF2m(n, find_irreducible(n), tables=True)
        curve = Koblitz(field)
        order = curve_order(n)
        all_pts = [P for x in range(1 << n) for P in curve.points_over(x)]
        assert len(all_pts) + 1 == order, (len(all_pts) + 1, order)
        cofactor = 1 << ((order & -order).bit_length() - 1)
        _TOY_CACHE[n] = (curve, order, cofactor, all_pts)
    return _TOY_CACHE[n]


def measure_toy_yield(n: int, m: int, l: int, targets: int, rng: random.Random):
    """Measure the fraction of targets in `<G>` that decompose, against the law."""
    curve, order, cofactor, all_pts = toy_instance(n)
    base = factor_base(curve, l)
    assert all(curve.on_curve(P) for P in base)

    pool = [curve.mul(P, cofactor) for P in rng.sample(all_pts, min(len(all_pts), targets))]
    pool = [P for P in pool if P is not None]

    base_set = set(base)
    pair_sums = {}
    if m >= 3:
        for i in range(len(base)):
            for j in range(i + 1, len(base)):
                s = curve.add(base[i], base[j])
                if s is not None:
                    pair_sums.setdefault(s, (i, j))

    hits, witness = 0, None
    for R in pool:
        found = None
        if m == 2:
            for P in base:
                Q = curve.add(R, curve.neg(P))
                if Q in base_set:
                    found = (P, Q)
                    break
        elif m == 3:
            for P in base:
                rest = curve.add(R, curve.neg(P))
                if rest in pair_sums:
                    i, j = pair_sums[rest]
                    found = (P, base[i], base[j])
                    break
        else:                                        # m == 4, pairs against pairs
            for s, (i, j) in pair_sums.items():
                rest = curve.add(R, curve.neg(s))
                if rest in pair_sums:
                    k, t = pair_sums[rest]
                    found = (base[i], base[j], base[k], base[t])
                    break
        if found is not None:
            hits += 1
            if witness is None:
                acc = None
                for P in found:
                    acc = curve.add(acc, P)
                assert acc == R, "witness does not sum to the target"
                witness = [hex(P[0]) for P in found]

    lam = yield_law(len(base), m, order)
    return {
        "n": n, "m": m, "l": l,
        "two_is_primitive_mod_n": n % 2 == 1 and is_prime(n) and mult_order(2, n) == n - 1,
        "curve_order": order, "cofactor": cofactor,
        "factor_base_size": len(base),
        "targets": len(pool),
        "decomposed": hits,
        "rate_measured": round(hits / len(pool), 4),
        "rate_predicted": round(1 - math.exp(-lam), 4),
        "expected_decompositions_per_target": round(lam, 4),
        "ratio_measured_over_predicted": round((hits / len(pool)) / (1 - math.exp(-lam)), 3),
        "witness_abscissae": witness,
    }


def mult_order(a: int, mod: int) -> int:
    k, x = 1, a % mod
    while x != 1:
        x = x * a % mod
        k += 1
    return k


# ── 2. admissibility and an existence witness at n = 131 ─────────────────

def sample_factor_base(curve: Koblitz, l: int, shift: int, count: int,
                       rng: random.Random):
    """An unbiased sample of `F = {P : x(P) in V}` for `V = span{z^shift, ...}`.

    The whole base has `~2^l` points and one inversion apiece, which at `n = 131`
    is minutes for the dimensions that matter, so draw abscissae instead of
    enumerating them.  The hit rate returned is the measured density of the base
    in the subspace and must come out near one half.
    """
    pts, tried = [], 0
    while len(pts) < count:
        tried += 1
        pts.extend(curve.points_over(rng.randrange(1, 1 << l) << shift))
    return pts[:count], len(pts) / (2 * tried)


def classes_at_131(curve: Koblitz, r: int, pts, rng: random.Random):
    """Histogram of the class of each sampled base point in `E / <G>` (order 4)."""
    cls = curve.batch_mul(pts, r)
    gen = next((c for c in cls if c is not None and curve.add(c, c) is not None), None)
    if gen is None:                                   # only 0 and the 2-torsion occur
        two = next((c for c in cls if c is not None), None)
        table = {None: 0} if two is None else {None: 0, two: 2}
    else:
        table = {None: 0, gen: 1, curve.add(gen, gen): 2,
                 curve.add(curve.add(gen, gen), gen): 3}
    idx = [table[c] for c in cls]
    parity_is_trace = all((idx[i] & 1) == curve.F.trace(pts[i][0]) for i in range(len(pts)))
    hist = {str(c): idx.count(c) for c in range(4)}
    present = sorted({c for c in range(4) if idx.count(c)})
    return hist, present, parity_is_trace, idx


def admissible_counts(present) -> list:
    """Which summand counts `m` can cancel the cofactor class, given the classes present."""
    out = []
    for m in range(2, 9):
        reach = {0}
        for _ in range(m):
            reach = {(a + c) % 4 for a in reach for c in present}
        if 0 in reach:
            out.append(m)
    return out


def witness_at_131(curve: Koblitz, r: int, l: int, m: int, rng: random.Random, tries: int = 400):
    """Plant a decomposition: `m` base points whose sum lies in `<G>`, certified."""
    F = curve.F
    for _ in range(tries):
        pts = []
        while len(pts) < m:
            cands = curve.points_over(rng.randrange(1, 1 << l))
            if cands:
                pts.append(cands[rng.randrange(len(cands))])
        R = None
        for P in pts:
            R = curve.add(R, P)
        if R is None or curve.mul(R, r) is not None:
            continue                                  # sum left <G>: classes did not cancel
        rec = {
            "m": m, "subspace_dim": l,
            "summand_abscissae": [hex(P[0]) for P in pts],
            "summand_ordinates": [hex(P[1]) for P in pts],
            "sum_x": hex(R[0]), "sum_y": hex(R[1]),
            "sum_is_on_curve": curve.on_curve(R),
            "sum_lies_in_subgroup_of_order_r": curve.mul(R, r) is None,
            "every_summand_abscissa_in_subspace": all(P[0] < (1 << l) for P in pts),
        }
        if m == 3:
            # Semaev's S_4 must vanish on the three abscissae and the target's
            rec["semaev_s4_at_summands_and_target"] = hex(
                curve.s4(pts[0][0], pts[1][0], pts[2][0], R[0]))
            rec["semaev_s4_vanishes"] = curve.s4(
                pts[0][0], pts[1][0], pts[2][0], R[0]) == 0
        return rec
    return None


# ── 3. can a decomposition be found?  the cost surface ───────────────────

def log2_comb(n_log2: float, k: int) -> float:
    """`log2 C(2^n, k)` to leading order: `k n - log2 k!`."""
    return k * n_log2 - math.log2(math.factorial(k))


def log2_add(*terms: float) -> float:
    """`log2` of a sum of powers of two."""
    hi = max(terms)
    return hi + math.log2(sum(2.0 ** (t - hi) for t in terms))


def cost_cell(m: int, l: float, n: int, split: int = 1, frobenius: bool = False):
    """One relation-collection-plus-linear-algebra run, in log2 operations.

    `split = v` tabulates the sums of every `v`-subset of the base once and
    enumerates the remaining `m - v` per target.  `v = 1` is the oracle this
    repository has actually built and measured: walk `C(|F|, m-1)` sub-tuples and
    root-find the last summand inside the subspace (`RESEARCH_SEMAEV_DECOMPOSITION.md`,
    "Don't enumerate triples"), with only `O(|F|)` memory.  `v >= 2` buys speed
    with a table, and `v = m` tabulates every `m`-subset sum -- at which point the
    table is a baby-step table and the method is a generic algorithm in costume.

    `frobenius` makes the base a union of `pi`-orbits.  GGMP's collapse needs only
    `pi(F) = F`, not that `F` be a subspace, so it is available at any `n` and any
    size: `|F|/n` unknowns, `|F|/n` relations, an `|F|/n`-square matrix.  What it
    costs is the subspace root-find -- an orbit union has no low-degree membership
    polynomial -- so the base has to be materialised and the last summand looked
    up instead, `2^l` stored points for the same `C(|F|, m-1)` per target.
    """
    assert 1 <= split <= m
    log_subsets = log2_comb(l, m)                      # C(|F|, m)
    log_yield = log_subsets - n                        # decompositions per target
    log_targets_per_relation = max(0.0, -log_yield)
    log_setup = log2_comb(l, split) if split > 1 else -math.inf
    log_oracle = log2_comb(l, m - split) if m > split else 0.0
    collapse = math.log2(n) if frobenius else 0.0      # |F| -> |F|/n unknowns

    log_collection = (l - collapse) + log_targets_per_relation + log_oracle
    log_linalg = math.log2(m) + 2 * (l - collapse)     # sparse Wiedemann, m per row
    store = max(log_setup, l) if frobenius else log_setup
    terms = [log_collection, log_linalg] + ([store] if store > -math.inf else [])
    return {
        "m": m, "l": round(l, 2), "split": split, "frobenius_stable": frobenius,
        "log2_factor_base": round(l, 2),
        "log2_decompositions_per_target": round(log_yield, 2),
        "log2_targets_per_relation": round(log_targets_per_relation, 2),
        "log2_oracle_per_target": round(log_oracle, 2),
        "log2_table_entries": None if store == -math.inf else round(store, 2),
        "log2_collection": round(log_collection, 2),
        "log2_linear_algebra": round(log_linalg, 2),
        "log2_total": round(log2_add(*terms), 2),
    }


def saturating_l(m: int, n: int) -> float:
    """The smallest subspace dimension at which a random target decomposes.

    `C(2^l, m) = 2^n`, i.e. `l = (n + log2 m!) / m`: below it targets have to be
    resampled, above it the factor base is larger than the method needs.
    """
    return (n + math.log2(math.factorial(m))) / m


def cost_row(m: int, n: int, split: int = 1, l_min: float = 8.0):
    """The best cell for one `(m, split)`, searched over `l`."""
    best = None
    for step in range(int(l_min * 4), 4 * n + 1):
        cell = cost_cell(m, step / 4.0, n, split)
        if best is None or cell["log2_total"] < best["log2_total"]:
            best = cell
    return best


def memory_capped_row(m: int, n: int, split: int, log2_entries: float, l_min: float = 8.0):
    """The best cell for one `(m, split)` whose table fits in `2^entries` slots."""
    best = None
    for step in range(int(l_min * 4), 4 * n + 1):
        cell = cost_cell(m, step / 4.0, n, split)
        if cell["log2_table_entries"] is not None and cell["log2_table_entries"] > log2_entries:
            continue
        if best is None or cell["log2_total"] < best["log2_total"]:
            best = cell
    return best


def free_oracle_floor(m: int, n: int):
    """Lowest cost any `m`-decomposition method can reach if the oracle were free.

    Charges one operation per target tried and `m 2^(2l)` for the linear algebra,
    and nothing at all for deciding a decomposition.  Nothing below this line is
    reachable however good the algebra gets.
    """
    best = None
    for step in range(1, 40001):
        l = step / 100.0
        if l > n:
            break
        log_targets = l + max(0.0, n - log2_comb(l, m))
        log_linalg = math.log2(m) + 2 * l
        tot = log2_add(log_targets, log_linalg)
        if best is None or tot < best[0]:
            best = (tot, l, log_targets, log_linalg)
    tot, l, lt, la = best
    return {"m": m, "l_star": round(l, 2), "log2_floor": round(tot, 2),
            "log2_targets_at_l_star": round(lt, 2),
            "log2_linear_algebra_at_l_star": round(la, 2)}


def oracle_budget(m: int, n: int, log2_rho: float):
    """How cheap the decomposition oracle would have to be to reach the rho line.

    At the `l` that minimises the free-oracle floor, an oracle costing `2^w` per
    target makes the whole run cost `2^(floor + w)` at best; solve for `w`.
    """
    fl = free_oracle_floor(m, n)
    l = fl["l_star"]
    log_targets = l + max(0.0, n - log2_comb(l, m))
    budget = log2_rho - log_targets
    search = log2_comb(l, m - 1)
    return {"m": m, "l_star": l,
            "log2_search_space_per_target": round(search, 2),
            "log2_oracle_budget_to_match_rho": round(budget, 2),
            "oracle_budget_operations": (None if budget > 64 else
                                         round(2 ** budget, 2) if budget > -20 else 0.0),
            "log2_speedup_over_search_required": round(budget - search, 2),
            "reachable": budget >= search}


def generic_with_memory(log2_entries: float, log2_r: float):
    """What a plain generic algorithm does with the same table.

    Baby-step giant-step on a cyclic group of order `r` with `2^mu` stored steps
    costs `max(2^mu, r/2^mu)` group operations, bottoming out at `sqrt(r)`.  Any
    decomposition method that only beats rho by spending memory has to beat this
    too, and it is the line the tabulated oracles are really competing against.
    """
    return {"log2_table_entries": round(log2_entries, 2),
            "log2_bsgs": round(max(log2_entries, log2_r - log2_entries), 2)}


def cyclotomic_coset_sizes(n: int):
    seen, sizes = set(), []
    for a in range(n):
        if a in seen:
            continue
        coset, x = set(), a
        while x not in coset:
            coset.add(x)
            x = x * 2 % n
        seen |= coset
        sizes.append(len(coset))
    return sorted(sizes)


def available_subspace_dimensions(n: int):
    """Dimensions of the Frobenius-stable `F_2`-subspaces of `F_2^n`: subset sums
    of the 2-cyclotomic coset sizes mod `n`."""
    dims = {0}
    for s in cyclotomic_coset_sizes(n):
        dims |= {d + s for d in dims}
    return sorted(dims)


# ── report ───────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=int, default=256,
                    help="targets per toy rung in the yield measurement")
    ap.add_argument("--sample", type=int, default=256,
                    help="factor-base points sampled for the n = 131 class histogram")
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    # ---- the target, and the boundary every cost is quoted against ---------
    assert poly_is_irreducible(IRR131, N131), "the challenge reduction polynomial"
    order = curve_order(N131)
    r = order // 4
    assert order % 4 == 0 and is_prime(r)
    log2r = math.log2(r)
    aut = 2 * N131                                     # <-1> x <pi> acting on <G>
    log2_rho = math.log2(math.sqrt(math.pi * r / (2 * aut)))

    field = GF2m(N131, IRR131)
    curve = Koblitz(field)
    target = {
        "curve": "K_0 : y^2 + xy = x^3 + 1 over F_2, used over F_2^131",
        "reduction_polynomial": "x^131 + x^13 + x^2 + x + 1",
        "order": str(order), "cofactor": 4, "r": str(r), "r_is_prime": True,
        "log2_r": round(log2r, 4),
        "group_structure": "cyclic Z/4r: the only rational 2-torsion point is (0, sqrt(b))",
        "rho_automorphism_order": aut,
        "log2_rho_reference": round(log2_rho, 4),
        "S_rho_reference": round(2 ** (log2_rho - log2r / 2), 6),
    }

    # ---- 1. existence: the yield law, measured on toy rungs ----------------
    toy_cells = []
    for n, plan in [(11, [(2, 6), (3, 4), (4, 4)]),
                    (13, [(2, 7), (3, 5), (4, 4)]),
                    (17, [(2, 9), (3, 6), (4, 5)]),
                    (19, [(2, 10), (3, 7), (4, 5)])]:
        for m, l in plan:
            toy_cells.append(measure_toy_yield(n, m, l, args.targets, rng))
            c = toy_cells[-1]
            print(f"  yield  n={n:3d} m={m} l={l:2d}  |F|={c['factor_base_size']:5d}  "
                  f"measured {c['rate_measured']:.3f}  predicted {c['rate_predicted']:.3f}  "
                  f"ratio {c['ratio_measured_over_predicted']:.2f}", flush=True)
    ratios = [c["ratio_measured_over_predicted"] for c in toy_cells]
    existence = {
        "law": "E[# m-subsets of F summing to a fixed R] = C(|F|, m) / #E "
               "~ 2^(m l) / (m! 2^n)",
        "consequence_at_131": "a decomposition exists for a typical target as soon as "
                              "m l >= 131 + log2 m!",
        "toy_rungs": toy_cells,
        "ratio_measured_over_predicted_min": round(min(ratios), 3),
        "ratio_measured_over_predicted_max": round(max(ratios), 3),
        "saturating_dimension_at_131": {str(m): round(saturating_l(m, N131), 2)
                                        for m in range(2, 9)},
    }

    # ---- 2. admissibility at n = 131, on the challenge curve ---------------
    admissibility = []
    for name, l, shift in [("span{1, z, ..., z^44}", 45, 0),
                           ("span{z, z^2, ..., z^45}, inside ker Tr", 45, 1)]:
        pts, density = sample_factor_base(curve, l, shift, args.sample, rng)
        hist, present, parity, idx = classes_at_131(curve, r, pts, rng)
        adm = admissible_counts(present)
        admissibility.append({
            "subspace": name, "dim": l,
            "measured_point_density_in_subspace": round(density, 4),
            "sampled": len(pts),
            "class_histogram_in_E_over_G": hist,
            "classes_present": present,
            "class_parity_equals_trace_of_abscissa": parity,
            "admissible_summand_counts": adm,
        })
        print(f"  classes  {name:42s} {hist}  admissible m = {adm}", flush=True)

    # ---- the existence witness, on the challenge curve ---------------------
    witnesses = [w for w in (witness_at_131(curve, r, 45, m, rng) for m in (2, 3, 4))
                 if w is not None]
    for w in witnesses:
        assert w["sum_is_on_curve"] and w["sum_lies_in_subgroup_of_order_r"]
        assert w["every_summand_abscissa_in_subspace"]
        if "semaev_s4_vanishes" in w:
            assert w["semaev_s4_vanishes"], "S_4 must vanish on a real 3-decomposition"
    print(f"  witnesses at n = 131: m = {[w['m'] for w in witnesses]}, "
          f"all certified in <G>", flush=True)

    # ---- 3. the cost surface ----------------------------------------------
    rows = []

    def add(label, cell, cls_):
        rows.append({"variant": label, **cell,
                     "log2_ratio_to_rho": round(cell["log2_total"] - log2_rho, 2),
                     "class": cls_})

    for m in range(2, 7):
        add(f"m = {m}, enumerate m-1 and root-find the last (the built oracle)",
            cost_cell(m, saturating_l(m, N131), N131, 1), "derived")
    for m in range(2, 7):
        add(f"m = {m}, Frobenius-stable orbit-union base, materialised",
            cost_cell(m, saturating_l(m, N131), N131, 1, frobenius=True), "derived")
    for m in (2, 3, 4):
        for split in range(2, m + 1):
            add(f"m = {m}, tabulate {split}-subset sums, unbounded memory",
                cost_row(m, N131, split),
                "generic in costume" if split == m else "derived")
    memory_rows = []
    for cap in (30, 40, 50, 60, 70):
        best, tag = None, None
        for m in range(2, 9):
            for split in range(1, m + 1):
                cell = memory_capped_row(m, N131, split, cap)
                if cell and (best is None or cell["log2_total"] < best["log2_total"]):
                    best, tag = cell, (m, split)
        gen = generic_with_memory(cap, log2r)
        add(f"best decomposition cell with at most 2^{cap} table entries "
            f"(m = {tag[0]}, tabulating {tag[1]}-subset sums)", best,
            "generic in costume" if tag[0] == tag[1] else "derived")
        memory_rows.append({
            "log2_entries_cap": cap, "m": tag[0], "split": tag[1],
            "log2_decomposition_total": best["log2_total"],
            "log2_bsgs_with_the_same_memory": gen["log2_bsgs"],
            "decomposition_beats_bsgs": best["log2_total"] < gen["log2_bsgs"],
            "log2_rho_needs_no_memory": round(log2_rho, 2),
        })

    flatness = [cost_cell(3, l, N131, 1) for l in (12, 20, 28, 36, 44)]
    floors = [free_oracle_floor(m, N131) for m in range(2, 9)]
    budgets = [oracle_budget(m, N131, log2_rho) for m in range(2, 9)]

    cost = {
        "unit": "log2 group operations, everything inside: targets tried, oracle work, "
                "table setup and linear algebra",
        "model": {
            "relations_needed": "|F| = 2^l (no Frobenius orbit collapse is available, see below)",
            "targets_per_relation": "2^n / C(|F|, m)",
            "oracle_per_target": "C(|F|, m - split) after tabulating C(|F|, split) sums",
            "linear_algebra": "m 2^(2l), sparse Wiedemann on an |F| x |F| matrix mod r",
        },
        "the_product_law": "2^l relations x 2^n/C(|F|,m) targets x C(|F|,m-1) oracle "
                           "= m 2^n, independent of l",
        "flatness_check_m3": flatness,
        "rows": rows,
        "free_oracle_floor": floors,
        "oracle_budget_to_match_rho": budgets,
        "memory_budgets_against_a_generic_algorithm": memory_rows,
    }
    for row in rows:
        print(f"  cost  {row['variant'][:58]:58s} l={row['l']:6.2f}  "
              f"2^{row['log2_total']:7.2f}  = 2^{row['log2_ratio_to_rho']:+.2f} x rho",
              flush=True)
    for mr in memory_rows:
        print(f"  memory 2^{mr['log2_entries_cap']:2d} entries: decomposition "
              f"2^{mr['log2_decomposition_total']:.2f}  vs BSGS "
              f"2^{mr['log2_bsgs_with_the_same_memory']:.2f}  vs rho "
              f"2^{log2_rho:.2f}", flush=True)

    # ---- 4. the Frobenius saving is not available at n = 131 --------------
    dims = available_subspace_dimensions(N131)
    assert dims == [0, 1, 130, 131], dims
    orbit_rows = [cost_cell(m, saturating_l(m, N131), N131, 1, frobenius=True)
                  for m in range(2, 9)]
    frobenius = {
        "statement": "no Frobenius-stable *subspace* of usable dimension exists at n = 131, "
                     "so the subspace root-finding oracle and the GGMP collapse cannot be had "
                     "at once; the collapse alone needs only pi(F) = F and is available from "
                     "any union of orbits, at the price of materialising the base",
        "closed_form_with_collapse": "collection = m 2^n / n",
        "orbit_union_rows": orbit_rows,
        "cyclotomic_coset_sizes_mod_131": cyclotomic_coset_sizes(N131),
        "available_invariant_dimensions": dims,
        "ord_131_of_2": mult_order(2, N131),
        "two_is_primitive_root_mod_131": mult_order(2, N131) == 130,
        "saving_forgone_relations": N131,
        "saving_forgone_linear_algebra": N131 ** 2,
        "log2_saving_forgone_linear_algebra": round(2 * math.log2(N131), 2),
        "note": "the relation collapse alone gives m 2^131 / 131 = 2^"
                f"{round(math.log2(3) + 131 - math.log2(131), 2)}, still 2^"
                f"{round(math.log2(3) + 131 - math.log2(131) - log2_rho, 2)} x rho; granting "
                "the linear-algebra saving on top changes nothing, because collection "
                "dominates at every m",
    }

    best_row = min(rows, key=lambda x: x["log2_total"])
    report = {
        "schema": "ecc2k130_point_decomposition/v1",
        "question": "can a target in <G> <= E(F_2^131) be written as a sum of m points "
                    "whose abscissae lie in an F_2-subspace, and can such a sum be found?",
        "verdict": "EXISTS AND IS ADMISSIBLE; NOT FINDABLE BELOW 2^131 WITH ANY ORACLE BUILT HERE",
        "unit": "log2 group operations; S = ops / sqrt(r)",
        "target": target,
        "existence": existence,
        "admissibility_at_131": admissibility,
        "witnesses_at_131": witnesses,
        "cost": cost,
        "frobenius_unavailable": frobenius,
        "headline": {
            "cheapest_implicit_base_variant_log2": min(
                x["log2_total"] for x in rows if x["log2_table_entries"] is None),
            "cheapest_implicit_base_ratio_log2": round(min(
                x["log2_ratio_to_rho"] for x in rows if x["log2_table_entries"] is None), 2),
            "cheapest_orbit_union_log2": min(
                x["log2_total"] for x in rows if x.get("frobenius_stable")),
            "cheapest_orbit_union_ratio_log2": round(min(
                x["log2_ratio_to_rho"] for x in rows if x.get("frobenius_stable")), 2),
            "cheapest_any_memory_log2": best_row["log2_total"],
            "cheapest_any_memory_ratio_log2": best_row["log2_ratio_to_rho"],
            "cheapest_any_memory_entries_log2": best_row["log2_table_entries"],
        },
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\n  rho reference 2^{log2_rho:.4f}   cheapest implicit-base variant "
          f"2^{report['headline']['cheapest_implicit_base_variant_log2']:.2f} "
          f"= 2^{report['headline']['cheapest_implicit_base_ratio_log2']:.2f} x rho; "
          f"cheapest orbit-union 2^{report['headline']['cheapest_orbit_union_log2']:.2f} "
          f"= 2^{report['headline']['cheapest_orbit_union_ratio_log2']:.2f} x rho")
    print(f"  wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
