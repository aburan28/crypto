#!/usr/bin/env python3
"""Normal bases of `F_2^m`, and the Frobenius-stable weight-two supports.

`run_four_point.py` rules out normalising a four-point relation by
Frobenius, on the grounds that the weight-two support is not a union of
`sigma`-orbits.  That is true of the challenge's *polynomial* basis and
false in general.  Weight-two is a property of a basis: in the polynomial
basis `z^i + z^j` squares to `z^{2i} + z^{2j}` and the reduction fires once
`2i >= 131`, but in a **normal** basis squaring is a cyclic shift of
coordinates and cannot leave the weight-two set at all.

Writing `alpha` for a normal element, the support

    x_{i,j} = alpha^(2^i) + alpha^(2^j),      0 <= i < j < m

satisfies `sigma(x_{i,j}) = x_{i+1, j+1}` with indices mod `m`.  For prime
`m` every orbit has exactly `m` members and there are `C(m,2)/m = (m-1)/2`
of them, one per difference class `d = j - i` in `1 .. (m-1)/2`.

Two consequences this module exists to provide:

  * lifting is constant on an orbit -- `sigma` preserves the curve, so an
    orbit lifts to curve points entirely or not at all, and a support is
    always `m * (number of lifting orbits)` points;
  * a pair sum's `sigma`-orbit can be keyed by a single canonical
    representative, so a collision table stores one entry per `m` pairs
    and every entry found stands for `m` relations.
"""

from __future__ import annotations

import random

from fastfield import FastGF2m


# ── normal elements ───────────────────────────────────────────────────────

def conjugates(F: FastGF2m, a: int):
    """`[a, a^2, a^4, ..., a^(2^(m-1))]`."""
    out, v = [], a
    for _ in range(F.deg):
        out.append(v)
        v = F.sqr(v)
    return out


def is_normal(F: FastGF2m, a: int) -> bool:
    """True when the conjugates of `a` are an `F_2`-basis of the field.

    Gaussian elimination over `F_2` on the `m x m` matrix whose rows are the
    conjugates written in the ambient (polynomial) basis.  A row that
    reduces to zero means the conjugates are dependent, which is exactly
    failure to be normal.
    """
    if a == 0:
        return False
    pivots: dict[int, int] = {}
    for v in conjugates(F, a):
        cur = v
        while cur:
            h = cur.bit_length() - 1
            if h not in pivots:
                pivots[h] = cur
                break
            cur ^= pivots[h]
        if cur == 0:
            return False
    return len(pivots) == F.deg


def find_normal_elements(F: FastGF2m, count: int, rng: random.Random,
                         max_tries: int = 10000):
    """`count` distinct normal elements, drawn uniformly at random.

    The density of normal elements is high -- about 1 in 3 for `m = 131` in
    practice -- so rejection sampling is the whole algorithm.
    """
    found, seen, tries = [], set(), 0
    while len(found) < count and tries < max_tries:
        tries += 1
        a = rng.getrandbits(F.deg) & F.mask
        if a and a not in seen and is_normal(F, a):
            seen.add(a)
            found.append(a)
    if len(found) < count:
        raise RuntimeError(
            f"only {len(found)} normal elements in {tries} tries for m={F.deg}")
    return found


# ── the support ───────────────────────────────────────────────────────────

class NormalSupport:
    """The `sigma`-stable weight-two support of a normal element.

    `orbits` holds the lifting difference-classes; `abscissae` is the flat
    support in orbit-major order, so `abscissae[t * m + k]` is
    `sigma^k` applied to the `t`-th orbit representative.  That layout is
    what makes the canonical anchoring in `anchored_pairs` a slice rather
    than a search.
    """

    def __init__(self, F: FastGF2m, E, alpha: int):
        if not is_normal(F, alpha):
            raise ValueError("alpha is not a normal element")
        m = F.deg
        self.F, self.E, self.alpha, self.m = F, E, alpha, m
        conj = conjugates(F, alpha)

        # One representative per difference class, then its whole orbit.
        self.orbits, flat, mixed = [], [], []
        for d in range(1, (m + 1) // 2 + (0 if m % 2 else 1)):
            if d > m // 2:
                break
            orbit = [conj[k] ^ conj[(k + d) % m] for k in range(m)]
            lifts = [E.has_point(x) for x in orbit]
            if all(lifts):
                self.orbits.append(d)
                flat.extend(orbit)
            elif any(lifts):
                mixed.append(d)
        self.mixed_orbits = mixed
        self.abscissae = flat
        self.orbit_count = len(self.orbits)

    # -- structure ------------------------------------------------------

    @property
    def size(self) -> int:
        return len(self.abscissae)

    def total_orbits(self) -> int:
        return self.m // 2

    def is_sigma_stable(self) -> bool:
        """Check closure under `sigma` directly, rather than by argument."""
        S = set(self.abscissae)
        return all(self.F.sqr(x) in S for x in self.abscissae)


# ── the support, indexed ──────────────────────────────────────────────────

    def index(self, orbit_pos: int, k: int) -> int:
        """Flat index of `sigma^k` applied to the `orbit_pos`-th orbit rep."""
        return orbit_pos * self.m + (k % self.m)

    def canonical_pair_orbits(self):
        """One representative signed pair per `sigma`-orbit of signed pairs.

        A signed pair is an unordered pair of distinct support abscissae
        together with a *relative* sign: `{a, b}` with `eps` standing for
        `R_a + eps R_b`.  (A global sign is not a degree of freedom: a sum
        and its negative share an abscissa, which is the only thing the
        collision search looks at.)  There are `B(B-1)` of these.

        `sigma` acts by shifting both positions, and no signed pair is fixed
        by a non-trivial rotation, so every orbit has exactly `m` members
        and there are `B(B-1)/m` of them.

        Enumerating one member per orbit is a matter of anchoring the first
        element at position 0 -- *and of not doing that twice*.  A pair with
        its two elements in different orbits can be rotated to put either
        element at position 0, and both results are anchored; taking both is
        the double count that makes a run report every relation twice.  The
        enumeration below is by `(t1, t2, delta)` instead, which names each
        orbit exactly once by construction:

          * `t1 < t2`: one orbit per `delta = k2 - k1` in `0 .. m-1`;
          * `t1 == t2`: `delta` and `m - delta` give the same unordered
            pair, so `delta` runs over `1 .. (m-1)/2` only.

        Yields `(ia, ib, eps)` with `ia`, `ib` flat support indices.
        """
        m, T = self.m, self.orbit_count
        for t1 in range(T):
            for t2 in range(t1, T):
                deltas = range(m) if t2 > t1 else range(1, (m + 1) // 2)
                for delta in deltas:
                    ia = self.index(t1, 0)
                    ib = self.index(t2, delta)
                    yield ia, ib, 1
                    yield ia, ib, -1

    def pair_orbit_count(self) -> int:
        m, T = self.m, self.orbit_count
        return 2 * (T * (T - 1) // 2 * m + T * ((m - 1) // 2))

    def sigma_canonical_x(self, x: int) -> int:
        """The least abscissa in the `sigma`-orbit of `x`.

        Keys a pair sum by its whole Frobenius orbit, so one stored entry
        stands for `m` pairs.  Negation needs no handling: `S` and `-S`
        already share an abscissa.
        """
        best, cur = x, x
        for _ in range(self.m - 1):
            cur = self.F.sqr(cur)
            if cur < best:
                best = cur
        return best

    def orbit_points(self, E):
        """Support points with signs consistent along each `sigma`-orbit.

        `sigma^k` applied to an orbit representative is a *specific* point,
        but an abscissa carries two, and picking one per abscissa
        independently -- by least `y`, say -- flips the sign at about half
        the positions.  The collision search does not notice, because it
        only ever looks at abscissae.  The *solve* does: it relies on

            log(sigma^k R) = s^k log(R)

        and a flipped sign negates that term, so a system built on
        independently chosen points is wrong in about half its entries
        while still looking perfectly well-formed.  Generating each orbit by
        applying `sigma` removes the choice.
        """
        pts = []
        for t in range(self.orbit_count):
            x = self.abscissae[t * self.m]
            cands = E.points_over(x)
            if not cands:
                raise RuntimeError("a lifting orbit has an abscissa with no point")
            P = min(cands, key=lambda p: p[1])
            for _ in range(self.m):
                pts.append(P)
                P = E.frobenius(P)
            if P != pts[t * self.m]:
                raise RuntimeError("sigma^m did not return the orbit rep")
        return pts
