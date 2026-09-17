#!/usr/bin/env python3
"""Homogeneous point relations on Koblitz curves, and whether they say anything.

Two things get kept apart here, because conflating them is the mistake this
whole study exists to record.

**A relation** is a multiset of factor-base points summing to the identity.
For three points that is exactly collinearity, and collinearity is exactly
Semaev's `S3`: substituting `y = lam x + mu` into `y^2 + xy = x^3 + a x^2 + b`
gives

    x^3 + (a + lam^2 + lam) x^2 + mu x + (b + mu^2) = 0

so the elementary symmetric functions of a collinear triple obey
`e2 = mu` and `e3 = b + mu^2`, i.e.

    x1 x2 x3 = b + (x1 x2 + x1 x3 + x2 x3)^2                        (*)

which rearranges to `(x1+x2)^2 x3^2 + (x1 x2) x3 + (x1 x2)^2 + b = 0`, the
`s3_in_last` quadratic already in `scripts/ecc2k130_point_decomposition.py`.
That is what makes the sweep affordable: a pair fixes a quadratic, so the
third abscissa is two candidates to look up, not a base to scan.  The
curve's `a` does not enter (*) at all; it enters only the solvability of
`lam^2 + lam = a + e1`, which is the separate trace condition checked here.

**Information about the discrete logarithm** is a strictly smaller thing.
If every base point is written `R_i = [u_i] P + [v_i] Q`, a relation
`sum R_i = O` gives

    sum u_i + log_P(Q) * sum v_i = 0   (mod r)

which determines `log_P(Q)` only when `sum v_i` is invertible mod `r`.  A
base *constructed* from `P` and `Q` comes with relations that hold by
construction -- `R_{1,0} + R_{0,1} - R_{1,1} = O` and its kin -- and those
have `sum u_i = sum v_i = 0`.  They are real relations, they contribute
real rank to a relation matrix, and they carry exactly zero information.
`useful_rank` below therefore quotients by the construction span before
reporting anything, and `solve_log` reports the recovered scalar only from
a relation whose `sum v_i` is a unit.
"""

from __future__ import annotations

import itertools
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fastfield import FastGF2m, IRR131, N131   # noqa: E402

# ── the challenge instance ────────────────────────────────────────────────
# Certicom ECC2K-130: y^2 + xy = x^3 + 1 over F_2^131, cofactor 4.
CHALLENGE_PX = 0x051C99BFA6F18DE467C80C23B98C7994AA
CHALLENGE_PY = 0x042EA2D112ECEC71FCF7E000D7EFC978BD
CHALLENGE_QX = 0x06C997F3E7F2C66A4A5D2FDA13756A37B1
CHALLENGE_QY = 0x04A38D11829D32D347BD0C0F584D546E9A
CHALLENGE_ELL = 680564733841876926932320129493409985129


class Koblitz:
    """`y^2 + xy = x^3 + a x^2 + b` over `F_2^m`, on `FastGF2m`."""

    def __init__(self, field: FastGF2m, a: int = 0, b: int = 1):
        self.F, self.a, self.b = field, a, b

    # -- points ----------------------------------------------------------

    def points_over(self, x: int):
        """The zero, one or two affine points with abscissa `x`."""
        F = self.F
        if x == 0:
            return [(0, F.frobenius(self.b, F.deg - 1))]      # y = sqrt(b)
        c = x ^ self.a ^ F.mul(self.b, F.sqr(F.inv(x)))       # x + a + b/x^2
        z = F.solve_artin_schreier(c)
        if z is None:
            return []
        return [(x, F.mul(x, z)), (x, F.mul(x, z ^ 1))]

    def has_point(self, x: int) -> bool:
        F = self.F
        if x == 0:
            return True
        return F.trace(x ^ self.a ^ F.mul(self.b, F.sqr(F.inv(x)))) == 0

    def on_curve(self, P) -> bool:
        if P is None:
            return True
        F = self.F
        x, y = P
        lhs = F.sqr(y) ^ F.mul(x, y)
        rhs = F.mul(x, F.sqr(x)) ^ F.mul(self.a, F.sqr(x)) ^ self.b
        return lhs == rhs

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
            x3 = F.sqr(lam) ^ lam ^ self.a
            return (x3, F.sqr(x1) ^ F.mul(lam ^ 1, x3))
        lam = F.mul(y1 ^ y2, F.inv(x1 ^ x2))
        x3 = F.sqr(lam) ^ lam ^ x1 ^ x2 ^ self.a
        return (x3, F.mul(lam, x1 ^ x3) ^ x3 ^ y1)

    def mul(self, P, k: int):
        if k < 0:
            P, k = self.neg(P), -k
        R, D = None, P
        while k:
            if k & 1:
                R = self.add(R, D)
            D = self.add(D, D)
            k >>= 1
        return R

    def frobenius(self, P):
        """`sigma(x, y) = (x^2, y^2)`, an endomorphism of a Koblitz curve."""
        if P is None:
            return None
        F = self.F
        return (F.sqr(P[0]), F.sqr(P[1]))

    def sum_points(self, pts):
        acc = None
        for P in pts:
            acc = self.add(acc, P)
        return acc

    # -- Semaev's third summation polynomial, as a quadratic in the last -
    def s3_in_last(self, x1: int, x2: int):
        """`(A, B, C)` with `A X^2 + B X + C = S3(x1, x2, X)`."""
        F = self.F
        return (F.sqr(x1 ^ x2), F.mul(x1, x2), F.sqr(F.mul(x1, x2)) ^ self.b)

    def s3(self, x1: int, x2: int, x3: int) -> int:
        A, B, C = self.s3_in_last(x1, x2)
        F = self.F
        return F.mul(A, F.sqr(x3)) ^ F.mul(B, x3) ^ C


def challenge_curve():
    F = FastGF2m(N131, IRR131)
    E = Koblitz(F, a=0, b=1)
    P = (CHALLENGE_PX, CHALLENGE_PY)
    Q = (CHALLENGE_QX, CHALLENGE_QY)
    return F, E, P, Q, CHALLENGE_ELL


# ── the three-point sweep ─────────────────────────────────────────────────

def solve_s3_last(E: Koblitz, pairs):
    """For each `(x1, x2)`, the abscissae `x3` with `S3(x1, x2, x3) = 0`.

    Returns one list of `(index, x3)` per input pair position.  Uses one
    batched inversion for the whole call, so the per-pair cost is a handful
    of multiplications rather than a field inversion.

    The quadratic `A X^2 + B X + C` in char 2 with `A, B != 0` has roots
    `X = (B/A) w` where `w^2 + w = A C / B^2`, which is solvable exactly
    when `Tr(A C / B^2) = 0`; the two roots are `w` and `w + 1`.  The
    degenerate cases `A = 0` (x1 = x2) and `B = 0` (x1 or x2 zero) are
    handled separately because they are linear, not quadratic.
    """
    F = E.F
    coeffs = [E.s3_in_last(x1, x2) for x1, x2 in pairs]
    # invert A * B^2 once per pair, in one batch
    prods, bsq = [], []
    for A, B, C in coeffs:
        Bs = F.sqr(B)
        bsq.append(Bs)
        prods.append(F.mul(A, Bs) if (A and Bs) else 0)
    inv = F.batch_inv(prods)

    out = []
    for i, (A, B, C) in enumerate(coeffs):
        Bs, iv = bsq[i], inv[i]
        if A == 0:
            # linear: B X + C = 0
            out.append([] if B == 0 else [F.mul(C, F.inv(B))])
            continue
        if B == 0:
            # A X^2 + C = 0 -> X = sqrt(C/A), always exactly one root
            out.append([F.frobenius(F.mul(C, F.inv(A)), F.deg - 1)])
            continue
        t = F.mul(F.mul(F.sqr(A), C), iv)          # A C / B^2 = A^2 C / (A B^2)
        if F.trace(t):
            out.append([])
            continue
        w = F.half_trace(t)
        scale = F.mul(F.mul(B, Bs), iv)            # B / A = B^3 / (A B^2)
        out.append([F.mul(scale, w), F.mul(scale, w ^ 1)])
    return out


def three_point_relations(E: Koblitz, base_x, *, allow_repeats=True,
                          chunk=20000, progress=None, budget_pairs=None):
    """Every homogeneous three-point relation with all abscissae in `base_x`.

    Enumerates unordered pairs (with repetition when `allow_repeats`), solves
    `S3` for the third abscissa, and keeps the hits that lie in the base.
    Returned triples are sorted abscissa triples, deduplicated.

    This is a *complete* search over the base, not a sample: if it returns
    nothing, the base has no three-point relation, full stop.
    """
    xs = sorted(set(base_x))
    member = set(xs)
    n = len(xs)
    found = set()
    seen_pairs = 0

    def pair_stream():
        for i in range(n):
            start = i if allow_repeats else i + 1
            for j in range(start, n):
                yield (xs[i], xs[j])

    buf = []
    for pair in pair_stream():
        buf.append(pair)
        if len(buf) >= chunk:
            seen_pairs += len(buf)
            _drain(E, buf, member, found, allow_repeats)
            buf = []
            if progress:
                progress(seen_pairs, len(found))
            if budget_pairs and seen_pairs >= budget_pairs:
                return sorted(found), seen_pairs, False
    if buf:
        seen_pairs += len(buf)
        _drain(E, buf, member, found, allow_repeats)
        if progress:
            progress(seen_pairs, len(found))
    return sorted(found), seen_pairs, True


def _drain(E, buf, member, found, allow_repeats):
    for (x1, x2), roots in zip(buf, solve_s3_last(E, buf)):
        for x3 in roots:
            if x3 not in member:
                continue
            if not allow_repeats and (x3 == x1 or x3 == x2):
                continue
            found.add(tuple(sorted((x1, x2, x3))))


def realise_triple(E: Koblitz, x1: int, x2: int, x3: int):
    """Signed points `(P1, P2, P3)` over the abscissae with `P1+P2+P3 = O`.

    `S3(x1,x2,x3) = 0` says a collinear triple exists over these abscissae;
    it does not say which of the two points above each abscissa to take.
    This pins that down by trying the (at most eight) sign choices and
    returning the first that actually sums to the identity, so every
    reported relation is checked on the curve rather than inferred from
    the polynomial.
    """
    cands = [E.points_over(x) for x in (x1, x2, x3)]
    if any(not c for c in cands):
        return None
    for P1 in cands[0]:
        for P2 in cands[1]:
            for P3 in cands[2]:
                if E.sum_points([P1, P2, P3]) is None:
                    return (P1, P2, P3)
    return None


# ── factor bases ──────────────────────────────────────────────────────────

def weight_two_abscissae(F: FastGF2m):
    """Every `z^i + z^j`, `i < j` -- the complete weight-two set."""
    return [(1 << i) ^ (1 << j)
            for i in range(F.deg) for j in range(i + 1, F.deg)]


def weight_at_most_two_abscissae(F: FastGF2m):
    out = [0] + [1 << i for i in range(F.deg)]
    out.extend(weight_two_abscissae(F))
    return out


def base_from_abscissae(E: Koblitz, xs):
    """Keep the abscissae that actually carry a point of the curve."""
    return [x for x in xs if E.has_point(x)]


def random_matched_base(E: Koblitz, size: int, rng: random.Random):
    """A uniformly random base of the same size, for a matched control.

    Matched on size only.  The weight-two base is not uniform -- it is a
    structured, Frobenius-stable set -- so this control isolates what the
    structure buys and nothing else.
    """
    F = E.F
    out = []
    seen = set()
    while len(out) < size:
        x = rng.getrandbits(F.deg)
        if x in seen or not E.has_point(x):
            continue
        seen.add(x)
        out.append(x)
    return out


def constructed_base(E: Koblitz, P, Q, r: int, size: int, rng: random.Random,
                     small: int = 0):
    """A base of points `[u]P + [v]Q`, with the coefficients recorded.

    Returns `(abscissae, coeffs)` where `coeffs[x] = (u, v, R)` and `R` is
    the point actually stored over that abscissa.  `R` is kept because an
    abscissa carries both `R` and `-R`, so a relation's coefficients depend
    on which of the two a triple used; `relation_vectors` resolves that
    sign against `R` rather than assuming it.

    When `small` is set the first entries use small coefficient pairs, which
    is what manufactures construction relations: `R_{1,0} + R_{0,1} = R_{1,1}`
    holds for free and contributes rank that means nothing.
    """
    pts, coeffs = [], {}
    combos = []
    if small:
        for u in range(small + 1):
            for v in range(small + 1):
                if u or v:
                    combos.append((u, v))
    while len(combos) < size:
        combos.append((rng.randrange(1, r), rng.randrange(1, r)))
    for u, v in combos[:size]:
        R = E.add(E.mul(P, u), E.mul(Q, v))
        if R is None or R[0] in coeffs:
            continue
        coeffs[R[0]] = (u % r, v % r, R)
        pts.append(R[0])
    return pts, coeffs


# ── rank, and whether any of it is useful ─────────────────────────────────

def _inv_mod(a: int, m: int):
    g, x = _egcd(a % m, m)
    return None if g != 1 else x % m


def _egcd(a: int, b: int):
    old_r, r = a, b
    old_s, s = 1, 0
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
    return old_r, old_s


def relation_vectors(E: Koblitz, triples, coeffs, r: int):
    """Each verified relation as `(sum u, sum v)` mod `r`.

    A relation says `sum_i eps_i ([u_i] P + [v_i] Q) = O`, i.e.
    `(sum eps_i u_i) + log_P(Q) (sum eps_i v_i) = 0 mod r`.  The sign
    `eps_i` is read off the realised point -- an abscissa carries both `R`
    and `-R` -- and a triple that does not realise as an on-curve relation
    is dropped rather than counted.
    """
    rows = []
    for t in triples:
        realised = realise_triple(E, *t)
        if realised is None:
            continue
        su = sv = 0
        ok = True
        for pt in realised:
            entry = coeffs.get(pt[0])
            if entry is None:
                ok = False
                break
            u, v, R = entry
            eps = 1 if pt == R else -1
            if eps == -1 and pt != E.neg(R):
                ok = False
                break
            su += eps * u
            sv += eps * v
        if ok:
            rows.append((su % r, sv % r))
    return rows


def useful_rank(rows, r: int):
    """Split relation rank into the part that constrains `log_P(Q)` and the rest.

    A row `(su, sv)` with `su = sv = 0` is a construction equation: it holds
    whatever `log_P(Q)` is, so it constrains nothing.  A row with `sv`
    invertible mod `r` determines `log_P(Q)` outright.  Anything else is
    reported as it is rather than counted as progress.

    Returns a dict with the raw row count, the count of construction rows,
    the count of rows that determine the logarithm, and the recovered
    logarithm if the rows agree on one.
    """
    trivial = [row for row in rows if row[0] == 0 and row[1] == 0]
    determining = []
    for su, sv in rows:
        if sv == 0:
            continue
        iv = _inv_mod(sv, r)
        if iv is None:
            continue
        determining.append((-su * iv) % r)
    consistent = len(set(determining)) <= 1
    return {
        "rows": len(rows),
        "construction_rows": len(trivial),
        "determining_rows": len(determining),
        "useful_rank": 0 if not determining else 1,
        "log": determining[0] if determining and consistent else None,
        "consistent": consistent,
    }


def verify_log(E: Koblitz, P, Q, k: int) -> bool:
    return k is not None and E.mul(P, k) == Q
