#!/usr/bin/env python3
"""E1-E6: running the six experiments `research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION_TARGETS.md`
pre-registered, against the boundaries frozen before any of them was run.

`experiments/ecc2k130_decomposition_targets.json` holds the predictions.  This
script holds the runs.  Nothing here re-derives a boundary; every number below
is produced by executing an algorithm and counting what it did.

    python3 scripts/ecc2k130_decomposition_experiments.py [--max-rung 19]

Writes `experiments/ecc2k130_decomposition_runs.json`.

The unit, fixed by the design: `Lambda = oracle operations / 2^n`, where an
*oracle operation* is one `(m-1)`-subset enumerated and completed -- the unit the
product law `2^l relations x 2^n/C(|F|,m) targets x C(|F|,m-1) oracle = m 2^n`
is written in.  Actual curve additions are counted alongside and reported, but
they are not the metric: two of them happen per subset.

Three counting rules are separated, because they are three different algorithms
and only the first is the one the law prices:

  full        enumerate every (m-1)-subset of the base for each target and
              harvest every decomposition found.  Predicted Lambda = m.
  first_hit   stop at the first decomposition of each target.  At an
              oversaturating `l` this is cheaper by about the yield, and it is a
              DIFFERENT algorithm -- reported in its own column, never mixed in.
  folded      full enumeration, but exploiting log(-P) = -log(P) so the base
              contributes |F|/2 unknowns instead of |F|.  A factor of two.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import random
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "experiments/ecc2k130_decomposition_runs.json"

_SPEC = importlib.util.spec_from_file_location(
    "ecc2k130_decomp", Path(__file__).with_name("ecc2k130_point_decomposition.py"))
assert _SPEC and _SPEC.loader
DEC = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(DEC)

PREDICTIONS = json.loads(
    (REPO / "experiments/ecc2k130_decomposition_targets.json").read_text())
LOG2_RHO_131 = PREDICTIONS["log2_rho_reference"]


# ── the rung: a curve, its largest prime subgroup, and a generator ───────

class Rung:
    """`K_0 : y^2 + xy = x^3 + 1` over `F_2^n`, posed in its largest prime
    subgroup so the relation matrix lives over a field.

    `E(F_2^n)` is not cyclic at every rung -- at `n = 11` it is `Z/23 x Z/92`
    (measured: order 2116, exponent 92), so `#E/p = 92` is the exponent itself
    and annihilates the whole group, and cannot be the projection.  The
    scalar that works is `exponent/p`: it kills everything of order coprime to
    `p` and lands every point in the order-`p` subgroup, cyclic or not.
    """

    #: above this degree the discrete-log tables that make field
    #: multiplication a lookup stop fitting in memory -- `2^29` entries is the
    #: allocation that killed the first attempt at the `n = 29` rung -- and
    #: without them a curve addition costs an inversion by exponentiation,
    #: roughly three orders of magnitude slower.  That, not the oracle work, is
    #: what keeps the ladder's upper rungs out of reach of this harness.
    MAX_TABLE_DEGREE = 20

    def __init__(self, n: int, tables: bool | None = None):
        self.n = n
        if tables is None:
            tables = n <= self.MAX_TABLE_DEGREE
        self.tables = tables
        self.field = DEC.GF2m(n, DEC.find_irreducible(n), tables=tables)
        self.curve = DEC.Koblitz(self.field)
        self.order = DEC.curve_order(n)
        self.p = max(DEC.prime_factors(self.order))
        assert DEC.is_prime(self.p)
        self.exponent = self._exponent()
        assert self.order % self.exponent == 0
        assert self.exponent % self.p == 0, "p must divide the group exponent"
        self.group_shape = ("cyclic" if self.exponent == self.order
                            else f"Z/{self.order // self.exponent} x Z/{self.exponent}")
        # `E = Z/d1 x Z/d2` with `d1 = order/exponent`.  The `p`-torsion has rank
        # two exactly when `p | d1`, and then there is no single cyclic subgroup
        # to pose the logarithm in: projecting lands in a two-dimensional group
        # and `log_G` is not defined on it.  `n = 11` is this case (`#E = 23^2`),
        # which is the distortion the design flagged in advance.
        self.p_torsion_rank = 2 if (self.order // self.exponent) % self.p == 0 else 1
        self.usable_end_to_end = self.p_torsion_rank == 1
        self.proj_scalar = self.exponent // self.p
        self.G = self._generator() if self.usable_end_to_end else None

    def _point_order(self, P):
        c, o = self.curve, self.order
        for q in sorted(set(DEC.prime_factors(self.order))):
            while o % q == 0 and c.mul(P, o // q) is None:
                o //= q
        return o

    def _exponent(self) -> int:
        """lcm of the orders of the first points found; stable long before it
        runs out of points, and asserted to divide the group order."""
        c, e, seen = self.curve, 1, 0
        for x in range(1, 1 << self.n):
            for P in c.points_over(x):
                e = e * self._point_order(P) // math.gcd(e, self._point_order(P))
                seen += 1
            if seen >= 40:
                break
        return e

    def project(self, P):
        """`[exponent/p] P`, the component of `P` in the order-`p` subgroup."""
        return self.curve.mul(P, self.proj_scalar)

    def _generator(self):
        c = self.curve
        for x in range(1, 1 << self.n):
            for P in c.points_over(x):
                H = self.project(P)
                if H is not None:
                    assert c.mul(H, self.p) is None, "projection missed the subgroup"
                    return H
        raise AssertionError("no generator of the prime subgroup")


# ── the factor base, and the oracle the product law prices ──────────────

def build_base(rung: Rung, l: int):
    """Every point whose abscissa lies in `span{1, z, ..., z^(l-1)}`.

    Returns the point list, a membership map to the index of the point's
    abscissa, and the abscissa count.  Both points over an abscissa are in the
    base; which of the two is the `+` representative is fixed by the y-order.
    """
    c = rung.curve
    points, index, abscissae = [], {}, []
    for x in range(1 << l):
        pts = c.points_over(x)
        if not pts:
            continue
        a = len(abscissae)
        abscissae.append(x)
        for P in sorted(pts):
            points.append(P)
            index[P] = a
    return points, index, abscissae


def collect(rung: Rung, l: int, m: int, rule: str, rng: random.Random,
            want: int, cap_targets: int = 1 << 30):
    """Collect `want` relations, counting what the oracle actually did.

    Each target is `R = [a] G`.  A decomposition `R = P_1 + ... + P_m` gives
    `sum_i log([c] P_i) = c a (mod p)` with `c = exponent/p`, because projecting
    is a homomorphism and lands everything in the order-`p` subgroup.  No class
    filtering is needed: `R` lies in the subgroup, so any decomposition of it
    already has summand classes cancelling.
    """
    assert rule in ("full", "first_hit", "folded")
    assert m == 3, "the pre-registered ladder is m = 3"
    c = rung.curve
    points, index, abscissae = build_base(rung, l)
    member = set(points)
    npts, nabs = len(points), len(abscissae)
    unknowns = nabs if rule == "folded" else npts

    # the + representative of each abscissa, for the folded accounting
    plus = {}
    for P in points:
        plus.setdefault(index[P], P)
    pos = {P: i for i, P in enumerate(points)}

    rows, rhs = [], []
    oracle_ops = group_adds = targets = 0
    while len(rows) < want and targets < cap_targets:
        a = rng.randrange(1, rung.p)
        R = c.mul(rung.G, a)
        if R is None:
            continue
        targets += 1
        negR = c.neg(R)
        found_this_target = []
        # `first_hit` takes one relation per target, so it must not scan the base
        # in a fixed order: doing so makes the low-indexed points appear in
        # almost every relation and the high-indexed ones in none, and the matrix
        # comes out rank-deficient no matter how many targets are thrown at it.
        # (Measured: at n = 19 a fixed-order scan never reached rank 279.)
        scan = list(range(npts))
        if rule == "first_hit":
            rng.shuffle(scan)
        for ii in range(npts):
            i = scan[ii]
            Pi = points[i]
            for jj in range(ii + 1, npts):
                j = scan[jj]
                oracle_ops += 1
                S = c.add(Pi, points[j])
                group_adds += 1
                if S is None:
                    continue
                T = c.add(negR, S)          # T = S - R, so R = Pi + Pj + (-T)
                group_adds += 1
                if T is None:
                    continue
                U = c.neg(T)
                # canonical order: a triple is recorded once, from its two
                # lowest-indexed members.  Without this every decomposition is
                # found m times -- once per (m-1)-subset of it -- and the yield
                # reads three times its true value.
                # canonicalisation is what stops `full` recording each triple
                # once per (m-1)-subset of it; `first_hit` keeps one relation per
                # target anyway, so it does not need -- and must not impose -- it
                if U in member and (rule == "first_hit" or pos[U] > j):
                    found_this_target.append((Pi, points[j], U))
                    if rule == "first_hit":
                        break
            if rule == "first_hit" and found_this_target:
                break
        for trio in found_this_target:
            acc = None                      # correctness gate: re-add in the group
            for P in trio:
                acc = c.add(acc, P)
            assert acc == R, "a collected relation does not sum to its target"
            row: dict = {}
            for P in trio:
                if rule == "folded":
                    k = index[P]
                    sign = 1 if P == plus[k] else -1
                else:
                    k = pos[P]
                    sign = 1
                row[k] = row.get(k, 0) + sign
            cols = {k: v % rung.p for k, v in row.items() if v % rung.p}
            if not cols:
                continue                     # trivial relation: never counted
            rows.append(cols)
            rhs.append((rung.proj_scalar % rung.p) * a % rung.p)
            if len(rows) >= want:
                break
    return {
        "points": points, "index": index, "abscissae": abscissae,
        "unknowns": unknowns, "plus": plus,
        "rows": rows, "rhs": rhs,
        "oracle_ops": oracle_ops, "group_adds": group_adds, "targets": targets,
    }


# ── linear algebra over F_p ──────────────────────────────────────────────

def solve_mod_p(rows, rhs, unknowns: int, p: int):
    """Gaussian elimination on the sparse relation matrix, densified.

    Returns the solution vector, or None if the relations do not determine it.
    """
    mat = [[0] * (unknowns + 1) for _ in range(len(rows))]
    for r, (cols, b) in enumerate(zip(rows, rhs)):
        for k, v in cols.items():
            mat[r][k] = v % p
        mat[r][unknowns] = b % p
    rank = 0
    where = [-1] * unknowns
    for col in range(unknowns):
        piv = None
        for r in range(rank, len(mat)):
            if mat[r][col]:
                piv = r
                break
        if piv is None:
            continue
        mat[rank], mat[piv] = mat[piv], mat[rank]
        inv = pow(mat[rank][col], p - 2, p)
        mat[rank] = [(v * inv) % p for v in mat[rank]]
        for r in range(len(mat)):
            if r != rank and mat[r][col]:
                f = mat[r][col]
                mat[r] = [(x - f * y) % p for x, y in zip(mat[r], mat[rank])]
        where[col] = rank
        rank += 1
        if rank == len(mat):
            break
    for r in range(rank, len(mat)):
        if not any(mat[r][:unknowns]) and mat[r][unknowns]:
            return None, rank                       # inconsistent
    if rank < unknowns:
        return None, rank
    sol = [0] * unknowns
    for col in range(unknowns):
        sol[col] = mat[where[col]][unknowns]
    return sol, rank


def matrix_rank_mod_p(rows, unknowns: int, p: int) -> int:
    mat = [[0] * unknowns for _ in range(len(rows))]
    for r, cols in enumerate(rows):
        for k, v in cols.items():
            mat[r][k] = v % p
    rank = 0
    for col in range(unknowns):
        piv = None
        for r in range(rank, len(mat)):
            if mat[r][col]:
                piv = r
                break
        if piv is None:
            continue
        mat[rank], mat[piv] = mat[piv], mat[rank]
        inv = pow(mat[rank][col], p - 2, p)
        mat[rank] = [(v * inv) % p for v in mat[rank]]
        for r in range(rank + 1, len(mat)):
            if mat[r][col]:
                f = mat[r][col]
                mat[r] = [(x - f * y) % p for x, y in zip(mat[r], mat[rank])]
        rank += 1
        if rank == len(mat):
            break
    return rank


# ── what a relation budget actually determines ──────────────────────────

def determined_columns(rows, unknowns: int, p: int):
    """The columns whose value every solution of the system agrees on.

    A pivot column is determined only if its echelon row is free of the
    non-pivot (free) columns; otherwise the kernel moves it.  Columns no
    relation ever touches are of course undetermined.  This is the set whose
    logarithms the descent may actually use.
    """
    mat = [dict((k, v % p) for k, v in r.items() if v % p) for r in rows]
    pivots: dict = {}                      # col -> reduced row
    for row in mat:
        cur = dict(row)
        while cur:
            c = min(cur)
            if c not in pivots:
                inv = pow(cur[c], p - 2, p)
                pivots[c] = {k: (v * inv) % p for k, v in cur.items()}
                break
            f, pr = cur[c], pivots[c]
            cur = {k: v for k in set(cur) | set(pr)
                   if (v := (cur.get(k, 0) - f * pr.get(k, 0)) % p)}
    # Back-substitute so each pivot row carries only free columns beside its own.
    # Every pivot is the least column of its row, so working from the highest
    # pivot column down means the row being eliminated is already reduced and
    # brings in free columns only -- which is why this terminates.
    for c in sorted(pivots, reverse=True):
        row = pivots[c]
        while True:
            k = next((k for k in row if k != c and k in pivots), None)
            if k is None:
                break
            f, pr = row[k], pivots[k]
            row = {j: v for j in set(row) | set(pr)
                   if (v := (row.get(j, 0) - f * pr.get(j, 0)) % p)}
        pivots[c] = row
    free = {k for k in range(unknowns) if k not in pivots}
    return {c for c, row in pivots.items() if not (set(row) - {c}) & free}


def budget_run(rung: Rung, l: int, m: int, rng: random.Random,
               slack_fraction: float = 0.0, descent_tries: int = 200):
    """Collect the `|F|` relations the product law budgets -- no more -- and see
    whether that is enough to finish.

    The harness that produced the first round of E1/E2 insisted on full rank
    over every column, which is the coupon-collector problem and costs
    `ln|F|/m` times the budget.  Index calculus does not need that: it drops the
    base points no relation determined and retries the descent when one turns
    up in a decomposition.  This measures both halves -- how much of the base
    `|F|` relations determine, and how many descents that costs.
    """
    points, index, abscissae = build_base(rung, l)
    unknowns = len(points)
    want = int(unknowns * (1.0 + slack_fraction))
    d = collect(rung, l, m, 'full', rng, want=want)
    rows = d["rows"][:want]
    known = determined_columns(rows, unknowns, rung.p)
    frac = len(known) / unknowns
    # `e^{-m}` predicts the columns no relation *touches*.  Determination is a
    # strictly stronger property -- a column can be hit and still be free, when
    # the only rows on it also touch a column nothing pins down -- so the two
    # fractions are reported against their own predictions and never conflated.
    hit = set().union(*(set(r) for r in rows)) if rows else set()
    frac_hit = len(hit & set(range(unknowns))) / unknowns
    # A base point in the torsion that projection kills carries the unknown
    # `log(O) = 0`, which no relation can pin and no descent can use.  Small
    # bases are where this bites: at `l = 5` three of 29 points are such, and
    # that -- not the stopping rule -- is why that cell behaves differently.
    dead = sum(1 for P in points if rung.project(P) is None)

    pos = {P: i for i, P in enumerate(points)}
    member = set(points)
    c = rung.curve
    secret = rng.randrange(2, rung.p)
    Q = c.mul(rung.G, secret)
    attempts, descent_ops, landed = 0, 0, False
    for _ in range(descent_tries):
        attempts += 1
        cc, dd = rng.randrange(1, rung.p), rng.randrange(1, rung.p)
        T = c.add(c.mul(Q, cc), c.mul(rung.G, dd))
        if T is None:
            continue
        negT = c.neg(T)
        hit = None
        for i in range(len(points)):
            for j in range(i + 1, len(points)):
                descent_ops += 1
                S = c.add(points[i], points[j])
                if S is None:
                    continue
                W = c.add(negT, S)
                if W is None:
                    continue
                U = c.neg(W)
                if U in member and pos[U] > j:
                    trio = (points[i], points[j], U)
                    if all(pos[P] in known for P in trio):
                        hit = trio
                        break
            if hit:
                break
        if hit:
            landed = True
            break
    return {
        "n": rung.n, "m": m, "l": l,
        "factor_base_size": unknowns,
        "relations_budgeted": want,
        "relations_per_unknown": round(want / unknowns, 3),
        "columns_determined": len(known),
        "fraction_determined": round(frac, 4),
        "fraction_undetermined": round(1 - frac, 4),
        "fraction_unhit": round(1 - frac_hit, 4),
        "poisson_prediction_unhit": round(math.exp(-m), 4),
        "base_points_killed_by_projection": dead,
        "descent_success_per_try_predicted": round(frac ** m, 4),
        "descent_attempts": attempts,
        "descent_landed": landed,
        "collection_oracle_ops": d["oracle_ops"],
        "descent_oracle_ops": descent_ops,
        "total_oracle_ops": d["oracle_ops"] + descent_ops,
        "lambda": round((d["oracle_ops"] + descent_ops) / 2 ** rung.n, 3),
    }


def full_rank_crossing(rung: Rung, l: int, m: int, rng: random.Random,
                       rule: str = "full", over: int = 4):
    """The relation count at which the rank first covers every column, against
    the coupon-collector law `|F| ln|F| / m` it obeys."""
    points, _, _ = build_base(rung, l)
    unknowns = len(points)
    d = collect(rung, l, m, rule, rng, want=unknowns * over)
    pivots: dict = {}
    rank, at = 0, None
    for i, row in enumerate(d["rows"], 1):
        cur = {k: v % rung.p for k, v in row.items() if v % rung.p}
        while cur:
            col = min(cur)
            if col not in pivots:
                inv = pow(cur[col], rung.p - 2, rung.p)
                pivots[col] = {k: (v * inv) % rung.p for k, v in cur.items()}
                rank += 1
                break
            f, pr = cur[col], pivots[col]
            cur = {k: v for k in set(cur) | set(pr)
                   if (v := (cur.get(k, 0) - f * pr.get(k, 0)) % rung.p)}
        if rank >= unknowns:
            at = i
            break
    coupon = unknowns * math.log(unknowns) / m
    return {
        "n": rung.n, "m": m, "l": l, "rule": rule,
        "factor_base_size": unknowns,
        "full_rank_at": at,
        "rho_measured": round(at / unknowns, 3) if at else None,
        "coupon_collector_relations": round(coupon, 1),
        "coupon_collector_rho": round(coupon / unknowns, 3),
        "measured_over_coupon": round(at / coupon, 3) if at else None,
    }


# ── E2: is the total really flat in the factor-base dimension? ──────────

def e2_rank_ceiling(n: int, l: int, m: int = 3):
    """The rank the relation matrix can reach at all, over every decomposition
    the base admits -- not over the ones a sample happened to find.

    E2's `l = 5` cell ran 15 455 targets and stopped at rank `28/29`, and the
    first round read that as the same relation-correlation effect it read into
    E1.  It is not.  Enumerating every triple of base points whose sum lands in
    the order-`p` subgroup -- 810 of them, exhaustively -- the row space has
    rank `28`, so no relation count whatsoever reaches full rank there.  The one
    kernel direction is constant on each abscissa, and the reason is the census
    below: at `l = 5` every reachable triple has exactly one even-abscissa
    summand, so `v = +1` on the odd abscissae and `v = -2` on the even ones
    annihilates every row.  From `l = 6` the base admits all-even triples too,
    that vector stops working, and the rank is full.
    """
    assert m == 3, "the pre-registered ladder is m = 3"
    rung = Rung(n)
    c, p = rung.curve, rung.p
    points, _, _ = build_base(rung, l)
    size = len(points)
    rows, census = [], {}
    for i in range(size):
        for j in range(i + 1, size):
            S = c.add(points[i], points[j])
            for k in range(j + 1, size):
                R = c.add(S, points[k]) if S is not None else points[k]
                if R is None or c.mul(R, p) is not None:
                    continue          # not a decomposition of a subgroup element
                rows.append({i: 1, j: 1, k: 1})
                odd = sum(points[t][0] & 1 for t in (i, j, k))
                census[odd] = census.get(odd, 0) + 1
    rank = matrix_rank_mod_p(rows, size, p)
    return {
        "n": n, "l": l, "m": m,
        "factor_base_size": size,
        "reachable_decompositions": len(rows),
        "rank_ceiling": rank,
        "deficiency": size - rank,
        "full_rank_reachable": rank == size,
        "odd_abscissa_census": {str(k): v for k, v in sorted(census.items())},
        "sampling": "exhaustive over every triple of base points",
    }


def e2_sweep(n: int, l: int, m: int = 3, rule: str = "full", seed: int = 4):
    """One rung, one dimension: collect to full rank and report `Lambda`."""
    rng = random.Random(seed * 1000 + l)
    rung = Rung(n)
    points, _, abscissae = build_base(rung, l)
    unknowns = len(abscissae) if rule == "folded" else len(points)
    want = unknowns + 8
    for _ in range(8):
        d = collect(rung, l, m, rule, rng, want)
        rank = matrix_rank_mod_p(d["rows"], d["unknowns"], rung.p)
        if rank >= d["unknowns"]:
            break
        want += max(8, unknowns // 3)
    lam = math.comb(len(points), m) / rung.order
    return {
        "n": n, "m": m, "l": l, "rule": rule,
        "factor_base_size": len(points), "unknowns": d["unknowns"],
        "saturating_l": round(DEC.saturating_l(m, n), 2),
        "oversaturated": lam > 1.0,
        "yield_predicted": round(lam, 4),
        "yield_measured": round(len(d["rows"]) / d["targets"], 4),
        "targets": d["targets"], "relations": len(d["rows"]),
        "matrix_rank": rank, "full_rank": rank >= d["unknowns"],
        "relations_per_unknown": round(len(d["rows"]) / d["unknowns"], 3),
        "oracle_ops": d["oracle_ops"],
        "log2_oracle_ops": round(math.log2(d["oracle_ops"]), 3),
        "lambda_metric": round(d["oracle_ops"] / 2 ** n, 3),
    }


# ── E4: large primes -- a real gain, or a relabelling? ──────────────────

def e4_large_prime(n: int, l: int, lp: int, m: int = 2, seed: int = 11,
                   targets: int = 400):
    """Allow the last summand's abscissa anywhere in `V' ⊃ V`, key the partial
    relation by that large prime, and pair partials that share one.

    The pairing cancels the large prime and leaves a relation over the base
    alone.  The design's caveat is that these paired relations may not be
    independent; this measures the **rank**, it does not assume it.
    """
    assert lp > l
    rng = random.Random(seed + n)
    rung = Rung(n)
    c = rung.curve
    points, index, abscissae = build_base(rung, l)
    pos = {P: i for i, P in enumerate(points)}
    member = set(points)
    big = 1 << lp
    npts = len(points)

    partials: dict = {}
    full_rows, full_rhs = [], []
    oracle_ops = 0
    guard_partials = 0
    for _ in range(targets):
        a = rng.randrange(1, rung.p)
        R = c.mul(rung.G, a)
        if R is None:
            continue
        negR = c.neg(R)
        for i in range(npts):
            oracle_ops += 1
            T = c.add(negR, points[i])        # T = P_i - R, so R = P_i + (-T)
            if T is None:
                continue
            U = c.neg(T)
            if U[0] >= big:
                continue                      # outside V', not a partial at all
            if U in member:
                continue                      # a plain relation, not a large prime
            guard_partials += 1
            key = U
            prev = partials.get(key)
            if prev is None:
                partials[key] = (pos[points[i]], a)
                continue
            j, b = prev                       # R - R' = P_i - P_j
            row = {}
            row[pos[points[i]]] = row.get(pos[points[i]], 0) + 1
            row[j] = row.get(j, 0) - 1
            cols = {k: v % rung.p for k, v in row.items() if v % rung.p}
            if not cols:
                continue
            full_rows.append(cols)
            full_rhs.append((rung.proj_scalar % rung.p) * ((a - b) % rung.p)
                            % rung.p)
    rank = matrix_rank_mod_p(full_rows, npts, rung.p) if full_rows else 0
    per_target = guard_partials / max(1, targets)
    return {
        "n": n, "m": m, "l": l, "large_prime_dim": lp,
        "factor_base_size": npts, "unknowns": npts,
        "targets": targets,
        "partial_relations": guard_partials,
        "partial_yield_per_target": round(per_target, 4),
        "guard_satisfied": per_target <= 1.0,
        "distinct_large_primes": len(partials),
        "paired_full_relations": len(full_rows),
        "matrix_rank": rank,
        "rank_over_relations": round(rank / len(full_rows), 4) if full_rows else None,
        "rank_over_unknowns": round(rank / npts, 4),
        "oracle_ops": oracle_ops,
    }


# ── E6: Frobenius-stable orbit-union bases, and relation independence ───

def frobenius_eigenvalue(rung: Rung) -> int:
    """`lambda` with `pi(P) = [lambda] P` on the order-`p` subgroup.

    For `K_0` the Frobenius satisfies `pi^2 + pi + 2 = 0` (trace `-1`), so
    `lambda` is a root of `T^2 + T + 2 mod p`; which of the two roots it is has
    to be settled by applying `pi` to a point, not guessed.
    """
    c, p_ = rung.curve, rung.p
    disc = (1 - 4 * 2) % p_                       # 1 - 4*2 = -7
    # Tonelli-Shanks for sqrt(disc) mod p_
    sq = _sqrt_mod(disc, p_)
    assert sq * sq % p_ == disc % p_, "T^2+T+2 has no root mod p"
    inv2 = pow(2, p_ - 2, p_)
    piG = _frob_point(rung, rung.G)
    for root in (((-1 + sq) * inv2) % p_, ((-1 - sq) * inv2) % p_):
        if c.mul(rung.G, root) == piG:
            return root
    raise AssertionError("neither root of T^2+T+2 acts as Frobenius")


def _sqrt_mod(a: int, p_: int) -> int:
    a %= p_
    if a == 0:
        return 0
    if p_ % 4 == 3:
        return pow(a, (p_ + 1) // 4, p_)
    q, sft = p_ - 1, 0
    while q % 2 == 0:
        q //= 2
        sft += 1
    z = 2
    while pow(z, (p_ - 1) // 2, p_) != p_ - 1:
        z += 1
    mt, cc = sft, pow(z, q, p_)
    t, r = pow(a, q, p_), pow(a, (q + 1) // 2, p_)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p_
            i += 1
        b = pow(cc, 1 << (mt - i - 1), p_)
        mt, cc = i, b * b % p_
        t, r = t * cc % p_, r * b % p_
    return r


def _frob_point(rung: Rung, P):
    if P is None:
        return None
    F = rung.field
    return (F.sqr(P[0]), F.sqr(P[1]))


def e6_orbit_union(n: int, orbits: int, m: int = 3, seed: int = 5,
                   slack: int = 10, relations_multiple: int = 1):
    """Build a base from `k` whole Frobenius orbits, collect relations, rewrite
    them onto orbit representatives with their `lambda^j` weights, and measure
    the **rank** against the `|F|/n` unknowns the collapse assumes.

    The saving the GGMP argument claims is exactly `n`; what it needs, and what
    nobody has checked, is that `|F|/n` orbit relations are independent.
    """
    rng = random.Random(seed + n)
    rung = Rung(n)
    c = rung.curve
    lam = frobenius_eigenvalue(rung)

    # whole orbits of the Frobenius, each of size exactly n (n prime, x not in F_2)
    seen, base, rep_of, power_of, reps = set(), [], {}, {}, []
    x = 2
    while len(reps) < orbits and x < (1 << n):
        pts = c.points_over(x)
        if not pts or pts[0] in seen:
            x += 1
            continue
        P = sorted(pts)[0]
        orbit, Q = [], P
        for j in range(n):
            orbit.append(Q)
            Q = _frob_point(rung, Q)
        if Q != P or len(set(orbit)) != n:
            x += 1
            continue                          # a short orbit: skip, not our case
        k = len(reps)
        reps.append(P)
        for j, Qj in enumerate(orbit):
            for R in (Qj, c.neg(Qj)):
                if R in seen:
                    continue
                seen.add(R)
                base.append(R)
                rep_of[R] = k
                power_of[R] = (j, 1 if R == Qj else -1)
        x += 1
    assert len(reps) == orbits, f"only {len(reps)} full orbits found"

    member = set(base)
    npts = len(base)
    unknowns = orbits
    lam_pow = [pow(lam, j, rung.p) for j in range(n)]

    rows, rhs = [], []
    targets = oracle_ops = 0
    want = max(unknowns * relations_multiple, unknowns + slack)
    while len(rows) < want and targets < 20000:
        a = rng.randrange(1, rung.p)
        R = c.mul(rung.G, a)
        if R is None:
            continue
        targets += 1
        negR = c.neg(R)
        for i in range(npts):
            for j in range(i + 1, npts):
                oracle_ops += 1
                S = c.add(base[i], base[j])
                if S is None:
                    continue
                T = c.add(negR, S)
                if T is None:
                    continue
                U = c.neg(T)
                if U in member:
                    trio = (base[i], base[j], U)
                    acc = None
                    for P in trio:
                        acc = c.add(acc, P)
                    assert acc == R
                    row: dict = {}
                    for P in trio:
                        k = rep_of[P]
                        jj, sgn = power_of[P]
                        w = lam_pow[jj] * sgn
                        row[k] = (row.get(k, 0) + w) % rung.p
                    cols = {k: v for k, v in row.items() if v % rung.p}
                    if not cols:
                        continue
                    rows.append(cols)
                    rhs.append((rung.proj_scalar % rung.p) * a % rung.p)
                    if len(rows) >= want:
                        break
            if len(rows) >= want:
                break
    rank = matrix_rank_mod_p(rows, unknowns, rung.p)
    sol, _ = solve_mod_p(rows, rhs, unknowns, rung.p)
    return {
        "n": n, "m": m, "orbits": orbits,
        "frobenius_eigenvalue_verified": True,
        "factor_base_size": npts,
        "orbit_unknowns": unknowns,
        "points_per_unknown": round(npts / unknowns, 2),
        "collapse_claimed": n,
        "targets": targets, "relations": len(rows),
        "matrix_rank": rank,
        "full_rank": rank >= unknowns,
        "solved": sol is not None,
        # `npts` counts a point and its negative separately and the sign fold is
        # applied alongside the orbit fold, so `npts/rank` measures BOTH
        # collapses.  The Frobenius part alone is half of it.
        "realised_saving_both_folds": round(npts / rank, 2) if rank else None,
        "realised_frobenius_saving": round(npts / 2 / rank, 2) if rank else None,
        "frobenius_saving_over_n": round(npts / 2 / rank / n, 3) if rank else None,
        "oracle_ops": oracle_ops,
    }


# ── E5: the yield distribution, not just its mean ───────────────────────

def odd_subgroup(rung: Rung):
    """The odd-order subgroup: every point of cofactor class zero, which is
    exactly the population an index-calculus target is drawn from."""
    c, o = rung.curve, rung.order
    odd = o
    while odd % 2 == 0:
        odd //= 2
    cof = o // odd
    H = None
    for x in range(1, 1 << rung.n):
        for P in c.points_over(x):
            Q = c.mul(P, cof)
            if Q is not None and rung._point_order(Q) == odd:
                H = Q
                break
        if H:
            break
    assert H is not None, "no generator of the odd part"
    elts, acc = {None}, None
    while True:
        acc = c.add(acc, H)
        if acc is None:
            break
        elts.add(acc)
    assert len(elts) == odd
    return elts, odd


def e5_dispersion(n: int, m: int, l: int):
    """Exact decomposition count for EVERY target, and its index of dispersion.

    The design asked for 450 targets a cell, enough to resolve a 20% effect at
    3 sigma.  This enumerates the whole odd subgroup instead -- 2,003 to 130,873
    targets -- so there is no sampling error to quote: `Var/mean` below is the
    population value, not an estimate of it.
    """
    rung = Rung(n)
    points, _, _ = build_base(rung, l)
    c = rung.curve
    subgroup, odd = odd_subgroup(rung)

    counts: dict = {}          # every m-subset
    useful: dict = {}          # m-subsets carrying a usable relation
    npts = len(points)
    neg = {P: c.neg(P) for P in points}

    def bump(S, degenerate):
        if S in subgroup:
            counts[S] = counts.get(S, 0) + 1
            if not degenerate:
                useful[S] = useful.get(S, 0) + 1

    if m == 2:
        for i in range(npts):
            Pi, Ni = points[i], neg[points[i]]
            for j in range(i + 1, npts):
                bump(c.add(Pi, points[j]), points[j] == Ni)
    elif m == 3:
        for i in range(npts):
            Pi, Ni = points[i], neg[points[i]]
            for j in range(i + 1, npts):
                Pj = points[j]
                Sij = c.add(Pi, Pj)
                pair_ij = Pj == Ni
                Nj = neg[Pj]
                for k in range(j + 1, npts):
                    Pk = points[k]
                    bump(c.add(Sij, Pk), pair_ij or Pk == Ni or Pk == Nj)
    else:
        raise ValueError("E5 runs at m = 2 and m = 3")

    def moments(tbl):
        t = sum(tbl.values())
        sq = sum(v * v for v in tbl.values())
        mu = t / odd
        return t, mu, sq / odd - mu * mu

    total, mean, var = moments(counts)
    utotal, umean, uvar = moments(useful)
    predicted = math.comb(npts, m) / rung.order
    return {
        "n": n, "m": m, "l": l,
        "factor_base_size": npts,
        "targets": odd,
        "sampling": "exhaustive over the odd-order subgroup",
        "decompositions_total": total,
        "mean_measured": round(mean, 6),
        "mean_predicted_C_over_order": round(predicted, 6),
        "mean_ratio": round(mean / predicted, 4) if predicted else None,
        "variance": round(var, 6),
        "index_of_dispersion_all_subsets": round(var / mean, 4) if mean else None,
        "degenerate_subsets": total - utotal,
        "useful_mean": round(umean, 6),
        "useful_variance": round(uvar, 6),
        "index_of_dispersion": round(uvar / umean, 4) if umean else None,
        "poisson_null": 1.0,
        "deviation_from_null": round(uvar / umean - 1.0, 4) if umean else None,
        "design_sample_standard_error": round(math.sqrt(2 / odd), 5),
        "note": "a subset containing a point and its negative sums to a target "
                "that is forced (the identity at m = 2, a base point at m = 3) "
                "and carries a trivial relation the matrix never sees.  Those "
                "subsets are a vanishing share of the mean and almost all of the "
                "variance, so the dispersion that matters is the useful one.",
    }


# ── E1: the scale-model ladder, end to end ──────────────────────────────

def ladder_dimension(n: int, m: int = 3) -> int:
    """`l = ceil((n + log2 m!)/m)`, the design's own choice.  It rounds up, so
    the base oversaturates and every target decomposes several times over."""
    return math.ceil((n + math.log2(math.factorial(m))) / m)


def individual_log(rung: Rung, l: int, logs, points, index, plus, rule: str,
                   pos, Q, rng: random.Random):
    """Descent: find `[c]Q + [d]G` over the base, then read `x` off the logs."""
    c = rung.curve
    member = set(points)
    npts = len(points)
    ops = 0
    for _ in range(4096):
        cc = rng.randrange(1, rung.p)
        dd = rng.randrange(1, rung.p)
        T = c.add(c.mul(Q, cc), c.mul(rung.G, dd))
        if T is None:
            continue
        negT = c.neg(T)
        for i in range(npts):
            hit = None
            for j in range(i + 1, npts):
                ops += 1
                S = c.add(points[i], points[j])
                if S is None:
                    continue
                W = c.add(negT, S)
                if W is None:
                    continue
                U = c.neg(W)
                # canonicalisation is what stops `full` recording each triple
                # once per (m-1)-subset of it; `first_hit` keeps one relation per
                # target anyway, so it does not need -- and must not impose -- it
                if U in member and (rule == "first_hit" or pos[U] > j):
                    hit = (points[i], points[j], U)
                    break
            if hit:
                total = 0
                for P in hit:
                    if rule == "folded":
                        k = index[P]
                        total += logs[k] if P == plus[k] else -logs[k]
                    else:
                        total += logs[pos[P]]
                total %= rung.p
                inv_proj = pow(rung.proj_scalar % rung.p, rung.p - 2, rung.p)
                x = (total * inv_proj - dd) % rung.p * pow(cc, rung.p - 2, rung.p)
                return x % rung.p, ops
    return None, ops


def run_e1_rung(n: int, m: int = 3, rule: str = "full", seed: int = 20260913,
                slack: int = 8, retries: int = 6, l_override: int | None = None):
    """One end-to-end rung: collect, solve, plant a logarithm, recover it.

    `l_override` runs the rung at a smaller factor-base dimension than the
    design's `ceil((n + log2 m!)/m)`.  That is a different cell of the same law,
    not a different law -- the product is flat in `l`, which is exactly what E2
    measures -- and at the larger rungs it is what makes the linear algebra fit:
    the oracle work is unchanged, but the matrix has `2^l` unknowns instead of
    `2^11`.  Rungs run this way are labelled with `l_is_design_value = False`.
    """
    rng = random.Random(seed + n)
    rung = Rung(n)
    if not rung.usable_end_to_end:
        return {
            "n": n, "m": m, "rule": rule, "l": ladder_dimension(n, m),
            "curve_order": rung.order, "group_shape": rung.group_shape,
            "log2_largest_prime": round(math.log2(rung.p), 2),
            "excluded": "p-torsion has rank 2, so no cyclic subgroup of order p "
                        "exists to pose the logarithm in",
            "p_torsion_rank": 2,
        }
    l = l_override if l_override is not None else ladder_dimension(n, m)
    points, index, abscissae = build_base(rung, l)
    pos = {P: i for i, P in enumerate(points)}
    plus = {}
    for P in points:
        plus.setdefault(index[P], P)
    unknowns = len(abscissae) if rule == "folded" else len(points)

    want = unknowns + slack
    for attempt in range(retries):
        d = collect(rung, l, m, rule, rng, want)
        sol, rank = solve_mod_p(d["rows"], d["rhs"], d["unknowns"], rung.p)
        if sol is not None:
            break
        want += max(8, unknowns // 4)
    assert sol is not None, f"n={n}: relations never determined the logs"

    # the planted logarithm, recovered by descent
    secret = rng.randrange(2, rung.p)
    Q = rung.curve.mul(rung.G, secret)
    got, descent_ops = individual_log(rung, l, sol, points, index, plus, rule,
                                      pos, Q, rng)
    assert got == secret, f"n={n}: planted log {secret} came back as {got}"

    ops = d["oracle_ops"] + descent_ops
    lam = ops / 2 ** n
    return {
        "n": n, "m": m, "rule": rule, "l": l,
        "l_design_value": ladder_dimension(n, m),
        "l_is_design_value": l == ladder_dimension(n, m),
        "curve_order": rung.order, "group_shape": rung.group_shape,
        "log2_largest_prime": round(math.log2(rung.p), 2),
        "subgroup_is_faithful": math.log2(rung.p) >= n - 4,
        "factor_base_size": len(points),
        "unknowns": d["unknowns"],
        "targets": d["targets"],
        "relations": len(d["rows"]),
        "matrix_rank": rank,
        "yield_measured": round(len(d["rows"]) / d["targets"], 3),
        "yield_predicted": round(math.comb(len(points), m) / rung.order, 3),
        "collection_oracle_ops": d["oracle_ops"],
        "descent_oracle_ops": descent_ops,
        "total_oracle_ops": ops,
        "group_additions": d["group_adds"],
        "log2_total_oracle_ops": round(math.log2(ops), 3),
        "lambda": round(lam, 3),
        "S": round(ops / math.sqrt(rung.p), 4),
        "planted_log_recovered": True,
        "relations_reverified_in_group": True,
        "trivial_relations_counted": 0,
    }


def least_squares_slope(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    den = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den if den else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-rung", type=int, default=19)
    args = ap.parse_args()

    # ---- the relation budget, and the stopping rule that decides it ----
    # The first round of E1/E2 stopped only when the relation matrix had full
    # rank over EVERY column.  That is the coupon-collector problem -- it costs
    # `|F| ln|F| / m` relations, not `|F|` -- and it is not what index calculus
    # needs.  Measured both ways here, because the difference is the whole of
    # the "2x" that round reported.
    crossings = [full_rank_crossing(Rung(n), l, 3, random.Random(9 + seed), rule)
                 for n, l in ((13, 6), (19, 7), (19, 8), (19, 9))
                 for rule in ("full", "first_hit") for seed in (0,)]
    for c in crossings:
        print(f"RHO n={c['n']:3d} l={c['l']} {c['rule']:10s} |F|={c['factor_base_size']:5d} "
              f"full rank at {c['full_rank_at']} (rho {c['rho_measured']}), coupon collector "
              f"predicts {c['coupon_collector_relations']} (rho {c['coupon_collector_rho']}), "
              f"ratio {c['measured_over_coupon']}", flush=True)
    budgets = [budget_run(Rung(n), l, 3, random.Random(100 + seed))
               for n, l in ((13, 6), (19, 7), (19, 8)) for seed in (1, 2, 3)]
    for b in budgets:
        print(f"BUDGET n={b['n']:3d} l={b['l']} |F|={b['factor_base_size']:5d} "
              f"determined {b['fraction_determined']:.3f} "
              f"(undetermined {b['fraction_undetermined']:.4f}) "
              f"descent attempts {b['descent_attempts']} landed {b['descent_landed']} "
              f"Lambda {b['lambda']:.3f}", flush=True)
    budget = {
        "question": "how many relations does the product law actually owe, and "
                    "what does its budget of |F| buy?",
        "superseded_claim": "the first round reported 'relations needed to reach "
                            "full rank run 2.0x |F|' and attributed it to "
                            "decompositions harvested from one target being "
                            "correlated.  Both halves are withdrawn: the excess "
                            "is the coupon collector, it is not a constant, and "
                            "the model does not owe it.",
        "full_rank_crossings": crossings,
        "budget_runs": budgets,
        "lambda_at_budget": [b["lambda"] for b in budgets],
        "mean_lambda_at_budget": round(sum(b["lambda"] for b in budgets)
                                       / len(budgets), 3),
        "predicted_lambda": 3,
        "every_descent_landed": all(b["descent_landed"] for b in budgets),
        "class": "accounting: the algorithm did not change, the harness's "
                 "stopping rule did, and the correction is to a number this "
                 "thread published",
    }

    # ---- E1 -----------------------------------------------------------
    e1_rows = []
    for n in (11, 13, 19, 29, 37):
        if n > args.max_rung:
            continue
        for rule in ("full", "first_hit", "folded"):
            row = run_e1_rung(n, rule=rule)
            e1_rows.append(row)
            if "excluded" in row:
                print(f"E1 n={n:3d} {rule:10s} EXCLUDED ({row['group_shape']})",
                      flush=True)
            else:
                print(f"E1 n={n:3d} {rule:10s} l={row['l']} |F|={row['factor_base_size']:5d} "
                      f"targets={row['targets']:5d} rels={row['relations']:5d} "
                      f"rank={row['matrix_rank']:5d} ops=2^{row['log2_total_oracle_ops']:6.2f} "
                      f"Lambda={row['lambda']:6.3f}  log recovered", flush=True)
    done = [r for r in e1_rows if "excluded" not in r]
    fits = {}
    for rule in ("full", "first_hit", "folded"):
        pts = [(r["n"], r["log2_total_oracle_ops"]) for r in done if r["rule"] == rule]
        if len(pts) >= 2:
            fits[rule] = {
                "rungs": [a for a, _ in pts],
                "slope": round(least_squares_slope([a for a, _ in pts],
                                                   [b for _, b in pts]), 4),
                "rungs_used": len(pts),
                "design_wanted_at_least": 4,
            }
    e1 = {
        "question": "does the product law actually hold on curves shaped like "
                    "ECC2K-130, or is the n = 131 derivation unsupported?",
        "stopping_rule_note": "the `rungs` below stop at full rank over every "
                              "column, which the relation-budget section shows "
                              "is the coupon collector and costs ln|F|/m times "
                              "the budget.  `lambda` there is inflated by that "
                              "factor; `budget_runs` carries the corrected one.",
        "primary_metric": "Lambda = oracle operations / 2^n, predicted flat at m = 3",
        "rungs": e1_rows,
        "slope_fits": fits,
        "falsifier": "least-squares slope of log2(total) against n below 0.95 "
                     "over four or more rungs, or any rung with Lambda < 0.5 m",
        "falsified": any(f["slope"] < 0.95 for f in fits.values())
                     or any(r["lambda"] < 1.5 for r in done),
        "correctness_gate": {
            "every_end_to_end_rung_recovered_its_planted_log":
                all(r["planted_log_recovered"] for r in done),
            "every_relation_reverified_in_the_group":
                all(r["relations_reverified_in_group"] for r in done),
            "trivial_relations_counted": 0,
        },
    }

    # ---- E2 -----------------------------------------------------------
    e2_rows = []
    for l in (5, 6, 7, 8, 9):
        for rule in ("full", "first_hit"):
            row = e2_sweep(19, l, rule=rule)
            e2_rows.append(row)
            print(f"E2 n=19 l={l} {rule:10s} |F|={row['factor_base_size']:4d} "
                  f"lam={row['yield_predicted']:9.3f} targets={row['targets']:6d} "
                  f"rank={row['matrix_rank']}/{row['unknowns']} "
                  f"Lambda={row['lambda_metric']:9.3f}", flush=True)
    flat = [r["lambda_metric"] for r in e2_rows if r["rule"] == "full"]
    # The full-rank sweep above tilts with `ln|F|/m`, which runs from 1.12 at
    # `|F| = 29` to 2.09 at `|F| = 527` -- a 1.87x drift that is the harness's,
    # not the dimension's, and the same size as the swing the flatness verdict
    # reads.  Re-swept at the budget the product law actually owes.
    e2_budget = [budget_run(Rung(19), l, 3, random.Random(200 + l))
                 for l in (5, 6, 7, 8, 9)]
    for b in e2_budget:
        print(f"E2-BUDGET n=19 l={b['l']} |F|={b['factor_base_size']:5d} "
              f"determined {b['fraction_determined']:.3f} "
              f"descent attempts {b['descent_attempts']} landed {b['descent_landed']} "
              f"Lambda {b['lambda']:.3f}", flush=True)
    # Why `l = 5` stalled at rank 28/29: exhaustively, not by sampling.
    ceilings = [e2_rank_ceiling(19, l) for l in (5, 6, 7)]
    for cc in ceilings:
        print(f"E2-CEILING n=19 l={cc['l']} |F|={cc['factor_base_size']:5d} "
              f"reachable={cc['reachable_decompositions']:7d} "
              f"rank ceiling {cc['rank_ceiling']}/{cc['factor_base_size']} "
              f"odd-abscissa census {cc['odd_abscissa_census']}", flush=True)
    flat_budget = [b["lambda"] for b in e2_budget]
    flat_budget_line = sum(flat_budget) / len(flat_budget)
    e2 = {
        "question": "is the total really flat in the factor-base dimension?",
        "rungs": e2_rows,
        "flat_line_lambda": round(sum(flat) / len(flat), 3) if flat else None,
        "flat_line_lambda_stopping_rule": "full rank over every column -- "
                                          "superseded, kept as the before mark",
        "budget_rule_rungs": e2_budget,
        "flat_line_lambda_at_budget": round(flat_budget_line, 3),
        "rank_ceilings": ceilings,
        "l5_outlier_explained": "not relation correlation: at l = 5 the row "
                                "space reachable over every one of the 810 "
                                "decompositions the base admits has rank 28 of "
                                "29, because every such triple uses exactly one "
                                "even-abscissa point.  No relation count "
                                "reaches full rank there.  From l = 6 the base "
                                "admits all-even triples and the rank is full.",
        "falsifier": "any dimension whose measured total is below half the flat line",
        "falsified": min(flat_budget) < 0.5 * flat_budget_line,
        "falsified_under_superseded_rule": bool(flat)
            and min(flat) < 0.5 * (sum(flat) / len(flat)),
    }

    # ---- E4 -----------------------------------------------------------
    e4_rows = []
    for l, lp, t in ((7, 10, 6000), (7, 11, 6000), (8, 11, 4000)):
        row = e4_large_prime(19, l, lp, targets=t)
        e4_rows.append(row)
        print(f"E4 n=19 l={l} l'={lp} partials={row['partial_relations']:6d} "
              f"yield/t={row['partial_yield_per_target']:6.3f} "
              f"guard={row['guard_satisfied']} paired={row['paired_full_relations']:5d} "
              f"rank={row['matrix_rank']:4d} rank/rels={row['rank_over_relations']} "
              f"rank/unk={row['rank_over_unknowns']}", flush=True)
    e4 = {
        "question": "do large primes move the product law, or only relabel it?",
        "rungs": e4_rows,
        "falsifier": "a guarded cell below the BSGS line at the same memory",
        "measured_caveat": "the design flagged relation independence as the most "
                           "likely place for the model to be wrong, and required "
                           "the rank to be measured rather than assumed",
    }

    # ---- E5 -----------------------------------------------------------
    e5_rows = []
    for n, m, l in ((13, 2, 7), (13, 3, 5), (13, 3, 6), (17, 2, 8), (17, 3, 6),
                    (19, 2, 9), (19, 3, 6), (19, 3, 7)):
        row = e5_dispersion(n, m, l)
        e5_rows.append(row)
        print(f"E5 n={n} m={m} l={l} T={row['targets']:6d} "
              f"mean={row['useful_mean']:.4f} (pred {row['mean_predicted_C_over_order']:.4f}) "
              f"D_all={row['index_of_dispersion_all_subsets']:.3f} "
              f"D={row['index_of_dispersion']:.4f}", flush=True)
    off = [r for r in e5_rows if abs(r["deviation_from_null"]) > 0.2]
    e5 = {
        "question": "is the yield Poisson, or is the note's central table resting "
                    "on the wrong tail?",
        "cells": e5_rows,
        "cells_off_null_by_more_than_0_2": [
            {"n": r["n"], "m": r["m"], "l": r["l"],
             "index_of_dispersion": r["index_of_dispersion"]} for r in off],
        "falsifier": "|Var/mean - 1| > 0.2 on three or more cells with a "
                     "consistent sign",
        "falsified": len(off) >= 3 and (
            all(r["deviation_from_null"] > 0 for r in off)
            or all(r["deviation_from_null"] < 0 for r in off)),
    }

    # ---- E6 -----------------------------------------------------------
    e6_rows = []
    for n, k, mult in ((13, 6, 20), (19, 5, 40), (19, 8, 30), (19, 12, 20)):
        row = e6_orbit_union(n, k, relations_multiple=mult)
        e6_rows.append(row)
        print(f"E6 n={n} orbits={k:2d} |F|={row['factor_base_size']:4d} "
              f"unknowns={row['orbit_unknowns']:3d} rels={row['relations']:4d} "
              f"rank={row['matrix_rank']:3d} solved={row['solved']} "
              f"frobenius_saving={row['realised_frobenius_saving']} "
              f"/n={row['frobenius_saving_over_n']}", flush=True)
    bad = [r for r in e6_rows if r["frobenius_saving_over_n"] is None
           or abs(r["frobenius_saving_over_n"] - 1.0) > 0.2]
    e6 = {
        "question": "do |F|/n orbit relations really carry n independent "
                    "unknowns, or does the collapse lose rank?",
        "rungs": e6_rows,
        "falsifier": "a realised saving differing from n by more than 20%",
        "falsified": bool(bad),
    }

    e3 = {
        "question": "how much of the difficulty is deciding, and how much is "
                    "producing the witness?",
        "status": "run in the previous round and recorded in section 3.2 of "
                  "research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md: a whole-base detector "
                  "localises its own witness by swapping, measured over eight "
                  "rungs with zero sub-base queries",
        "artefact": "experiments/ecc2k130_point_decomposition.json "
                    "-> swap_localisation",
        "item_1_now_run": "the design's first item -- does a real solver charge "
                          "the same for R - P + Q as for R -- is measured by "
                          "scripts/ecc2k130_e3_solver_panel.py and recorded in "
                          "experiments/ecc2k130_e3_solver_panel.json.  It is a "
                          "stage diagnostic under AGENTS.md section 8: one "
                          "oracle call on one rung, priced, with nothing "
                          "inferred about a full ECDLP and no speedup claimed.",
    }

    report = {
        "schema": "ecc2k130_decomposition_runs/v1",
        "purpose": "the six pre-registered experiments of "
                   "research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION_TARGETS.md, run",
        "unit": "Lambda = oracle operations / 2^n; an oracle operation is one "
                "(m-1)-subset enumerated",
        "boundaries_frozen_before_the_runs":
            "experiments/ecc2k130_decomposition_targets.json",
        "log2_rho_reference_at_131": LOG2_RHO_131,
        "relation_budget_and_stopping_rule": budget,
        "E1_scale_model_ladder": e1,
        "E2_dimension_flatness": e2,
        "E3_detector_target_agnosticism": e3,
        "E4_large_primes": e4,
        "E5_yield_distribution": e5,
        "E6_orbit_union_bases": e6,
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
