#!/usr/bin/env python3
"""Correctness for the relation machinery, and the small controls.

The controls are the point of this file.  A sweep that finds no relations
is only evidence if the same code finds them when they are there, so every
negative result in this study is backed by a positive control on a curve
small enough to check exhaustively.
"""

import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "scripts"))

from fastfield import FastGF2m                                      # noqa: E402
import relations as R                                               # noqa: E402
from ecc2k130_point_decomposition import find_irreducible, curve_order, is_prime  # noqa: E402


def small_curve(m: int):
    """Koblitz `y^2 + xy = x^3 + 1` over `F_2^m`, with its prime subgroup."""
    F = FastGF2m(m, find_irreducible(m))
    E = R.Koblitz(F, a=0, b=1)
    order = curve_order(m, trace=-1)
    assert order % 4 == 0
    r = order // 4
    return F, E, order, r


def subgroup_generator(E, order, r, rng):
    for _ in range(500):
        x = rng.getrandbits(E.F.deg)
        pts = E.points_over(x)
        if not pts:
            continue
        G = E.mul(pts[0], order // r)
        if G is not None and E.mul(G, r) is None:
            return G
    raise AssertionError("no generator found")


# ── the field/curve layer ────────────────────────────────────────────────

def test_challenge_points_are_on_the_curve_and_in_the_subgroup():
    F, E, P, Q, r = R.challenge_curve()
    assert E.on_curve(P) and E.on_curve(Q)
    assert E.mul(P, r) is None, "P must have order dividing r"
    assert E.mul(Q, r) is None, "Q must have order dividing r"
    assert is_prime(r)
    assert curve_order(131, trace=-1) == 4 * r


def test_points_over_agrees_with_has_point():
    F, E, _, _ = small_curve(13)[0], None, None, None
    F, E, order, r = small_curve(13)
    for x in range(1 << 13):
        pts = E.points_over(x)
        assert E.has_point(x) == bool(pts)
        for p in pts:
            assert E.on_curve(p) and p[0] == x


def test_curve_order_matches_an_exhaustive_count():
    F, E, order, r = small_curve(13)
    count = 1                                     # the point at infinity
    for x in range(1 << 13):
        count += len(E.points_over(x))
    assert count == order


def test_frobenius_is_an_endomorphism():
    F, E, order, r = small_curve(13)
    rng = random.Random(1)
    G = subgroup_generator(E, order, r, rng)
    for _ in range(20):
        A = E.mul(G, rng.randrange(1, r))
        B = E.mul(G, rng.randrange(1, r))
        assert E.frobenius(E.add(A, B)) == E.add(E.frobenius(A), E.frobenius(B))
        assert E.on_curve(E.frobenius(A))


# ── S3 and the three-point sweep ─────────────────────────────────────────

def test_s3_vanishes_exactly_on_collinear_triples():
    F, E, order, r = small_curve(11)
    rng = random.Random(2)
    xs = [x for x in range(1 << 11) if E.has_point(x)]
    checked_zero = checked_nonzero = 0
    for _ in range(300):
        x1, x2 = rng.choice(xs), rng.choice(xs)
        for x3 in rng.sample(xs, 6):
            val = E.s3(x1, x2, x3)
            realised = R.realise_triple(E, x1, x2, x3)
            if val == 0:
                # S3 = 0 must mean some sign choice sums to O
                assert realised is not None, (x1, x2, x3)
                checked_zero += 1
            else:
                assert realised is None, (x1, x2, x3)
                checked_nonzero += 1
    assert checked_zero > 0 and checked_nonzero > 0


def test_solve_s3_last_returns_genuine_roots():
    F, E, order, r = small_curve(11)
    rng = random.Random(3)
    xs = [x for x in range(1 << 11) if E.has_point(x)]
    pairs = [(rng.choice(xs), rng.choice(xs)) for _ in range(400)]
    got = R.solve_s3_last(E, pairs)
    for (x1, x2), roots in zip(pairs, got):
        for x3 in roots:
            assert E.s3(x1, x2, x3) == 0, (x1, x2, x3)


def test_solve_s3_last_is_complete_against_brute_force():
    """Every root the batched solver reports, and no root it misses."""
    F, E, order, r = small_curve(9)
    allx = list(range(1 << 9))
    rng = random.Random(4)
    pairs = [(rng.randrange(1 << 9), rng.randrange(1 << 9)) for _ in range(120)]
    got = R.solve_s3_last(E, pairs)
    for (x1, x2), roots in zip(pairs, got):
        brute = {x for x in allx if E.s3(x1, x2, x) == 0}
        assert set(roots) == brute, ((x1, x2), sorted(roots), sorted(brute))


def test_three_point_sweep_matches_brute_force_on_a_small_base():
    F, E, order, r = small_curve(11)
    rng = random.Random(5)
    xs = [x for x in range(1 << 11) if E.has_point(x)]
    base = sorted(rng.sample(xs, 60))
    found, pairs, complete = R.three_point_relations(E, base, chunk=97)
    assert complete
    brute = set()
    n = len(base)
    for i in range(n):
        for j in range(i, n):
            for k in range(j, n):
                if E.s3(base[i], base[j], base[k]) == 0:
                    brute.add((base[i], base[j], base[k]))
    assert set(found) == brute, (sorted(set(found) ^ brute))


def test_every_reported_relation_realises_on_the_curve():
    F, E, order, r = small_curve(13)
    rng = random.Random(6)
    xs = [x for x in range(1 << 13) if E.has_point(x)]
    base = sorted(rng.sample(xs, 200))
    found, _, complete = R.three_point_relations(E, base)
    assert complete and found, "control base must produce relations"
    for t in found:
        pts = R.realise_triple(E, *t)
        assert pts is not None
        assert E.sum_points(list(pts)) is None
        for p in pts:
            assert E.on_curve(p)


# ── batched group operations ─────────────────────────────────────

def test_batch_add_matches_scalar_add_including_degenerate_cases():
    F, E, order, r = small_curve(17)
    rng = random.Random(21)
    G = subgroup_generator(E, order, r, rng)
    As = [E.mul(G, rng.randrange(1, r)) for _ in range(60)]
    Bs = [E.mul(G, rng.randrange(1, r)) for _ in range(60)]
    # plant every degenerate case the batch path branches on
    As[0], Bs[0] = None, Bs[0]
    As[1], Bs[1] = As[1], None
    As[2], Bs[2] = As[2], E.neg(As[2])          # sums to the identity
    As[3], Bs[3] = As[3], As[3]                 # doubling
    got = R.batch_add(E, As, Bs)
    for A, B, g in zip(As, Bs, got):
        assert g == E.add(A, B), (A, B, g, E.add(A, B))
        assert E.on_curve(g)


def test_batch_combos_matches_one_at_a_time():
    F, E, order, r = small_curve(17)
    rng = random.Random(22)
    G = subgroup_generator(E, order, r, rng)
    k = rng.randrange(2, r)
    Q = E.mul(G, k)
    pairs = [(rng.randrange(1, r), rng.randrange(1, r)) for _ in range(40)]
    pairs += [(1, 0), (0, 1), (1, 1), (2, 3)]
    got = R.batch_combos(E, G, Q, pairs)
    for (u, v), g in zip(pairs, got):
        assert g == E.add(E.mul(G, u), E.mul(Q, v)), (u, v)


# ── useful rank ──────────────────────────────────────────────────────────

def test_construction_rows_are_recognised_as_carrying_nothing():
    r = 1009
    rows = [(0, 0), (0, 0), (0, 0)]
    out = R.useful_rank(rows, r)
    assert out["rows"] == 3
    assert out["construction_rows"] == 3
    assert out["determining_rows"] == 0
    assert out["useful_rank"] == 0
    assert out["log"] is None


def test_a_determining_row_recovers_the_scalar():
    r = 1009
    k = 421
    # su + k*sv = 0 mod r
    sv = 7
    su = (-k * sv) % r
    out = R.useful_rank([(su, sv)], r)
    assert out["useful_rank"] == 1
    assert out["log"] == k


def test_mixed_rows_report_only_the_determining_ones():
    r = 1009
    k = 421
    rows = [(0, 0), ((-k * 3) % r, 3), (0, 0)]
    out = R.useful_rank(rows, r)
    assert out["construction_rows"] == 2
    assert out["determining_rows"] == 1
    assert out["log"] == k


# ── the small controls: recover a planted scalar end to end ──────────────

def _control(m: int, size: int, seed: int):
    F, E, order, r = small_curve(m)
    rng = random.Random(seed)
    G = subgroup_generator(E, order, r, rng)
    k = rng.randrange(2, r)
    Q = E.mul(G, k)
    base, coeffs = R.constructed_base(E, G, Q, r, size, rng)
    found, _, complete = R.three_point_relations(E, base)
    assert complete
    rows = R.relation_vectors(E, found, coeffs, r)
    out = R.useful_rank(rows, r)
    return E, G, Q, k, r, found, rows, out


def test_small_control_m17_recovers_the_planted_scalar():
    E, G, Q, k, r, found, rows, out = _control(17, 90, seed=11)
    assert found, "control found no relations at all"
    assert out["determining_rows"] > 0, "relations found but none determining"
    assert out["log"] == k, (out["log"], k)
    assert R.verify_log(E, G, Q, out["log"])


def test_small_control_m19_recovers_the_planted_scalar():
    E, G, Q, k, r, found, rows, out = _control(19, 150, seed=12)
    assert found
    assert out["log"] == k
    assert R.verify_log(E, G, Q, out["log"])


def test_control_relations_are_consistent_with_each_other():
    E, G, Q, k, r, found, rows, out = _control(17, 90, seed=13)
    assert out["consistent"], "two relations disagreed on the logarithm"


def test_a_base_built_from_small_coefficients_produces_construction_rows():
    """The accounting result, on a curve where it can be checked by hand.

    Small-coefficient bases carry relations that hold whatever the logarithm
    is.  They must show up as construction rows, not as useful rank.
    """
    F, E, order, r = small_curve(17)
    rng = random.Random(14)
    G = subgroup_generator(E, order, r, rng)
    k = rng.randrange(2, r)
    Q = E.mul(G, k)
    base, coeffs = R.constructed_base(E, G, Q, r, 120, rng, small=10)
    found, _, complete = R.three_point_relations(E, base)
    assert complete and found
    rows = R.relation_vectors(E, found, coeffs, r)
    out = R.useful_rank(rows, r)
    assert out["construction_rows"] > 0, "expected free relations, found none"
    # and whatever is left must still be honest about the logarithm
    if out["log"] is not None:
        assert out["log"] == k


# ── base builders ────────────────────────────────────────────────────────

def test_weight_two_set_is_complete_and_distinct():
    F = FastGF2m(131, R.IRR131)
    xs = R.weight_two_abscissae(F)
    assert len(xs) == 131 * 130 // 2
    assert len(set(xs)) == len(xs)
    assert all(bin(x).count("1") == 2 for x in xs)


def test_weight_two_set_is_not_frobenius_stable_in_the_polynomial_basis():
    """A correction to the obvious guess, and it constrains the sweeps.

    `z^i + z^j` squares to `z^{2i} + z^{2j}`, which stays weight two only
    while both exponents stay below the degree.  Once `2i >= 131` the
    reduction `z^131 = z^13 + z^2 + z + 1` fires and the image has weight
    four or more.  So the weight-two set is stable exactly on the pairs
    with `i < j <= 65`, which is `C(66, 2) = 2145` of its 8515 elements.

    This matters beyond bookkeeping: it rules out normalising a relation
    by Frobenius to put one summand at an orbit representative *of the
    base*, because the base is not a union of orbits.  A normal basis
    would make the low-weight set stable; the challenge's polynomial basis
    does not.
    """
    F = FastGF2m(131, R.IRR131)
    xs = set(R.weight_two_abscissae(F))
    stable = {x for x in xs if F.sqr(x) in xs}
    assert len(stable) == 2145, len(stable)
    expected = {(1 << i) ^ (1 << j) for i in range(66) for j in range(i + 1, 66)}
    assert stable == expected
    assert len(xs) - len(stable) == 6370


def test_frobenius_orbit_of_a_weight_two_element_leaves_the_set():
    """Concretely: iterating sigma from inside the stable part still escapes."""
    F = FastGF2m(131, R.IRR131)
    xs = set(R.weight_two_abscissae(F))
    x = (1 << 60) ^ (1 << 65)
    assert x in xs
    seen_outside = False
    cur = x
    for _ in range(131):
        cur = F.sqr(cur)
        if cur not in xs:
            seen_outside = True
    assert seen_outside
    assert cur == x, "sigma^131 must be the identity"


def test_random_matched_base_is_the_requested_size():
    F, E, order, r = small_curve(17)
    rng = random.Random(15)
    b = R.random_matched_base(E, 50, rng)
    assert len(b) == 50 and len(set(b)) == 50
    assert all(E.has_point(x) for x in b)


if __name__ == "__main__":
    import traceback
    fails = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
            print(f"ok   {name}")
        except Exception as exc:
            fails += 1
            print(f"FAIL {name}: {exc}")
            traceback.print_exc()
    print("failures:", fails)
    sys.exit(1 if fails else 0)
