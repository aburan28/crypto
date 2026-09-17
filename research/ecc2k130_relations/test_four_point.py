#!/usr/bin/env python3
"""The four-point collision sweep, checked where four-point relations exist.

The ECC2K-130 run is expected to find nothing -- the counting floor says so
in advance -- so the only thing that makes it evidence is that the same code
finds relations on a curve small enough to carry them, and finds all of them.
"""

import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "scripts"))

import relations as R                                   # noqa: E402
import run_four_point as F4                             # noqa: E402
from test_relations import small_curve, subgroup_generator   # noqa: E402


def test_decode_pair_inverts_the_enumeration():
    n = 9
    seen = []
    for i in range(n):
        for j in range(i + 1, n):
            seen.append((i, j, 1))
            seen.append((i, j, -1))
    for idx, expect in enumerate(seen):
        assert F4.decode_pair(idx, n) == expect, (idx, expect)


def test_representatives_are_one_point_per_abscissa():
    F, E, order, r = small_curve(13)
    xs = [x for x in range(1 << 13) if E.has_point(x)][:300]
    reps, kept = F4.representative_points(E, xs)
    assert len(reps) == len(kept) == len(xs)
    assert len({p[0] for p in reps}) == len(reps)
    for p in reps:
        assert E.on_curve(p)


def test_sweep_finds_every_four_point_relation_on_a_small_base():
    """Complete against brute force over all sign-canonical quadruples."""
    F, E, order, r = small_curve(13)
    rng = random.Random(31)
    xs = [x for x in range(1 << 13) if E.has_point(x)]
    base = sorted(rng.sample(xs, 70))
    reps, kept = F4.representative_points(E, base)
    n = len(reps)

    keys, consumed, complete, _ = F4.pair_sum_sweep(E, reps, log=lambda *a: None)
    assert complete
    groups = F4.find_collisions(keys)
    cls = F4.verify_quadruples(E, reps, n, groups)
    found = cls["distinct"]
    got = {tuple(sorted(tuple(p) for p in q)) for q in
           [[(int(a, 16), int(b, 16)) for a, b in quad] for quad in found]}

    # brute force: all four distinct representatives, all sign patterns
    brute = set()
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for l in range(k + 1, n):
                    for s in range(8):
                        quad = [reps[i],
                                reps[j] if s & 1 else E.neg(reps[j]),
                                reps[k] if s & 2 else E.neg(reps[k]),
                                reps[l] if s & 4 else E.neg(reps[l])]
                        if E.sum_points(quad) is None:
                            brute.add(tuple(sorted(tuple(p) for p in quad)))
    # a relation and its global negation are the same relation
    def canon(q):
        neg = tuple(sorted(tuple(E.neg(p)) for p in q))
        return min(q, neg)
    assert {canon(q) for q in got} == {canon(q) for q in brute}, (
        len(got), len(brute))


def test_every_reported_relation_actually_sums_to_the_identity():
    F, E, order, r = small_curve(13)
    rng = random.Random(32)
    xs = [x for x in range(1 << 13) if E.has_point(x)]
    base = sorted(rng.sample(xs, 120))
    reps, kept = F4.representative_points(E, base)
    keys, _, complete, _ = F4.pair_sum_sweep(E, reps, log=lambda *a: None)
    groups = F4.find_collisions(keys)
    cls = F4.verify_quadruples(E, reps, len(reps), groups)
    found = cls["distinct"] + cls["repeated"] + cls["frobenius"]
    assert found, "control base carried no four-point relations"
    for quad in found:
        pts = [(int(a, 16), int(b, 16)) for a, b in quad]
        for p in pts:
            assert E.on_curve(p)
        assert E.sum_points(pts) is None
    # only the `distinct` class promises four distinct summands; `repeated`
    # is the class that exists precisely because some relations repeat one
    for quad in cls["distinct"]:
        pts = [(int(a, 16), int(b, 16)) for a, b in quad]
        assert len({p[0] for p in pts}) == 4, "degenerate quad in distinct"


def test_the_negation_twin_is_not_reported_as_a_relation():
    """The trap the sign-canonical enumeration exists to avoid.

    A full negation-closed base makes `{P, Q}` and `{-P, -Q}` collide on the
    abscissa of their sums for every pair. The representative enumeration
    must not produce those at all.
    """
    F, E, order, r = small_curve(13)
    rng = random.Random(33)
    xs = [x for x in range(1 << 13) if E.has_point(x)]
    base = sorted(rng.sample(xs, 80))
    reps, _ = F4.representative_points(E, base)
    keys, _, _, _ = F4.pair_sum_sweep(E, reps, log=lambda *a: None)
    groups = F4.find_collisions(keys)
    cls = F4.verify_quadruples(E, reps, len(reps), groups)
    for quad in cls["distinct"]:
        pts = [(int(a, 16), int(b, 16)) for a, b in quad]
        # no reported relation may be {P, Q, -P, -Q}
        negs = {tuple(E.neg(p)) for p in pts}
        overlap = len(negs & {tuple(p) for p in pts})
        assert overlap == 0, f"trivial negation quad reported: {quad}"


def test_frobenius_characteristic_equation_holds_on_every_point():
    """sigma^2 + sigma + 2 = 0, the identity that makes the repeated-summand
    relations free. #E(F_2) = 4 gives trace t = -1, so the Frobenius
    characteristic equation sigma^2 - t sigma + 2 is sigma^2 + sigma + 2."""
    F, E, P, Q, r = R.challenge_curve()
    rng = random.Random(41)
    checked = 0
    while checked < 12:
        pts = E.points_over(rng.getrandbits(131))
        if not pts:
            continue
        for X in pts:
            s1 = E.frobenius(X)
            s2 = E.frobenius(s1)
            assert E.sum_points([s2, s1, X, X]) is None
        checked += 1
    for X in (P, Q):
        s1 = E.frobenius(X)
        s2 = E.frobenius(s1)
        assert E.sum_points([s2, s1, X, X]) is None


def _any_point(E, start=2):
    """The first abscissa from `start` that carries a point."""
    for x in range(start, 100000):
        pts = E.points_over(x)
        if pts:
            return pts[0]
    raise AssertionError("no point found")


def test_frobenius_pattern_reads_off_signed_exponents():
    F, E, P, Q, r = R.challenge_curve()
    R0 = _any_point(E)
    s1, s2 = E.frobenius(R0), E.frobenius(E.frobenius(R0))
    pat = F4.frobenius_pattern(E, [R0, R0, s1, s2])
    assert pat == [(0, 1), (0, 1), (1, 1), (2, 1)], pat
    pat = F4.frobenius_pattern(E, [R0, E.neg(s1)])
    assert pat == [(0, 1), (1, -1)], pat
    # a point outside the orbit is rejected rather than forced
    other = None
    for x in range(3, 4000):
        cand = E.points_over(x)
        if cand and F4.frobenius_pattern(E, [R0, cand[0]]) is None:
            other = cand[0]
            break
    assert other is not None
    assert F4.frobenius_pattern(E, [R0, other]) is None


def test_is_frobenius_identity_accepts_the_characteristic_equation():
    """Both shapes it generates, and nothing that merely looks like them."""
    F, E, P, Q, r = R.challenge_curve()
    R0 = _any_point(E)
    s1 = E.frobenius(R0)
    s2 = E.frobenius(s1)
    s3 = E.frobenius(s2)
    # sigma^2 + sigma + 2 = 0
    assert E.sum_points([R0, R0, s1, s2]) is None
    assert F4.is_frobenius_identity(E, [R0, R0, s1, s2])
    # sigma^3 = 2 - sigma, so 2R - sigma(R) - sigma^3(R) = O
    assert E.sum_points([R0, R0, E.neg(s1), E.neg(s3)]) is None
    assert F4.is_frobenius_identity(E, [R0, R0, E.neg(s1), E.neg(s3)])
    # a pattern that does not annihilate every point is not an identity
    assert not F4.is_frobenius_identity(E, [R0, R0, s1, s1])
    assert not F4.is_frobenius_identity(E, [R0, s1])


def test_the_stable_support_yields_frobenius_relations_and_no_others():
    """The structural finding, on the support built to carry it.

    A sigma-stable support contains x, x^2 and x^4 together for some of its
    elements, and each such element gives 2R + sigma(R) + sigma^2(R) = O for
    free. These are genuine four-point relations occurring astronomically
    above the random rate -- and every one of them is an instance of the
    characteristic equation, so none is evidence about a logarithm.
    """
    F, E, P, Q, r = R.challenge_curve()
    stable = [(1 << i) ^ (1 << j) for i in range(66) for j in range(i + 1, 66)]
    xs = R.base_from_abscissae(E, stable)
    reps, kept = F4.representative_points(E, xs)
    keys, _, complete, _ = F4.pair_sum_sweep(E, reps, log=lambda *a: None)
    assert complete
    cls = F4.verify_quadruples(E, reps, len(reps), F4.find_collisions(keys))
    assert cls["repeated"], "expected the free Frobenius relations, found none"
    assert len(cls["repeated"]) > 50, len(cls["repeated"])
    assert len(cls["frobenius"]) == len(cls["repeated"]), (
        "a repeated-summand relation appeared that is not an identity of the "
        "endomorphism ring -- that would be a genuinely new structure and "
        "must not be silently folded in")
    assert not cls["distinct"], (
        "a four-distinct relation at this base size would contradict the "
        "counting floor by ~2^89 and needs investigating, not asserting")


def test_digest_is_not_the_low_bits_of_x():
    """The mix must actually move the structured low bits."""
    assert F4.digest_of(0) == 0
    low_equal = sum(1 for x in range(1, 4000)
                    if F4.digest_of(x) == (x & F4.DIGEST_MASK))
    assert low_equal == 0


def test_budget_stops_the_sweep_and_marks_it_incomplete():
    F, E, order, r = small_curve(17)
    rng = random.Random(34)
    xs = [x for x in range(1 << 17) if E.has_point(x)][:1200]
    reps, _ = F4.representative_points(E, xs)
    keys, consumed, complete, elapsed = F4.pair_sum_sweep(
        E, reps, budget_seconds=1.0, chunk=2000, log=lambda *a: None)
    assert not complete
    assert consumed < len(reps) * (len(reps) - 1)


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
