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
    found, fp, degen = F4.verify_quadruples(E, reps, n, groups)
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
    found, fp, degen = F4.verify_quadruples(E, reps, len(reps), groups)
    assert found, "control base carried no four-point relations"
    for quad in found:
        pts = [(int(a, 16), int(b, 16)) for a, b in quad]
        assert len({p[0] for p in pts}) == 4, "degenerate quad reported"
        for p in pts:
            assert E.on_curve(p)
        assert E.sum_points(pts) is None


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
    found, fp, degen = F4.verify_quadruples(E, reps, len(reps), groups)
    for quad in found:
        pts = [(int(a, 16), int(b, 16)) for a, b in quad]
        # no reported relation may be {P, Q, -P, -Q}
        negs = {tuple(E.neg(p)) for p in pts}
        overlap = len(negs & {tuple(p) for p in pts})
        assert overlap == 0, f"trivial negation quad reported: {quad}"


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
