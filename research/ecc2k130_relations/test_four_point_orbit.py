#!/usr/bin/env python3
"""Controls for the orbit-grouped four-point sweep.

The ECC2K-130 run finds nothing, and the counting floor says so in advance,
so the run is evidence only insofar as the pipeline is checked where the
answers are known.  Each defect that once produced plausible output has a
test here that fails if it comes back.
"""

import random
import sys
from itertools import combinations, product
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import relations as R                                        # noqa: E402
import smallcurve as SC                                      # noqa: E402
from fourpoint import (configurations, exhaustive_four_point,  # noqa: E402
                       supply_prediction)
from normalbasis import (NormalSupport, conjugates,          # noqa: E402
                         find_normal_elements, is_normal)
from planted import frobenius_scalar, recover_planted        # noqa: E402

SEED = 20260917


def _poly_is_irreducible(p: int, deg: int) -> bool:
    """`x^(2^deg) = x mod p` and no proper subfield factor."""
    def mulmod(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
        while r.bit_length() - 1 >= deg:
            r ^= p << (r.bit_length() - 1 - deg)
        return r

    def frob(e):
        v = 2
        for _ in range(e):
            v = mulmod(v, v)
        return v

    if frob(deg) != 2:
        return False
    for q in (f for f in range(2, deg) if deg % f == 0):
        g, h = frob(deg // q) ^ 2, p
        while g:
            if g.bit_length() > h.bit_length():
                g, h = h, g
            if not g:
                break
            h ^= g << (h.bit_length() - g.bit_length())
            if h.bit_length() < g.bit_length():
                g, h = h, g
        if h != 1:
            return False
    return True


def test_small_field_polynomials_are_irreducible():
    for deg, poly in SC.IRR.items():
        assert _poly_is_irreducible(poly, deg), f"m={deg} polynomial reducible"


def test_curve_order_recursion_matches_an_exhaustive_count():
    for m in (11, 13):
        _, E, order = SC.small_curve(m)
        assert order == SC.count_points(E), f"order mismatch at m={m}"


def test_four_torsion_is_cyclic_of_order_four():
    for m in (11, 13, 17, 19):
        _, E, order = SC.small_curve(m)
        e4 = SC.four_torsion(E, order)
        assert len(e4) == 4 and len(set(map(str, e4))) == 4
        for T in e4:
            assert E.mul(T, 4) is None


def test_normal_elements_have_independent_conjugates():
    for m in (11, 13):
        F, _, _ = SC.small_curve(m)
        rng = random.Random(SEED)
        for a in find_normal_elements(F, 5, rng):
            assert is_normal(F, a)
            assert len(set(conjugates(F, a))) == m
        assert not is_normal(F, 0)
        assert not is_normal(F, 1)          # 1 is fixed by sigma


def test_weight_two_normal_support_is_frobenius_stable():
    """The premise #404 got backwards, checked rather than argued."""
    for m in (11, 13, 19):
        F, E, _ = SC.small_curve(m)
        rng = random.Random(SEED)
        for a in find_normal_elements(F, 4, rng):
            s = NormalSupport(F, E, a)
            assert s.is_sigma_stable()
            assert not s.mixed_orbits, "lifting must be constant on an orbit"
            assert s.size == s.m * s.orbit_count
            assert s.total_orbits() == (m - 1) // 2


def test_challenge_support_is_frobenius_stable():
    F, E, _, _, _ = R.challenge_curve()
    rng = random.Random(SEED)
    s = NormalSupport(F, E, find_normal_elements(F, 1, rng)[0])
    assert s.is_sigma_stable()
    assert s.total_orbits() == 65, "C(131,2)/131 = 65 orbits, not 66"
    assert not s.mixed_orbits
    assert s.size % 131 == 0


def test_pair_orbits_are_enumerated_exactly_once():
    """Regression for the double count that reported 125,953 'relations'."""
    for m in (11, 13):
        F, E, _ = SC.small_curve(m)
        rng = random.Random(SEED)
        for a in find_normal_elements(F, 3, rng):
            s = NormalSupport(F, E, a)
            if s.size == 0:
                continue
            pairs = list(s.canonical_pair_orbits())
            assert len(pairs) == s.pair_orbit_count()
            assert len(pairs) == s.size * (s.size - 1) // s.m
            keys = set()
            for ia, ib, eps in pairs:
                t1, k1 = divmod(ia, s.m)
                t2, k2 = divmod(ib, s.m)
                keys.add(min(
                    tuple(sorted(((t1, (k1 + d) % s.m), (t2, (k2 + d) % s.m))))
                    + (eps,) for d in range(s.m)))
            assert len(keys) == len(pairs), "a pair-orbit was named twice"


def test_orbit_points_respect_the_frobenius_log_law():
    """Regression for independently signed points breaking the solve."""
    for m in (11, 13):
        F, E, order = SC.small_curve(m)
        r = order
        while r % 2 == 0:
            r //= 2
        rng = random.Random(SEED)
        s = NormalSupport(F, E, find_normal_elements(F, 1, rng)[0])
        if s.size == 0:
            continue
        pts = s.orbit_points(E)
        for t in range(s.orbit_count):
            base = pts[t * s.m]
            cur = base
            for k in range(s.m):
                assert pts[t * s.m + k] == cur
                cur = E.frobenius(cur)
            assert cur == base


def test_configurations_excludes_cancelling_multisets():
    """`8 C(B,4)`: four *distinct* abscissae, never `R + (-R) + T + (-T)`."""
    assert configurations(3) == 0
    assert configurations(4) == 8
    assert configurations(5) == 40


def test_exhaustive_counter_matches_brute_force_quadruples():
    F, E, order = SC.small_curve(11)
    e4 = SC.four_torsion(E, order)
    rng = random.Random(SEED)
    s = NormalSupport(F, E, find_normal_elements(F, 1, rng)[0])
    xs = s.abscissae[:22]
    fast = exhaustive_four_point(E, xs, e4)

    reps = [min(E.points_over(x), key=lambda p: p[1])
            for x in xs if E.points_over(x)]
    e4set = set(map(str, e4))
    usable = exact = 0
    for comb in combinations(range(len(reps)), 4):
        for signs in product((1, -1), repeat=3):
            pts = [reps[comb[0]]] + [
                reps[i] if sg > 0 else E.neg(reps[i])
                for i, sg in zip(comb[1:], signs)]
            S = E.sum_points(pts)
            if str(S) in e4set:
                usable += 1
                if S is None:
                    exact += 1
    assert fast["usable_relations"] == usable
    assert fast["exact_relations"] == exact


def test_support_lies_in_the_index_two_subgroup():
    """Why the coset factor is 2 and not 4."""
    for m in (11, 13):
        F, E, order = SC.small_curve(m)
        r = order
        while r % 2 == 0:
            r //= 2
        rng = random.Random(SEED)
        s = NormalSupport(F, E, find_normal_elements(F, 1, rng)[0])
        for P in s.orbit_points(E):
            assert E.mul(P, 2 * r) is None, "support point outside H"


def test_no_relation_reaches_an_order_four_point():
    """The two order-4 elements of `E[4]` are unreachable, measured."""
    F, E, order = SC.small_curve(13)
    e4 = SC.four_torsion(E, order)
    order4 = [T for T in e4 if T is not None and E.mul(T, 2) is not None]
    assert len(order4) == 2
    rng = random.Random(SEED)
    s = NormalSupport(F, E, find_normal_elements(F, 1, rng)[0])
    assert exhaustive_four_point(E, s.abscissae, order4)["usable_relations"] == 0


def test_supply_model_matches_ground_truth():
    """The model that once missed by 2x-43x, against exhaustive counts."""
    ratios = []
    for m in (11, 13):
        F, E, order = SC.small_curve(m)
        e4 = SC.four_torsion(E, order)
        rng = random.Random(SEED)
        seen = set()
        for a in find_normal_elements(F, 20, rng):
            s = NormalSupport(F, E, a)
            if s.size in seen or s.size < 4 * m:
                continue
            seen.add(s.size)
            pred = supply_prediction(s.size, order, len(e4))
            meas = exhaustive_four_point(E, s.abscissae, e4)
            if pred["expected_usable"] >= 100:
                ratios.append(meas["usable_relations"] / pred["expected_usable"])
    assert len(ratios) >= 3
    mean = sum(ratios) / len(ratios)
    assert 0.8 < mean < 1.25, f"supply model off by {mean:.2f}x"
    assert all(0.5 < r < 2.0 for r in ratios), ratios


def test_frobenius_scalar_is_found_without_a_logarithm():
    F, E, P, Q, r = R.challenge_curve()
    s = frobenius_scalar(E, P, r)
    assert E.mul(P, s) == E.frobenius(P)
    assert E.mul(Q, s) == E.frobenius(Q)
    assert pow(s, 131, r) == 1
    assert (s * s + s + 2) % r == 0


def test_planted_logarithm_is_recovered_end_to_end():
    for seed in (3, 5):
        out = recover_planted(13, seed=seed, log=lambda *a: None)
        assert out["recovered"], out
        assert out["recovered_d"] == out["planted_d"]
        assert out["unknowns"] == out["orbits"] + 1


if __name__ == "__main__":
    import traceback
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"ok   {name}")
            except Exception:
                fails += 1
                print(f"FAIL {name}")
                traceback.print_exc()
    raise SystemExit(1 if fails else 0)
