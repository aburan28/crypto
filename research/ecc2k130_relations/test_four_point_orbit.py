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
from run_four_point_orbit import is_identity_relation         # noqa: E402
from target_boundary import cost as logarithm_cost            # noqa: E402
from validate_amortised_attack import (aggregate as amortised_aggregate,  # noqa: E402
                                       attack as amortised_attack)
from validate_target_model import (_signed_subset_sums,       # noqa: E402
                                   _support)

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


def test_characteristic_equation_patterns_are_identities():
    """The only collisions the ECC2K-130 run finds, identified rather than binned.

    `sigma^2 + sigma + 2 = 0` makes `2A + sigma A + sigma^2 A` vanish on every
    point of the curve, and `sigma^3 = 2 - sigma` makes `2A - sigma A -
    sigma^3 A` do the same.  Both are construction equations: real relations
    carrying no information at all.
    """
    _, E, _, _, _ = R.challenge_curve()
    assert is_identity_relation(E, 1, 2, +1)
    assert is_identity_relation(E, 1, 3, -1)
    # neighbouring patterns must not be mistaken for identities
    assert not is_identity_relation(E, 2, 3, +1)
    assert not is_identity_relation(E, 1, 4, -1)
    assert not is_identity_relation(E, 1, 2, -1)


def test_planted_logarithm_is_recovered_end_to_end():
    for seed in (3, 5):
        out = recover_planted(13, seed=seed, log=lambda *a: None)
        assert out["recovered"], out
        assert out["recovered_d"] == out["planted_d"]
        assert out["unknowns"] == out["orbits"] + 1


def test_subset_sum_set_is_sigma_closed():
    """Why sigma buys no hit rate: the reachable sums are already sigma-closed."""
    for mdeg, orbits, n in ((17, 2, 2), (19, 2, 2)):
        F, E, _ = SC.small_curve(mdeg)
        rng = random.Random(5)
        sup = _support(F, E, orbits, rng)
        assert sup is not None
        sums = _signed_subset_sums(E, sup.orbit_points(E), n)
        assert all(E.frobenius(S) in sums for S in sums if S is not None)


def test_homogeneous_relations_are_not_priced_as_a_logarithm():
    """The m-sweep crosses rho only when relations carry no information.

    Guards the section 5 conclusion: whatever the homogeneous table says, the
    cost of a logarithm must stay above the rho reference.
    """
    import json
    from pathlib import Path
    data = json.loads((HERE / "results" / "target_boundary.json").read_text())
    homog = data["homogeneous_extension_past_m8"]
    assert any(r["log2_total_vs_rho"] < 0 for r in homog), (
        "the homogeneous accounting is supposed to cross; section 5 exists to "
        "explain why that is not a result")
    best = data["best_logarithm_cost"]
    assert best["log2_total_vs_rho"] > 0, (
        "a logarithm priced against known targets must not beat rho here")
    assert data["memory_charged"] is False


def test_amortising_the_table_is_a_real_saving_and_correctly_signed():
    """Building once must be cheaper than rebuilding per attempt, never dearer.

    This caught the counterfactual charging `attempts * max(build, stream)`,
    which drops the build entirely whenever streaming dominates and so made
    rebuilding look *cheaper* than amortising at T=10, n=15. A faithful
    rebuild pays both terms on every attempt.
    """
    for T, n, s in ((1, 14, 13), (10, 15, 8), (3, 8, 4), (2, 12, 11)):
        for quot in (True, False):
            once = logarithm_cost(T, n, s, amortise=True, quotient_table=quot)
            each = logarithm_cost(T, n, s, amortise=False, quotient_table=quot)
            assert once["log2_total_cost"] <= each["log2_total_cost"] + 1e-9


def test_probe_cost_and_e4_translates_are_priced():
    """Regression: probes were counted but not costed, and E[4] was ignored."""
    q = logarithm_cost(2, 12, 11, quotient_table=True)
    f = logarithm_cost(2, 12, 11, quotient_table=False)
    # a quotiented table must charge canonicalisation per probe
    assert q["log2_probe_cost"] > 0
    # a full table must not, and must pay 131x the build instead
    assert f["log2_probe_cost"] == 0
    assert f["log2_build_once"] > q["log2_build_once"]
    # both E[4] translates are probed, so probes exceed the streamed count
    from target_boundary import REACHABLE_E4
    assert REACHABLE_E4 == 2
    assert q["log2_probes_per_attempt"] > 1


def test_a_skipped_or_failed_cell_cannot_be_reported_as_success():
    """Regression for an eighth defect, again in validation code.

    Filtering to cells that carry a recovery and then reporting "N of N"
    over the survivors lets a skipped support or a cell that ran out of
    targets vanish from the denominator, so the artifact claims success for
    validation that never ran.
    """
    good = {"verified_by_point_identity": True, "recovered_d": 1, "planted_d": 1}
    skipped = {"skipped": "no support of that size"}
    failed = {"recovered": False, "reason": "ran out of targets"}

    assert amortised_aggregate([good, good], 2)["all_recovered"]
    for bad in (skipped, failed):
        agg = amortised_aggregate([good, bad], 2)
        assert not agg["all_recovered"], bad
        assert agg["cells_recovered"] == 1
        assert agg["cells_expected"] == 2
    # a cell that never ran at all must not shrink the denominator either
    assert not amortised_aggregate([good], 2)["all_recovered"]


def test_amortised_attack_recovers_a_planted_logarithm():
    """The lopsided single-orbit structure the optimum uses, run for real."""
    for spec in ((13, 3, 1), (19, 3, 1)):
        out = amortised_attack(*spec)
        assert out.get("verified_by_point_identity"), out
        assert out["recovered_d"] == out["planted_d"]
        assert out["unknowns"] == out["orbits"] + 1
        assert out["relations"] == out["unknowns"]


def test_logarithm_cost_is_monotone_in_the_obvious_places():
    a = logarithm_cost(10, 15, 8)
    assert a is not None
    # more unknowns cannot be cheaper at fixed relation shape
    assert logarithm_cost(200, 15, 8)["log2_total_cost"] > a["log2_total_cost"]
    # an impossible split returns nothing rather than a number
    assert logarithm_cost(1, 500, 250) is None


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
