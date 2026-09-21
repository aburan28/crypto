#!/usr/bin/env python3
"""Tests for the six experiment runs.

These are runs, not derivations, so the tests check the machinery the runs rest
on: that the group is posed correctly, that the Frobenius eigenvalue is the one
that actually acts, that a collected relation really holds in the group, and
that the linear algebra recovers a logarithm it was given.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import random
import unittest

SCRIPT = Path(__file__).with_name("ecc2k130_decomposition_experiments.py")
SPEC = importlib.util.spec_from_file_location("ecc2k130_experiments", SCRIPT)
assert SPEC and SPEC.loader
X = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(X)
DEC = X.DEC


class RungTests(unittest.TestCase):
    def test_the_exponent_divides_the_order_and_carries_the_prime(self) -> None:
        for n in (11, 13, 17, 19):
            r = X.Rung(n)
            self.assertEqual(r.order % r.exponent, 0, f"n={n}")
            self.assertEqual(r.exponent % r.p, 0, f"n={n}")
            self.assertTrue(DEC.is_prime(r.p))

    def test_rank_two_p_torsion_is_detected_and_excluded(self) -> None:
        """`n = 11` has `#E = 2116 = 23^2` with `E = Z/23 x Z/92`, so the
        23-torsion is two-dimensional and no single cyclic subgroup exists to
        pose the logarithm in.  Every other ladder rung here is usable."""
        eleven = X.Rung(11)
        self.assertEqual(eleven.p_torsion_rank, 2)
        self.assertFalse(eleven.usable_end_to_end)
        self.assertIsNone(eleven.G)
        for n in (13, 17, 19):
            r = X.Rung(n)
            self.assertEqual(r.p_torsion_rank, 1, f"n={n}")
            self.assertTrue(r.usable_end_to_end)

    def test_the_projection_lands_every_point_in_the_order_p_subgroup(self) -> None:
        r = X.Rung(13)
        c = r.curve
        for x in range(1, 40):
            for P in c.points_over(x):
                Q = r.project(P)
                self.assertIsNone(c.mul(Q, r.p), hex(x))

    def test_the_generator_really_has_order_p(self) -> None:
        for n in (13, 19):
            r = X.Rung(n)
            self.assertIsNone(r.curve.mul(r.G, r.p))
            self.assertIsNotNone(r.G)


class FrobeniusTests(unittest.TestCase):
    def test_the_eigenvalue_is_a_root_of_its_characteristic_polynomial(self) -> None:
        """`pi^2 + pi + 2 = 0` for `K_0`, and the root that acts is settled by
        applying Frobenius to a point, not by picking one."""
        for n in (13, 19):
            r = X.Rung(n)
            lam = X.frobenius_eigenvalue(r)
            self.assertEqual((lam * lam + lam + 2) % r.p, 0, f"n={n}")
            self.assertEqual(r.curve.mul(r.G, lam), X._frob_point(r, r.G), f"n={n}")

    def test_frobenius_fixes_the_curve_and_has_order_n_on_a_generic_point(self) -> None:
        r = X.Rung(19)
        c = r.curve
        P = next(Q for x in range(2, 40) for Q in c.points_over(x))
        self.assertTrue(c.on_curve(X._frob_point(r, P)))
        Q = P
        for _ in range(19):
            Q = X._frob_point(r, Q)
        self.assertEqual(Q, P, "Frobenius must have order dividing n")


class CollectionTests(unittest.TestCase):
    def test_the_base_is_exactly_the_points_over_the_subspace(self) -> None:
        r = X.Rung(13)
        points, index, abscissae = X.build_base(r, 6)
        self.assertTrue(all(r.curve.on_curve(P) for P in points))
        self.assertTrue(all(P[0] < (1 << 6) for P in points))
        self.assertEqual(sorted(set(P[0] for P in points)), sorted(abscissae))
        self.assertEqual(len(points), 2 * len(abscissae) - (1 if 0 in abscissae else 0))

    def test_canonicalisation_makes_the_measured_yield_match_the_law(self) -> None:
        """Without it every decomposition is found once per `(m-1)`-subset of
        itself and the yield reads three times its true value."""
        r = X.Rung(13)
        points, _, _ = X.build_base(r, 6)
        predicted = math.comb(len(points), 3) / r.order
        d = X.collect(r, 6, 3, "full", random.Random(3), want=60)
        measured = len(d["rows"]) / d["targets"]
        self.assertAlmostEqual(measured / predicted, 1.0, delta=0.35)

    def test_every_collected_relation_holds_in_the_group(self) -> None:
        """`collect` asserts this inline; this checks the assertion can fire."""
        r = X.Rung(13)
        c = r.curve
        P = next(Q for x in range(2, 20) for Q in c.points_over(x))
        self.assertNotEqual(c.add(P, P), P, "a broken group law must be visible")


class LinearAlgebraTests(unittest.TestCase):
    def test_it_solves_a_system_it_was_handed(self) -> None:
        p = 2003
        secret = [11, 202, 1999, 7]
        rng = random.Random(1)
        rows, rhs = [], []
        for _ in range(24):
            cols = {k: rng.randrange(1, p) for k in range(4)}
            rows.append(cols)
            rhs.append(sum(v * secret[k] for k, v in cols.items()) % p)
        sol, rank = X.solve_mod_p(rows, rhs, 4, p)
        self.assertEqual(rank, 4)
        self.assertEqual(sol, secret)

    def test_it_refuses_an_underdetermined_system(self) -> None:
        p = 2003
        rows = [{0: 1, 1: 1}, {0: 2, 1: 2}]
        sol, rank = X.solve_mod_p(rows, [5, 10], 2, p)
        self.assertIsNone(sol)
        self.assertEqual(rank, 1)
        self.assertEqual(X.matrix_rank_mod_p(rows, 2, p), 1)


class EndToEndTests(unittest.TestCase):
    def test_the_planted_logarithm_comes_back(self) -> None:
        row = X.run_e1_rung(13, rule="full")
        self.assertTrue(row["planted_log_recovered"])
        self.assertEqual(row["matrix_rank"], row["unknowns"])
        self.assertGreater(row["lambda"], 0.5 * row["m"])

    def test_the_folded_rule_halves_the_unknowns(self) -> None:
        full = X.run_e1_rung(13, rule="full")
        fold = X.run_e1_rung(13, rule="folded")
        self.assertEqual(fold["unknowns"] * 2, full["unknowns"] + 1)
        self.assertTrue(fold["planted_log_recovered"])


class StoppingRuleTests(unittest.TestCase):
    """The first round stopped only at full rank over every column and read the
    cost of doing so as a property of the method.  These pin the two halves of
    the correction: what `determined_columns` means, and that the budget the
    product law prices is enough to finish."""

    def test_a_column_is_determined_only_when_no_free_column_reaches_it(self) -> None:
        p = 2003
        # Columns 0 and 1 are pivoted alone; 2 and 3 appear only together, so
        # neither is pinned even though both are hit.
        rows = [{0: 1}, {1: 1}, {2: 1, 3: 1}]
        self.assertEqual(X.determined_columns(rows, 4, p), {0, 1})
        # One more independent row on the pair determines both.
        rows.append({2: 1, 3: 2})
        self.assertEqual(X.determined_columns(rows, 4, p), {0, 1, 2, 3})

    def test_hit_is_weaker_than_determined(self) -> None:
        """Every column above is touched by a relation; half of them are still
        free.  Poisson predicts the first fraction, never the second."""
        p = 2003
        rows = [{0: 1}, {1: 1}, {2: 1, 3: 1}]
        hit = set().union(*(set(r) for r in rows))
        self.assertEqual(hit, {0, 1, 2, 3})
        self.assertLess(len(X.determined_columns(rows, 4, p)), len(hit))

    def test_the_budgeted_relation_count_finishes_a_descent(self) -> None:
        """`|F|` relations -- what the product law prices, and no more."""
        b = X.budget_run(X.Rung(13), 6, 3, random.Random(100))
        self.assertEqual(b["relations_budgeted"], b["factor_base_size"])
        self.assertTrue(b["descent_landed"])
        self.assertGreater(b["fraction_determined"], 0.8)
        self.assertLess(b["lambda"], 3 * b["m"])

    def test_full_rank_costs_the_coupon_collector_not_a_constant(self) -> None:
        """The superseded claim was a constant `2.0x`.  It tracks `ln|F|/m`,
        which is a different number at every factor-base size."""
        small = X.full_rank_crossing(X.Rung(13), 6, 3, random.Random(9))
        large = X.full_rank_crossing(X.Rung(19), 9, 3, random.Random(9))
        self.assertLess(small["coupon_collector_rho"], large["coupon_collector_rho"])
        for c in (small, large):
            self.assertGreater(c["measured_over_coupon"], 0.5, c)
            self.assertLess(c["measured_over_coupon"], 2.0, c)


class RankCeilingTests(unittest.TestCase):
    def test_the_l5_cell_cannot_reach_full_rank_at_any_relation_count(self) -> None:
        """E2's `l = 5` stall is the base, not the sample: exhaustively over
        every decomposition it admits, the row space is one short."""
        c = X.e2_rank_ceiling(19, 5)
        self.assertEqual(c["factor_base_size"], 29)
        self.assertEqual(c["rank_ceiling"], 28)
        self.assertFalse(c["full_rank_reachable"])
        # and the reason: every reachable triple has exactly one even abscissa
        self.assertEqual(set(c["odd_abscissa_census"]), {"2"})

    def test_the_next_dimension_up_admits_the_triples_that_break_it(self) -> None:
        c = X.e2_rank_ceiling(19, 6)
        self.assertTrue(c["full_rank_reachable"])
        self.assertGreater(int(c["odd_abscissa_census"]["0"]), 0)


class DispersionTests(unittest.TestCase):
    def test_the_mean_matches_the_yield_law_and_the_tail_is_poisson(self) -> None:
        d = X.e5_dispersion(13, 2, 7)
        self.assertEqual(d["sampling"], "exhaustive over the odd-order subgroup")
        self.assertAlmostEqual(d["mean_ratio"], 1.0, delta=0.1)
        self.assertAlmostEqual(d["index_of_dispersion"], 1.0, delta=0.2)

    def test_degenerate_subsets_are_what_makes_the_raw_dispersion_large(self) -> None:
        """A point and its negative sum to a forced target; those subsets carry
        no relation and almost all of the raw variance."""
        d = X.e5_dispersion(13, 2, 7)
        self.assertGreater(d["index_of_dispersion_all_subsets"], 2.0)
        self.assertGreater(d["degenerate_subsets"], 0)
        self.assertLess(d["index_of_dispersion"], d["index_of_dispersion_all_subsets"])


class OrbitTests(unittest.TestCase):
    def test_orbit_relations_are_independent_and_the_saving_is_n(self) -> None:
        d = X.e6_orbit_union(13, 6, relations_multiple=20)
        self.assertTrue(d["full_rank"])
        self.assertTrue(d["solved"])
        self.assertAlmostEqual(d["frobenius_saving_over_n"], 1.0, delta=0.2)

    def test_too_few_relations_look_like_a_dependency_and_are_not_one(self) -> None:
        """The first pass at this read rank-deficient at n = 19 on 15 relations;
        with enough of them the rank is full.  Starved is not dependent."""
        starved = X.e6_orbit_union(19, 5, relations_multiple=1, slack=10)
        fed = X.e6_orbit_union(19, 5, relations_multiple=40)
        self.assertLessEqual(starved["matrix_rank"], fed["matrix_rank"])
        self.assertTrue(fed["full_rank"])


if __name__ == "__main__":
    unittest.main()
