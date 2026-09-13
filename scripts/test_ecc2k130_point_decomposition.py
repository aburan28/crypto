#!/usr/bin/env python3
"""Tests for the ECC2K-130 point-decomposition boundary.

The cheap parts are checked against independent computations here; the
expensive ones (the toy yield sweep, the class histogram at `n = 131`) are the
script's own run and are asserted inside it.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import random
import unittest

SCRIPT = Path(__file__).with_name("ecc2k130_point_decomposition.py")
SPEC = importlib.util.spec_from_file_location("ecc2k130_decomp", SCRIPT)
assert SPEC and SPEC.loader
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)

ORDER = 2722258935367507707729280517973639940516
R = ORDER // 4


class FieldTests(unittest.TestCase):
    """The challenge field, against textbook identities rather than a table."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.F = mod.GF2m(mod.N131, mod.IRR131)

    def test_the_challenge_reduction_polynomial_is_irreducible(self) -> None:
        self.assertTrue(mod.poly_is_irreducible(mod.IRR131, 131))
        # and a near miss is not, so the test is not vacuous
        self.assertFalse(mod.poly_is_irreducible((1 << 131) | (1 << 13) | 1, 131))

    def test_multiplication_squaring_and_inversion_agree(self) -> None:
        rng = random.Random(7)
        for _ in range(40):
            a = rng.randrange(1, 1 << 131)
            b = rng.randrange(1, 1 << 131)
            self.assertEqual(self.F.sqr(a), self.F.mul(a, a))
            self.assertEqual(self.F.mul(a, self.F.inv(a)), 1)
            self.assertEqual(self.F.mul(a, b), self.F.mul(b, a))
            self.assertEqual(self.F.sqrt(self.F.sqr(a)), a)

    def test_trace_is_linear_onto_F2_and_gates_artin_schreier(self) -> None:
        rng = random.Random(8)
        hits = 0
        for _ in range(60):
            a = rng.randrange(1, 1 << 131)
            b = rng.randrange(1, 1 << 131)
            self.assertEqual(self.F.trace(a ^ b), self.F.trace(a) ^ self.F.trace(b))
            z = self.F.solve_artin_schreier(a)
            if self.F.trace(a):
                self.assertIsNone(z)
            else:
                hits += 1
                self.assertEqual(self.F.sqr(z) ^ z, a)
        self.assertGreater(hits, 0, "both branches must be exercised")


class CurveTests(unittest.TestCase):
    """Group law and order on the challenge curve."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.K = mod.Koblitz(mod.GF2m(mod.N131, mod.IRR131))

    def test_order_matches_the_koblitz_recurrence(self) -> None:
        self.assertEqual(mod.curve_order(131), ORDER)
        self.assertTrue(mod.is_prime(R))
        self.assertEqual(R.bit_length(), 130)
        self.assertAlmostEqual(math.log2(R), 129.0, places=4)

    def test_points_are_on_the_curve_and_killed_by_the_order(self) -> None:
        pts = [P for x in range(1, 12) for P in self.K.points_over(x)]
        self.assertTrue(pts)
        for P in pts:
            self.assertTrue(self.K.on_curve(P))
            self.assertIsNone(self.K.mul(P, ORDER))
            self.assertIsNone(self.K.add(P, self.K.neg(P)))

    def test_batch_scalar_multiplication_matches_the_plain_one(self) -> None:
        pts = [P for x in range(1, 9) for P in self.K.points_over(x)][:8]
        self.assertEqual(self.K.batch_mul(pts, R), [self.K.mul(P, R) for P in pts])

    def test_class_parity_is_the_trace_of_the_abscissa(self) -> None:
        """`P` is in `2E` exactly when `Tr(x(P)) = Tr(a) = 0`, and `2E` is the
        even half of `E / <G> = Z/4`."""
        pts = [P for x in range(1, 24) for P in self.K.points_over(x)]
        cls = self.K.batch_mul(pts, R)
        two = next(c for c in cls if c is not None and self.K.add(c, c) is None)
        for P, c in zip(pts, cls):
            even = c is None or c == two
            self.assertEqual(even, self.K.F.trace(P[0]) == 0, hex(P[0]))

    def test_semaev_s4_vanishes_exactly_on_four_summing_abscissae(self) -> None:
        pts = [P for x in range(1, 40) for P in self.K.points_over(x)]
        P1, P2, P3 = pts[0], pts[2], pts[4]
        R4 = self.K.add(self.K.add(P1, P2), P3)
        self.assertEqual(self.K.s4(P1[0], P2[0], P3[0], R4[0]), 0)
        # a perturbed fourth abscissa must not vanish
        other = next(P for P in pts if P[0] not in (P1[0], P2[0], P3[0], R4[0]))
        self.assertNotEqual(self.K.s4(P1[0], P2[0], P3[0], other[0]), 0)


class ToyTests(unittest.TestCase):
    """The yield law, on a field small enough to enumerate completely."""

    def test_yield_law_tracks_the_measurement_at_n_11(self) -> None:
        cell = mod.measure_toy_yield(11, 3, 5, 96, random.Random(3))
        self.assertEqual(cell["curve_order"], mod.curve_order(11))
        self.assertAlmostEqual(cell["ratio_measured_over_predicted"], 1.0, delta=0.35)
        self.assertIsNotNone(cell["witness_abscissae"])

    def test_a_trace_zero_base_spreads_its_sums_over_half_the_curve(self) -> None:
        """§3.1: `ker Tr` lies in the index-2 subgroup `2E`, so pair sums cannot
        leave it and a target in `<G>` is hit twice as often as `C(|F|,2)/#E`."""
        for n, l in ((13, 7), (17, 8)):
            b = mod.trace_zero_bonus(n, l)
            self.assertEqual(b["targets"], b["subgroup_order"])   # exhaustive
            self.assertAlmostEqual(b["generic"]["implied_spread_over_curve_order"],
                                   1.0, delta=0.05, msg=f"n={n}")
            self.assertAlmostEqual(b["trace_zero"]["implied_spread_over_curve_order"],
                                   0.5, delta=0.03, msg=f"n={n}")
            self.assertAlmostEqual(b["measured_yield_factor"], 2.0, delta=0.1,
                                   msg=f"n={n}")

    def test_a_whole_base_detector_localises_its_own_witness(self) -> None:
        """§3.2: swapping a candidate summand for a class-matched base point
        turns a yes/no detector into a witness finder, with no sub-base query."""
        rng = random.Random(19)
        for n, l, k in ((17, 6, 1), (17, 6, 2), (19, 7, 2)):
            d = mod.swap_localisation(n, l, 2, k, 24, rng)
            self.assertEqual(d["sub_base_queries"], 0)
            self.assertEqual(d["queries_per_target"], k * d["factor_base_size"])
            # both failure modes are collisions, so both stay within an O(1)
            # factor of the collision scale they are predicted from
            self.assertLess(d["miss_rate_over_k_times_collision_scale"], 2.0, d)
            self.assertLess(d["false_positive_rate_over_collision_scale_to_k"],
                            2.0, d)

    def test_a_second_swap_query_squares_the_false_positive_rate(self) -> None:
        rng = random.Random(23)
        one = mod.swap_localisation(19, 7, 2, 1, 24, rng)
        two = mod.swap_localisation(19, 7, 2, 2, 24, rng)
        self.assertLess(two["false_positive_rate_per_non_summand"],
                        one["false_positive_rate_per_non_summand"] / 3)
        self.assertGreater(two["miss_rate_per_summand"],
                           one["miss_rate_per_summand"])       # and costs misses

    def test_the_trace_zero_span_really_lies_in_the_kernel_of_the_trace(self) -> None:
        for n, l in ((13, 7), (17, 9), (19, 10)):
            field = mod.GF2m(n, mod.find_irreducible(n))
            span = mod.trace_zero_subspace(field, l)
            self.assertEqual(len(span), 1 << l)
            self.assertEqual(len(set(span)), 1 << l)      # an honest subspace
            self.assertTrue(all(field.trace(x) == 0 for x in span))


class CostTests(unittest.TestCase):
    """The derived cost surface: the identities it rests on, checked directly."""

    def test_the_product_law_is_flat_in_the_subspace_dimension(self) -> None:
        totals = {mod.cost_cell(3, l, 131, 1)["log2_collection"]
                  for l in (12, 20, 28, 36, 44)}
        self.assertEqual(len(totals), 1, totals)
        self.assertAlmostEqual(totals.pop(), 131 + math.log2(3), places=2)

    def test_collection_is_m_times_the_field_size_for_every_m(self) -> None:
        for m in range(2, 8):
            l = mod.saturating_l(m, 131) - 4
            cell = mod.cost_cell(m, l, 131, 1)
            self.assertAlmostEqual(cell["log2_collection"], 131 + math.log2(m), places=2)

    def test_saturating_dimension_makes_a_target_decompose_once(self) -> None:
        for m in range(2, 9):
            cell = mod.cost_cell(m, mod.saturating_l(m, 131), 131, 1)
            self.assertAlmostEqual(cell["log2_decompositions_per_target"], 0.0, places=2)

    def test_required_oracle_speedup_is_the_product_law_read_backwards(self) -> None:
        rho = math.log2(math.sqrt(math.pi * R / (2 * 2 * 131)))
        for m in range(2, 9):
            got = mod.oracle_budget(m, 131, rho)["log2_speedup_over_search_required"]
            self.assertAlmostEqual(got, rho - 131 - math.log2(m), places=1)
            self.assertFalse(mod.oracle_budget(m, 131, rho)["reachable"])

    def test_a_generic_algorithm_with_the_same_memory_is_cheaper(self) -> None:
        log2r = math.log2(R)
        for cap in (30, 40, 50, 60, 70):
            best = min(
                (c for m in range(2, 9) for split in range(1, m + 1)
                 if (c := mod.memory_capped_row(m, 131, split, cap)) is not None),
                key=lambda c: c["log2_total"])
            bsgs = mod.generic_with_memory(cap, log2r)["log2_bsgs"]
            self.assertLess(bsgs, best["log2_total"],
                            f"at 2^{cap} entries BSGS must still win")

    def test_no_invariant_factor_base_of_usable_size_exists_at_131(self) -> None:
        self.assertEqual(mod.available_subspace_dimensions(131), [0, 1, 130, 131])
        self.assertEqual(mod.mult_order(2, 131), 130)
        # the construction is not empty in general -- a control
        self.assertIn(3, mod.available_subspace_dimensions(7))


if __name__ == "__main__":
    unittest.main()
