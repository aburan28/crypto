#!/usr/bin/env python3
"""Tests for the follow-on experiment boundaries.

These are derivations, so the tests check the identities they rest on and the
guards that keep each model from degenerating into a generic algorithm.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import unittest

SCRIPT = Path(__file__).with_name("ecc2k130_decomposition_targets.py")
SPEC = importlib.util.spec_from_file_location("ecc2k130_targets", SCRIPT)
assert SPEC and SPEC.loader
T = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(T)
DEC = T.DEC

R = DEC.curve_order(131) // 4
RHO = math.log2(math.sqrt(math.pi * R / (2 * 2 * 131)))


class LadderTests(unittest.TestCase):
    def test_every_rung_has_the_ecc2k130_subspace_obstruction(self) -> None:
        rungs = T.ladder()
        self.assertEqual([c["n"] for c in rungs],
                         [11, 13, 19, 29, 37, 53, 59, 61, 67])
        for c in rungs:
            n = c["n"]
            self.assertEqual(c["available_invariant_dimensions"], [0, 1, n - 1, n])
            self.assertEqual(DEC.mult_order(2, n), n - 1)

    def test_lambda_is_flat_at_the_summand_count(self) -> None:
        for c in T.ladder():
            self.assertAlmostEqual(c["predicted_lambda_m3"], 3.0, delta=0.8)

    def test_the_predicted_slope_against_n_is_one(self) -> None:
        rungs = T.ladder()
        xs = [c["n"] for c in rungs]
        ys = [c["predicted_log2_total_m3"] for c in rungs]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        slope = (sum((x - mx) * (y - my) for x, y in zip(xs, ys))
                 / sum((x - mx) ** 2 for x in xs))
        self.assertAlmostEqual(slope, 1.0, places=2)

    def test_end_to_end_rungs_are_the_ones_that_fit_a_budget(self) -> None:
        for c in T.ladder():
            self.assertEqual(c["end_to_end_feasible"],
                             c["predicted_log2_total_m3"] <= 40.0)


class DetectorTests(unittest.TestCase):
    def test_localising_is_never_dearer_than_detector_only(self) -> None:
        for m in range(2, 9):
            a = T.detector_floor(m, 131, True)["log2_floor"]
            b = T.detector_floor(m, 131, False)["log2_floor"]
            self.assertLessEqual(a, b + 1e-9, f"m={m}")

    def test_the_gap_widens_with_the_summand_count(self) -> None:
        gaps = [T.detector_floor(m, 131, False)["log2_floor"]
                - T.detector_floor(m, 131, True)["log2_floor"] for m in range(2, 9)]
        self.assertEqual(gaps, sorted(gaps))
        self.assertLess(gaps[0], 1.0)        # at m = 2 the two coincide
        self.assertGreater(gaps[-1], 30.0)   # by m = 8 they are far apart

    def test_a_free_non_localising_detector_still_loses_to_rho(self) -> None:
        """The headline of E3: deciding for free is not enough, localising is."""
        for m in range(2, 9):
            self.assertGreater(T.detector_floor(m, 131, False)["log2_floor"], RHO)
        beats = [m for m in range(2, 9)
                 if T.detector_floor(m, 131, True)["log2_floor"] < RHO]
        self.assertEqual(beats, [4, 5, 6, 7, 8])

    def test_the_frobenius_collapse_lowers_both_floors(self) -> None:
        for m in range(2, 9):
            for loc in (True, False):
                plain = T.detector_floor(m, 131, loc)["log2_floor"]
                orbit = T.detector_floor(m, 131, loc, frobenius=True)["log2_floor"]
                self.assertLess(orbit, plain, f"m={m} localising={loc}")


class LargePrimeTests(unittest.TestCase):
    def test_the_guard_keeps_the_oracle_filtering(self) -> None:
        for m in (2, 3, 4):
            c = T.large_prime(m, 131)
            y = (DEC.log2_comb(c["l_star"], m - 1) + c["large_prime_dim"]
                 - math.log2(m) - 131)
            self.assertLessEqual(y, 1e-6, "partial yield must stay under one a target")
            self.assertGreater(c["large_prime_dim"], c["l_star"])

    def test_it_beats_the_plain_law_and_still_loses_to_a_generic_algorithm(self) -> None:
        plain = DEC.cost_cell(3, DEC.saturating_l(3, 131), 131, 1)["log2_total"]
        for mu in (30, 40, 50, 60):
            c = T.large_prime(3, 131, mu)
            self.assertLess(c["log2_total"], plain)
            bsgs = DEC.generic_with_memory(mu, math.log2(R))["log2_bsgs"]
            self.assertGreater(c["log2_total"], bsgs,
                               f"at 2^{mu} entries BSGS must still win")
            self.assertGreater(c["log2_total"], RHO)


class OrbitTests(unittest.TestCase):
    def test_the_collapse_divides_collection_by_n(self) -> None:
        for m in range(2, 9):
            l = DEC.saturating_l(m, 131)
            plain = DEC.cost_cell(m, l, 131, 1)["log2_collection"]
            orbit = DEC.cost_cell(m, l, 131, 1, frobenius=True)["log2_collection"]
            # cost_cell rounds to two decimals, so compare within that
            self.assertAlmostEqual(plain - orbit, math.log2(131), delta=0.02)
            self.assertAlmostEqual(orbit, 131 + math.log2(m) - math.log2(131),
                                   delta=0.02)

    def test_the_orbit_base_has_to_be_materialised(self) -> None:
        cell = DEC.cost_cell(3, DEC.saturating_l(3, 131), 131, 1, frobenius=True)
        self.assertIsNotNone(cell["log2_table_entries"])
        self.assertAlmostEqual(cell["log2_table_entries"], cell["l"], places=2)

    def test_the_collapse_never_reaches_rho(self) -> None:
        for m in range(2, 9):
            cell = DEC.cost_cell(m, DEC.saturating_l(m, 131), 131, 1, frobenius=True)
            self.assertGreater(cell["log2_total"] - RHO, 60.0)


class DispersionTests(unittest.TestCase):
    def test_sample_size_is_the_dispersion_standard_error(self) -> None:
        d = T.poisson_test_design(effect=0.20, z=3.0)
        self.assertEqual(d["targets_per_cell"], math.ceil(2 * (3.0 / 0.20) ** 2))
        self.assertEqual(d["null_value"], 1.0)
        tighter = T.poisson_test_design(effect=0.10, z=3.0)["targets_per_cell"]
        self.assertEqual(tighter, 4 * d["targets_per_cell"])


if __name__ == "__main__":
    unittest.main()
