"""Independent integer and ordinal checks for the v2 size frontier."""

import importlib.util
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("size_frontier_v2", HERE / "size_frontier_v2.py")
frontier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(frontier)


class SizeFrontierV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.design = frontier.build()
        cls.frozen = json.loads((HERE / "design.json").read_text())

    def test_exact_primary_thresholds_are_minimal(self):
        order = 2417851639230796216685689
        expected = {2: 13247128046, 3: 1469216, 4: 16627, 5: 1182, 6: 209}
        for summands, columns in expected.items():
            with self.subTest(summands=summands):
                below = math.comb(166 * (columns - 1) + summands - 1, summands)
                at = math.comb(166 * columns + summands - 1, summands)
                self.assertLess(below, order - 1)
                self.assertGreaterEqual(at, order - 1)
                self.assertEqual(frontier.threshold(order, summands), columns)

    def test_v1_prefix_and_each_new_size_are_addressable(self):
        design = self.design
        self.assertEqual(design["base_specs"][:1050], self.frozen["base_specs"])
        self.assertEqual(design["base_spec_count"], 1320)
        self.assertEqual(design["new_base_spec_count"], 270)
        self.assertEqual(design["total_tuple_count"], 57024000)
        self.assertEqual(design["v1_tuple_count_preserved"], 45360000)
        for ordinal in (0, 43199, 45359999, 45360000, 57023999):
            with self.subTest(ordinal=ordinal):
                row = frontier.case(design, ordinal)
                self.assertEqual(row["ordinal"], ordinal)
                self.assertEqual(row["base"], design["base_specs"][ordinal // 43200])
                self.assertEqual(row["solver_index"], ordinal % 43200)
                self.assertFalse(row["executed"])
                self.assertEqual(row["disposition"],
                                 frontier.disposition(row["base"], row["solver"]))
        new_sizes = {base["parameter"] for base in design["base_specs"][1050:]}
        self.assertEqual(new_sizes, {1182, 2048, 4096, 8192, 16627})
        with self.assertRaises(ValueError):
            frontier.case(design, 57024000)
        self.assertEqual(sum(design["disposition_counts"].values()), 57024000)
        self.assertEqual(frontier.case(design, 45360000)["disposition"],
                         "incompatible_compact_four_sum_domain")
        selected = design["base_specs"][1052]
        self.assertEqual(selected["closure"], "signed_frobenius")
        solver = {"summands": 5, "split_bits": 0, "symmetry": "none",
                  "backend": "f4", "enumeration": "binary", "max_large_primes": 2,
                  "threads": 1, "linear_algebra": "modular_biguint"}
        self.assertEqual(frontier.disposition(selected, solver),
                         "wide_large_prime_partial_producer_required")

    def test_aggregate_dispositions_reproduce_frozen_rust_counts(self):
        counts = frontier.count_dispositions(
            self.frozen["base_specs"], self.frozen["solver_axes"])
        counts["large_prime_adapter_is_single_word"] = counts.pop(
            "wide_large_prime_partial_producer_required")
        self.assertEqual(counts, self.frozen["disposition_counts"])

    def test_every_solver_ordinal_decodes_bijectively(self):
        radices = self.design["solver_radices"]
        for ordinal in range(self.design["solver_tuple_count"]):
            digits = frontier.solver_digits(ordinal, radices)
            encoded = 0
            for digit, radix in zip(digits, radices):
                self.assertGreaterEqual(digit, 0)
                self.assertLess(digit, radix)
                encoded = encoded * radix + digit
            self.assertEqual(encoded, ordinal)

    def test_support_ceiling_and_rank_query_screen_are_exact(self):
        order = 2417851639230796216685689
        rows = {(row["curve_a"], row["orbit_columns"], row["summands"]): row
                for row in self.design["size_cases"]}
        self.assertEqual(len(rows), 70)
        for columns in (600, 900, 1182, 2048, 4096, 8192, 16627):
            for summands in (4, 5, 6):
                with self.subTest(columns=columns, summands=summands):
                    row = rows[0, columns, summands]
                    points = 166 * columns
                    multisets = math.comb(points + summands - 1, summands)
                    ceiling = Fraction(min(multisets, order - 1), order - 1)
                    self.assertEqual(row["conditional_signed_point_count"], points)
                    self.assertEqual(row["unordered_multisets_with_repetition"], str(multisets))
                    self.assertEqual(
                        Fraction(int(row["uniform_nonidentity_hit_ceiling"]["numerator"]),
                                 int(row["uniform_nonidentity_hit_ceiling"]["denominator"])),
                        ceiling,
                    )
                    self.assertEqual(row["full_pair_multisets_if_materialized"],
                                     math.comb(points + 1, 2))
                    necessary = max(columns, math.ceil(Fraction(columns, 2 * ceiling)))
                    self.assertEqual(row["minimum_queries_for_half_rank_ceiling"], necessary)
        self.assertLess(Fraction(
            int(rows[0, 900, 5]["uniform_nonidentity_hit_ceiling"]["numerator"]),
            order - 1), 1)
        self.assertEqual(rows[0, 1182, 5]["uniform_nonidentity_hit_ceiling"]["numerator"],
                         str(order - 1))
        self.assertEqual(rows[0, 16627, 4]["uniform_nonidentity_hit_ceiling"]["numerator"],
                         str(order - 1))


if __name__ == "__main__":
    unittest.main()
