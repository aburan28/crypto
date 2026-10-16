#!/usr/bin/env python3
"""Small exact controls for the mixed-prefix support arithmetic."""

from itertools import combinations_with_replacement, product
import unittest

import screen


class PrefixScreenTests(unittest.TestCase):
    def test_two_disjoint_domains_match_explicit_multisets(self):
        small = ("s0", "s1")
        residual = ("r0", "r1", "r2")
        explicit = list(product(combinations_with_replacement(small, 3),
                                combinations_with_replacement(residual, 2)))
        self.assertEqual(len(explicit), 24)
        rows = screen.count_support(len(small), len(residual))
        row = next(row for row in rows if (row["small_summands"], row["residual_summands"]) == (3, 2))
        self.assertEqual(int(row["unordered_point_multisets_with_repetition"]), len(explicit))

    def test_k600_two_residuals_crosses_vacuous_ceiling(self):
        rows = screen.count_support(99600, 240368)
        one = next(row for row in rows if (row["small_summands"], row["residual_summands"]) == (4, 1))
        two = next(row for row in rows if (row["small_summands"], row["residual_summands"]) == (3, 2))
        self.assertEqual(int(one["unordered_point_multisets_with_repetition"]),
                         985664065933326097483200)
        self.assertEqual(int(two["unordered_point_multisets_with_repetition"]),
                         4757348390396427029827200)
        self.assertLess(int(one["unordered_point_multisets_with_repetition"]), screen.ORDER - 1)
        self.assertGreater(int(two["unordered_point_multisets_with_repetition"]), screen.ORDER - 1)
        self.assertEqual(two["uniform_nonidentity_hit_probability_ceiling"]["numerator"],
                         str(screen.ORDER - 1))

    def test_manifest_selection_requires_unique_exact_label(self):
        row = {"a": 0, "policy": "public_x_hash", "seed": 2026100801, "columns": 600}
        self.assertEqual(screen.sole_base({"bases": [row]}, 600), row)
        with self.assertRaisesRegex(ValueError, "unique"):
            screen.sole_base({"bases": [row, row]}, 600)
        with self.assertRaisesRegex(ValueError, "unique"):
            screen.sole_base({"bases": [row]}, 2048)


if __name__ == "__main__":
    unittest.main()
