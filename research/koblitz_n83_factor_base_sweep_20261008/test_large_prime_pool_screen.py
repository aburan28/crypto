"""Independent small-count controls for the retained large-prime pool screen."""

import itertools
import math
import unittest

import large_prime_pool_screen as screen


class LargePrimePoolScreenTests(unittest.TestCase):
    def test_exact_residual_counts_and_full_partition(self):
        for base in range(1, 4):
            for residual in range(1, 4):
                for summands in range(2, 5):
                    total = 0
                    for count in range(summands + 1):
                        direct = sum(
                            1
                            for _ in itertools.product(
                                itertools.combinations_with_replacement(
                                    range(base), summands - count
                                ),
                                itertools.combinations_with_replacement(
                                    range(residual), count
                                ),
                            )
                        )
                        self.assertEqual(
                            screen.partial_support(base, residual, summands, count),
                            direct,
                        )
                        total += direct
                    self.assertEqual(
                        total, math.comb(base + residual + summands - 1, summands)
                    )

    def test_prefix_requires_both_orbit_and_candidate_index_prefixes(self):
        def header(columns):
            return {
                "curve_a": 0,
                "policy": "public_x_hash",
                "seed": 7,
                "subgroup_order": "2417851639230796216685689",
                "closure": "signed_frobenius",
                "orbit_columns": columns,
                "point_count": 166 * columns,
                "representatives": [[str(i), str(i + 1)] for i in range(columns)],
                "accepted_candidate_indices": list(range(1, columns + 1)),
            }

        smaller, larger = header(2), header(3)
        screen.assert_prefix(smaller, larger)
        changed = header(3)
        changed["representatives"][0] = ["wrong", "point"]
        with self.assertRaises(AssertionError):
            screen.assert_prefix(smaller, changed)
        changed = header(3)
        changed["accepted_candidate_indices"][1] = 99
        with self.assertRaises(AssertionError):
            screen.assert_prefix(smaller, changed)

    def test_invalid_counts_fail_closed(self):
        with self.assertRaises(ValueError):
            screen.partial_support(10, 20, 3, 4)
        with self.assertRaises(ValueError):
            screen.multiset_count(-1, 2)


if __name__ == "__main__":
    unittest.main()
