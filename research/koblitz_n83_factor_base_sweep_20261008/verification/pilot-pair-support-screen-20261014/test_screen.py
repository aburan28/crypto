"""Independent small-group and archive-integrity checks for the support bound."""

from itertools import combinations_with_replacement
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import screen


class SupportScreenTests(unittest.TestCase):
    def test_multiset_bound_on_an_exhaustive_small_group(self):
        order = 11
        small = (1, 2)
        residual = (3, 4)
        tuples = [(a + b) % order for a in small for b in residual]
        self.assertEqual(tuples, [4, 5, 5, 6])
        self.assertEqual(screen.multiset_count(2, 2, 2, 1), 4)
        self.assertEqual(len(set(tuples) - {0}), 3)
        self.assertLessEqual(len(set(tuples) - {0}),
                             screen.multiset_count(2, 2, 2, 1))
        three_small = list(combinations_with_replacement(small, 3))
        two_residual = list(combinations_with_replacement(residual, 2))
        self.assertEqual(len(three_small) * len(two_residual),
                         screen.multiset_count(2, 2, 5, 2))
        sums = {(sum(a) + sum(b)) % order
                for a in three_small for b in two_residual}
        self.assertLessEqual(len(sums - {0}),
                             screen.multiset_count(2, 2, 5, 2))

    def test_all_retained_five_summand_shapes_have_nonvacuous_bounds(self):
        result = screen.calculate()
        self.assertEqual((result["ordered_pair_count"],
                          result["worker_admissible_pair_count"]), (648, 558))
        self.assertEqual([row["nonvacuous_hit_ceiling_count"]
                          for row in result["shape_summaries"]],
                         [558, 558, 209, 198])
        self.assertFalse(result["new_n83_solver_search_executed"])
        self.assertIsNone(result["selected_best_total_runtime"])

    def test_same_base_pools_have_exactly_distinct_screened_support(self):
        rows = screen.calculate()["ordered_pairs"]
        arms = {row["pool_policy"]: row for row in rows
                if row["base_policy"] == "public_x_sequential"
                and row["base_seed"] == 2026100802
                and row["base_orbits"] == 600
                and row["pool_seed"] == 2026100801
                and row["pool_orbits"] == 600
                and row["pool_policy"] in
                    ("public_x_gray_prefix", "public_x_hash")}
        self.assertEqual(len(arms), 2)
        gray, hashed = arms["public_x_gray_prefix"], arms["public_x_hash"]
        self.assertEqual((gray["residual_signed_orbits"],
                          hashed["residual_signed_orbits"]), (97, 600))
        gray_bound = gray["shape_bounds"][1]
        hash_bound = hashed["shape_bounds"][1]
        self.assertEqual((gray_bound["summands"],
                          gray_bound["residual_summands"]), (5, 2))
        self.assertEqual((gray_bound["minimum_queries_before_half_hit_ceiling_is_possible"],
                          hash_bound["minimum_queries_before_half_hit_ceiling_is_possible"]),
                         (57, 2))
        self.assertGreater(int(hash_bound["point_multisets_including_repetition"]),
                           30 * int(gray_bound["point_multisets_including_repetition"]))

    def test_altered_reciprocal_overlap_is_rejected(self):
        receipt = json.loads(screen.OVERLAP.read_text())
        pair = next(row for row in receipt["ordered_pairs"]
                    if row["base_orbits"] == row["pool_orbits"] == 600
                    and 0 < row["shared_signed_orbits"] < 598)
        reverse = next(row for row in receipt["ordered_pairs"]
                       if row["base_object_blake3"] == pair["pool_object_blake3"]
                       and row["pool_object_blake3"] == pair["base_object_blake3"])
        reverse["shared_signed_orbits"] += 1
        reverse["shared_x_values"] += 83
        reverse["residual_signed_orbits"] -= 1
        with TemporaryDirectory() as directory:
            path = Path(directory) / "altered.json"
            path.write_text(json.dumps(receipt))
            original = screen.OVERLAP
            try:
                screen.OVERLAP = path
                with self.assertRaisesRegex(ValueError, "reciprocal overlap differs"):
                    screen.calculate()
            finally:
                screen.OVERLAP = original


if __name__ == "__main__":
    unittest.main()
