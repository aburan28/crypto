"""Check the exact half-ceiling threshold and the archived panel result."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import rank_floor


class RankFloorTests(unittest.TestCase):
    def test_half_ceiling_floor_is_tight_on_small_group(self):
        # Two candidate multisets over ten nonidentity query targets,
        # with a goal of two decompositions across the queries.
        floor = rank_floor.query_floor(order=11, multiset_count=2,
                                       desired_rows=2)
        self.assertEqual(floor, 5)
        self.assertLess((floor - 1) * 2 * 2, 2 * (11 - 1))
        self.assertGreaterEqual(floor * 2 * 2, 2 * (11 - 1))
        self.assertEqual(rank_floor.query_floor(11, 20, 2), 1)
        with self.assertRaisesRegex(ValueError, "positive order"):
            rank_floor.query_floor(11, 0, 2)

    def test_all_archived_pair_shapes_and_K600_floor(self):
        result = rank_floor.calculate()
        self.assertEqual((result["admissible_pair_count"],
                          result["pair_shape_count"]), (558, 2232))
        self.assertEqual(len(result["shape_and_K_summaries"]), 12)
        five_two = next(row for row in result["shape_and_K_summaries"]
                        if (row["summands"], row["residual_summands"],
                            row["desired_decompositions_K"]) == (5, 2, 600))
        self.assertEqual((five_two["pair_count"],
                          five_two["minimum_query_floor"],
                          five_two["lower_median_query_floor"],
                          five_two["maximum_query_floor"]),
                         (166, 889, 4878, 79681727))
        self.assertEqual({(row["pool_policy"],
                           row["half_K_decomposition_query_floor"])
                          for row in result["selected_same_base_contrast"]},
                         {("public_x_gray_prefix", 33975),
                          ("public_x_hash", 889)})
        self.assertFalse(result["new_n83_solver_search_executed"])
        self.assertIsNone(result["selected_best_total_runtime"])

    def test_changed_parent_receipt_is_rejected_before_analysis(self):
        parent = json.loads(rank_floor.SUPPORT.read_text())
        parent["ordered_pairs"][0]["base_orbits"] += 1
        with TemporaryDirectory() as directory:
            path = Path(directory) / "altered.json"
            path.write_text(json.dumps(parent))
            original = rank_floor.SUPPORT
            try:
                rank_floor.SUPPORT = path
                with self.assertRaisesRegex(ValueError, "receipt changed"):
                    rank_floor.calculate()
            finally:
                rank_floor.SUPPORT = original


if __name__ == "__main__":
    unittest.main()
