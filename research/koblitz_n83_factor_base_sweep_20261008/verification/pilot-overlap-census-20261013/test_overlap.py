"""Independent arithmetic and population checks for the pilot overlap receipt."""

import json
import unittest

import overlap


class PilotOverlapTests(unittest.TestCase):
    def test_binary_field_and_frobenius_cycle(self):
        x, y = map(int, overlap.GENERATOR)
        self.assertTrue(overlap.on_curve(x, y))
        self.assertEqual(overlap.square(x), overlap.multiply(x, x))
        self.assertEqual(overlap.square(y), overlap.multiply(y, y))
        phase_x = set()
        for _ in range(overlap.N):
            self.assertNotIn(x, phase_x)
            phase_x.add(x)
            x, y = overlap.square(x), overlap.square(y)
        self.assertEqual([str(x), str(y)], overlap.GENERATOR)

    def test_saved_receipt_has_exact_reciprocal_overlap(self):
        row = json.loads(overlap.OUTPUT.read_text())
        pairs = row["ordered_pairs"]
        self.assertEqual((row["primary_archive_count"], len(pairs)), (27, 648))
        self.assertEqual((row["disjoint_pair_count"],
                          row["fully_absorbed_pool_count"],
                          row["singleton_residual_count"],
                          row["passes_current_min_two_residual_screen_count"]),
                         (378, 33, 57, 558))
        self.assertEqual(sum(row["overlap_histogram"].values()), 648)
        by_direction = {(p["base_object_blake3"], p["pool_object_blake3"]): p
                        for p in pairs}
        self.assertEqual(len(by_direction), 648)
        for pair in pairs:
            reverse = by_direction[(pair["pool_object_blake3"],
                                    pair["base_object_blake3"])]
            self.assertEqual(pair["shared_signed_orbits"],
                             reverse["shared_signed_orbits"])
            self.assertEqual(pair["shared_x_values"],
                             overlap.N * pair["shared_signed_orbits"])
            self.assertEqual(pair["residual_signed_orbits"],
                             pair["pool_orbits"] - pair["shared_signed_orbits"])
            self.assertEqual(pair["passes_current_min_two_residual_screen"],
                             pair["residual_signed_orbits"] >= 2)
            self.assertEqual(len(pair["source_derived_exact_model_shapes"]),
                             4 if pair["residual_signed_orbits"] >= 2 else 0)
        shapes = row["source_derived_exact_model_shape_summaries"]
        self.assertEqual(len(shapes), 4)
        self.assertTrue(all(shape["pair_count"] == 558 for shape in shapes))
        five_two = next(shape for shape in shapes
                        if (shape["summands"], shape["residual_summands"])
                           == (5, 2))
        self.assertEqual((five_two["minimum_first_pattern_total_variables"],
                          five_two["maximum_first_pattern_total_variables"]),
                         (91472, 97080))
        self.assertEqual(row["first_retained_pair"]["shared_signed_orbits"], 0)
        self.assertIsNone(row["selected_best_total_runtime"])
        self.assertFalse(row["new_sat_search_executed"])

    def test_hash_streams_are_empirically_disjoint(self):
        row = json.loads(overlap.OUTPUT.read_text())
        for name, counts in row["policy_pair_summaries"].items():
            if "public_x_hash" in name:
                self.assertEqual(counts["disjoint"], counts["ordered_pairs"])


if __name__ == "__main__":
    unittest.main()
