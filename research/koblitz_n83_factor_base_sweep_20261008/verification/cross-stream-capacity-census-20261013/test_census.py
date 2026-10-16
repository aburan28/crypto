"""Population, boundary and retained-receipt checks for the source census."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import census


class CrossStreamCapacityCensusTests(unittest.TestCase):
    def test_saved_census_covers_every_v4_recipe_and_shape(self):
        result = census.build()
        self.assertEqual(census.canonical(result), census.OUTPUT.read_bytes())
        self.assertEqual(result["ordered_pair_count"], 169884)
        self.assertEqual(result["shape_count_per_pair"], 4)
        self.assertEqual(len(result["size_shape_bounds"]), 484)
        for row in result["shape_cap_totals"]:
            self.assertEqual(sum(row["cap_classes"].values()), 169884)
            self.assertEqual(row["cap_classes"][census.NOT_APPLICABLE],
                             152460)
            self.assertEqual(sum(row["cap_classes"][name]
                                 for name in census.CAP_CLASSES), 17424)
        pilot = [row for row in result["source_class_shape_counts"]
                 if row["source_class"]
                 == "pilot_exact_import_gate_available_unexecuted"]
        self.assertEqual(len(pilot), 4)
        self.assertTrue(all(row["cap_classes"][census.CAP_CLASSES[0]] == 648
                            and sum(row["cap_classes"].values()) == 648
                            for row in pilot))
        self.assertIsNone(result["selected_best_total_runtime"])
        self.assertFalse(result["sat_search_executed_by_census"])

    def test_bounds_are_monotone_and_classes_have_witnesses(self):
        rows = census.build()["size_shape_bounds"]
        self.assertEqual({row["cap_class"] for row in rows},
                         set(census.CAP_CLASSES))
        for row in rows:
            self.assertLessEqual(row["minimum_total_variables"],
                                 row["maximum_total_variables"])
            self.assertLessEqual(row["minimum_domain_clauses"],
                                 row["maximum_domain_clauses"])
        known = census.bounds(600, 600, 5, 2)
        self.assertEqual((known["maximum_total_variables"],
                          known["maximum_domain_clauses"], known["xor_rows"]),
                         (97080, 269830, 4980))
        self.assertEqual(known["cap_class"], census.CAP_CLASSES[0])

    def test_retained_pair_and_raw_receipt_are_bound_to_audit(self):
        row = census.build()["retained_first_pair_replay"]
        self.assertEqual(row["pair_index"], 36953)
        self.assertEqual(row["shared_signed_orbits"], 0)
        self.assertEqual(row["search_status"], "UNKNOWN_search_cap")
        with TemporaryDirectory() as directory:
            changed = Path(directory) / "result.json"
            changed.write_bytes(census.RETAINED.read_bytes() + b" ")
            with patch.object(census, "RETAINED", changed):
                with self.assertRaisesRegex(ValueError,
                                            "frozen input changed"):
                    census.build()


if __name__ == "__main__":
    unittest.main()
