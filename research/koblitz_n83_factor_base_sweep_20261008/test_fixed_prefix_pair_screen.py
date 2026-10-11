"""Exact small-group and frozen-input checks for the fixed-prefix screen."""

import json
from math import comb
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import fixed_prefix_pair_screen as screen


class FixedPrefixPairScreenTests(unittest.TestCase):
    def test_small_group_counts_actual_fixed_prefix_coverage(self):
        order = 17
        base = (1, 2, 4)
        pairs = {(a + b) % order for a in base for b in base}
        self.assertEqual(len(pairs), comb(len(base) + 1, 2))
        for prefixes in ((1,), (1, 5), (1, 5, 10)):
            covered = {(prefix + pair) % order for prefix in prefixes for pair in pairs}
            self.assertLessEqual(len(covered), len(prefixes) * len(pairs))
        bounds = screen.thresholds(order, len(base), 2)
        self.assertEqual(bounds["minimum_prefix_checks_for_half_single_hit_ceiling"], 2)
        self.assertEqual(bounds["minimum_prefix_checks_for_half_rank_ceiling"], 3)
        self.assertLess((bounds["minimum_prefix_checks_for_half_single_hit_ceiling"] - 1) *
                        len(pairs) * 2, order)
        self.assertGreaterEqual(bounds["minimum_prefix_checks_for_half_single_hit_ceiling"] *
                                len(pairs) * 2, order)

    def test_retained_k2048_is_bound_to_s3_receipt_and_exact_threshold(self):
        result = screen.calculate()
        self.assertEqual(len(result["cases"]), 21)
        self.assertEqual(result["subgroup_order"], str(screen.PRIMARY_ORDER))
        self.assertIsNone(result["selected_best_total_runtime"])
        row = next(row for row in result["cases"] if row["orbit_columns"] == 2048)
        self.assertEqual(row["point_count_basis"], "retained_k2048_manifest")
        self.assertEqual(row["point_count"], 339968)
        self.assertEqual(row["pair_multisets_with_repetition"], "57789290496")
        self.assertEqual(row["minimum_prefix_checks_for_half_single_hit_ceiling"], 20919547709261)
        self.assertEqual(row["minimum_prefix_checks_for_half_rank_ceiling"], 42843233708565920)
        for row in result["cases"]:
            pairs = int(row["pair_multisets_with_repetition"])
            single = row["minimum_prefix_checks_for_half_single_hit_ceiling"]
            rank = row["minimum_prefix_checks_for_half_rank_ceiling"]
            self.assertLess((single - 1) * 2 * pairs, screen.PRIMARY_ORDER)
            self.assertGreaterEqual(single * 2 * pairs, screen.PRIMARY_ORDER)
            self.assertLess((rank - 1) * 2 * pairs,
                            row["orbit_columns"] * screen.PRIMARY_ORDER)
            self.assertGreaterEqual(rank * 2 * pairs,
                                    row["orbit_columns"] * screen.PRIMARY_ORDER)

    def test_changed_s3_round_trip_fails_closed(self):
        with tempfile.TemporaryDirectory(prefix="n83-prefix-bound-") as temporary:
            upload = json.loads(screen.INPUTS["k2048_upload"].read_text())
            upload["objects"][0]["downloaded_hash_matches"] = False
            replacement = Path(temporary) / "upload.json"
            replacement.write_text(json.dumps(upload))
            inputs = dict(screen.INPUTS, k2048_upload=replacement)
            with patch.object(screen, "INPUTS", inputs):
                with self.assertRaisesRegex(ValueError, "S3 round trip"):
                    screen.calculate()


if __name__ == "__main__":
    unittest.main()
