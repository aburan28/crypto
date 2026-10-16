"""Controls for the complete retained-pair source frontier."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest
from unittest import mock


sys.path.insert(0, str(Path(__file__).resolve().parent))
import pareto


class ParetoControls(unittest.TestCase):
    def test_strict_product_order_preserves_ties_and_tradeoffs(self):
        self.assertFalse(pareto.dominates((1, 2, 3), (1, 2, 3)))
        self.assertTrue(pareto.dominates((1, 2, 3), (1, 3, 3)))
        self.assertFalse(pareto.dominates((1, 3, 2), (2, 2, 3)))

    def test_saved_census_replays_and_frontier_is_nondominated(self):
        saved = json.loads(pareto.RECEIPT.read_text())
        self.assertEqual(pareto.make_receipt(), saved)
        self.assertEqual((saved["ordered_pair_count"],
                          saved["worker_admissible_pair_count"],
                          saved["unavailable_pair_count"]), (648, 558, 90))
        self.assertEqual((saved["source_pareto_pair_count"],
                          saved["source_pareto_unique_shape_count"]), (311, 24))
        frontier = [pareto.metrics(row) for row in saved["frontier_shapes"]]
        self.assertEqual(len(frontier), len(set(frontier)))
        for left in frontier:
            self.assertFalse(any(pareto.dominates(right, left) for right in frontier))
        selected = saved["selected_paired_arms"]
        self.assertFalse(selected["gray_pool"]["source_pareto"])
        self.assertEqual(selected["gray_pool"]["strict_dominator_pair_count"], 114)
        self.assertTrue(selected["hash_pool"]["source_pareto"])
        self.assertEqual(selected["hash_pool"]["strict_dominator_pair_count"], 0)
        self.assertIsNone(saved["selected_best_total_runtime"])
        self.assertFalse(saved["new_n83_solver_search_executed"])

    def test_changed_parent_bytes_are_rejected(self):
        with mock.patch.dict(pareto.PINS, {pareto.SCREEN: "0" * 64}):
            with self.assertRaisesRegex(ValueError, "pinned source bytes changed"):
                pareto.make_receipt()


if __name__ == "__main__":
    unittest.main()
