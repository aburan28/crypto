#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parent))
import compose_koblitz_stage26_stage32 as composition


class Stage26Stage32CompositionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.score = composition.compose_result()

    def test_only_two_expected_rows_are_corrected(self) -> None:
        correction = self.score["stage32_wdsat_correction"]
        self.assertEqual(
            [row["blind_instance_id"] for row in correction["replacements"]],
            composition.CORRECTED_IDS,
        )
        self.assertEqual(correction["solver_error_rows_before"], 2)
        self.assertEqual(correction["solver_error_rows_after"], 0)
        self.assertFalse(correction["classification_counts_changed"])
        self.assertTrue(
            all(row["corrected_status"] == "timeout_inconclusive" for row in correction["replacements"])
        )

    def test_classification_and_resource_totals(self) -> None:
        self.assertEqual(
            self.score["classification_counts"],
            {"inconclusive": 219, "true_negative": 120, "true_positive": 141},
        )
        self.assertAlmostEqual(self.score["backend_resources"]["wdsat"]["total_core_seconds"], 7404.749336)
        self.assertAlmostEqual(
            self.score["stage32_wdsat_correction"]["charged_increment_total_core_seconds"],
            426.718011,
        )
        self.assertAlmostEqual(self.score["charged_total_core_seconds_available"], 21038.596137)
        self.assertEqual(self.score["maximum_recorded_individual_or_tree_rss_bytes"], 1526329344)

    def test_open_gates_remain_explicit(self) -> None:
        self.assertFalse(self.score["licensed_magma_f4_same_instance_panel_complete"])
        self.assertFalse(self.score["independent_external_reproduction_satisfied"])
        self.assertFalse(self.score["full_cost_gate_passed"])
        self.assertFalse(self.score["koblitz_index_calculus_sota"])
        self.assertIn("unaffiliated reproduction and novelty review remain absent", self.score["remaining_blockers"])


if __name__ == "__main__":
    unittest.main()
