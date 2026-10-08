#!/usr/bin/env python3
"""Small-field point-law controls for the variable-b support and [4] formulas."""
from __future__ import annotations

import unittest
from pathlib import Path
import tempfile

import produce
import verify


class SupportFormulaTest(unittest.TestCase):
    def test_all_x_and_b_against_two_independent_fields_and_group_law(self) -> None:
        polynomial = (1 << 5) | (1 << 2) | 1
        producer_field = produce.gate.Field(5, [0, 2])
        independent_field = verify.FastGF2m(5, polynomial)
        third_field = verify.ref.Field(5, polynomial)
        trace_mask = produce.gate.trace_mask(producer_field)
        for a2 in (0, 1):
            for b in (1, 2, 3):
                curve = verify.Koblitz(independent_field, a=a2, b=b)
                for x in range(1 << 5):
                    actual = produce.lift_and_project(
                        producer_field, trace_mask, a2, b, x)
                    self.assertEqual(actual, verify.row_result(
                        independent_field, a2, b, x))
                    if x:
                        inv = third_field.inverse(x)
                        rhs = x ^ a2 ^ third_field.mul(b, third_field.square(inv))
                        self.assertEqual(third_field.trace(rhs), actual[0] == 0)
                    points = curve.points_over(x)
                    self.assertEqual(len(points), actual[0])
                    for point in points:
                        self.assertTrue(curve.on_curve(point))
                        image = curve.mul(point, 4)
                        self.assertEqual(-1 if image is None else image[0], actual[1])

    def test_gray_update_matches_independent_mask_evaluation(self) -> None:
        basis = [3, 5, 9, 17]
        x = 0
        for ordinal in range(1 << len(basis)):
            if ordinal:
                x ^= basis[(ordinal & -ordinal).bit_length() - 1]
            mask = ordinal ^ (ordinal >> 1)
            self.assertEqual(x, verify.x_from_mask(basis, mask))

    def test_source_low_slot_replays_the_prior_frozen_census(self) -> None:
        config = produce.checked_inputs(require_lock=False)
        field = produce.gate.Field(131, [0, 1, 2, 13])
        bases, _ = produce.slot_bases(field, config["normal_beta"])
        trace_mask = produce.gate.trace_mask(field)
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "source_low_0"
            summary, columns = produce.scan(
                field, trace_mask, config["curves"][0], "low_0", bases["low_0"],
                output, config)
            self.assertEqual((summary["physical_points"],
                              summary["projected_sign_classes"]), (7977, 3988))
            self.assertEqual(summary["rows_sha256"],
                             "a3b7b7d01e8575b3a8950e9727705ec6993ead41f1b0a00f89e9524661244cbc")
            alternate = verify.ref.Field(131, int(config["field_modulus_hex"], 16))
            fast = verify.FastGF2m(131, int(config["field_modulus_hex"], 16))
            checked, replay_columns, point_images = verify.replay_scan(
                fast, alternate, config["curves"][0], "low_0", bases["low_0"],
                output, config)
            self.assertEqual(checked, summary)
            self.assertEqual(replay_columns, columns)
            self.assertEqual(point_images, 8)


if __name__ == "__main__":
    unittest.main()
