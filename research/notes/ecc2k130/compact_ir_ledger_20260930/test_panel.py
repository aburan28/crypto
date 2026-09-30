"""Fail-closed source/input freeze and Callgrind-total parser tests."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from run_panel import load_cell, parse_ir
from verify_panel import replay_ir


class FrozenPanelTest(unittest.TestCase):
    def test_all_six_cells_resolve_hashed_public_points_and_labels(self) -> None:
        for cell_id in ("n37_L1", "n37_L1024", "n41_L1", "n41_L1024",
                        "n53_L1", "n53_L1024"):
            with self.subTest(cell_id=cell_id):
                _config, cell, spec, points, fixture = load_cell(cell_id)
                self.assertEqual(cell["L"], spec["L"])
                self.assertTrue(points.is_file() and fixture.is_file())

    def test_callgrind_total_uses_ir_column_and_rejects_partial_profiles(self) -> None:
        cases = (
            ("events: Dr Ir Dw\nsummary: 2 101 3\ntotals: 4 103 5\n", 103),
            ("events: Ir\nsummary: 53\n", 53),
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "callgrind.out"
            for content, expected in cases:
                path.write_text(content)
                self.assertEqual(parse_ir(path), expected)
                self.assertEqual(replay_ir(path), expected)
            for content in ("events: Dr\nsummary: 4\n",
                            "events: Ir\n", "events: Ir Ir\nsummary: 1 2\n",
                            "events: Ir Dr\nsummary: 0 1\n"):
                path.write_text(content)
                with self.assertRaises(AssertionError):
                    parse_ir(path)
                with self.assertRaises(AssertionError):
                    replay_ir(path)


if __name__ == "__main__":
    unittest.main()
