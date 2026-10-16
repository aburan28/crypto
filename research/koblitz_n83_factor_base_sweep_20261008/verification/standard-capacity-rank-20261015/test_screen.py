"""Boundary and frozen-input checks for the conditional N83 capacity screen."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parent / "screen.py"
SPEC = importlib.util.spec_from_file_location("n83_standard_capacity", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SCREEN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SCREEN)
AUDIT_SPEC = importlib.util.spec_from_file_location("n83_standard_capacity_audit",
                                                SCRIPT.with_name("audit.py"))
assert AUDIT_SPEC is not None and AUDIT_SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(AUDIT_SPEC)
AUDIT_SPEC.loader.exec_module(AUDIT)


class ScreenTests(unittest.TestCase):
    def test_half_rank_ceiling_integer_boundary(self) -> None:
        floor = SCREEN.rank_floor(columns=3, multisets=2, order=7)
        self.assertEqual(floor, 5)
        self.assertLess((floor - 1) * 2 * 2, 3 * (7 - 1))
        self.assertGreaterEqual(floor * 2 * 2, 3 * (7 - 1))

    def test_changed_parent_receipt_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.json"
            path.write_text('{"status":"PASS"}\n')
            digest = SCREEN.sha256(path.read_bytes())
            self.assertEqual(SCREEN.load_pinned(path, digest)["status"], "PASS")
            path.write_text('{"status":"PASS"} \n')
            with self.assertRaisesRegex(ValueError, "pinned input changed"):
                SCREEN.load_pinned(path, digest)

    def test_archived_frontier_and_reference_cells(self) -> None:
        receipt = SCREEN.generate()
        self.assertEqual(len(receipt["rows"]), 85)
        self.assertEqual([m["source_dimension"] for m in SCREEN.admitted_manifests()],
                         list(range(2, 19)))
        cells = {(row["dimension"], row["arity"]): row for row in receipt["rows"]}
        self.assertEqual(cells[15, 6]["queries_before_half_column_row_ceiling"], 11718)
        self.assertEqual(cells[18, 5]["queries_before_half_column_row_ceiling"], 15525)
        self.assertEqual(cells[18, 5]["nominal_all_affine_chain_variables"], 339)
        self.assertEqual(cells[16, 6]["nominal_all_affine_chain_equations"], 415)

    def test_rehashed_semantic_multiset_corruption_is_rejected(self) -> None:
        receipt = SCREEN.generate()
        row = next(row for row in receipt["rows"]
                   if row["dimension"] == 16 and row["arity"] == 6)
        row["unordered_multisets"] = str(int(row["unordered_multisets"]) + 1)
        receipt["full_table_sha256"] = SCREEN.sha256(json.dumps(
            receipt["rows"], sort_keys=True, separators=(",", ":")).encode())
        with self.assertRaisesRegex(ValueError, "multiset count"):
            AUDIT.audit(receipt)


if __name__ == "__main__":
    unittest.main()
