"""Synthetic admission-gate checks; the Rust test constructs the actual base."""
from __future__ import annotations

import copy
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import cold_base

HERE = Path(__file__).resolve().parent
FIXTURE = (HERE.parent / "autolab_orbit_extract_20260924" /
           "independent_replay_20260924_codex/base_header.jsonl.gz")


class CertifiedBaseAdmission(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certified = json.loads(gzip.decompress(FIXTURE.read_bytes()))

    def generated(self):
        result = {key: copy.deepcopy(self.certified[key])
                  for key in cold_base.MATCH_KEYS}
        result["kind"] = "point_defined_factor_base"
        result["cofactor"] = int(self.certified["cofactor"])
        result["scanned_x"] = self.certified["field_x_values_scanned"]
        result["selection_mode"] = "legacy_rank_fixture_lcg_v1"
        return result

    def run_gate(self, generated):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "producer.json"
            raw.write_text(json.dumps({
                "factor_base_input_path": None,
                "factor_base_input_hash": None,
                "factor_base_input_blake3": None,
                "compact_orbit_base_header": generated,
            }))
            receipt = cold_base.materialize(raw, FIXTURE,
                                            root / "fresh.jsonl", root / "receipt.json")
            header = json.loads((root / "fresh.jsonl").read_text())
            return receipt, header

    def test_exact_order_and_cofactor_type_normalization(self):
        receipt, header = self.run_gate(self.generated())
        self.assertEqual(receipt["classification"], "FRESH_BASE_MATCHES_CERTIFIED_ARM")
        self.assertTrue(receipt["point_order_identical"])
        self.assertEqual(header["cofactor"], self.certified["cofactor"])
        for key in cold_base.POINT_KEYS:
            self.assertEqual(header[key], self.certified[key])

    def test_point_index_change_is_rejected(self):
        altered = self.generated()
        altered["factor_base_point_coordinates"][0], altered["factor_base_point_coordinates"][1] = (
            altered["factor_base_point_coordinates"][1], altered["factor_base_point_coordinates"][0])
        with self.assertRaisesRegex(AssertionError, "factor_base_point_coordinates"):
            self.run_gate(altered)

    def test_selection_and_cofactor_changes_are_rejected(self):
        for field, replacement in (("selection_mode", "ascending_x_v1"),
                                   ("cofactor", 2)):
            altered = self.generated()
            altered[field] = replacement
            with self.subTest(field=field), self.assertRaises(AssertionError):
                self.run_gate(altered)


if __name__ == "__main__":
    unittest.main()
