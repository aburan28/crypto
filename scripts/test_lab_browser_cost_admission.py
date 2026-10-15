"""UI admission checks using frozen native receipts; no solver or arithmetic."""
import copy
import json
from pathlib import Path
import unittest
from unittest.mock import patch
import build_lab_browser as browser

SESSION = "research/conductor_gap_cm_orbits_20261015/sessions/primary-0"


class NativeCostAdmission(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(browser.ROOT) / SESSION / "records.jsonl"
        cls.records = [json.loads(line) for line in path.read_text().splitlines()]

    def test_unpriced_native_receipts_are_searchable_without_becoming_s_zero(self):
        curves = {r["workload"]["curve"]["slug"]: {"ecbench": [], "factor_bases": []} for r in self.records}
        with patch.object(browser, "ECBENCH_SESSIONS_GLOBS", (SESSION,)):
            sessions, methods, _, yields = browser.ecbench_sessions(curves)
        self.assertEqual(len(sessions), 1)
        self.assertTrue(any(m["method"] == "ic.cm_folded_native_v1" for m in methods))
        for curve in curves.values():
            for entry in curve["ecbench"]:
                if entry["arm"] == "rho":
                    self.assertGreater(entry["mean_s"], 0)
                else:
                    self.assertIsNone(entry["mean_s"])
                    self.assertIsNone(entry["mean_ratio_to_floor"])
                    self.assertGreater(entry["unpriced_zero_runs"], 0)
        self.assertEqual(len(yields), 50)
        for entry in yields:
            self.assertIsNone(entry["s"])
            self.assertIsNone(entry["total_gae"])
            self.assertEqual(entry["charged_gae_lower_bound"], 0)

    def test_zero_without_unpriced_work_is_preserved_as_a_count(self):
        record = copy.deepcopy(next(r for r in self.records if r["arm"] == "velu"))
        record["cost"]["unpriced"] = []
        row = browser.yield_row(record, record["session_id"])
        self.assertEqual(row["s"], 0)
        self.assertEqual(row["total_gae"], 0)


if __name__ == "__main__":
    unittest.main()
