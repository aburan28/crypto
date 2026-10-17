"""Focused guards for the frozen hosted six-summand comparison."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("m6_plan", HERE / "plan.py")
assert SPEC and SPEC.loader
plan = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(plan)


class FrozenPlanTests(unittest.TestCase):
    def test_saved_plan_replays(self) -> None:
        self.assertEqual(json.loads(plan.SAVED.read_text()), plan.build())

    def test_matched_trial_and_shape(self) -> None:
        row = plan.build()
        self.assertEqual([a["role"] for a in row["arms"]], ["gray_pool", "hash_pool"])
        self.assertEqual(len({a["base_object_blake3"] for a in row["arms"]}), 1)
        for arm in row["arms"]:
            argv = arm["worker_search_argv"]
            self.assertEqual(argv[1], "run-cross-pilot")
            self.assertEqual(argv[8], "/corpus.json")
            self.assertEqual(argv[9:14], ["2300", "20000", "64", "4", "2"])
            self.assertTrue(arm["source_capacity"]["current_worker_cap_pass"])
        self.assertEqual(row["selected_best_total_runtime"], None)

    def test_hash_pin_rejects_changed_receipt(self) -> None:
        original = plan.HASHES["support"]
        try:
            plan.HASHES["support"] = "0" * 64
            with self.assertRaisesRegex(ValueError, "frozen prior plan or support receipt"):
                plan.build()
        finally:
            plan.HASHES["support"] = original


if __name__ == "__main__":
    unittest.main()
