"""Fail-closed qualification controls for the short-reservation replay."""

import importlib.util
import unittest
from pathlib import Path


SOURCE = Path(__file__).with_name("run_segments.py")
SPEC = importlib.util.spec_from_file_location("run_segments", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class QualificationTest(unittest.TestCase):
    def test_only_complete_exact_uncontended_block_qualifies(self):
        isolation = {
            "exit_status": 0,
            "left_on_reserved": {"user_threads": []},
            "contended_samples": 0,
        }
        paired = {"status": "complete", "runs": [{"status": "ok"} for _ in range(22)]}
        self.assertTrue(MODULE.qualifies(0, isolation, paired))
        self.assertFalse(MODULE.qualifies(1, isolation, paired))
        self.assertFalse(MODULE.qualifies(0, None, paired))
        self.assertFalse(MODULE.qualifies(0, isolation, None))
        self.assertFalse(MODULE.qualifies(0, {**isolation, "contended_samples": 1}, paired))
        self.assertFalse(MODULE.qualifies(0, {**isolation, "left_on_reserved": {}}, paired))
        self.assertFalse(MODULE.qualifies(0, {**isolation, "left_on_reserved": {"user_threads": ["other"]}}, paired))
        self.assertFalse(MODULE.qualifies(0, isolation, {**paired, "runs": paired["runs"][:-1]}))
        self.assertFalse(MODULE.qualifies(0, isolation, {**paired, "runs": paired["runs"][:-1] + [{"status": "failure"}]}))


if __name__ == "__main__":
    unittest.main()
