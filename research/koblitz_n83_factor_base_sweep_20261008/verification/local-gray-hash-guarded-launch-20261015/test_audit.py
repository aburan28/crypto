#!/usr/bin/env python3
"""Fail-closed controls for the retained quiet-host watch audit."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import tempfile
import unittest

import audit


class ReceiptAuditTest(unittest.TestCase):
    def setUp(self) -> None:
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.path = Path(self.scratch.name)
        for name in ("config.json", "samples.jsonl", "result.json"):
            shutil.copyfile(audit.HERE / name, self.path / name)

    def test_replays_retained_receipt(self) -> None:
        receipt = audit.audit(self.path, check_live_output=False)
        self.assertEqual(receipt["status"], "PASS_BOUNDED_NO_LAUNCH_ONLY")

    def test_rejects_tampered_raw_cpu(self) -> None:
        path = self.path / "samples.jsonl"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        rows[0]["top"]["stdout"] = rows[0]["top"]["stdout"].replace(
            "0.0% idle", "50.0% idle")
        path.write_text("".join(json.dumps(row) + "\n" for row in rows))
        with self.assertRaises(AssertionError):
            audit.audit(self.path, check_live_output=False)

    def test_rejects_false_launch_result(self) -> None:
        path = self.path / "result.json"
        result = json.loads(path.read_text())
        result["launched"] = True
        path.write_text(json.dumps(result))
        with self.assertRaises(AssertionError):
            audit.audit(self.path, check_live_output=False)


if __name__ == "__main__":
    unittest.main()
