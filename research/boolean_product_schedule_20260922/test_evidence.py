#!/usr/bin/env python3
"""Small integrity regressions for the completed experiment's verifier."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("schedule_analysis", HERE / "analyze.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="boolean-evidence-test-")
        self.root = Path(self.temp.name) / "run"
        shutil.copytree(HERE / "run_01", self.root)

    def tearDown(self):
        self.temp.cleanup()

    def rejected(self):
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(AssertionError):
            analysis.main(self.root)

    def test_complete_evidence_replays(self):
        expected = (self.root / "results.json").read_bytes()
        with contextlib.redirect_stdout(io.StringIO()):
            analysis.main(self.root)
        self.assertEqual(expected, (self.root / "results.json").read_bytes())

    def test_changed_source_is_rejected(self):
        with (self.root / "worker.rs").open("a") as out:
            out.write("// changed after measurement\n")
        self.rejected()

    def test_missing_raw_sample_is_rejected(self):
        path = self.root / "raw-n6.jsonl"
        lines = path.read_text().splitlines(keepends=True)
        path.write_text("".join(lines[:-1]))
        self.rejected()

    def test_duplicate_receipt_is_rejected(self):
        path = self.root / "receipts.json"
        rows = json.loads(path.read_text())
        rows[-1] = rows[0]
        path.write_text(json.dumps(rows))
        self.rejected()

    def test_failed_worker_is_rejected(self):
        path = self.root / "receipts.json"
        rows = json.loads(path.read_text())
        rows[0]["exit_code"] = 1
        path.write_text(json.dumps(rows))
        self.rejected()

    def test_incomplete_campaign_is_rejected(self):
        path = self.root / "metadata.json"
        metadata = json.loads(path.read_text())
        metadata["complete"] = False
        path.write_text(json.dumps(metadata))
        self.rejected()


if __name__ == "__main__":
    unittest.main()
