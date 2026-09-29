#!/usr/bin/env python3
"""Regression checks for immutable evidence and the coefficient-reuse contract."""
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("envelope_analysis", HERE / "analyze.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="envelope-evidence-")
        self.root = Path(self.temp.name) / "run"
        shutil.copytree(HERE / "run_01", self.root)

    def tearDown(self):
        self.temp.cleanup()

    def rejected(self):
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(AssertionError):
            analysis.main(self.root)

    def rewrite_cell(self, predicate, mutate):
        # Update the receipt hash as well: the semantic contract must catch the
        # false claim even if the raw transport bytes are internally consistent.
        path = self.root / "raw-n6.jsonl"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        row = next(r for r in rows if predicate(r))
        cell = row["cell"]
        mutate(row)
        row["raw_line"] = json.dumps({k: v for k, v in row.items()
                                      if k not in ("cell", "split", "raw_line")}) + "\n"
        path.write_text("".join(json.dumps(r) + "\n" for r in rows))
        receipts_path = self.root / "receipts.json"
        receipts = json.loads(receipts_path.read_text())
        receipt = next(r for r in receipts if r["cell"] == cell)
        receipt["stdout_sha256"] = hashlib.sha256("".join(
            r["raw_line"] for r in rows if r["cell"] == cell).encode()).hexdigest()
        receipts_path.write_text(json.dumps(receipts))

    def test_manifest_and_complete_replay(self):
        manifest = json.loads((self.root / "manifest.json").read_text())
        for name, expected in manifest["files"].items():
            self.assertEqual(analysis.sha(self.root / name), expected, name)
        expected = (self.root / "results.json").read_bytes()
        report = (self.root / "RESULT.md").read_bytes()
        with contextlib.redirect_stdout(io.StringIO()):
            analysis.main(self.root)
        self.assertEqual(expected, (self.root / "results.json").read_bytes())
        self.assertEqual(report, (self.root / "RESULT.md").read_bytes())

    def test_changed_source_rejected(self):
        with (self.root / "worker.rs").open("a") as out:
            out.write("// post-measurement mutation\n")
        self.rejected()

    def test_confirmation_uses_unchanged_worker_and_replays(self):
        confirmation = Path(self.temp.name) / "confirmation"
        shutil.copytree(HERE / "run_02", confirmation)
        self.assertEqual((self.root / "worker.rs").read_bytes(), (confirmation / "worker.rs").read_bytes())
        protocol = json.loads((confirmation / "protocol.json").read_text())
        self.assertEqual(protocol["primary_worker_sha256"], analysis.sha(self.root / "worker.rs"))
        primary = json.loads((self.root / "protocol.json").read_text())
        self.assertFalse(set(primary["discovery_seeds"] + primary["holdout_seeds"])
                         & set(protocol["discovery_seeds"] + protocol["holdout_seeds"]))
        for field in ("performance_gate", "dramatic_gate", "variants", "families", "batches", "repetitions", "limits"):
            self.assertEqual(primary[field], protocol[field], field)
        for name, expected in json.loads((confirmation / "manifest.json").read_text())["files"].items():
            self.assertEqual(analysis.sha(confirmation / name), expected, name)
        expected = (confirmation / "results.json").read_bytes()
        with contextlib.redirect_stdout(io.StringIO()):
            analysis.main(confirmation)
        self.assertEqual(expected, (confirmation / "results.json").read_bytes())

    def test_packed_direct_context_reuse_is_checked(self):
        self.rewrite_cell(
            lambda r: r["type"] == "sample" and r["variant"] == "packed_direct" and "escape-b4" in r["cell"],
            lambda r: r.update(hits=3, fallbacks=1, changed_hits=2))
        self.rejected()

    def test_packed_envelope_missing_hit_rejected(self):
        self.rewrite_cell(
            lambda r: r["type"] == "sample" and r["variant"] == "packed_envelope" and "coefficients-b4" in r["cell"],
            lambda r: r.update(changed_hits=0))
        self.rejected()

    def test_missing_raw_sample_rejected(self):
        path = self.root / "raw-n6.jsonl"
        path.write_text("".join(path.read_text().splitlines(keepends=True)[:-1]))
        self.rejected()

    def test_duplicate_receipt_rejected(self):
        path = self.root / "receipts.json"
        rows = json.loads(path.read_text())
        rows[-1] = rows[0]
        path.write_text(json.dumps(rows))
        self.rejected()

    def test_failed_worker_rejected(self):
        path = self.root / "receipts.json"
        rows = json.loads(path.read_text())
        rows[0]["exit_code"] = 1
        path.write_text(json.dumps(rows))
        self.rejected()

    def test_incomplete_campaign_rejected(self):
        path = self.root / "metadata.json"
        metadata = json.loads(path.read_text())
        metadata["complete"] = False
        path.write_text(json.dumps(metadata))
        self.rejected()

    def test_missing_changed_hit_rejected(self):
        self.rewrite_cell(
            lambda r: r["type"] == "sample" and r["variant"] == "envelope" and "coefficients-b4" in r["cell"],
            lambda r: r.update(changed_hits=0))
        self.rejected()

    def test_unreported_escape_fallback_rejected(self):
        self.rewrite_cell(
            lambda r: r["type"] == "sample" and r["variant"] == "envelope" and "escape-b4" in r["cell"],
            lambda r: r.update(hits=4, fallbacks=0, changed_hits=3))
        self.rejected()

    def test_missing_degree_drop_multipliers_rejected(self):
        self.rewrite_cell(
            lambda r: r["type"] == "fixture" and "degree_cycle-b4" in r["cell"],
            lambda r: r.update(newly_required_multipliers=0))
        self.rejected()

    def test_envelope_expansion_rejected(self):
        self.rewrite_cell(
            lambda r: r["type"] == "fixture" and "escape-b4" in r["cell"],
            lambda r: r["envelope"][0].append(7))
        self.rejected()


if __name__ == "__main__":
    unittest.main()
