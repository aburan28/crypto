#!/usr/bin/env python3
"""Integrity and scientific-contract checks for the fixed scaling experiment."""
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
spec = importlib.util.spec_from_file_location("rank_analysis", HERE / "analyze.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="boolean-rank-evidence-")
        self.root = Path(self.temp.name) / "run"
        shutil.copytree(HERE / "run_01", self.root)

    def tearDown(self):
        self.temp.cleanup()

    def rejected(self):
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(AssertionError):
            analysis.main(self.root)

    def mutate_record(self, n, predicate, mutate):
        path = self.root / f"raw-n{n}.jsonl"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        row = next(r for r in rows if predicate(r))
        cell = row["cell"]
        mutate(row)
        row["raw_line"] = json.dumps({k:v for k,v in row.items() if k not in ("cell","split","raw_line")}) + "\n"
        path.write_text("".join(json.dumps(r) + "\n" for r in rows))
        receipts_path = self.root / "receipts.json"
        receipts = json.loads(receipts_path.read_text())
        receipt = next(r for r in receipts if r["cell"] == cell)
        receipt["stdout_sha256"] = hashlib.sha256("".join(r["raw_line"] for r in rows if r["cell"] == cell).encode()).hexdigest()
        receipts_path.write_text(json.dumps(receipts))

    def test_complete_evidence_and_manifest_replay(self):
        for name, expected in json.loads((self.root/"manifest.json").read_text())["files"].items():
            self.assertEqual(analysis.sha(self.root/name), expected, name)
        expected = (self.root/"results.json").read_bytes()
        report = (self.root/"RESULT.md").read_bytes()
        with contextlib.redirect_stdout(io.StringIO()):
            analysis.main(self.root)
        self.assertEqual(expected, (self.root/"results.json").read_bytes())
        self.assertEqual(report, (self.root/"RESULT.md").read_bytes())
        result = json.loads(expected)
        self.assertEqual([x["n"] for x in result["dense_not_executed"]], [20,28,36])
        self.assertTrue(all(x["total_ns"] is None for x in result["dense_not_executed"]))

    def test_changed_source_rejected(self):
        with (self.root/"worker.rs").open("a") as out:
            out.write("// post-measurement change\n")
        self.rejected()

    def test_sparse_successor_replays_with_additive_legend_correction(self):
        successor=Path(self.temp.name)/"successor"
        shutil.copytree(HERE/"run_02",successor)
        for name,expected in json.loads((successor/"manifest.json").read_text())["files"].items():
            self.assertEqual(analysis.sha(successor/name),expected,name)
        expected=(successor/"results.json").read_bytes()
        original=(successor/"RESULT.md").read_bytes()
        primary=json.loads((self.root/"protocol.json").read_text())
        protocol=json.loads((successor/"protocol.json").read_text())
        self.assertFalse(set(primary['holdout_seeds']) & set(protocol['holdout_seeds']))
        self.assertEqual(protocol['predecessor_worker_sha256'],analysis.sha(self.root/'worker.rs'))
        self.assertEqual(protocol['repetitions'],20)
        with contextlib.redirect_stdout(io.StringIO()):
            analysis.main(successor)
        self.assertEqual(expected,(successor/"results.json").read_bytes())
        self.assertNotEqual(original,(successor/"RESULT.md").read_bytes())
        self.assertEqual((HERE/"run_02_REPORT_CORRECTED.md").read_bytes(),(successor/"RESULT.md").read_bytes())

    def test_correction_is_bound_to_original_and_changes_only_legend(self):
        correction=json.loads((HERE/'REPORT_CORRECTION.json').read_text())
        original=HERE/correction['original_report']
        corrected=HERE/correction['corrected_report']
        self.assertEqual(analysis.sha(original),correction['original_sha256'])
        self.assertEqual(analysis.sha(corrected),correction['corrected_sha256'])
        for entry in [correction['authority'],correction['unchanged_numerical_results']]:
            self.assertEqual(analysis.sha(HERE/entry['path']),entry['sha256'])
        self.assertEqual(original.read_text().replace('twelve balanced repetitions','20 balanced repetitions'),corrected.read_text())

    def test_missing_sample_rejected(self):
        path = self.root/"raw-n36.jsonl"
        path.write_text("".join(path.read_text().splitlines(keepends=True)[:-1]))
        self.rejected()

    def test_duplicate_receipt_rejected(self):
        path = self.root/"receipts.json"
        rows = json.loads(path.read_text())
        rows[-1] = rows[0]
        path.write_text(json.dumps(rows))
        self.rejected()

    def test_failed_worker_rejected(self):
        path = self.root/"receipts.json"
        rows = json.loads(path.read_text())
        rows[0]["exit_code"] = 1
        path.write_text(json.dumps(rows))
        self.rejected()

    def test_incomplete_campaign_rejected(self):
        path = self.root/"metadata.json"
        row = json.loads(path.read_text())
        row["complete"] = False
        path.write_text(json.dumps(row))
        self.rejected()

    def test_exponential_lookup_claim_rejected(self):
        self.mutate_record(36, lambda r:r["type"]=="sample" and r["variant"]=="ranked",
                           lambda r:r.update(lookup_entries=1<<36))
        self.rejected()

    def test_unexecuted_dense_arm_cannot_be_substituted(self):
        self.mutate_record(36, lambda r:r["type"]=="sample" and r["variant"]=="binary",
                           lambda r:r.update(variant="dense"))
        self.rejected()

    def test_restricted_mask_contract_checked(self):
        self.mutate_record(36, lambda r:r["type"]=="fixture" and r["family"]=="restricted_cycle",
                           lambda r:r.update(active=(1<<36)-1))
        self.rejected()

    def test_dropped_output_rows_rejected(self):
        self.mutate_record(36, lambda r:r["type"]=="sample" and r["variant"]=="ranked",
                           lambda r:r.update(total_rows=0))
        self.rejected()


if __name__ == "__main__":
    unittest.main()
