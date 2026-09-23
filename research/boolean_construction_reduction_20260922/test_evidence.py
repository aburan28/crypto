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
spec = importlib.util.spec_from_file_location("reduction_analysis", HERE / "analyze.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="boolean-reduction-evidence-")
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
        self.assertIsNone(result["full_polynomial_solver_cost"])
        self.assertIsNone(result["rho_ratio"])

    def test_changed_source_rejected(self):
        with (self.root/"worker.rs").open("a") as out:
            out.write("// post-measurement change\n")
        self.rejected()

    def test_incidence_successor_replays_and_keeps_fresh_holdouts(self):
        successor=Path(self.temp.name)/"successor"
        shutil.copytree(HERE/"run_02",successor)
        for name,expected in json.loads((successor/"manifest.json").read_text())["files"].items():
            self.assertEqual(analysis.sha(successor/name),expected,name)
        protocol=json.loads((successor/"protocol.json").read_text())
        primary=json.loads((self.root/"protocol.json").read_text())
        self.assertFalse(set(protocol['holdout_seeds']) & set(primary['holdout_seeds']))
        self.assertEqual(protocol['predecessor_worker_sha256'],analysis.sha(self.root/'worker.rs'))
        self.assertEqual((successor/'worker.rs').read_bytes(),(HERE/'worker.rs').read_bytes())
        expected=(successor/'results.json').read_bytes()
        report=(successor/'RESULT.md').read_bytes()
        with contextlib.redirect_stdout(io.StringIO()): analysis.main(successor)
        self.assertEqual(expected,(successor/'results.json').read_bytes())
        self.assertEqual(report,(successor/'RESULT.md').read_bytes())

    def test_incidence_auxiliary_storage_cannot_disappear(self):
        shutil.rmtree(self.root)
        shutil.copytree(HERE/'run_02',self.root)
        self.mutate_record(36,lambda r:r['type']=='sample' and r['variant']=='incidence_reduce',
                           lambda r:r.update(auxiliary_max_bytes=0))
        self.rejected()

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

    def test_fused_phase_cannot_be_reported_as_zero(self):
        self.mutate_record(36, lambda r:r["type"]=="sample" and r["variant"]=="stream_reduce",
                           lambda r:r.update(construction_ns=0))
        self.rejected()

    def test_missing_reduction_phase_rejected(self):
        self.mutate_record(36, lambda r:r["type"]=="sample" and r["variant"]=="ranked_reduce",
                           lambda r:r.update(reduction_ns=None))
        self.rejected()

    def test_dropped_source_rows_rejected(self):
        self.mutate_record(36, lambda r:r["type"]=="sample" and r["variant"]=="stream_reduce",
                           lambda r:r.update(source_rows=0))
        self.rejected()

    def test_changed_logical_reduction_work_rejected(self):
        self.mutate_record(36, lambda r:r["type"]=="sample" and r["variant"]=="stream_reduce",
                           lambda r:r.update(row_xors=r["row_xors"]+1))
        self.rejected()

    def test_mismatched_canonical_rank_rejected(self):
        self.mutate_record(36, lambda r:r["type"]=="sample" and r["variant"]=="stream_reduce",
                           lambda r:r.update(total_rank=r["total_rank"]-1))
        self.rejected()

    def test_fake_dense_control_above_bridge_rejected(self):
        self.mutate_record(36, lambda r:r["type"]=="sample" and r["variant"]=="sparse_reduce",
                           lambda r:r.update(variant="dense_reduce"))
        self.rejected()

    def test_restricted_mask_contract_checked(self):
        self.mutate_record(36, lambda r:r["type"]=="fixture" and r["family"]=="restricted_cycle",
                           lambda r:r.update(active=(1<<36)-1))
        self.rejected()


if __name__ == "__main__":
    unittest.main()
