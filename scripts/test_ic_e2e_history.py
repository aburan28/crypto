"""Tests for scripts/ic_e2e_history.py."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ic_e2e_history as hist  # noqa: E402
from test_ic_e2e_benchmark import synthetic_report, PARAMS  # noqa: E402


def record(run_id, ratios, host="A", counters=None, sha="s"):
    rungs = {}
    for p, ratio in ratios.items():
        rungs[p] = {"degree": 31, "subgroup_bits": 20.0, "targets": 4, "ic_verified": 4, "rho_verified": 4,
                    "counters": counters if counters is not None else {"collection_trials": 128},
                    "ic_whole_process_seconds": 1.0, "rho_seconds_total": ratio, "whole_process_ratio": ratio,
                    "charged_ratio": 10.0, "whole_process_crossover": ratio > 1}
    return {"schema_version": 1, "run_id": run_id, "url": None, "started_at": run_id, "host": {"platform": host, "cpu_count": 4},
            "ic_binary_sha256": sha, "git": {}, "reference": "ref", "rungs": rungs}


class IntervalTests(unittest.TestCase):
    def test_t_interval_matches_hand_calculation(self):
        r = hist.interval([2.0, 2.2, 2.1, 2.3])
        self.assertEqual(r["n"], 4)
        self.assertAlmostEqual(r["mean"], 2.15)
        # sd = sqrt(((.15)^2+(.05)^2+(.05)^2+(.15)^2)/3) = 0.1291; half = 3.182*0.1291/2 = 0.2054
        self.assertAlmostEqual(r["ci95"][0], 2.15 - 0.2054, places=3)
        self.assertAlmostEqual(r["ci95"][1], 2.15 + 0.2054, places=3)

    def test_single_sample_has_no_interval(self):
        self.assertIsNone(hist.interval([2.0])["ci95"])


class ReportTests(unittest.TestCase):
    def test_claim_needs_three_runs_and_interval_above_one(self):
        recs = [record("1", {PARAMS: 2.2}), record("2", {PARAMS: 2.3}), record("3", {PARAMS: 2.25})]
        rep = hist.build_report(recs, None)
        r = rep["rungs"][PARAMS]
        self.assertTrue(r["runtime_claim_faster_than_rho_end_to_end"])
        self.assertEqual(rep["problems"], [])
        two = hist.build_report(recs[:2], None)["rungs"][PARAMS]
        self.assertFalse(two["runtime_claim_faster_than_rho_end_to_end"])
        noisy = [record("1", {PARAMS: 1.2}), record("2", {PARAMS: 0.9}), record("3", {PARAMS: 1.3})]
        self.assertFalse(hist.build_report(noisy, None)["rungs"][PARAMS]["runtime_claim_faster_than_rho_end_to_end"])

    def test_counter_drift_between_runs_voids_the_claim(self):
        recs = [record("1", {PARAMS: 2.2}), record("2", {PARAMS: 2.3}), record("3", {PARAMS: 2.25}, counters={"collection_trials": 64})]
        rep = hist.build_report(recs, None)
        self.assertFalse(rep["rungs"][PARAMS]["runtime_claim_faster_than_rho_end_to_end"])
        self.assertTrue(any("differ between runs" in p for p in rep["problems"]))

    def test_reference_mismatch_is_a_problem(self):
        recs = [record("1", {PARAMS: 2.2}), record("2", {PARAMS: 2.3}), record("3", {PARAMS: 2.25})]
        ref = {"rungs": {PARAMS: {"counters": {"collection_trials": 999}}}}
        rep = hist.build_report(recs, ref)
        self.assertTrue(any("frozen reference" in p for p in rep["problems"]))
        self.assertFalse(rep["rungs"][PARAMS]["runtime_claim_faster_than_rho_end_to_end"])

    def test_two_hosts_are_not_pooled_silently(self):
        recs = [record("1", {PARAMS: 2.2}), record("2", {PARAMS: 2.3}, host="B"), record("3", {PARAMS: 2.25})]
        rep = hist.build_report(recs, None)
        self.assertEqual(rep["distinct_hosts"], 2)
        self.assertTrue(any("more than one host" in p for p in rep["problems"]))

    def test_markdown_has_a_row_per_rung(self):
        recs = [record("1", {PARAMS: 2.2, "b.json": 0.5}), record("2", {PARAMS: 2.3, "b.json": 0.6}), record("3", {PARAMS: 2.25, "b.json": 0.55})]
        text = hist.render_markdown(hist.build_report(recs, None))
        self.assertEqual(text.count("| `k0n-test.json` |"), 1)
        self.assertEqual(text.count("| `b.json` |"), 1)
        self.assertIn("faster than rho, end to end", text)
        self.assertIn("not established", text)


class AddTests(unittest.TestCase):
    def test_extract_from_an_artifact_and_refuse_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run = root / "art"
            (run / "k0n-test").mkdir(parents=True)
            (run / "k0n-test" / "report.json").write_text(json.dumps(synthetic_report()))
            (run / "manifest.json").write_text(json.dumps({"started_at": "t", "host": {"platform": "P", "cpu_count": 4},
                                                            "ic_binary": {"sha256": "abc"}, "git": {"commit": "c"}}))
            (run / "summary.json").write_text(json.dumps({"ok": True, "reference": {"path": "ref"}, "rungs": [{
                "params": PARAMS, "degree": 31, "subgroup_bits": 20.0, "targets": 4, "ic_verified": 4, "rho_verified": 4,
                "wall": {"now": {"ic_whole_process_seconds": 0.6, "rho_seconds_total": 2.0, "whole_process_ratio": 2.0 / 0.6,
                                 "charged_ratio": 20.0, "whole_process_crossover": True}}}]}))
            out = root / "runs"
            self.assertEqual(hist.main(["add", "--run-dir", str(run), "--run-id", "r1", "--out-dir", str(out)]), 0)
            rec = json.loads((out / "r1.json").read_text())
            self.assertEqual(rec["rungs"][PARAMS]["counters"]["collection_trials"], 128)
            self.assertEqual(rec["ic_binary_sha256"], "abc")
            self.assertEqual(hist.main(["add", "--run-dir", str(run), "--run-id", "r1", "--out-dir", str(out)]), 2)
            code = hist.main(["report", "--runs", str(out / "r1.json"), "--summary-json", str(root / "s.json")])
            self.assertEqual(code, 0)
            s = json.loads((root / "s.json").read_text())
            self.assertIsNone(s["rungs"][PARAMS]["ratio"]["ci95"])


if __name__ == "__main__":
    unittest.main()
