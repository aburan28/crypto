import contextlib
import io
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

import study
import report


class ArithmeticTests(unittest.TestCase):
    def test_exhaustive_field_against_independent_integer_reference(self):
        rng = random.Random(912)
        polynomials = [[], [0], [256], [0] * 25, [256] * 25]
        polynomials += [[rng.randrange(257) for _ in range(n)] for n in (2, 7, 25)]
        for coefficients in polynomials:
            for x in range(257):
                expected = study.reference(coefficients, x)
                self.assertEqual(study.power_sum(coefficients, x), expected)
                self.assertEqual(study.horner(coefficients, x), expected)

    def test_wire_format_preserves_boundary_values(self):
        values = (0, 1, 255, 256)
        self.assertEqual(study.pack(values), b"\x00\x00\x01\x00\xff\x00\x00\x01")
        self.assertEqual(study.unpack(study.pack(values)), values)

    def test_worker_accounts_for_costs_and_matches_inputs(self):
        receipts = [study.worker({"variant": v, "seed": 1, "count": 128}) for v in study.VARIANTS]
        for r in receipts:
            self.assertTrue(r["verified"])
            self.assertEqual(set(r["phases_ns"]), set(study.PHASES))
            self.assertGreaterEqual(r["worker_ns"], sum(r["phases_ns"].values()))
            self.assertEqual(r["input_bytes"], 2 * (25 + 128))
            self.assertEqual(r["output_bytes"], 2 * 128)
            self.assertIsNone(r["memory"]["python_traced_peak_bytes"])
        for key in ("input_sha256", "output_sha256"):
            self.assertEqual(receipts[0][key], receipts[1][key])

    def test_separate_memory_pass(self):
        r = study.worker({"variant": "horner", "seed": 1, "count": 64, "memory_pass": True})
        self.assertGreater(r["memory"]["python_traced_peak_bytes"], 0)
        self.assertTrue(r["verified"])

    def test_incorrect_arithmetic_is_detected(self):
        with patch.object(study, "horner", return_value=999):
            r = study.worker({"variant": "horner", "seed": 1, "count": 64})
        self.assertFalse(r["verified"])
        self.assertEqual(r["status"], "incorrect")

    def test_real_child_timeout_retained(self):
        r = study.launch({"variant": "power_sum", "seed": 1, "count": 4096, "memory_pass": False}, 0.000001)
        self.assertEqual(r["status"], "timeout")
        self.assertFalse(r["verified"])
        self.assertGreater(r["parent_wall_ns"], 0)

    def test_real_child_failure_retained(self):
        r = study.launch({"variant": "invalid", "seed": 1, "count": 64, "memory_pass": False}, 5)
        self.assertEqual(r["status"], "failed")
        self.assertFalse(r["verified"])


class StatisticsTests(unittest.TestCase):
    def test_censoring_is_not_dropped(self):
        rows = [{"status": "completed", "operations": 10},
                {"status": "censored", "operations": 100},
                {"status": "censored", "operations": 100}]
        result = study.coverage_summary(rows, 100)
        self.assertEqual(result["restricted_mean_operations"], 70)
        self.assertIsNone(result["completion_median"])
        self.assertIsNone(result["completion_p90"])
        self.assertEqual(result["censored"], 2)

    def test_null_control_has_no_speedup(self):
        result = study.bootstrap_ratio([(10, 10), (100, 100), (7, 7)], 77)
        self.assertEqual(result["ratio"], 1)
        self.assertEqual(result["ci95"], [1, 1])

    def test_constant_speedup_and_invalid_ratios(self):
        result = study.bootstrap_ratio([(20, 10), (200, 100)], 77)
        self.assertEqual(result["ci95"], [2, 2])
        with self.assertRaises(ValueError):
            study.bootstrap_ratio([(0, 2)], 1)

    def test_reproducible_independent_split_seeds(self):
        seeds = [study.seed_for(split, trial) for split in ("discovery", "holdout") for trial in range(128)]
        self.assertEqual(len(set(seeds)), len(seeds))
        self.assertEqual(study.seed_for("discovery", 0), seeds[0])

    def test_coverage_cap_and_determinism(self):
        for variant in ("nearest", "refresh"):
            self.assertEqual(study.cover_trial(variant, 32, 4, 1)["status"], "censored")
            self.assertEqual(study.cover_trial(variant, 32, 4, 8192), study.cover_trial(variant, 32, 4, 8192))
            self.assertEqual(study.cover_trial(variant, 32, 4, 8192)["covered"], 32)


class IntegrationTests(unittest.TestCase):
    def test_archived_evidence_hashes_summaries_and_trial_counts(self):
        for name in ("local-20260911", "final-20260911"):
            directory = Path(__file__).parent / "results" / name
            manifest, summary, rows = report.load_verified(directory)
            expected = 8 * manifest["trials_per_split_and_size"] + 4 * manifest["arithmetic_pairs_per_split"] + 4
            self.assertEqual(len(rows), expected)
            self.assertTrue(all(r["all_verified_and_paired"] for r in summary["arithmetic"].values()))

    def test_complete_receipts_no_overwrite_and_failure_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "result"
            with contextlib.redirect_stdout(io.StringIO()):
                summary = study.run_studies(directory, trials=2, arithmetic_trials=2, count=64, budget=32)
            rows = [json.loads(line) for line in (directory / "trials.jsonl").read_text().splitlines()]
            self.assertEqual(len(rows), 28)
            report.render(directory)
            self.assertTrue((directory / "REPORT.md").is_file())
            for split in ("discovery", "holdout"):
                self.assertTrue(summary["arithmetic"][split]["all_verified_and_paired"])
                self.assertEqual(summary["coverage"][f"{split}/128"]["nearest"]["censored"], 2)
            with self.assertRaises(FileExistsError):
                study.run_studies(directory)
            # A failed trial must invalidate the comparison, not be filtered out.
            failure = next(r for r in rows if r["study"] == "arithmetic" and not r["memory_pass"])
            failure["status"] = "timeout"
            failure["verified"] = False
            failed = study.summarize(rows, 32)["arithmetic"][failure["split"]]
            self.assertFalse(failed["all_verified_and_paired"])
            self.assertIsNone(failed["kernel_ratio"])
            self.assertEqual(failed["decision"], "inconclusive")
            with (directory / "trials.jsonl").open("a") as ledger:
                ledger.write("{}\n")
            with self.assertRaisesRegex(ValueError, "ledger hash mismatch"):
                report.load_verified(directory)

    def test_different_fixture_digests_block_speedup_claim(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "result"
            with contextlib.redirect_stdout(io.StringIO()):
                study.run_studies(directory, trials=2, arithmetic_trials=2, count=64, budget=32)
            rows = [json.loads(line) for line in (directory / "trials.jsonl").read_text().splitlines()]
            row = next(r for r in rows if r["study"] == "arithmetic" and not r["memory_pass"])
            row["input_sha256"] = "different"
            result = study.summarize(rows, 32)["arithmetic"][row["split"]]
            self.assertFalse(result["all_verified_and_paired"])
            self.assertIsNone(result["full_wall_ratio"])


if __name__ == "__main__":
    unittest.main()
