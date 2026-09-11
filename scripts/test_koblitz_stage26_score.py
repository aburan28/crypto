#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_koblitz_stage26_affinity as score


CELLS = {
    "n31-l5-m3-standard-a1-f0",
    "n31-l5-m3-ggmp-a0-f0",
    "n41-l5-m3-standard-a1-f0",
    "n59-l9-m3-standard-a1-f0",
}


class Stage26ScoreTests(unittest.TestCase):
    def test_completion_gate_audit_preserves_open_gates(self) -> None:
        value = score.completion_gate_audit()
        self.assertEqual(value["2_same_instance_backend_matrix"]["status"], "partial")
        self.assertIn("licensed Magma F4 execution on all 160 inputs", value["2_same_instance_backend_matrix"]["missing"])
        self.assertEqual(value["5_unknown_scalar"]["status"], "finite_degree23_complete")
        self.assertEqual(value["7_external_review"]["status"], "missing")
        self.assertIn("not a Koblitz index-calculus SOTA", value["overall"])

    def test_matched_direct_mitm_binding(self) -> None:
        value = score.matched_direct_mitm()
        self.assertEqual(value["instances"], 160)
        self.assertEqual(value["outcomes"], 160)
        self.assertEqual(value["classification_counts"], {"true_negative": 80, "true_positive": 80})
        self.assertIsNone(value["conflicts"])

    def test_unknown_scalar_control_binding(self) -> None:
        value = score.unknown_scalar_control()
        self.assertEqual(value["completed_unknown_scalar_targets"], 5)
        self.assertFalse(value["factor_base_logs_known_by_construction"])
        self.assertEqual(value["attempt_totals"]["relation_found"], 252)
        self.assertGreater(value["ratios"]["setup_charged_ic_over_rho_core"], 1600)

    def test_natural_relation_yield_binding(self) -> None:
        value = score.natural_relation_yield_control()
        self.assertEqual(value["natural"]["targets"], 256)
        self.assertEqual(value["natural"]["hits"], 163)
        self.assertFalse(value["factor_base"]["factor_base_discrete_log_labels_constructed"])
        self.assertEqual(value["verification_and_admission"]["measurement_admission_status"], "pending_independent_payload_replay")
        self.assertEqual(value["current_measurement_admission_status"], "finite_public_synthetic_internal_replay_complete")
        self.assertTrue(value["factor_base_materialization_replayed"])
        self.assertEqual(value["independent_internal_payload_replay"]["canonical_pairs"], 9165621)
        self.assertFalse(value["independent_external_reproduction_satisfied"])

    def test_elapsed_time(self) -> None:
        self.assertEqual(score.elapsed_time("1:02.50"), 62.5)
        self.assertEqual(score.elapsed_time("1:02:03"), 3723.0)
        with self.assertRaises(score.Stage26ScoreError):
            score.elapsed_time("7")

    def test_gnu_time_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.time"
            path.write_text(
                "\n".join(
                    (
                        "User time (seconds): 1.25",
                        "System time (seconds): 0.50",
                        "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.10",
                        "Maximum resident set size (kbytes): 1024",
                        "Exit status: 0",
                    )
                )
                + "\n"
            )
            receipt = score.parse_gnu_time(path)
            self.assertEqual(receipt["total_core_seconds"], 1.75)
            self.assertEqual(receipt["wall_seconds"], 2.1)
            self.assertEqual(receipt["peak_rss_bytes"], 1024 * 1024)

    def test_workflow_parallel_span(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "workflow.json"
            jobs = []
            for index, cell in enumerate(sorted(CELLS)):
                jobs.append(
                    {
                        "name": f"one-cpu-cell ({cell})",
                        "databaseId": 100 + index,
                        "status": "completed",
                        "conclusion": "success",
                        "startedAt": f"2026-09-11T18:2{index}:00Z",
                        "completedAt": f"2026-09-11T18:3{index}:00Z",
                        "url": f"https://example.test/{index}",
                    }
                )
            path.write_text(
                json.dumps(
                    {
                        "databaseId": score.EXPECTED_RUN,
                        "headSha": score.EXPECTED_COMMIT,
                        "status": "completed",
                        "conclusion": "success",
                        "createdAt": "2026-09-11T18:00:00Z",
                        "updatedAt": "2026-09-11T18:40:00Z",
                        "url": "https://example.test/run",
                        "jobs": jobs,
                    }
                )
            )
            result = score.workflow_accounting(path, CELLS)
            self.assertEqual(result["workflow_wall_seconds"], 2400)
            self.assertEqual(result["parallel_cell_job_span_seconds"], 780)
            self.assertEqual(set(result["cell_jobs"]), CELLS)


if __name__ == "__main__":
    unittest.main()
