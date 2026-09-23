#!/usr/bin/env python3
"""Focused custody, accounting and scoring tests for the native-F4 arm."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import unittest

import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4


class SelectionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.rows = [
            {"blind_instance_id": "b-one", "cell_id": "n31"},
            {"blind_instance_id": "b-two", "cell_id": "n59"},
            {"blind_instance_id": "b-three", "cell_id": "n59"},
        ]
        self.bundle = {"instances": self.rows}

    def test_cell_selection_preserves_authenticated_order_and_cap(self) -> None:
        args = SimpleNamespace(
            blind_instance_id=None, cell="n59", all=False, max_instances=1
        )
        self.assertEqual(
            native_f4.select_instances(self.bundle, args), [self.rows[1]]
        )

    def test_named_selection_rejects_missing_ids(self) -> None:
        args = SimpleNamespace(
            blind_instance_id=["b-one", "b-missing"],
            cell=None,
            all=False,
            max_instances=None,
        )
        with self.assertRaises(phase_b.PhaseBError):
            native_f4.select_instances(self.bundle, args)


class AccountingTests(unittest.TestCase):
    def test_process_resources_keep_wall_cpu_and_peak_separate(self) -> None:
        records = [
            {
                "metrics": {
                    "wall_seconds": 2.0,
                    "total_core_seconds": 1.5,
                    "single_core_seconds": 1.5,
                    "peak_rss_bytes": 100,
                }
            },
            {
                "metrics": {
                    "wall_seconds": 3.0,
                    "total_core_seconds": 2.5,
                    "single_core_seconds": 2.5,
                    "peak_rss_bytes": 250,
                }
            },
        ]
        got = native_f4.process_resources(records, single_thread_requested=True)
        self.assertEqual(got["summed_process_wall_seconds"], 5.0)
        self.assertEqual(got["total_core_seconds"], 4.0)
        self.assertEqual(got["single_core_seconds"], 4.0)
        self.assertEqual(got["maximum_individual_process_rss_bytes"], 250)
        parallel = native_f4.process_resources(
            records, single_thread_requested=False
        )
        self.assertEqual(parallel["total_core_seconds"], 4.0)
        self.assertIsNone(parallel["single_core_seconds"])

    def test_watchdog_timeout_needs_no_backend_json_identity(self) -> None:
        native_f4.validate_f4_result(
            {
                "solver": "native-f4",
                "status": "timeout_inconclusive",
                "timed_out": True,
                "source_instance_id": None,
                "backend_report": None,
            },
            {"source_instance": {"id_blake3": "source-id"}},
        )


class ScoringTests(unittest.TestCase):
    def make_run(self, root: Path, solver_status: str) -> tuple[Path, Path]:
        run = root / "run"
        run.mkdir()
        (run / "tasks").mkdir()
        blind_id = "b-" + "1" * 64
        source_id = "x-" + "2" * 64
        plan = {
            "schema": native_f4.RUN_PLAN_SCHEMA,
            "selected_blind_instance_ids": [blind_id],
        }
        phase_b.write_json_new(run / "execution-plan.json", plan)
        phase_b.write_json_new(
            run / "run-summary.json",
            {"schema": native_f4.RUN_SUMMARY_SCHEMA, "status": "complete"},
        )
        task_root = phase_b.task_directory(run, 0, blind_id)
        task_root.mkdir()
        phase_b.write_json_new(
            task_root / "task-result.json",
            {
                "schema": native_f4.TASK_SCHEMA,
                "blind_instance_id": blind_id,
                "source_system_id": source_id,
                "cell_id": "n59-l9-m3-standard-a1-f0",
                "result": {"status": solver_status},
            },
        )
        inventory = phase_b.all_regular_inventory(run, {"run-seal.json"})
        seal = {
            "schema": native_f4.RUN_SEAL_SCHEMA,
            "status": "native_f4_run_frozen_before_truth_scoring",
            "inventory": inventory,
            "inventory_sha256": phase_b.canonical_sha256(inventory),
        }
        seal["seal_payload_sha256"] = phase_b.canonical_sha256(seal)
        phase_b.write_json_new(run / "run-seal.json", seal)

        truth = root / "truth.json"
        rows = [
            {
                "backend": backend,
                "blind_instance_id": blind_id,
                "cell_id": "n59-l9-m3-standard-a1-f0",
                "source_system_id": source_id,
                "target_class": "decomposable",
            }
            for backend in ("native-xor", "wdsat", "cryptominisat")
        ]
        truth.write_text(
            json.dumps(
                {
                    "schema": "koblitz_pdp_phase_b_score.v1",
                    "full_panel_scored": True,
                    "rows": rows,
                }
            )
        )
        return run, truth

    def test_scoring_opens_truth_only_after_valid_run_seal(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            run, truth = self.make_run(root, "sat")
            output = root / "score"
            score = native_f4.score_panel(
                SimpleNamespace(run=run, truth=truth, output=output)
            )
            self.assertEqual(score["classification_counts"], {"true_positive": 1})
            self.assertTrue((output / "score-seal.json").is_file())

    def test_false_classification_fails_closed(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            run, truth = self.make_run(root, "unsat")
            with self.assertRaises(phase_b.PhaseBError):
                native_f4.score_panel(
                    SimpleNamespace(run=run, truth=truth, output=root / "score")
                )


if __name__ == "__main__":
    unittest.main()
