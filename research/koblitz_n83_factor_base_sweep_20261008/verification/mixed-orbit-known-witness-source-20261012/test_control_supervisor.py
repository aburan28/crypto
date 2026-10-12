#!/usr/bin/env python3
"""Fail-closed checks that do not launch the N83 control."""

import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import control_supervisor as gate


class ControlSupervisorTests(unittest.TestCase):
    def authorization(self, build: Path, output: Path) -> dict:
        return {
            "schema": "n83.additional-local-budget/v3",
            "study": gate.STUDY,
            "purpose": "one_mixed_orbit_known_witness_control",
            "prior_total_charged_seconds": "15000",
            "additional_authorized_seconds": "1500",
            "build_dir": str(build.resolve()),
            "output_dir": str(output.resolve()),
            "user_approval_reference": "explicit future approval fixture",
        }

    def test_authorization_rejects_stale_or_underfunded_grant(self):
        with tempfile.TemporaryDirectory() as folder:
            build, output = Path(folder) / "build", Path(folder) / "run"
            good = self.authorization(build, output)
            self.assertIsNone(gate.authorization_error(good, build, output))
            for key, value in (
                ("schema", "n83.additional-local-budget/v2"),
                ("prior_total_charged_seconds", "12600"),
                ("additional_authorized_seconds", "1499"),
                ("build_dir", str(output)),
                ("user_approval_reference", ""),
            ):
                altered = {**good, key: value}
                self.assertIsNotNone(gate.authorization_error(altered, build, output))

    def test_stale_grant_never_reaches_container_launcher(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            build, output = root / "build", root / "run"
            build.mkdir()
            small, envelope = root / "small", root / "envelope"
            small.write_bytes(b"")
            envelope.write_bytes(b"")
            stale = self.authorization(build, output)
            stale["prior_total_charged_seconds"] = "12600"
            auth = root / "authorization.json"
            auth.write_text(json.dumps(stale))
            with mock.patch.object(gate, "stage", side_effect=AssertionError("launched")):
                with self.assertRaisesRegex(ValueError, "authorization"):
                    gate.run(auth, build, output, small, envelope)
            self.assertFalse(output.exists())

    def test_capped_search_is_not_success_and_fake_success_fails(self):
        fixture = json.loads((gate.HERE / "fixture.json").read_text())
        config = {
            "schema": "n83.mixed-orbit-known-witness-config/v1",
            "scope": "known_answer_control_only", "study": gate.STUDY,
            "small_summands": 3, "residual_summands": 2,
            "conflict_budget": gate.CONFLICTS, "max_models": gate.MODELS,
            "budget_seconds": gate.WORKER_WALL,
            "memory_cgroup_limit_bytes": 4 * 1024**3,
            "memory_cgroup_swap_limit_bytes": 0,
        }
        result = {
            "schema": "n83.mixed-orbit-known-witness-result/v1",
            "scope": "known_answer_control_only", "target": fixture["equation"]["target"],
            "natural_relations": None, "total_index_calculus_runtime_ms": None,
            "selected_best_total_runtime": None, "candidate": None,
            "status": "UNKNOWN_search_cap",
            "search": {"solver_calls": 1, "conflicts": 5000, "models": 0,
                       "invalid_models": 0, "patterns_attempted": 1,
                       "exhausted": True, "refuted": False},
            "base_import_elapsed_ms_L0": 8.0, "worker_wall_ms_L0": 434.0,
            "memory_cgroup_peak_bytes": 468664320,
        }
        self.assertIsNone(gate.control_receipt_error(config, result, fixture,
                                                      Path("unused"), Path("unused")))
        fake = copy.deepcopy(result)
        fake["status"] = "PASS_known_witness_solver_access"
        self.assertIsNotNone(gate.control_receipt_error(config, fake, fixture,
                                                        Path("unused"), Path("unused")))
        fake = copy.deepcopy(result)
        fake["status"] = "UNKNOWN_preflight_cap"
        self.assertIsNotNone(gate.control_receipt_error(config, fake, fixture,
                                                        Path("unused"), Path("unused")))
        fake = copy.deepcopy(result)
        fake["natural_relations"] = 0
        self.assertIsNotNone(gate.control_receipt_error(config, fake, fixture,
                                                        Path("unused"), Path("unused")))

    def test_worker_wall_cap_requires_a_matching_receipt(self):
        cap = {"schema": "n83.mixed-orbit-known-witness-cap/v1",
               "status": "UNKNOWN_budget", "budget_seconds": gate.WORKER_WALL,
               "elapsed_process_wall_ms": 840001.0}
        self.assertTrue(gate.valid_worker_cap(cap))
        self.assertFalse(gate.valid_worker_cap({**cap, "budget_seconds": 839}))
        self.assertFalse(gate.valid_worker_cap({**cap, "elapsed_process_wall_ms": float("nan")}))


if __name__ == "__main__":
    unittest.main()
