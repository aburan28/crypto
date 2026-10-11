"""Fail-closed controls for the proposed N83 structured-orbit gate."""

import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import structured_orbit_search_supervisor as gate


def verification() -> dict:
    return {
        "schema": "n83.structured-orbit-import/v1", "status": "PASS",
        "curve_slug": gate.CURVE, "object_blake3": gate.OBJECT_BLAKE3,
        "points": 332166, "signed_orbits": 2001, "x_coordinates": 166083,
        "solver_search_executed": False,
    }


def search_records() -> tuple[dict, dict]:
    config = {
        "schema": "n83.structured-orbit-cold-config/v1", "study": gate.STUDY,
        "object_blake3": gate.OBJECT_BLAKE3,
        "manifest_blake3": gate.MANIFEST_BLAKE3,
        "public_corpus_path": "/corpus.json", "validation_path": "/validation.json",
        "budget_seconds": gate.WORKER_WALL, "max_trials": gate.TRIALS,
        "conflict_budget": gate.CONFLICTS, "max_models": gate.MODELS,
        "summands": 5, "strategy": "ChainedS3FrobeniusOrbit",
        "memory_cgroup_limit_bytes": 4 * 1024**3,
        "memory_cgroup_swap_limit_bytes": 0,
        "source_blake3": "a" * 64,
    }
    result = {
        "schema": "n83.structured-orbit-cold-result/v1", "study": gate.STUDY,
        "status": "UNKNOWN_trial_cap", "object_blake3": gate.OBJECT_BLAKE3,
        "manifest_blake3": gate.MANIFEST_BLAKE3,
        "public_corpus_canonical_json_blake3": gate.CORPUS_BLAKE3,
        "factor_base_size": 332166, "orbit_count": 2001,
        "trials": 1, "relations": 0, "independent_relations": 0,
        "dependent_relations": 0, "inconsistent_relations": 0,
        "verification_failures": 0, "sat_calls": 1, "sat_unknowns": 0,
        "sat_invalid_models": 0, "sat_models": 0, "sat_conflicts": 31,
        "matrix_rows": 0, "matrix_columns": 2002,
        "terminal_matrix_rank": 0, "linear_solve_attempts": 0,
        "attempt_records": [{}], "relation_matrix": [], "rank_history": [],
        "candidate_log": None, "candidate_verified": False,
        "memory_cgroup_peak_bytes": 190_000_000,
        "base_import_elapsed_ms_L0": 12.0, "worker_wall_ms_L0": 30.0,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    return config, result


class StructuredOrbitGateTests(unittest.TestCase):
    def test_default_search_image_matches_pinned_build_image(self):
        script = (Path(gate.__file__).parent / "verification" /
                  "structured-orbit-search-gate-20261010" / "build-linux.sh")
        self.assertIn(f'image = "{gate.IMAGE}"', script.read_text())
        self.assertEqual(gate.run.__defaults__[0], gate.IMAGE)

    def test_archive_and_search_receipts_reject_false_claims(self):
        self.assertIsNone(gate.verify_receipt_error(verification()))
        wrong = verification()
        wrong["x_coordinates"] -= 1
        self.assertIn("x_coordinates", gate.verify_receipt_error(wrong))
        config, result = search_records()
        self.assertIsNone(gate.search_receipt_error(config, result))
        result["status"] = "PASS_verified_target"
        result["candidate_verified"] = True
        result["candidate_log"] = "3"
        self.assertIn("status", gate.search_receipt_error(config, result))
        config, result = search_records()
        result["total_index_calculus_runtime_ms"] = 30.0
        self.assertIn("total_index_calculus_runtime_ms", gate.search_receipt_error(config, result))
        config, result = search_records()
        result["sat_unknowns"] = 1
        self.assertIn("trial-cap", gate.search_receipt_error(config, result))
        result["status"] = "UNKNOWN_solver_cap"
        result["sat_calls"] = 0  # construction can exhaust before search
        self.assertIsNone(gate.search_receipt_error(config, result))
        config, result = search_records()
        result["relations"] = result["independent_relations"] = 1
        result["relation_matrix"] = [{"row": ["0"], "summands": [], "summand_negated": []}]
        self.assertIn("column vector", gate.search_receipt_error(config, result))
        cap = {"schema": "n83.structured-orbit-cold-cap/v1", "status": "UNKNOWN_budget",
               "budget_seconds": gate.WORKER_WALL, "elapsed_process_wall_ms": 710_000.0}
        self.assertTrue(gate.valid_worker_cap(cap))
        cap["budget_seconds"] -= 1
        self.assertFalse(gate.valid_worker_cap(cap))

    def test_missing_additional_budget_cannot_launch(self):
        with tempfile.TemporaryDirectory(prefix="n83-orbit-budget-") as folder:
            root = Path(folder)
            auth = root / "authorization.json"
            auth.write_text(json.dumps({
                "schema": "n83.additional-local-budget/v1", "study": gate.STUDY,
                "purpose": "structured_orbit_import_and_one_m5_trial",
                "prior_pilot_charged_seconds": "7200",
                "output_dir": str(root / "output"),
                "build_dir": str(root),
                "additional_authorized_seconds": "1199",
                "user_approval_reference": "synthetic-test-only",
            }))
            binary, obj = root / "worker", root / "object"
            binary.write_bytes(b"\x7fELFsynthetic")
            obj.write_bytes(b"synthetic")
            with patch.object(gate.subprocess, "run", side_effect=AssertionError("Docker must not start")):
                with self.assertRaisesRegex(ValueError, "below the full conservative charge"):
                    gate.run(obj, binary, auth, root / "output")

    def test_offline_build_rejects_insufficient_budget_before_creating_output(self):
        script = (Path(gate.__file__).parent / "verification" /
                  "structured-orbit-search-gate-20261010" / "build-linux.sh")
        with tempfile.TemporaryDirectory(prefix="n83-orbit-build-guard-") as folder:
            root = Path(folder)
            output = root / "fresh-build"
            auth = root / "authorization.json"
            auth.write_text(json.dumps({
                "schema": "n83.additional-local-budget/v1", "study": gate.STUDY,
                "purpose": "structured_orbit_import_and_one_m5_trial",
                "prior_pilot_charged_seconds": "7200",
                "build_dir": str(output), "output_dir": str(root / "fresh-run"),
                "additional_authorized_seconds": "1799",
                "user_approval_reference": "synthetic-test-only",
            }))
            result = subprocess.run(["sh", str(script), str(output), str(auth)],
                                    cwd=gate.ROOT, capture_output=True, text=True, timeout=5)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("1800-second", result.stderr)
            self.assertFalse(output.exists())

    def test_timeout_kill_failure_still_returns_a_stage_receipt(self):
        with tempfile.TemporaryDirectory(prefix="n83-orbit-timeout-") as folder:
            output = Path(folder)

            class FakeChild:
                waits = 0

                def wait(self, timeout=None):
                    self.waits += 1
                    if self.waits == 1:
                        raise subprocess.TimeoutExpired("docker run", timeout)
                    return 137

                def kill(self):
                    pass

            def fake_run(command, **_kwargs):
                if command[:2] == ["docker", "kill"]:
                    raise subprocess.TimeoutExpired("docker kill", 15)
                if command[:2] == ["docker", "inspect"]:
                    return subprocess.CompletedProcess(command, 0,
                                                       json.dumps({"ExitCode": 137, "OOMKilled": False}), "")
                return subprocess.CompletedProcess(command, 0, "", "")

            with patch.object(gate.subprocess, "Popen", return_value=FakeChild()), \
                 patch.object(gate.subprocess, "run", side_effect=fake_run):
                receipt = gate.stage(name="verify", command=["docker", "run"], wall=1,
                                     output=output, container_name="synthetic")
            self.assertTrue(receipt["timed_out"])
            self.assertIn("docker kill failed", receipt["cleanup_error"])
            self.assertTrue((output / "verify.stdout").is_file())

    def test_import_failure_prevents_search_and_success_charges_both_caps(self):
        with tempfile.TemporaryDirectory(prefix="n83-orbit-gate-") as folder:
            root = Path(folder)
            binary, obj = root / "worker", root / "object"
            binary.write_bytes(b"\x7fELFsynthetic")
            binary.chmod(0o755)
            obj.write_bytes(b"synthetic object")
            actual_sha256 = gate.sha256_file

            def synthetic_hash(path):
                return gate.OBJECT_SHA256 if Path(path) == obj else actual_sha256(path)

            def image_inspect(command, **_kwargs):
                self.assertEqual(command[:3], ["docker", "image", "inspect"])
                return subprocess.CompletedProcess(command, 0, "sha256:" + "e" * 64 + " arm64 linux\n", "")

            for valid_import, expected_status in ((False, "PRODUCER_FAILURE_import_receipt"),
                                                  (True, "UNKNOWN_trial_cap")):
                with self.subTest(valid_import=valid_import):
                    output = root / ("valid" if valid_import else "invalid")
                    auth = root / ("auth-valid.json" if valid_import else "auth-invalid.json")
                    auth.write_text(json.dumps({
                        "schema": "n83.additional-local-budget/v1", "study": gate.STUDY,
                        "purpose": "structured_orbit_import_and_one_m5_trial",
                        "prior_pilot_charged_seconds": "7200",
                        "output_dir": str(output),
                        "build_dir": str(root),
                        "additional_authorized_seconds": "1800",
                        "user_approval_reference": "synthetic-test-only",
                    }))
                    (root / "build.json").write_text(json.dumps({
                        "schema": "n83.structured-orbit-linux-build/v1",
                        "status": "PASS_build_only", "source_commit": "a" * 40,
                        "binary_sha256": actual_sha256(binary),
                        "authorization_sha256": hashlib.sha256(auth.read_bytes()).hexdigest(),
                        "container_image_id": "sha256:" + "e" * 64,
                        "binary": str(binary), "conservative_charge_seconds": 600,
                        "container_cleanup_verified": True,
                        "source_clean": True, "source_same": True,
                        "solver_search_executed": False,
                        "exit_code": 0, "timed_out": False,
                        "launch_error": None, "cleanup_error": None,
                        "container_state": {"ExitCode": 0, "OOMKilled": False},
                        "elapsed_seconds_L0": 12.0,
                    }))
                    launched = []

                    def fake_stage(*, name, command, wall, output, container_name):
                        launched.append(name)
                        self.assertEqual(command[command.index("--network") + 1], "none")
                        self.assertEqual(command[command.index("--memory-swap") + 1], "4096m")
                        self.assertEqual(command[command.index("--cpus") + 1], "1")
                        self.assertEqual(wall, gate.VERIFY_WALL if name == "verify" else gate.SEARCH_WALL)
                        if name == "verify":
                            receipt = verification()
                            if not valid_import:
                                receipt["points"] -= 1
                            (output / "verify.stdout").write_text(json.dumps(receipt))
                        else:
                            config, result = search_records()
                            (output / "run").mkdir()
                            (output / "run" / "config.json").write_text(json.dumps(config))
                            (output / "run" / "result.json").write_text(json.dumps(result))
                        return {
                            "name": name, "exit_code": 0, "timed_out": False,
                            "cleanup_error": None, "launch_error": None,
                            "container_state": {"ExitCode": 0, "OOMKilled": False},
                        }

                    with patch.object(gate, "clean_commit", return_value="a" * 40), \
                         patch.object(gate, "sha256_file", side_effect=synthetic_hash), \
                         patch.object(gate.subprocess, "run", side_effect=image_inspect), \
                         patch.object(gate, "stage", side_effect=fake_stage):
                        outer = gate.run(obj, binary, auth, output)
                    self.assertEqual(outer["status"], expected_status)
                    self.assertEqual(launched, ["verify", "search"] if valid_import else ["verify"])
                    self.assertEqual(json.loads((output / "budget.json").read_text())
                                     ["conservative_charge_seconds"], "1800")


if __name__ == "__main__":
    unittest.main()
