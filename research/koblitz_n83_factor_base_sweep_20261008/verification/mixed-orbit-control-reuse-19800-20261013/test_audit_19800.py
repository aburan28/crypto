"""Synthetic controls for the independent v7 N83 receipt auditor."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import audit as checker


class IndependentAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = json.loads(checker.FIXTURE.read_text())

    def config(self):
        return {
            "schema": "n83.mixed-orbit-known-witness-config/v1",
            "scope": "known_answer_control_only", "study": checker.STUDY,
            "fixture_blake3": checker.FIXTURE_BLAKE3,
            "small_object_blake3": self.fixture["archive"]["small"]["s3_uri"].split("/")[-1][:-9],
            "envelope_object_blake3": self.fixture["archive"]["envelope"]["s3_uri"].split("/")[-1][:-9],
            "small_summands": 3, "residual_summands": 2,
            "budget_seconds": 840, "conflict_budget": 5000, "max_models": 32,
            "memory_cgroup_limit_bytes": 4294967296,
            "memory_cgroup_swap_limit_bytes": 0,
        }

    def capped_result(self):
        return {
            "schema": "n83.mixed-orbit-known-witness-result/v1",
            "scope": "known_answer_control_only", "status": "UNKNOWN_search_cap",
            "target": self.fixture["equation"]["target"], "candidate": None,
            "natural_relations": None, "total_index_calculus_runtime_ms": None,
            "selected_best_total_runtime": None,
            "base_import_elapsed_ms_L0": 10.0, "worker_wall_ms_L0": 100.0,
            "memory_cgroup_peak_bytes": 500000000,
            "search": {"conflicts": 5000, "models": 0, "invalid_models": 0,
                       "group_rejected_models": 0, "solver_calls": 1,
                       "patterns_attempted": 1, "domain_clauses": 465738,
                       "peak_variables": 101636, "peak_xors": 4980,
                       "peak_clauses": 500000, "exhausted": True,
                       "refuted": False},
        }

    def test_constructed_witness_replays_in_independent_group_arithmetic(self):
        points = [checker.parse_point(term["point"]) for term in self.fixture["terms"]]
        target = checker.parse_point(self.fixture["equation"]["target"])
        self.assertTrue(checker.candidate_sum_matches(points, target))
        self.assertFalse(checker.candidate_sum_matches(points[:-1], target))

    def test_conflict_cap_is_censored_and_tampering_is_rejected(self):
        config, result = self.config(), self.capped_result()
        missing = Path("/does-not-exist")
        self.assertEqual(checker.classify_worker(config, result, self.fixture,
                                                 missing, missing), "VALIDATED_CENSORED")
        result["search"]["conflicts"] = 5001
        with self.assertRaisesRegex(ValueError, "conflicts outside cap"):
            checker.classify_worker(config, result, self.fixture, missing, missing)

    def test_success_cannot_promote_without_archive_and_group_replay(self):
        config, result = self.config(), self.capped_result()
        result["status"] = "PASS_known_witness_solver_access"
        result["search"].update({"conflicts": 10, "models": 1, "exhausted": False})
        result["candidate"] = {
            "a": "0", "b": "1", "small_point_indices": [0, 166, 332],
            "residual_points": [term["point"] for term in self.fixture["terms"][3:]],
        }
        missing = Path("/does-not-exist")
        with self.assertRaisesRegex(ValueError, "archive/group replay"):
            checker.classify_worker(config, result, self.fixture, missing, missing)

    def test_stage_network_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            run = Path(folder)
            for suffix in ("stdout", "stderr"):
                (run / f"verify.{suffix}").write_bytes(b"")
            command = ["docker", "run", "--pull", "never", "--platform", "linux/arm64",
                       "--network", "bridge", "--cpus", "1", "--memory", "4096m",
                       "--memory-swap", "4096m", "--pids-limit", "64",
                       "--security-opt", "no-new-privileges", "--cap-drop", "ALL",
                       "--read-only"]
            command += [f"type=bind,src=/source,dst=/{item},readonly" for item in
                        ("worker", "small.jsonl.gz", "envelope.jsonl.gz",
                         "small-manifest.json", "envelope-manifest.json", "fixture.json")]
            command += [checker.IMAGE, "/worker", "verify-control",
                        "/small-manifest.json", "/small.jsonl.gz",
                        "/envelope-manifest.json", "/envelope.jsonl.gz",
                        "/fixture.json"]
            stage = {"name": "verify", "exit_code": 0, "timed_out": False,
                     "launch_error": None, "cleanup_error": None,
                     "cleanup_verified": True, "container_state": {
                         "ExitCode": 0, "OOMKilled": False,
                         "Running": False, "Status": "exited"},
                     "wall_elapsed_seconds_L0": 1.0, "command": command,
                     "stdout_sha256": checker.sha(run / "verify.stdout"),
                     "stderr_sha256": checker.sha(run / "verify.stderr")}
            with self.assertRaisesRegex(ValueError, "stage --network differs"):
                checker.check_stage(stage, "verify", run)
            command[command.index("--network") + 1] = "none"
            checker.check_stage(stage, "verify", run)

    def test_missing_raw_run_is_recorded_as_unvalidated(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            output = root / "audit.json"
            arguments = ["audit.py", "--run-dir", str(root / "missing-run"),
                         "--authorization", str(root / "missing-authorization"),
                         "--small-object", str(root / "missing-small"),
                         "--envelope-object", str(root / "missing-envelope"),
                         "--output", str(output)]
            with mock.patch("sys.argv", arguments):
                self.assertEqual(checker.main(), 1)
            row = json.loads(output.read_text())
            self.assertEqual(row["status"], "UNVALIDATED")
            self.assertEqual(row["raw_sha256"], {})

    def test_validated_cap_keeps_nonzero_process_disposition(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            output = root / "audit.json"
            arguments = ["audit.py", "--run-dir", str(root / "missing-run"),
                         "--authorization", str(root / "missing-authorization"),
                         "--small-object", str(root / "missing-small"),
                         "--envelope-object", str(root / "missing-envelope"),
                         "--output", str(output)]
            with (mock.patch("sys.argv", arguments),
                  mock.patch.object(checker, "audit", return_value={
                      "schema": "n83.known-witness-v7-independent-audit/v1",
                      "status": "VALIDATED_CENSORED"})):
                self.assertEqual(checker.main(), 1)
            self.assertEqual(json.loads(output.read_text())["status"],
                             "VALIDATED_CENSORED")


if __name__ == "__main__":
    unittest.main()
