"""Synthetic process guards for retained N83 chained-S3 search."""

import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import chain_search_supervisor as guard


def receipts(config: dict, status: str = "PREFLIGHT_ONLY") -> tuple[dict, dict]:
    limits = {
        "max_variables": config["max_variables"],
        "max_domain_clauses": config["max_domain_clauses"],
        "max_models": config["max_models"],
        "conflict_budget": config["conflict_budget"],
    }
    trials = 0 if status == "PREFLIGHT_ONLY" else 1
    calls = 0 if status == "PREFLIGHT_ONLY" else 1
    verified = status == "PASS_verified_target_only"
    unknown = status == "UNKNOWN_solver_cap"
    report = {key: 0 for key in guard.COUNTERS}
    report.update({
        "factor_base_size": 166 * config["columns"],
        "orbit_count": config["columns"], "trials": trials,
        "sat_calls": calls, "sat_unknowns": int(unknown),
        "relations": int(verified), "independent_relations": int(verified),
        "linear_solve_attempts": int(verified),
        "relation_collection_ns": "1000", "linear_algebra_ns": "0",
    })
    worker_config = {
        "schema": guard.WORKER_CONFIG_SCHEMA, "study": guard.STUDY,
        "curve_a": 0, "fixture": 0, "orbit_columns": config["columns"],
        "summands": config["summands"], "strategy": "chain-s3",
        "policy": config["policy"], "seed": config["seed"],
        "max_trials": config["max_trials"],
        "budget_seconds": config["worker_wall_seconds"],
        "allow_direct_relation": False,
        "source_commit": config["source_commit"],
        "source_attestation_mode": "supervised_snapshot",
        "fixture_dir": "/fixtures",
        "chain_s3_limits": limits,
        "memory_cgroup_limit_bytes": config["memory_cgroup_limit_bytes"],
        "memory_cgroup_swap_limit_bytes": 0,
        "source_adapter_blake3": "a" * 64,
    }
    summary = {
        **worker_config, "schema": guard.WORKER_SCHEMA,
        "status": status, "object": "objects/" + "b" * 64 + ".jsonl.gz",
        "point_set_blake3": "c" * 64,
        "panel_manifest_blake3": "d" * 64,
        "public_corpus_canonical_json_blake3": "e" * 64,
        "source_chain_blake3": "f" * 64,
        "source_index_calculus_blake3": "1" * 64,
        "memory_cgroup_peak_bytes": 12345,
        "base_import_ms": 1.0, "target_validation_ms": 2.0,
        "solver_ms": 3.0, "post_solver_validation_ms": 0.0,
        "pre_summary_process_wall_ms": 6.0,
        "column_log_verification": False,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
        "solver_stage_executed": calls > 0,
        "report": report,
        "verified_log": "3" if verified else None,
    }
    return worker_config, summary


class ChainSearchGuardTests(unittest.TestCase):
    def setUp(self):
        self.config = {
            "columns": 64, "policy": "public_x_hash", "seed": guard.SEEDS[0],
            "summands": 5, "max_trials": 0, "max_variables": 100000,
            "max_domain_clauses": 300000, "max_models": 1,
            "conflict_budget": 1000, "worker_wall_seconds": 1,
            "memory_cgroup_limit_bytes": 256 * 1024 * 1024,
            "source_commit": "a" * 40,
        }

    def test_receipt_checks_rank_solver_caps_and_false_completion(self):
        worker_config, worker = receipts(self.config)
        self.assertIsNone(guard.receipt_error(worker, worker_config, self.config))
        worker["report"]["relations"] = 1
        worker["report"]["independent_relations"] = 1
        self.assertIn("preflight", guard.receipt_error(worker, worker_config, self.config))
        config = dict(self.config, max_trials=1)
        worker_config, worker = receipts(config, "UNKNOWN_solver_cap")
        self.assertIsNone(guard.receipt_error(worker, worker_config, config))
        worker["report"]["sat_unknowns"] = 0
        self.assertIn("solver-cap", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        self.assertIsNone(guard.receipt_error(worker, worker_config, config))
        worker["report"]["independent_relations"] = 0
        self.assertIsNotNone(guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker["total_index_calculus_runtime_ms"] = 6.0
        self.assertIn("unmeasured", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker_config["source_commit"] = "0" * 40
        self.assertIn("source_commit", guard.receipt_error(worker, worker_config, config))

    def test_outer_retains_preflight_timeout_resource_and_bad_receipts(self):
        with tempfile.TemporaryDirectory(prefix="n83-chain-search-guard-") as temporary:
            root = Path(temporary)
            panel = root / "panel"
            panel.mkdir()
            (panel / "manifest.json").write_text("{}")
            (panel / "replay.json").write_text("{}")
            (panel / "upload-receipt.json").write_text("{}")
            fixtures = root / "fixtures"
            fixtures.mkdir()
            (fixtures / "probe-corpus.json").write_text("{}")
            (fixtures / "probe-validation.json").write_text("{}")
            binary = root / "worker"
            binary.write_bytes(b"\x7fELFsynthetic")
            binary.chmod(0o755)
            checkout = Path(guard.__file__).resolve().parents[2]
            expected = {
                "preflight": "PREFLIGHT_ONLY",
                "timeout": "UNKNOWN_wall_cap",
                "resource": "UNKNOWN_resource_or_worker_exit",
                "malformed": "PRODUCER_FAILURE_worker_receipt",
                "source_changed": "PRODUCER_FAILURE_frozen_input_changed",
                "input_changed": "PRODUCER_FAILURE_frozen_input_changed",
            }
            for mode, status in expected.items():
                with self.subTest(mode=mode):
                    output = root / mode

                    def fake_run(command, **_kwargs):
                        if command[:3] == ["docker", "image", "inspect"]:
                            return subprocess.CompletedProcess(command, 0, "sha256:" + "f" * 64 + "\n", "")
                        return subprocess.CompletedProcess(command, 0, "", "")

                    class FakeChild:
                        calls = 0

                        def wait(self, timeout=None):
                            self.calls += 1
                            if mode == "timeout" and self.calls == 1:
                                raise subprocess.TimeoutExpired("docker run", timeout)
                            return 137 if mode == "resource" else 0

                        def kill(self):
                            pass

                    def fake_popen(command, **_kwargs):
                        config = json.loads((output / "config.json").read_text())
                        self.assertIn("--network", command)
                        self.assertEqual(command[command.index("--network") + 1], "none")
                        self.assertEqual(command[command.index("--memory-swap") + 1], "256m")
                        self.assertEqual(command[-16:-14], ["/worker", "primary-chain-cold"])
                        if mode not in ("timeout", "resource"):
                            worker_config, worker = receipts(config)
                            if mode == "malformed":
                                worker["report"]["sat_invalid_models"] = 1
                            run_dir = output / "run"
                            run_dir.mkdir()
                            (run_dir / "config.json").write_text(json.dumps(worker_config))
                            (run_dir / "summary.json").write_text(json.dumps(worker))
                        if mode == "source_changed":
                            (output / "source" / guard.SOURCES[2]).write_bytes(b"changed")
                        if mode == "input_changed":
                            (fixtures / "probe-corpus.json").write_text('{"changed":true}')
                        return FakeChild()

                    with patch.object(guard, "clean_commit", return_value="a" * 40), \
                         patch.object(guard.subprocess, "run", side_effect=fake_run), \
                         patch.object(guard.subprocess, "Popen", side_effect=fake_popen):
                        outer = guard.run(
                            panel, fixtures, 64, "public_x_hash", guard.SEEDS[0], 5,
                            0, 100000, 300000, 1, 1000, 1.0, 256,
                            binary, output, checkout=checkout,
                        )
                    self.assertEqual(outer["status"], status)
                    self.assertEqual(json.loads((output / "outer.json").read_text())["status"], status)


if __name__ == "__main__":
    unittest.main()
