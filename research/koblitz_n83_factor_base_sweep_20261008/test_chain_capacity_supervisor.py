"""Synthetic process guards for exact retained-base chained-S3 construction."""

import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import chain_capacity_supervisor as guard


def worker_receipt(config: dict, status: str = "PASS_model_construction_only") -> dict:
    built = status == "PASS_model_construction_only"
    required_variables = 84743 if config["summands"] == 5 else 105908
    required_clauses = 200000
    return {
        "schema": guard.WORKER_SCHEMA, "study": guard.STUDY,
        "status": status, "curve_a": 0, "fixture": 0,
        "orbit_columns": config["columns"], "policy": config["policy"],
        "seed": config["seed"], "summands": config["summands"],
        "identity_mask": config["identity_mask"],
        "max_variables": config["max_variables"],
        "max_domain_clauses": config["max_domain_clauses"],
        "memory_cgroup_limit_bytes": config["memory_cgroup_limit_bytes"],
        "memory_cgroup_swap_limit_bytes": 0, "memory_cgroup_peak_bytes": 12345,
        "source_commit": config["source_commit"],
        "source_attestation_mode": "supervised_snapshot",
        "source_exporter_blake3": "a" * 64,
        "source_adapter_blake3": "b" * 64,
        "source_chain_blake3": "c" * 64,
        "source_index_calculus_blake3": "d" * 64,
        "panel_manifest_blake3": "e" * 64,
        "point_set_blake3": "f" * 64,
        "public_corpus_canonical_json_blake3": "1" * 64,
        "object": "factor-bases/icv1/sample.jsonl",
        "required_max_variables": required_variables,
        "required_domain_clauses": required_clauses,
        "legal_x_coordinates": 83 * config["columns"],
        "sat_variables": 84743 if built else None,
        "sat_clauses": 250000 if built else None,
        "sat_xor_rows": 1000 if built else None,
        "sat_and_gates": 82668 if built else None,
        "sat_s3_nodes": 4 if built else None,
        "base_import_ms": 1.0, "target_validation_ms": 2.0,
        "model_construction_ms": 3.0, "process_wall_ms": 6.0,
        "solver_search_executed": False,
        "relation_stage_executed": False,
        "rank_stage_executed": False,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }


class ChainCapacityGuardTests(unittest.TestCase):
    def test_v2_sizes_match_frozen_design(self):
        design = json.loads((Path(guard.__file__).parent / "size-frontier-v2.json").read_text())
        self.assertEqual(list(guard.V2_SIZES), design["new_sizes"])

    def setUp(self):
        self.config = {
            "columns": 64, "policy": "public_x_hash", "seed": guard.SEEDS[0],
            "summands": 5, "identity_mask": 0,
            "max_variables": 100000, "max_domain_clauses": 300000,
            "memory_cgroup_limit_bytes": 256 * 1024 * 1024,
            "source_commit": "a" * 40,
        }

    def test_receipt_rejects_false_stage_and_invalid_caps(self):
        self.assertIsNone(guard.receipt_error(worker_receipt(self.config), self.config))
        v2_config = dict(self.config, columns=1182)
        self.assertIsNone(guard.receipt_error(worker_receipt(v2_config), v2_config))
        receipt = worker_receipt(self.config)
        receipt["solver_search_executed"] = True
        self.assertIsNotNone(guard.receipt_error(receipt, self.config))
        receipt = worker_receipt(self.config)
        receipt["identity_mask"] = 1
        self.assertIsNotNone(guard.receipt_error(receipt, self.config))
        receipt = worker_receipt(self.config)
        receipt["sat_clauses"] = 1
        self.assertIsNotNone(guard.receipt_error(receipt, self.config))
        for status in ("UNKNOWN_variable_cap", "UNKNOWN_domain_clause_cap"):
            config = dict(self.config)
            if status == "UNKNOWN_variable_cap":
                config["max_variables"] = 80000
            else:
                config["max_domain_clauses"] = 100000
            receipt = worker_receipt(config, status)
            self.assertIsNone(guard.receipt_error(receipt, config))
            receipt["sat_variables"] = 5
            self.assertIsNotNone(guard.receipt_error(receipt, config))

    def test_outer_retains_success_caps_timeout_memory_and_mutation(self):
        with tempfile.TemporaryDirectory(prefix="n83-chain-guard-") as temporary:
            root = Path(temporary)
            panel = root / "panel"
            panel.mkdir()
            binary = root / "worker"
            binary.write_bytes(b"\x7fELFsynthetic")
            binary.chmod(0o755)
            checkout = Path(guard.__file__).resolve().parents[2]
            modes = (
                ("success", "PASS_model_construction_only"),
                ("v2_success", "PASS_model_construction_only"),
                ("variable_cap", "UNKNOWN_variable_cap"),
                ("domain_cap", "UNKNOWN_domain_clause_cap"),
                ("timeout", "UNKNOWN_wall_cap"),
                ("memory", "UNKNOWN_resource_or_worker_exit"),
                ("malformed", "PRODUCER_FAILURE_worker_receipt"),
                ("source_changed", "PRODUCER_FAILURE_frozen_input_changed"),
            )
            for mode, expected in modes:
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
                            return 137 if mode == "memory" else 0

                        def kill(self):
                            pass

                    def fake_popen(command, **_kwargs):
                        config = json.loads((output / "config.json").read_text())
                        if mode not in ("timeout", "memory"):
                            status = {
                                "variable_cap": "UNKNOWN_variable_cap",
                                "domain_cap": "UNKNOWN_domain_clause_cap",
                            }.get(mode, "PASS_model_construction_only")
                            receipt = worker_receipt(config, status)
                            if mode == "malformed":
                                receipt["rank_stage_executed"] = True
                            (output / "worker.json").write_text(json.dumps(receipt))
                        if mode == "source_changed":
                            (output / "source" / guard.SOURCES[2]).write_bytes(b"mutated")
                        return FakeChild()

                    with patch.object(guard, "clean_commit", return_value="a" * 40), \
                         patch.object(guard.subprocess, "run", side_effect=fake_run), \
                         patch.object(guard.subprocess, "Popen", side_effect=fake_popen):
                        outer = guard.run(
                            panel, 1182 if mode == "v2_success" else 64,
                            "public_x_hash", guard.SEEDS[0], 5, 0,
                            80000 if mode == "variable_cap" else 100000,
                            100000 if mode == "domain_cap" else 300000,
                            1.0, 256, binary, output,
                            checkout=checkout,
                        )
                    self.assertEqual(outer["status"], expected)
                    if mode == "v2_success":
                        self.assertEqual(json.loads((output / "config.json").read_text())["columns"], 1182)
                    self.assertEqual(json.loads((output / "outer.json").read_text())["status"], expected)


if __name__ == "__main__":
    unittest.main()
