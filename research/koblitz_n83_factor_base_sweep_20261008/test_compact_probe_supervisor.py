"""Synthetic Docker-process controls for the retained compact-probe guard."""

import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import compact_probe_supervisor as guard


def worker_receipt(config: dict, hit: bool = False) -> dict:
    status = "HIT" if hit else "MISS"
    row = None
    indices = None
    if hit:
        indices = [0, 2, 4, 6]
        row = {
            "schema": "n83.primary-wide-relation-row/v1",
            "modulus_decimal": guard.ORDER,
            "orbit_columns": config["columns"],
            "nonzero": [{"column": 0, "coefficient_decimal": "4"}],
            "dense_decimal_blake3": "a" * 64,
            "group_verified": True,
        }
    return {
        "schema": guard.WORKER_SCHEMA, "study": guard.STUDY,
        "status": status, "curve_a": 0, "fixture": 0,
        "orbit_columns": config["columns"], "policy": config["policy"],
        "seed": config["seed"], "max_candidate_states": config["max_candidate_states"],
        "memory_cgroup_limit_bytes": config["memory_cgroup_limit_bytes"],
        "memory_cgroup_swap_limit_bytes": 0,
        "memory_cgroup_peak_bytes": 12345,
        "source_commit": config["source_commit"],
        "source_attestation_mode": "supervised_snapshot",
        "source_exporter_blake3": "a" * 64,
        "source_adapter_blake3": "b" * 64,
        "source_compact_blake3": "c" * 64,
        "panel_manifest_blake3": "d" * 64,
        "point_set_blake3": "e" * 64,
        "public_corpus_canonical_json_blake3": "f" * 64,
        "base_import_ms": 1.0, "target_validation_ms": 2.0,
        "point_probe_ms": 3.0, "relation_row_ms": 4.0,
        "process_wall_ms": 10.0,
        "relation_row_stage_executed": hit, "relation_row": row,
        "rank_stage_executed": False, "column_log_verification": False,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
        "probe": {
            "schema": "compact-four-sum-point-probe/v1", "status": status,
            "n": 83, "curve_a": 0, "orbit_columns": config["columns"],
            "unordered_pairs": config["pair_mode"] == "unordered",
            "candidate_states": config["max_candidate_states"],
            "point_indices": indices, "group_verified": True if hit else None,
            "rank_stage_executed": False, "total_index_calculus_runtime_ms": None,
        },
    }


class CompactProbeSupervisorTests(unittest.TestCase):
    def test_receipt_gate_checks_full_width_row_and_configuration(self):
        config = {
            "columns": 64, "policy": "public_x_hash", "seed": guard.SEEDS[0],
            "pair_mode": "unordered", "max_candidate_states": guard.candidate_states(64, "unordered"),
            "memory_cgroup_limit_bytes": 256 * 1024 * 1024,
            "source_commit": "a" * 40,
        }
        for hit in (False, True):
            receipt = worker_receipt(config, hit)
            self.assertIsNone(guard.receipt_error(receipt, config))
            receipt["rank_stage_executed"] = True
            self.assertIsNotNone(guard.receipt_error(receipt, config))
        receipt = worker_receipt(config, True)
        receipt["relation_row"]["modulus_decimal"] = "7"
        self.assertIsNotNone(guard.receipt_error(receipt, config))
        receipt = worker_receipt(config, False)
        receipt["probe"]["candidate_states"] -= 1
        self.assertIsNotNone(guard.receipt_error(receipt, config))

    def test_outer_keeps_success_timeout_memory_and_malformed_outcomes(self):
        with tempfile.TemporaryDirectory(prefix="n83-compact-guard-") as temporary:
            root = Path(temporary)
            panel = root / "panel"
            panel.mkdir()
            binary = root / "worker"
            binary.write_bytes(b"\x7fELFsynthetic")
            binary.chmod(0o755)
            for mode, expected in (
                ("miss", "PASS_probe_miss"),
                ("hit", "PASS_probe_hit"),
                ("timeout", "UNKNOWN_wall_cap"),
                ("memory", "UNKNOWN_resource_or_worker_exit"),
                ("malformed", "PRODUCER_FAILURE_worker_receipt"),
                ("source_changed", "PRODUCER_FAILURE_frozen_input_changed"),
                ("config_changed", "PRODUCER_FAILURE_frozen_input_changed"),
                ("binary_changed", "PRODUCER_FAILURE_frozen_input_changed"),
            ):
                with self.subTest(mode=mode):
                    output = root / mode
                    commands = []

                    def fake_run(command, **_kwargs):
                        commands.append(command)
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
                        commands.append(command)
                        config = json.loads((output / "config.json").read_text())
                        if mode in ("miss", "hit"):
                            guard.write_new_json(output / "worker.json", worker_receipt(config, mode == "hit"))
                        elif mode == "malformed":
                            (output / "worker.json").write_text("[]\n")
                        elif mode == "source_changed":
                            (output / "source" / "attestation.json").write_text("changed\n")
                        elif mode == "config_changed":
                            (output / "config.json").write_text("changed\n")
                        elif mode == "binary_changed":
                            binary.write_bytes(b"\x7fELFchanged")
                        return FakeChild()

                    with patch.object(guard, "clean_commit", return_value="a" * 40), \
                            patch.object(guard.subprocess, "run", side_effect=fake_run), \
                            patch.object(guard.subprocess, "Popen", side_effect=fake_popen):
                        outer = guard.run(
                            panel, 64, "public_x_hash", guard.SEEDS[0], "unordered",
                            0.5, 256, binary, output, checkout=Path(__file__).resolve().parents[2],
                        )
                    self.assertEqual(outer["status"], expected)
                    self.assertEqual(json.loads((output / "outer.json").read_text()), outer)
                    self.assertTrue((output / "config.json").exists())
                    self.assertTrue((output / "stdout.log").exists())
                    self.assertTrue((output / "stderr.log").exists())
                    docker_run = next(command for command in commands if command[:2] == ["docker", "run"])
                    self.assertIn("--network", docker_run)
                    self.assertIn("none", docker_run)
                    self.assertIn("--memory-swap", docker_run)
                    self.assertIn("ICV1_FROZEN_SOURCE_DIR=/source", docker_run)
                    self.assertEqual(docker_run[-10:], [
                        "/worker", "primary-s3-probe", "/panel", "64", "public_x_hash",
                        str(guard.SEEDS[0]), "unordered", str(guard.candidate_states(64, "unordered")),
                        "256", "/out/worker.json",
                    ])
                    if mode == "timeout":
                        self.assertTrue(any(command[:2] == ["docker", "kill"] for command in commands))
                        self.assertTrue(any(command[:2] == ["docker", "ps"] for command in commands))


if __name__ == "__main__":
    unittest.main()
