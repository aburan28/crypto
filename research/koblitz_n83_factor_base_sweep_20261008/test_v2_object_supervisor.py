"""Synthetic guard tests; they do not construct an N83 factor base."""

import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import v2_object_supervisor as guard


def manifest(config: dict) -> dict:
    digest = "a" * 64
    object_name = f"objects/{digest}.jsonl.gz"
    return {
        "schema": "n83.factor-base-panel/v2-size-frontier",
        "study": guard.STUDY,
        "status": "completed_factor_base_object",
        "source_commit": config["source_commit"],
        "source_blake3": "b" * 64,
        "size_design_blake3": "c" * 64,
        "budget_seconds": config["worker_wall_seconds"],
        "expected_base_count": 1,
        "completed_base_count": 1,
        "selected_best_total_runtime": None,
        "bases": [{
            "a": config["curve_a"], "policy": config["policy"],
            "columns": config["columns"], "seed": config["seed"],
            "points": 166 * config["columns"],
            "runtime_rank_eligible": False,
            "total_index_calculus_runtime_ms": None,
            "relation_yield": None, "matrix_rank": None, "verified_dlp": None,
            "compressed_blake3": digest, "plain_blake3": "d" * 64,
            "point_set_blake3": "e" * 64,
            "object": object_name, "bytes": 3,
            "s3_uri": ("s3://crypto-autoresearcher/factor-bases/icv1/etc/"
                       f"{guard.STUDY}/v2-size-frontier/a{config['curve_a']}/{object_name}"),
        }],
    }


class V2ObjectSupervisorTests(unittest.TestCase):
    def test_manifest_rejects_wrong_selection_destination_and_missing_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = {
                "source_commit": "f" * 40, "worker_wall_seconds": 2,
                "curve_a": 0, "policy": "public_x_hash", "columns": 1182,
                "seed": guard.SEEDS[0],
            }
            row = manifest(config)
            self.assertIn("missing", guard.manifest_error(row, config, root))
            path = root / row["bases"][0]["object"]
            path.parent.mkdir()
            path.write_bytes(b"abc")
            self.assertIsNone(guard.manifest_error(row, config, root))
            row["bases"][0]["columns"] = 2048
            self.assertIn("columns", guard.manifest_error(row, config, root))
            row["bases"][0]["columns"] = 1182
            row["bases"][0]["s3_uri"] = "s3://wrong/prefix"
            self.assertIn("destination", guard.manifest_error(row, config, root))

    def test_docker_command_pins_memory_swap_network_image_and_paths(self):
        config = {
            "memory_mib": 512, "container_name": "fixed-test",
            "container_image_id": "sha256:" + "f" * 64,
            "curve_a": 0, "policy": "public_x_hash", "columns": 1182,
            "seed": guard.SEEDS[0], "worker_wall_seconds": 5,
        }
        command = guard.docker_command(
            config, Path("/repo/worktrees/case"), Path("/repo"),
            Path("/worker-binary"), Path("/out-dir"),
        )
        self.assertEqual(command[command.index("--memory") + 1], "512m")
        self.assertEqual(command[command.index("--memory-swap") + 1], "512m")
        self.assertEqual(command[command.index("--network") + 1], "none")
        self.assertIn("--read-only", command)
        self.assertIn("--pull", command)
        self.assertIn("never", command)
        self.assertEqual(command[-8:], [
            "/worker", "v2-construct-one", "/out/object", "0",
            "public_x_hash", "1182", str(guard.SEEDS[0]), "5",
        ])

    def test_outer_retains_success_timeout_resource_and_malformed_receipts(self):
        with tempfile.TemporaryDirectory(prefix="n83-v2-guard-") as temporary:
            root = Path(temporary)
            repository = root / "repo"
            checkout = repository / "worktrees" / "case"
            git_common = repository / ".git"
            checkout.mkdir(parents=True)
            git_common.mkdir()
            binary = root / "worker"
            binary.write_bytes(b"\x7fELFsynthetic")
            binary.chmod(0o755)
            expected = {
                "success": "PASS_construction_only",
                "timeout": "UNKNOWN_wall_cap",
                "resource": "UNKNOWN_resource_or_worker_exit",
                "malformed": "PRODUCER_FAILURE_manifest",
                "source_changed": "PRODUCER_FAILURE_frozen_input_changed",
            }
            for mode, status in expected.items():
                with self.subTest(mode=mode):
                    output = root / mode

                    def fake_run(command, **_kwargs):
                        if "--git-common-dir" in command:
                            return subprocess.CompletedProcess(command, 0, str(git_common) + "\n", "")
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
                        if mode not in ("timeout", "resource"):
                            value = manifest(config)
                            if mode == "malformed":
                                value["bases"][0]["points"] -= 1
                            object_dir = output / "object"
                            object_path = object_dir / value["bases"][0]["object"]
                            object_path.parent.mkdir(parents=True)
                            object_path.write_bytes(b"abc")
                            (object_dir / "manifest.json").write_text(json.dumps(value))
                        return FakeChild()

                    commits = ["f" * 40, "0" * 40] if mode == "source_changed" else ["f" * 40] * 2
                    with patch.object(guard, "clean_commit", side_effect=commits), \
                         patch.object(guard.subprocess, "run", side_effect=fake_run), \
                         patch.object(guard.subprocess, "Popen", side_effect=fake_popen):
                        receipt = guard.run(checkout, repository, binary, output,
                                            0, "public_x_hash", 1182, guard.SEEDS[0], 1, 512)
                    self.assertEqual(receipt["status"], status)
                    self.assertEqual(json.loads((output / "outer.json").read_text())["status"], status)
                    self.assertFalse(receipt["replay_executed"])
                    self.assertFalse(receipt["s3_upload_executed"])


if __name__ == "__main__":
    unittest.main()
