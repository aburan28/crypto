#!/usr/bin/env python3
"""Offline and synthetic tests for the Stage-22 licensed Magma replay."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import textwrap
import unittest

import run_koblitz_phase_b_magma_replay as replay


def executable(path: Path, text: str) -> Path:
    path.write_text(textwrap.dedent(text).lstrip())
    path.chmod(0o755)
    return path


def fake_tools(root: Path) -> tuple[Path, Path, Path]:
    magma = executable(
        root / "magma",
        r'''
        #!/usr/bin/env python3
        import pathlib
        import re
        import sys
        if "--version" in sys.argv:
            print("Magma V2.fake")
            raise SystemExit(0)
        script = pathlib.Path(sys.argv[-1]).read_text()
        for key, value in re.findall(r'printf "(KOBLITZ_REPLAY_[A-Z0-9_]+)=([^"\\]+)\\n";', script):
            print(key + "=" + value)
        print("KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1")
        print("KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse")
        print("KOBLITZ_MAGMA_STATUS=SAT")
        print("KOBLITZ_MAGMA_F4_DEGREES=[ 2 ]")
        print("KOBLITZ_MAGMA_BASIS_SIZE=2")
        print("KOBLITZ_MAGMA_CPU_SECONDS=0.001")
        print("KOBLITZ_MAGMA_WALL_SECONDS=0.002")
        if "KOBLITZ_MODEL_SCHEMA" in script:
            variables = int(re.search(r"KOBLITZ_MODEL_VARIABLES=(\d+)", script).group(1))
            print("KOBLITZ_MODEL_SCHEMA=koblitz_magma_f4_model.v1")
            print("KOBLITZ_MODEL_STATUS=SAT")
            print("KOBLITZ_MODEL_BITS=" + "0" * variables)
            print("KOBLITZ_MODEL_VARIABLES=" + str(variables))
            print("KOBLITZ_MODEL_EXCLUDED_COUNT=0")
            print("KOBLITZ_MODEL_CPU_SECONDS=0.001")
            print("KOBLITZ_MODEL_WALL_SECONDS=0.002")
        ''',
    )
    minisat = executable(
        root / "minisat",
        """
        #!/bin/sh
        echo 'minisat synthetic 1'
        """,
    )
    backend = executable(
        root / "backend",
        r'''
        #!/usr/bin/env python3
        import json
        import pathlib
        import sys
        if "--version" in sys.argv:
            print("synthetic backend 1")
            raise SystemExit(0)
        manifest = json.loads(pathlib.Path(sys.argv[2]).read_text())
        assignment = json.loads(pathlib.Path(sys.argv[3]).read_text())
        print(json.dumps({
            "schema": "koblitz_pdp_assignment_validation.v1",
            "status": "valid_point_witness",
            "source_instance_id": manifest["source_instance"]["id_blake3"],
            "source_instance_verified": True,
            "regenerated_source_exact": True,
            "assignment_values": len(assignment),
            "assignment_blake3": "a" * 64,
            "source_assignment": assignment,
            "source_model_valid": True,
            "source_witness_valid": True
        }))
        ''',
    )
    return magma, minisat, backend


def tiny_magma() -> str:
    return (
        'SetNthreads(1);\nSetGPU(false);\nSetSeed(1);\n'
        'R := BooleanPolynomialRing(1, "grevlex");\nX := [R.1];\n'
        'F := [X[1]];\nI := ideal<R | F>;\n' + replay.legacy.F4_SUFFIX
    )


def synthetic_packet(root: Path) -> Path:
    packet = root / "packet"
    packet.mkdir()
    protocol_bytes = replay.regular_bytes(replay.DEFAULT_PROTOCOL, "protocol")
    replay.write_new(packet / "protocol.json", protocol_bytes)
    blind, _ = replay.read_json(
        replay.DEFAULT_BUNDLE / "inputs/blind-bundle.json", "blind bundle"
    )
    tasks = []
    for ordinal, source in enumerate(blind["instances"]):
        task_id = f"{ordinal:06d}-{source['blind_instance_id']}"
        directory = packet / "instances" / task_id
        directory.mkdir(parents=True)
        magma = directory / "instance.magma"
        manifest = directory / "manifest.json"
        replay.write_new(magma, tiny_magma().encode())
        replay.write_json_new(manifest, {
            "blind_instance_id": source["blind_instance_id"],
            "source_variables": 1,
            "source_instance": {"id_blake3": "a" * 64},
        })
        tasks.append({
            "ordinal": ordinal,
            "id": task_id,
            "blind_instance_id": source["blind_instance_id"],
            "source_system_id": source["source_system_id"],
            "source_nonce": source["source_nonce"],
            "source_instance_blake3": "a" * 64,
            "source_variables": 1,
            "cell_id": source["cell_id"],
            "n": source["n"], "ell": source["ell"], "m": source["m"],
            "basis": source["basis"], "curve_a": source["curve_a"],
            "factor_index": source["factor_index"], "target": source["target"],
            "files": {
                "magma_boolean_f4": replay.file_record(magma, packet, "Magma input"),
                "manifest": replay.file_record(manifest, packet, "manifest"),
            },
        })
    protocol, _ = replay.read_json(replay.DEFAULT_PROTOCOL, "protocol")
    manifest = replay.with_self_hash({
        "schema": replay.PACKET_MANIFEST_SCHEMA,
        "status": "solver_inputs_frozen",
        "created_at": replay.now(),
        "protocol_path": "protocol.json",
        "protocol_sha256": replay.sha256_file(packet / "protocol.json"),
        "source_terminal_bundle_seal_sha256": protocol["source"]["terminal_bundle_seal_sha256"],
        "source_terminal_run_seal_sha256": protocol["source"]["terminal_run_seal_sha256"],
        "source_terminal_run_inventory_sha256": protocol["source"]["terminal_run_inventory_sha256"],
        "source_terminal_score_seal_sha256": protocol["source"]["terminal_score_seal_sha256"],
        "order": "exact blind-bundle order", "instance_count": 160,
        "tasks": tasks, "ground_truth_included": False,
        "execution_policy": protocol["execution"], "claim_boundary": protocol["claim_boundary"],
    }, "manifest_payload_sha256")
    replay.write_json_new(packet / "packet-manifest.json", manifest)
    files = replay.inventory(packet, {"packet-seal.json"})
    seal = replay.with_self_hash({
        "schema": replay.PACKET_SEAL_SCHEMA,
        "status": "external_solver_packet_frozen",
        "manifest_path": "packet-manifest.json",
        "manifest_sha256": replay.sha256_file(packet / "packet-manifest.json"),
        "instance_count": 160, "magma_input_count": 160,
        "inventory": files, "inventory_sha256": replay.canonical_sha256(files),
        "claim_boundary": protocol["claim_boundary"],
    }, "seal_payload_sha256")
    replay.write_json_new(packet / "packet-seal.json", seal)
    return packet


def reseal_return(root: Path) -> None:
    seal_path = root / "return-seal.json"
    seal = json.loads(seal_path.read_text())
    seal["preflight_sha256"] = replay.sha256_file(root / "preflight.json")
    files = replay.inventory(root, {"return-seal.json"})
    seal["inventory"] = files
    seal["inventory_sha256"] = replay.canonical_sha256(files)
    seal.pop("seal_payload_sha256", None)
    seal["seal_payload_sha256"] = replay.canonical_sha256(seal)
    seal_path.write_bytes(replay.pretty_bytes(seal))


def refresh_task_receipt(root: Path, task: dict, field: str, receipt_path: Path) -> None:
    receipt = json.loads(receipt_path.read_text())
    reference = replay.file_record(receipt_path, root, "test receipt")
    reference["process"] = receipt["process"]
    task[field]["receipt"] = reference


def rehash_start(path: Path, command: list[str]) -> None:
    start = json.loads(path.read_text())
    start["command"] = command
    start.pop("start_payload_sha256")
    start["start_payload_sha256"] = replay.canonical_sha256(start)
    path.write_bytes(replay.pretty_bytes(start))


class MagmaReplayTests(unittest.TestCase):
    def test_committed_dry_run_binds_all_160_exports(self) -> None:
        result = replay.dry_run(replay.DEFAULT_PROTOCOL, replay.DEFAULT_BUNDLE)
        self.assertEqual(result["blind_instance_count"], 160)
        self.assertEqual(result["sealed_magma_export_count"], 160)
        self.assertFalse(result["licensed_magma_executed"])

    def test_synthetic_one_task_run_is_sealed_then_scored(self) -> None:
        with tempfile.TemporaryDirectory(prefix="magma-replay-") as directory:
            root = Path(directory)
            packet = synthetic_packet(root)
            tools = root / "tools"
            tools.mkdir()
            magma, minisat, backend = fake_tools(tools)
            returned = root / "return"
            seal = replay.run(
                packet, returned, magma, minisat, backend, replay.DEFAULT_METER,
                "Synthetic licensed-host statement for tests only",
                synthetic_test_mode=True, max_tasks=1,
            )
            self.assertEqual(seal["status"], "external_results_frozen")
            checked = replay.validate_return(packet, returned)
            self.assertEqual(len(checked["results"]), 1)
            self.assertEqual(checked["results"][0]["final_status"], "sat")
            self.assertEqual(checked["results"][0]["attempts"], {
                "primary_f4": 1, "f4_plus_sat_g": 1, "point_validation": 1,
            })
            scored = root / "score"
            score_seal = replay.score(
                packet, returned, replay.DEFAULT_BUNDLE / "score/score.json",
                replay.DEFAULT_BUNDLE / "score/score-seal.json", scored,
            )
            self.assertEqual(score_seal["status"], "post_return_score_frozen")
            score, _ = replay.read_json(scored / "score.json", "score")
            self.assertEqual(score["classifications"], {"true_positive": 1})
            self.assertFalse(score["independent_external_reproduction_satisfied"])

    def test_packet_and_return_are_write_once_and_tamper_evident(self) -> None:
        with tempfile.TemporaryDirectory(prefix="magma-replay-") as directory:
            root = Path(directory)
            packet = synthetic_packet(root)
            with self.assertRaisesRegex(replay.ReplayError, "inventory"):
                (packet / "unexpected").write_text("tamper")
                replay.validate_packet(packet)
            (packet / "unexpected").unlink()
            tools = root / "tools"
            tools.mkdir()
            magma, minisat, backend = fake_tools(tools)
            returned = root / "return"
            replay.run(
                packet, returned, magma, minisat, backend, replay.DEFAULT_METER,
                "Synthetic licensed-host statement for tests only",
                synthetic_test_mode=True, max_tasks=1,
            )
            with self.assertRaisesRegex(replay.ReplayError, "new; retries and resume are forbidden"):
                replay.run(
                    packet, returned, magma, minisat, backend, replay.DEFAULT_METER,
                    "Synthetic licensed-host statement for tests only",
                    synthetic_test_mode=True, max_tasks=1,
                )
            summary = returned / "summary.json"
            summary.write_bytes(summary.read_bytes() + b" ")
            with self.assertRaisesRegex(replay.ReplayError, "inventory"):
                replay.validate_return(packet, returned)

    def test_resealed_status_selection_and_preflight_forgeries_fail(self) -> None:
        with tempfile.TemporaryDirectory(prefix="magma-replay-") as directory:
            root = Path(directory)
            packet = synthetic_packet(root)
            tools = root / "tools"
            tools.mkdir()
            magma, minisat, backend = fake_tools(tools)

            returned = root / "status-return"
            replay.run(
                packet, returned, magma, minisat, backend, replay.DEFAULT_METER,
                "Synthetic licensed-host statement for tests only",
                synthetic_test_mode=True, max_tasks=1,
            )
            task_id = json.loads((returned / "packet-binding.json").read_text())["selected_task_ids"][0]
            task_path = returned / "tasks" / task_id / "task-result.json"
            task = json.loads(task_path.read_text())
            task["final_status"] = "unsat"
            task_path.write_bytes(replay.pretty_bytes(task))
            reseal_return(returned)
            with self.assertRaisesRegex(replay.ReplayError, "final status"):
                replay.validate_return(packet, returned)

            returned = root / "selection-return"
            replay.run(
                packet, returned, magma, minisat, backend, replay.DEFAULT_METER,
                "Synthetic licensed-host statement for tests only",
                synthetic_test_mode=True, max_tasks=1,
            )
            binding_path = returned / "packet-binding.json"
            binding = json.loads(binding_path.read_text())
            binding["selected_task_ids"] *= 2
            binding.pop("binding_payload_sha256")
            binding["binding_payload_sha256"] = replay.canonical_sha256(binding)
            binding_path.write_bytes(replay.pretty_bytes(binding))
            seal_path = returned / "return-seal.json"
            seal = json.loads(seal_path.read_text())
            seal["selected_tasks"] = 2
            seal["terminal_task_records"] = 2
            seal_path.write_bytes(replay.pretty_bytes(seal))
            reseal_return(returned)
            with self.assertRaisesRegex(replay.ReplayError, "unique blind-order prefix"):
                replay.validate_return(packet, returned)

            returned = root / "preflight-return"
            replay.run(
                packet, returned, magma, minisat, backend, replay.DEFAULT_METER,
                "Synthetic licensed-host statement for tests only",
                synthetic_test_mode=True, max_tasks=1,
            )
            preflight_path = returned / "preflight.json"
            preflight = json.loads(preflight_path.read_text())
            preflight["synthetic_test_mode"] = False
            preflight["status"] = "licensed_host_ready"
            preflight["tools"]["magma"]["recognized_version_banner"] = True
            preflight["source"]["state"] = {
                "commit": preflight["source"]["state"]["commit"],
                "dirty": False,
                "porcelain": [],
            }
            preflight.pop("preflight_payload_sha256")
            preflight["preflight_payload_sha256"] = replay.canonical_sha256(preflight)
            preflight_path.write_bytes(replay.pretty_bytes(preflight))
            reseal_return(returned)
            with self.assertRaisesRegex(replay.ReplayError, "recognized Magma"):
                replay.validate_return(packet, returned)

    def test_resealed_command_assignment_and_metadata_swaps_fail(self) -> None:
        with tempfile.TemporaryDirectory(prefix="magma-replay-") as directory:
            root = Path(directory)
            packet = synthetic_packet(root)
            tools = root / "tools"
            tools.mkdir()
            magma, minisat, backend = fake_tools(tools)

            returned = root / "command-return"
            replay.run(
                packet, returned, magma, minisat, backend, replay.DEFAULT_METER,
                "Synthetic licensed-host statement for tests only",
                synthetic_test_mode=True, max_tasks=1,
            )
            task_id = json.loads((returned / "packet-binding.json").read_text())["selected_task_ids"][0]
            task_path = returned / "tasks" / task_id / "task-result.json"
            task = json.loads(task_path.read_text())
            receipt_path = returned / "tasks" / task_id / "primary-f4/receipt.json"
            metrics_path = receipt_path.parent / "metrics.json"
            receipt = json.loads(receipt_path.read_text())
            metrics = json.loads(metrics_path.read_text())
            forged_command = list(metrics["command"])
            forged_command[-1] = str((returned / "tasks" / task_id / "f4-plus-sat-g.magma").resolve())
            metrics["command"] = forged_command
            metrics_path.write_bytes(replay.pretty_bytes(metrics))
            receipt["process"] = metrics
            receipt["meter_launcher_command"][-len(forged_command):] = forged_command
            receipt_path.write_bytes(replay.pretty_bytes(receipt))
            rehash_start(returned / "tasks" / task_id / "primary-f4-start.json", forged_command)
            refresh_task_receipt(returned, task, "primary_f4", receipt_path)
            task_path.write_bytes(replay.pretty_bytes(task))
            reseal_return(returned)
            with self.assertRaisesRegex(replay.ReplayError, "different script"):
                replay.validate_return(packet, returned)

            returned = root / "assignment-return"
            replay.run(
                packet, returned, magma, minisat, backend, replay.DEFAULT_METER,
                "Synthetic licensed-host statement for tests only",
                synthetic_test_mode=True, max_tasks=1,
            )
            task_id = json.loads((returned / "packet-binding.json").read_text())["selected_task_ids"][0]
            task_path = returned / "tasks" / task_id / "task-result.json"
            task = json.loads(task_path.read_text())
            assignment_path = returned / "tasks" / task_id / "assignment.json"
            assignment_path.write_text("[true]\n")
            task["point_validation"]["assignment"] = replay.file_record(
                assignment_path, returned, "forged assignment"
            )
            receipt_path = returned / "tasks" / task_id / "point-validation/receipt.json"
            receipt = json.loads(receipt_path.read_text())
            stdout_path = receipt_path.parent / "stdout"
            report = json.loads(stdout_path.read_text())
            report["source_assignment"] = [True]
            stdout_path.write_text(json.dumps(report) + "\n")
            receipt["stdout"] = {
                "path": receipt["stdout"]["path"],
                "bytes": stdout_path.stat().st_size,
                "sha256": replay.sha256_file(stdout_path),
            }
            receipt_path.write_bytes(replay.pretty_bytes(receipt))
            task["point_validation"]["report"] = report
            refresh_task_receipt(returned, task, "point_validation", receipt_path)
            task_path.write_bytes(replay.pretty_bytes(task))
            reseal_return(returned)
            with self.assertRaisesRegex(replay.ReplayError, r"differs from SAT\(G\)"):
                replay.validate_return(packet, returned)

            returned = root / "metadata-return"
            replay.run(
                packet, returned, magma, minisat, backend, replay.DEFAULT_METER,
                "Synthetic licensed-host statement for tests only",
                synthetic_test_mode=True, max_tasks=1,
            )
            task_id = json.loads((returned / "packet-binding.json").read_text())["selected_task_ids"][0]
            task_path = returned / "tasks" / task_id / "task-result.json"
            task = json.loads(task_path.read_text())
            task["ordinal"] = 99
            task["cell_id"] = "forged-cell"
            task_path.write_bytes(replay.pretty_bytes(task))
            reseal_return(returned)
            with self.assertRaisesRegex(replay.ReplayError, "identity changed"):
                replay.validate_return(packet, returned)


if __name__ == "__main__":
    unittest.main()
