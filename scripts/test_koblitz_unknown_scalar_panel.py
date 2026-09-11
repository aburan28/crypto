#!/usr/bin/env python3
"""Focused control-plane tests for the Stage-23 unknown-scalar panel."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import types
import unittest
from unittest import mock


SCRIPT = Path(__file__).with_name("run_koblitz_unknown_scalar_panel.py")
SPEC = importlib.util.spec_from_file_location("stage23", SCRIPT)
assert SPEC and SPEC.loader
stage23 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage23)


def target(index: int = 0) -> dict:
    return {
        "ordinal": index,
        "target_id": f"target-{index}",
        "point": {"x": str(index + 1), "y": str(index + 2)},
        "packed_point": index + 10,
        "draw_counter": index,
        "candidate_attempts": 1,
        "ic_seed": 11 + index,
        "rho_seed": 21 + index,
    }


def ic_result(item: dict, *, complete: bool = True) -> dict:
    scalar = "3" if complete else None
    value = {
        "schema": "koblitz_unknown_scalar_ic_result.v1",
        "status": "complete_verified" if complete else "incomplete",
        "profile": "smoke",
        "target_id": item["target_id"],
        "target": item["point"],
        "seed": item["ic_seed"],
        "target_scalar_constructed_or_supplied": False,
        "factor_base_logs_constructed_or_supplied": False,
        "report": {
            "trials": 1,
            "relations": 1,
            "outcomes": {
                "relation_found": 1,
                "refuted": 0,
                "unknown": 0,
                "invalid_model": 0,
                "direct_skipped": 0,
                "direct_solved": 0,
            },
            "recovered_scalar": scalar,
            "recovered_scalar_point_verified": complete,
            "direct_relation": False,
            "sat_invalid_models": 0,
            "matrix_rows": 1,
            "matrix_columns": 1,
            "terminal_matrix_rank": 1 if complete else 0,
            "sat_conflicts": 0,
            "sat_calls": 1,
            "sat_models": 1,
            "rank_checks": 1,
            "linear_solve_attempts": 1,
            "rank_history": [
                {
                    "rows": 1,
                    "columns": 1,
                    "rank": 1 if complete else 0,
                    "candidate_produced": complete,
                    "candidate_verified": complete,
                }
            ],
        },
        "attempt_records": [{"trial": 1, "disposition": "relation_found", "conflicts": 0, "solver_calls": 1, "models": 1}],
        "relation_matrix": [{"row": ["1"]}],
        "progress": [
            {"event": "relation_attempt_finished", "trial": 1},
            {"event": "matrix_rank", "rows": 1, "columns": 1, "rank": 1 if complete else 0},
        ],
    }
    value["attempt_records_blake3"] = stage23.custody.json_blake3(value["attempt_records"])
    value["relation_matrix_blake3"] = stage23.custody.json_blake3(value["relation_matrix"])
    value["progress_blake3"] = stage23.custody.json_blake3(value["progress"])
    return value


def rho_result(item: dict, *, complete: bool = True) -> dict:
    value = {
        "schema": "koblitz_unknown_scalar_rho_result.v1",
        "status": "complete_verified" if complete else "incomplete",
        "profile": "smoke",
        "target_id": item["target_id"],
        "target": item["point"],
        "seed": item["rho_seed"],
        "target_scalar_constructed_or_supplied": False,
        "report": {
            "recovered_scalar": "3" if complete else None,
            "recovered_scalar_point_verified": complete,
            "exhausted": not complete,
            "iterations": 2,
            "jump_table_rebuilds": 1,
            "restarts_attempted": 1,
        },
        "charges": {
            "coefficient_draws": 34,
            "setup_scalar_multiplications": 34,
            "setup_group_additions": 17,
            "walk_group_additions": 6,
            "canonicalizations": 7,
            "frobenius_maps": 49,
            "negations_examined": 49,
            "partition_hashes": 6,
            "collisions": 1,
            "failed_collisions": 0,
        },
        "progress": [
            {"event": "rho_restart_started"},
            {"event": "rho_jump_table_ready"},
            {"event": "rho_collision"},
            {"event": "rho_finished"},
        ],
        "timing_ns": {
            "target_and_subgroup_validation": 10,
            "rho_setup": 20,
            "rho_walk": 30,
            "candidate_verification": 5,
            "end_to_end": 80,
        },
    }
    value["progress_blake3"] = stage23.custody.json_blake3(value["progress"])
    return value


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def process_record(
    command: list[str], watchdog: float, *, returncode: int, orphan: bool = False
) -> dict:
    return {
        "command": command,
        "returncode": returncode,
        "watchdog_seconds": watchdog,
        "timed_out": False,
        "orphan_group_terminated": orphan,
        "metrics": {
            "wall_seconds": 1.0,
            "user_seconds": 0.08,
            "system_seconds": 0.02,
            "total_core_seconds": 0.1,
            "single_core_seconds": 0.1,
            "peak_rss_bytes": 100,
            "meter": "fresh-process getrusage(RUSAGE_CHILDREN)",
        },
    }


def write_task(
    run_root: Path,
    name: str,
    command: list[str],
    watchdog: float,
    inputs: list[dict],
    environment: dict[str, str],
    *,
    returncode: int,
) -> dict:
    root = run_root / "tasks" / name
    root.mkdir(parents=True)
    intent = {
        "schema": "koblitz_unknown_scalar_process_intent.v1",
        "name": name,
        "command": command,
        "watchdog_seconds": watchdog,
        "environment": environment,
        "inputs": inputs,
        "meter": stage23.custody.file_identity(stage23.METER),
        "meter_exclusive_create": True,
        "stdin": "devnull",
        "close_fds": True,
    }
    write_json(root / "intent.json", intent)
    (root / "stdout").write_bytes(b"")
    (root / "stderr").write_bytes(b"synthetic failed build\n")
    process = process_record(command, watchdog, returncode=returncode)
    write_json(root / "metrics.json", process)
    receipt = {
        "schema": "koblitz_unknown_scalar_process_receipt.v1",
        "name": name,
        "process": process,
        "intent": stage23.custody.file_identity(root / "intent.json"),
        "stdout": stage23.custody.file_identity(root / "stdout"),
        "stderr": stage23.custody.file_identity(root / "stderr"),
        "metrics": stage23.custody.file_identity(root / "metrics.json"),
        "inputs": inputs,
        "result": None,
        "parse_error": None,
        "terminal_status": "complete" if returncode == 0 else "incomplete",
    }
    write_json(root / "receipt.json", receipt)
    return {"receipt": receipt, "result": None}


def make_incomplete_build_run(parent: Path) -> tuple[Path, Path]:
    run_root = (parent / "run").resolve()
    (run_root / "tasks").mkdir(parents=True)
    (run_root / "binaries").mkdir()
    (run_root / "inputs").mkdir()
    frozen = stage23.protocol("smoke")
    source = stage23.source_binding()
    host = stage23.host_binding()
    shutil.copyfile(stage23.PROTOCOL, run_root / "inputs/protocol.json")
    write_json(run_root / "inputs/source.json", source)
    write_json(run_root / "inputs/host.json", host)
    source_identity = stage23.custody.file_identity(run_root / "inputs/source.json")
    protocol_identity = stage23.custody.file_identity(run_root / "inputs/protocol.json")
    command = [
        source["cargo"]["path"],
        "build",
        "--release",
        "--locked",
        "--jobs",
        "1",
        "--target-dir",
        str(run_root / "build-target"),
        "--example",
        "koblitz_public_factor_base_discovery",
        "--example",
        "koblitz_unknown_scalar_panel",
    ]
    task = write_task(
        run_root,
        "00-build",
        command,
        float(frozen["execution"]["build_watchdog_seconds"]),
        [source_identity, protocol_identity],
        stage23.environment(),
        returncode=1,
    )
    inner = [
        source["python"]["path"],
        source["direct_sources"][str(stage23.RUNNER_SOURCE.relative_to(stage23.REPO))]["path"],
        "run",
        "--profile",
        "smoke",
        "--output",
        str(run_root),
        "--meter",
        source["direct_sources"][str(stage23.METER.relative_to(stage23.REPO))]["path"],
        "--allow-dirty",
    ]
    stage23.finish_run(
        output=run_root,
        profile="smoke",
        frozen=frozen,
        source=source,
        host=host,
        tasks=[task],
        rows=[],
        targets=None,
        binaries={},
        expected_inner=inner,
        critical_failure="build",
    )
    outer_path = parent / "outer.metrics.json"
    outer = process_record(
        inner,
        float(frozen["execution"]["whole_driver_watchdog_seconds"]),
        returncode=0,
    )
    outer["metrics"].update(
        {
            "wall_seconds": 2.0,
            "user_seconds": 0.16,
            "system_seconds": 0.04,
            "total_core_seconds": 0.2,
            "single_core_seconds": 0.2,
            "peak_rss_bytes": 200,
        }
    )
    write_json(outer_path, outer)
    return run_root, outer_path


def reseal(run_root: Path) -> None:
    seal_path = run_root / "run-seal.json"
    seal = json.loads(seal_path.read_text())
    seal["summary"] = stage23.custody.file_identity(run_root / "run-summary.json")
    inventory = stage23.custody.inventory(run_root, {"run-seal.json"})
    seal["inventory"] = inventory
    seal["inventory_sha256"] = stage23.custody.canonical_sha256(inventory)
    seal.pop("seal_payload_sha256", None)
    seal["seal_payload_sha256"] = stage23.custody.canonical_sha256(seal)
    write_json(seal_path, seal)


class Stage23Tests(unittest.TestCase):
    def test_protocol_and_self_test(self) -> None:
        result = stage23.self_test()
        self.assertEqual(result["self_test"], "pass")
        self.assertEqual(stage23.profile_values("production", stage23.protocol("production"))["targets"], 5)

    def test_target_validation_rejects_scalar_labels_and_duplicates(self) -> None:
        good = {
            "schema": "koblitz_unknown_scalar_target_panel.v1",
            "status": "complete",
            "target_scalar_constructed_or_recorded": False,
            "factor_base_log_labels_constructed_or_recorded": False,
            "identity": {"profile": "smoke", "targets": [target(0), target(1)]},
        }
        good["identity_blake3"] = stage23.custody.json_blake3(good["identity"])
        self.assertEqual(len(stage23.validate_targets(good, "smoke", 2)), 2)
        bad_scalar = dict(good, secret=7)
        with self.assertRaises(stage23.Stage23Error):
            stage23.validate_targets(bad_scalar, "smoke", 2)
        duplicate = dict(good)
        duplicate["identity"] = {"profile": "smoke", "targets": [target(0), target(0)]}
        with self.assertRaises(stage23.Stage23Error):
            stage23.validate_targets(duplicate, "smoke", 2)

    def test_ic_inventory_and_completion_checks(self) -> None:
        item = target()
        self.assertTrue(stage23.validate_ic(ic_result(item), item))
        incomplete = ic_result(item, complete=False)
        self.assertFalse(stage23.validate_ic(incomplete, item))
        rankless = ic_result(item, complete=False)
        rankless["report"]["rank_checks"] = 0
        rankless["report"]["linear_solve_attempts"] = 0
        rankless["report"]["rank_history"] = []
        rankless["progress"] = [
            row for row in rankless["progress"] if row["event"] != "matrix_rank"
        ]
        rankless["progress_blake3"] = stage23.custody.json_blake3(rankless["progress"])
        with self.assertRaises(stage23.Stage23Error):
            stage23.validate_ic(rankless, item)
        tampered = ic_result(item)
        tampered["attempt_records"] = []
        with self.assertRaises(stage23.Stage23Error):
            stage23.validate_ic(tampered, item)

    def test_rho_operation_ledger(self) -> None:
        item = target()
        self.assertTrue(stage23.validate_rho(rho_result(item), item))
        tampered = rho_result(item)
        tampered["charges"]["setup_group_additions"] = 33
        with self.assertRaises(stage23.Stage23Error):
            stage23.validate_rho(tampered, item)
        overlapping = rho_result(item)
        overlapping["timing_ns"]["rho_walk"] = 100
        with self.assertRaises(stage23.Stage23Error):
            stage23.validate_rho(overlapping, item)

    def test_orphaned_process_is_never_complete(self) -> None:
        with tempfile.TemporaryDirectory(prefix="stage23-orphan-") as temporary:
            root = Path(temporary)
            (root / "tasks").mkdir()

            def fake_run(*_args, **_kwargs):
                task_root = root / "tasks/orphan"
                (task_root / "stdout").write_bytes(b"")
                (task_root / "stderr").write_bytes(b"")
                write_json(
                    task_root / "metrics.json",
                    process_record(["/usr/bin/true"], 1.0, returncode=0, orphan=True),
                )
                return types.SimpleNamespace(returncode=0)

            with mock.patch.object(stage23.subprocess, "run", side_effect=fake_run):
                task = stage23.run_metered(
                    root=root,
                    name="orphan",
                    command=["/usr/bin/true"],
                    timeout=1.0,
                    meter=stage23.METER,
                    env=stage23.environment(),
                    inputs=[],
                    expect_json=False,
                )
            self.assertEqual(task["receipt"]["terminal_status"], "incomplete")
            self.assertTrue(task["receipt"]["process"]["orphan_group_terminated"])

    def test_full_verifier_reconstructs_and_rejects_coherent_tampering(self) -> None:
        with tempfile.TemporaryDirectory(prefix="stage23-verify-") as temporary:
            parent = Path(temporary)
            run_root, outer = make_incomplete_build_run(parent / "baseline")
            args = types.SimpleNamespace(
                run_root=run_root,
                outer_metrics=outer,
                output=parent / "baseline-verification",
            )
            verified = stage23.verify(args)
            self.assertEqual(verified["status"], "verification_frozen")

            forged_root, forged_outer = make_incomplete_build_run(parent / "forged-command")
            task_root = forged_root / "tasks/00-build"
            intent = json.loads((task_root / "intent.json").read_text())
            intent["command"] = [*intent["command"], "--forged"]
            write_json(task_root / "intent.json", intent)
            metrics = json.loads((task_root / "metrics.json").read_text())
            metrics["command"] = intent["command"]
            write_json(task_root / "metrics.json", metrics)
            receipt = json.loads((task_root / "receipt.json").read_text())
            receipt["process"] = metrics
            receipt["intent"] = stage23.custody.file_identity(task_root / "intent.json")
            receipt["metrics"] = stage23.custody.file_identity(task_root / "metrics.json")
            write_json(task_root / "receipt.json", receipt)
            reseal(forged_root)
            forged_args = types.SimpleNamespace(
                run_root=forged_root,
                outer_metrics=forged_outer,
                output=parent / "forged-command-verification",
            )
            with self.assertRaises(stage23.Stage23Error):
                stage23.verify(forged_args)

            state_root, state_outer = make_incomplete_build_run(parent / "forged-state")
            summary_path = state_root / "run-summary.json"
            summary = json.loads(summary_path.read_text())
            summary["status"] = "complete_verified_panel"
            summary["completed_rows"] = 2
            summary["ratios"]["verdict"] = "finite_complete_panel"
            write_json(summary_path, summary)
            state_seal_path = state_root / "run-seal.json"
            state_seal = json.loads(state_seal_path.read_text())
            state_seal["panel_complete"] = True
            write_json(state_seal_path, state_seal)
            reseal(state_root)
            state_args = types.SimpleNamespace(
                run_root=state_root,
                outer_metrics=state_outer,
                output=parent / "forged-state-verification",
            )
            with self.assertRaises(stage23.Stage23Error):
                stage23.verify(state_args)

    def test_plan_is_write_once_and_outside_checkout(self) -> None:
        parser_args = type("Args", (), {})()
        parser_args.profile = "smoke"
        parser_args.meter = stage23.METER
        parser_args.allow_dirty = True
        with tempfile.TemporaryDirectory(prefix="stage23-plan-") as temporary:
            parser_args.output = Path(temporary) / "new-output"
            rendered = stage23.plan(parser_args)
            self.assertEqual(rendered["expected_targets"], 2)
            self.assertEqual(rendered["expected_child_processes"], 8)
            parser_args.output.mkdir()
            with self.assertRaises(stage23.custody.Stage21Error):
                stage23.plan(parser_args)
        parser_args.output = stage23.REPO / "forbidden"
        with self.assertRaises(stage23.custody.Stage21Error):
            stage23.plan(parser_args)


if __name__ == "__main__":
    unittest.main()
