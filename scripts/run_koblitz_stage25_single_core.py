#!/usr/bin/env python3
"""Run and verify a Linux one-CPU Stage-23 production control."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace
from typing import Any

import run_koblitz_relation_yield_bridge as custody
import run_koblitz_unknown_scalar_panel as stage23


REPO = Path(__file__).resolve().parents[1]
PROTOCOL = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-25-single-core-protocol.json"
RUNNER = Path(__file__).resolve()
WORKFLOW = REPO / ".github/workflows/koblitz-stage25-single-core.yml"
RESULT_SCHEMA = "koblitz_stage25_single_core_result.v1"
AFFINITY_SCHEMA = "koblitz_stage25_affinity_receipt.v1"
SEAL_SCHEMA = "koblitz_stage25_single_core_seal.v1"


class Stage25Error(RuntimeError):
    pass


def load_json(path: Path, context: str) -> dict[str, Any]:
    value, _ = custody.read_json(path, context)
    return value


def require_exact_keys(value: Any, keys: set[str], context: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != keys:
        raise Stage25Error(f"{context} uses an unexpected schema")
    return value


def probe_child_affinity() -> list[int]:
    code = "import json,os; print(json.dumps(sorted(os.sched_getaffinity(0))))"
    completed = subprocess.run(
        [sys.executable, "-c", code],
        text=True,
        capture_output=True,
        check=True,
    )
    value = json.loads(completed.stdout)
    if not isinstance(value, list) or any(type(item) is not int for item in value):
        raise Stage25Error("child affinity probe returned an invalid CPU list")
    return value


def choose_and_set_cpu(requested: int | None) -> tuple[list[int], int, list[int]]:
    if platform.system() != "Linux" or not hasattr(os, "sched_getaffinity"):
        raise Stage25Error("Stage-25 requires Linux sched_getaffinity/sched_setaffinity")
    before = sorted(os.sched_getaffinity(0))
    if not before:
        raise Stage25Error("initial Linux affinity set is empty")
    selected = before[0] if requested is None else requested
    if selected not in before:
        raise Stage25Error("requested CPU is outside the initial allowed affinity set")
    os.sched_setaffinity(0, {selected})
    effective = sorted(os.sched_getaffinity(0))
    if effective != [selected]:
        raise Stage25Error("kernel did not restrict the Stage-25 parent to one CPU")
    return before, selected, effective


def validate_affinity(value: Any) -> dict[str, Any]:
    value = require_exact_keys(
        value,
        {
            "schema",
            "platform",
            "uname",
            "initial_allowed_cpus",
            "selected_cpu",
            "parent_effective_before",
            "child_effective_before",
            "parent_effective_after",
            "child_effective_after",
            "inheritance_contract",
        },
        "affinity receipt",
    )
    if value["schema"] != AFFINITY_SCHEMA or value["platform"] != "Linux":
        raise Stage25Error("affinity receipt has the wrong platform or schema")
    selected = value["selected_cpu"]
    if type(selected) is not int:
        raise Stage25Error("selected CPU is not an integer")
    initial = value["initial_allowed_cpus"]
    if (
        not isinstance(initial, list)
        or not initial
        or any(type(cpu) is not int or cpu < 0 for cpu in initial)
        or initial != sorted(set(initial))
    ):
        raise Stage25Error("initial affinity set is invalid")
    if selected not in initial:
        raise Stage25Error("selected CPU was not initially allowed")
    for field in (
        "parent_effective_before",
        "child_effective_before",
        "parent_effective_after",
        "child_effective_after",
    ):
        if value[field] != [selected]:
            raise Stage25Error(f"{field} is not the selected singleton CPU")
    if value["inheritance_contract"] != "Linux process affinity is inherited by fork/exec descendants unless explicitly changed":
        raise Stage25Error("affinity inheritance contract changed")
    uname = value["uname"]
    if not isinstance(uname, list) or len(uname) != 6 or any(not isinstance(item, str) for item in uname):
        raise Stage25Error("affinity uname record is invalid")
    return value


def affinity_reset_source_audit() -> dict[str, Any]:
    tokens = ("sched_setaffinity", "pthread_setaffinity", "taskset", "core_affinity")
    paths: list[Path] = []
    for root in (REPO / "src", REPO / "examples"):
        paths.extend(path for path in root.rglob("*") if path.is_file())
    paths.extend(
        (
            stage23.RUNNER_SOURCE,
            stage23.METER,
            REPO / "scripts/run_koblitz_relation_yield_bridge.py",
        )
    )
    matches = []
    for path in sorted(set(paths)):
        try:
            text = path.read_text()
        except (OSError, UnicodeDecodeError):
            continue
        for token in tokens:
            if token in text:
                matches.append({"path": str(path.relative_to(REPO)), "token": token})
    if matches:
        raise Stage25Error(f"Stage-23 source contains an affinity-reset primitive: {matches}")
    return {
        "roots": ["src", "examples", "Stage-23 runner/meter/custody scripts"],
        "forbidden_tokens": list(tokens),
        "matches": [],
    }


def source_binding() -> dict[str, Any]:
    upstream = stage23.source_binding()
    if upstream["git"]["dirty"]:
        raise Stage25Error("Stage-25 production requires a clean Git checkout")
    direct = {
        str(path.relative_to(REPO)): custody.file_identity(path)
        for path in (RUNNER, PROTOCOL, WORKFLOW)
    }
    return {
        "schema": "koblitz_stage25_source_binding.v1",
        "git": upstream["git"],
        "python": upstream["python"],
        "stage23_source": upstream,
        "stage25_direct_sources": direct,
        "affinity_reset_source_audit": affinity_reset_source_audit(),
    }


def complete_process(metrics: dict[str, Any]) -> bool:
    return (
        metrics["returncode"] == 0
        and metrics["timed_out"] is False
        and metrics["orphan_group_terminated"] is False
    )


def validate_stage23_artifacts(output: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Reconstruct the embedded Stage-23 run before using any summary fields."""
    run_root = output / "stage23-run"
    outer_root = output / "stage23-run-outer"
    verification_root = output / "stage23-verification"
    run_seal = load_json(run_root / "run-seal.json", "Stage-23 run seal")
    run_seal = stage23.require_exact_keys(
        run_seal,
        {
            "schema",
            "status",
            "profile",
            "protocol_sha256",
            "source_commit",
            "summary",
            "inventory",
            "inventory_sha256",
            "expected_inner_command",
            "panel_complete",
            "scientific_measurement_admitted",
            "seal_payload_sha256",
        },
        "Stage-23 run seal",
    )
    run_payload = dict(run_seal)
    claimed_run_payload = run_payload.pop("seal_payload_sha256")
    if claimed_run_payload != custody.canonical_sha256(run_payload):
        raise Stage25Error("embedded Stage-23 run seal hash is invalid")
    actual_inventory = custody.inventory(run_root, {"run-seal.json"})
    if (
        actual_inventory != run_seal["inventory"]
        or custody.canonical_sha256(actual_inventory) != run_seal["inventory_sha256"]
    ):
        raise Stage25Error("embedded Stage-23 run inventory changed")
    stage23.validate_bound_file(
        run_seal["summary"],
        run_root / "run-summary.json",
        run_root,
        "Stage-23 run summary",
    )
    summary = load_json(run_root / "run-summary.json", "Stage-23 run summary")
    reconstructed = stage23.reconstruct_run(run_root, summary, run_seal)
    outer = custody.process_metrics(
        outer_root / "driver.metrics.json", "Stage-23 outer metrics"
    )
    if (
        outer["command"] != run_seal["expected_inner_command"]
        or outer["watchdog_seconds"]
        != float(reconstructed["frozen"]["execution"]["whole_driver_watchdog_seconds"])
        or not stage23.process_completed_cleanly(outer)
    ):
        raise Stage25Error("embedded Stage-23 outer receipt is invalid")
    child = reconstructed["resources"]
    if (
        outer["metrics"]["total_core_seconds"] + 1e-9
        < child["total_core_seconds"]
        or outer["metrics"]["wall_seconds"] + 1e-9
        < child["summed_process_wall_seconds"]
        or outer["metrics"]["peak_rss_bytes"] < child["peak_process_rss_bytes"]
    ):
        raise Stage25Error("embedded Stage-23 outer receipt does not enclose its children")

    verification_seal = load_json(
        verification_root / "verification-seal.json", "Stage-23 verification seal"
    )
    verification_seal = stage23.require_exact_keys(
        verification_seal,
        {
            "schema",
            "status",
            "verification",
            "inventory",
            "inventory_sha256",
            "seal_payload_sha256",
        },
        "Stage-23 verification seal",
    )
    verification_payload = dict(verification_seal)
    claimed_verification_payload = verification_payload.pop("seal_payload_sha256")
    if claimed_verification_payload != custody.canonical_sha256(verification_payload):
        raise Stage25Error("embedded Stage-23 verification seal hash is invalid")
    verification_inventory = custody.inventory(
        verification_root, {"verification-seal.json"}
    )
    if (
        verification_inventory != verification_seal["inventory"]
        or custody.canonical_sha256(verification_inventory)
        != verification_seal["inventory_sha256"]
    ):
        raise Stage25Error("embedded Stage-23 verification inventory changed")
    stage23.validate_bound_file(
        verification_seal["verification"],
        verification_root / "verification.json",
        verification_root,
        "Stage-23 verification",
    )
    verification = load_json(
        verification_root / "verification.json", "Stage-23 verification"
    )
    if (
        verification.get("panel_status") != summary.get("status")
        or verification.get("completed_rows") != summary.get("completed_rows")
        or verification.get("expected_rows") != summary.get("expected_rows")
        or verification.get("resources", {}).get("charged_children") != child
        or verification.get("resources", {}).get("whole_driver") != outer["metrics"]
        or verification.get("scientific_measurement_admitted") is not False
        or verification.get("full_cost_gate_passed") is not False
        or verification.get("koblitz_index_calculus_sota") is not False
    ):
        raise Stage25Error("embedded Stage-23 project verification is not derivable")
    return summary, outer, verification


def build_result(
    *,
    profile: str,
    source: dict[str, Any],
    affinity: dict[str, Any],
    wrapper_wall: float,
    wrapper_child_user: float,
    wrapper_child_system: float,
    output: Path,
) -> dict[str, Any]:
    validate_affinity(affinity)
    run_root = output / "stage23-run"
    outer_root = output / "stage23-run-outer"
    verification_root = output / "stage23-verification"
    summary, outer, verification = validate_stage23_artifacts(output)
    run_seal = load_json(run_root / "run-seal.json", "Stage-23 run seal")
    if (
        summary.get("profile") != profile
        or summary.get("status") != "complete_verified_panel"
        or summary.get("completed_rows") != summary.get("expected_rows")
        or run_seal.get("panel_complete") is not True
        or verification.get("panel_status") != "complete_verified_panel"
        or verification.get("completed_rows") != verification.get("expected_rows")
        or not complete_process(outer)
    ):
        raise Stage25Error("Stage-23 production panel did not complete and verify")
    if profile != "production" or summary["expected_rows"] != 5:
        raise Stage25Error("Stage-25 admission requires the five-target production profile")
    child = summary["resources"]
    ratios = summary["ratios"]
    ic_core = sum(row["ic"]["resources"]["total_core_seconds"] for row in summary["rows"])
    rho_core = sum(row["rho"]["resources"]["total_core_seconds"] for row in summary["rows"])
    result = {
        "schema": RESULT_SCHEMA,
        "status": "complete_single_cpu_affinity_control",
        "profile": profile,
        "source_commit": source["git"]["commit"],
        "source_tree_sha256": source["git"]["committed_tree_sha256"],
        "affinity": affinity,
        "stage23": {
            "run_seal": custody.file_identity(run_root / "run-seal.json"),
            "run_summary": custody.file_identity(run_root / "run-summary.json"),
            "outer_metrics": custody.file_identity(outer_root / "driver.metrics.json"),
            "verification_seal": custody.file_identity(verification_root / "verification-seal.json"),
            "completed_rows": summary["completed_rows"],
            "expected_rows": summary["expected_rows"],
        },
        "measurements": {
            "single_core_elapsed_seconds": outer["metrics"]["wall_seconds"],
            "outer_total_core_seconds": outer["metrics"]["total_core_seconds"],
            "outer_user_seconds": outer["metrics"]["user_seconds"],
            "outer_system_seconds": outer["metrics"]["system_seconds"],
            "outer_peak_rss_bytes": outer["metrics"]["peak_rss_bytes"],
            "stage23_child_total_core_seconds": child["total_core_seconds"],
            "stage23_summed_process_wall_seconds": child["summed_process_wall_seconds"],
            "stage23_peak_process_rss_bytes": child["peak_process_rss_bytes"],
            "ic_total_core_seconds": ic_core,
            "rho_total_core_seconds": rho_core,
            "online_ic_over_rho_core": ratios["online_ic_over_rho_core"],
            "setup_charged_ic_over_rho_core": ratios["setup_charged_ic_over_rho_core"],
            "affinity_wrapper_wall_seconds": wrapper_wall,
            "affinity_wrapper_child_user_seconds": wrapper_child_user,
            "affinity_wrapper_child_system_seconds": wrapper_child_system,
        },
        "accounting_boundary": {
            "all_stage23_descendants_inherit_one_cpu_affinity": True,
            "single_core_elapsed_is_outer_wall_under_kernel_affinity": True,
            "simultaneous_process_tree_memory_measured": False,
            "source_dependency_toolchain_or_license_acquisition_charged": False,
            "covers_stage20_n31_n41_n59_backend_matrix": False,
            "licensed_magma_executed": False,
        },
        "retained_mathematical_witness_replay_completed": True,
        "independent_mathematical_payload_replay_completed": False,
        "scientific_measurement_admitted": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    return result


def validate_result(value: Any) -> dict[str, Any]:
    required = {
        "schema",
        "status",
        "profile",
        "source_commit",
        "source_tree_sha256",
        "affinity",
        "stage23",
        "measurements",
        "accounting_boundary",
        "retained_mathematical_witness_replay_completed",
        "independent_mathematical_payload_replay_completed",
        "scientific_measurement_admitted",
        "independent_external_reproduction_satisfied",
        "full_cost_gate_passed",
        "koblitz_index_calculus_sota",
    }
    value = require_exact_keys(value, required, "Stage-25 result")
    if value["schema"] != RESULT_SCHEMA or value["status"] != "complete_single_cpu_affinity_control":
        raise Stage25Error("Stage-25 result is not complete")
    validate_affinity(value["affinity"])
    measurements = value["measurements"]
    expected_measurements = {
        "single_core_elapsed_seconds",
        "outer_total_core_seconds",
        "outer_user_seconds",
        "outer_system_seconds",
        "outer_peak_rss_bytes",
        "stage23_child_total_core_seconds",
        "stage23_summed_process_wall_seconds",
        "stage23_peak_process_rss_bytes",
        "ic_total_core_seconds",
        "rho_total_core_seconds",
        "online_ic_over_rho_core",
        "setup_charged_ic_over_rho_core",
        "affinity_wrapper_wall_seconds",
        "affinity_wrapper_child_user_seconds",
        "affinity_wrapper_child_system_seconds",
    }
    if not isinstance(measurements, dict) or set(measurements) != expected_measurements:
        raise Stage25Error("Stage-25 measurement schema changed")
    for field, number in measurements.items():
        if (
            isinstance(number, bool)
            or not isinstance(number, (int, float))
            or not math.isfinite(number)
            or number < 0
        ):
            raise Stage25Error(f"Stage-25 measurement {field} is invalid")
    if measurements["single_core_elapsed_seconds"] <= 0:
        raise Stage25Error("single-core elapsed time is not positive")
    boundary = value["accounting_boundary"]
    if boundary != {
        "all_stage23_descendants_inherit_one_cpu_affinity": True,
        "single_core_elapsed_is_outer_wall_under_kernel_affinity": True,
        "simultaneous_process_tree_memory_measured": False,
        "source_dependency_toolchain_or_license_acquisition_charged": False,
        "covers_stage20_n31_n41_n59_backend_matrix": False,
        "licensed_magma_executed": False,
    }:
        raise Stage25Error("Stage-25 accounting boundary changed")
    if value["retained_mathematical_witness_replay_completed"] is not True:
        raise Stage25Error("retained witness replay status changed")
    for field in (
        "independent_mathematical_payload_replay_completed",
        "scientific_measurement_admitted",
        "independent_external_reproduction_satisfied",
        "full_cost_gate_passed",
        "koblitz_index_calculus_sota",
    ):
        if value[field] is not False:
            raise Stage25Error(f"Stage-25 widened {field}")
    return value


def run(args: argparse.Namespace) -> tuple[dict[str, Any], bool]:
    protocol_identity = custody.file_identity(PROTOCOL)
    protocol = load_json(PROTOCOL, "Stage-25 protocol")
    if protocol.get("schema") != "koblitz_stage25_single_core_protocol.v1" or protocol.get("status") != "frozen_before_execution":
        raise Stage25Error("Stage-25 protocol is not frozen")
    output = custody.safe_new_directory(args.output, "Stage-25 output")
    source = source_binding()
    initial, selected, effective = choose_and_set_cpu(args.cpu)
    child_before = probe_child_affinity()
    output.mkdir(parents=True)
    run_root = output / "stage23-run"
    outer_root = output / "stage23-run-outer"
    verification_root = output / "stage23-verification"
    outer_root.mkdir()
    shutil.copyfile(PROTOCOL, output / "protocol.json")
    archived_protocol_identity = custody.file_identity(output / "protocol.json")
    if (
        archived_protocol_identity["bytes"],
        archived_protocol_identity["sha256"],
    ) != (protocol_identity["bytes"], protocol_identity["sha256"]):
        raise Stage25Error("archived Stage-25 protocol differs from its source")
    custody.write_json_new(output / "source.json", source)
    plan_command = [
        str(Path(sys.executable).resolve()),
        str(stage23.RUNNER_SOURCE),
        "plan",
        "--profile",
        args.profile,
        "--output",
        str(run_root),
        "--meter",
        str(stage23.METER),
    ]
    plan_process = subprocess.run(plan_command, text=True, capture_output=True, check=True)
    plan = json.loads(plan_process.stdout)
    custody.write_json_new(output / "stage23-plan.json", plan)
    start_wall = time.monotonic()
    before_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    outer_process = subprocess.run(plan["outer_command"], text=True, capture_output=True)
    custody.write_new(output / "outer-launch.stdout", outer_process.stdout.encode())
    custody.write_new(output / "outer-launch.stderr", outer_process.stderr.encode())
    verify_process = subprocess.run(
        [
            str(Path(sys.executable).resolve()),
            str(stage23.RUNNER_SOURCE),
            "verify",
            "--run-root",
            str(run_root),
            "--outer-metrics",
            str(outer_root / "driver.metrics.json"),
            "--output",
            str(verification_root),
        ],
        text=True,
        capture_output=True,
    )
    custody.write_new(output / "stage23-verify.stdout", verify_process.stdout.encode())
    custody.write_new(output / "stage23-verify.stderr", verify_process.stderr.encode())
    after_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    wrapper_wall = time.monotonic() - start_wall
    parent_after = sorted(os.sched_getaffinity(0))
    child_after = probe_child_affinity()
    affinity = {
        "schema": AFFINITY_SCHEMA,
        "platform": platform.system(),
        "uname": list(platform.uname()),
        "initial_allowed_cpus": initial,
        "selected_cpu": selected,
        "parent_effective_before": effective,
        "child_effective_before": child_before,
        "parent_effective_after": parent_after,
        "child_effective_after": child_after,
        "inheritance_contract": "Linux process affinity is inherited by fork/exec descendants unless explicitly changed",
    }
    custody.write_json_new(output / "affinity.json", affinity)
    complete = outer_process.returncode == 0 and verify_process.returncode == 0
    if complete:
        result = build_result(
            profile=args.profile,
            source=source,
            affinity=affinity,
            wrapper_wall=wrapper_wall,
            wrapper_child_user=after_usage.ru_utime - before_usage.ru_utime,
            wrapper_child_system=after_usage.ru_stime - before_usage.ru_stime,
            output=output,
        )
    else:
        result = {
            "schema": RESULT_SCHEMA,
            "status": "incomplete",
            "profile": args.profile,
            "source_commit": source["git"]["commit"],
            "affinity": affinity,
            "outer_returncode": outer_process.returncode,
            "verify_returncode": verify_process.returncode,
            "full_cost_gate_passed": False,
            "koblitz_index_calculus_sota": False,
        }
    custody.write_json_new(output / "result.json", result)
    seal_payload = {
        "schema": SEAL_SCHEMA,
        "status": "complete" if complete else "incomplete",
        "protocol": archived_protocol_identity,
        "source": custody.file_identity(output / "source.json"),
        "affinity": custody.file_identity(output / "affinity.json"),
        "result": custody.file_identity(output / "result.json"),
        "stage23_plan": custody.file_identity(output / "stage23-plan.json"),
        "stage23_run_seal": custody.file_identity(run_root / "run-seal.json") if (run_root / "run-seal.json").exists() else None,
        "stage23_outer_metrics": custody.file_identity(outer_root / "driver.metrics.json") if (outer_root / "driver.metrics.json").exists() else None,
        "stage23_verification_seal": custody.file_identity(verification_root / "verification-seal.json") if (verification_root / "verification-seal.json").exists() else None,
    }
    seal = dict(seal_payload)
    seal["seal_payload_sha256"] = custody.canonical_sha256(seal_payload)
    custody.write_json_new(output / "result-seal.json", seal)
    return seal, complete


def verify(args: argparse.Namespace) -> dict[str, Any]:
    output = args.output.resolve()
    seal = load_json(output / "result-seal.json", "Stage-25 seal")
    expected_keys = {
        "schema",
        "status",
        "protocol",
        "source",
        "affinity",
        "result",
        "stage23_plan",
        "stage23_run_seal",
        "stage23_outer_metrics",
        "stage23_verification_seal",
        "seal_payload_sha256",
    }
    seal = require_exact_keys(seal, expected_keys, "Stage-25 seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256")
    if claimed != custody.canonical_sha256(payload):
        raise Stage25Error("Stage-25 seal hash is invalid")
    if seal["schema"] != SEAL_SCHEMA or seal["status"] != "complete":
        raise Stage25Error("Stage-25 seal is not complete")
    for field in (
        "protocol",
        "source",
        "affinity",
        "result",
        "stage23_plan",
        "stage23_run_seal",
        "stage23_outer_metrics",
        "stage23_verification_seal",
    ):
        identity = seal[field]
        if identity is None:
            raise Stage25Error(f"Stage-25 seal omits {field}")
        custody.validate_file_identity(identity, field)
    current_protocol = custody.file_identity(PROTOCOL)
    if (
        seal["protocol"]["bytes"],
        seal["protocol"]["sha256"],
    ) != (current_protocol["bytes"], current_protocol["sha256"]):
        raise Stage25Error("archived Stage-25 protocol differs from the current source")
    source = load_json(output / "source.json", "Stage-25 source binding")
    if source != source_binding():
        raise Stage25Error("Stage-25 source binding differs from the current clean checkout")
    result = validate_result(load_json(output / "result.json", "Stage-25 result"))
    affinity = validate_affinity(load_json(output / "affinity.json", "Stage-25 affinity"))
    if result["affinity"] != affinity:
        raise Stage25Error("Stage-25 result and affinity receipt differ")
    measurements = result["measurements"]
    reconstructed = build_result(
        profile=result["profile"],
        source=source,
        affinity=affinity,
        wrapper_wall=measurements["affinity_wrapper_wall_seconds"],
        wrapper_child_user=measurements["affinity_wrapper_child_user_seconds"],
        wrapper_child_system=measurements["affinity_wrapper_child_system_seconds"],
        output=output,
    )
    if result != reconstructed:
        raise Stage25Error("Stage-25 result differs from reconstructed Stage-23 evidence")
    outer = custody.process_metrics(output / "stage23-run-outer/driver.metrics.json", "Stage-23 outer metrics")
    summary = load_json(output / "stage23-run/run-summary.json", "Stage-23 summary")
    if result["measurements"]["single_core_elapsed_seconds"] != outer["metrics"]["wall_seconds"]:
        raise Stage25Error("single-core elapsed time differs from the outer wall receipt")
    if result["measurements"]["stage23_child_total_core_seconds"] != summary["resources"]["total_core_seconds"]:
        raise Stage25Error("child core-seconds differ from the Stage-23 summary")
    if result["measurements"]["affinity_wrapper_wall_seconds"] < outer["metrics"]["wall_seconds"]:
        raise Stage25Error("affinity wrapper wall does not enclose the Stage-23 outer wall")
    return {
        "schema": "koblitz_stage25_single_core_verification.v1",
        "status": "single_cpu_affinity_receipt_verified",
        "selected_cpu": affinity["selected_cpu"],
        "single_core_elapsed_seconds": result["measurements"]["single_core_elapsed_seconds"],
        "outer_total_core_seconds": result["measurements"]["outer_total_core_seconds"],
        "outer_peak_rss_bytes": result["measurements"]["outer_peak_rss_bytes"],
        "completed_rows": result["stage23"]["completed_rows"],
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def self_test() -> dict[str, Any]:
    sample = {
        "schema": AFFINITY_SCHEMA,
        "platform": "Linux",
        "uname": ["Linux", "host", "kernel", "version", "x86_64", ""],
        "initial_allowed_cpus": [2, 3],
        "selected_cpu": 2,
        "parent_effective_before": [2],
        "child_effective_before": [2],
        "parent_effective_after": [2],
        "child_effective_after": [2],
        "inheritance_contract": "Linux process affinity is inherited by fork/exec descendants unless explicitly changed",
    }
    validate_affinity(sample)
    broken = dict(sample)
    broken["child_effective_after"] = [2, 3]
    try:
        validate_affinity(broken)
    except Stage25Error:
        pass
    else:
        raise AssertionError("multi-CPU child affinity was accepted")
    return {
        "schema": "koblitz_stage25_self_test.v1",
        "status": "PASS",
        "checks": 2,
        "protocol_sha256": custody.file_identity(PROTOCOL)["sha256"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--profile", choices=("production",), default="production")
    run_parser.add_argument("--cpu", type=int)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--output", type=Path, required=True)
    subparsers.add_parser("self-test")
    args = parser.parse_args()
    try:
        if args.command == "run":
            value, complete = run(args)
            print(json.dumps(value, indent=2, sort_keys=True))
            if not complete:
                raise Stage25Error("Stage-25 run retained an incomplete terminal")
        elif args.command == "verify":
            print(json.dumps(verify(args), indent=2, sort_keys=True))
        else:
            print(json.dumps(self_test(), indent=2, sort_keys=True))
    except (
        Stage25Error,
        stage23.Stage23Error,
        custody.Stage21Error,
        OSError,
        ValueError,
        KeyError,
        subprocess.SubprocessError,
    ) as error:
        parser.exit(1, f"stage25: {error}\n")


if __name__ == "__main__":
    main()
