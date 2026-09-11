#!/usr/bin/env python3
"""Plan, run, seal, and verify the Stage-23 public unknown-scalar panel."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from typing import Any

import run_koblitz_relation_yield_bridge as custody


REPO = Path(__file__).resolve().parents[1]
PROTOCOL = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-23-unknown-scalar-protocol.json"
METER = REPO / "scripts/process_meter.py"
DISCOVERY_SOURCE = REPO / "examples/koblitz_public_factor_base_discovery.rs"
PANEL_SOURCE = REPO / "examples/koblitz_unknown_scalar_panel.rs"
LOCK = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock"
RUNNER_SOURCE = Path(__file__).resolve()
RUN_SEAL_SCHEMA = "koblitz_unknown_scalar_run_seal.v1"
RUN_SUMMARY_SCHEMA = "koblitz_unknown_scalar_run_summary.v1"
VERIFICATION_SCHEMA = "koblitz_unknown_scalar_verification.v1"
VERIFICATION_SEAL_SCHEMA = "koblitz_unknown_scalar_verification_seal.v1"
SOURCE_RELATIVE_PATHS = (
    "Cargo.toml",
    "Cargo.lock",
    "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock",
    "examples/koblitz_public_factor_base_discovery.rs",
    "examples/koblitz_unknown_scalar_panel.rs",
    "scripts/run_koblitz_unknown_scalar_panel.py",
    "scripts/process_meter.py",
    "scripts/run_koblitz_relation_yield_bridge.py",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/koblitz_pdp_phase_a.rs",
    "src/cryptanalysis/sat.rs",
    "src/cryptanalysis/semaev_sat.rs",
)


class Stage23Error(RuntimeError):
    pass


def load_json(path: Path, context: str) -> dict[str, Any]:
    value, _ = custody.read_json(path, context)
    return value


def validate_protocol_value(value: dict[str, Any], profile: str) -> dict[str, Any]:
    if value.get("schema") != "koblitz_unknown_scalar_panel_protocol.v1":
        raise Stage23Error("wrong Stage-23 protocol schema")
    if value.get("status") != "frozen_before_production_execution":
        raise Stage23Error("Stage-23 protocol is not frozen")
    if profile not in {"production", "smoke"}:
        raise Stage23Error("profile must be production or smoke")
    if value["targets"]["count"] != 5 or value["factor_base"]["divisor_indices"] != [0, 2]:
        raise Stage23Error("production panel or factor base changed")
    if value["index_calculus"]["parallel_workers"] != 1 or value["execution"]["parallel_workers"] != 1:
        raise Stage23Error("Stage-23 must remain single-worker")
    return value


def protocol(profile: str) -> dict[str, Any]:
    return validate_protocol_value(load_json(PROTOCOL, "Stage-23 protocol"), profile)


def profile_values(profile: str, frozen: dict[str, Any]) -> dict[str, Any]:
    if profile == "production":
        return {
            "targets": frozen["targets"]["count"],
            "n": frozen["curve"]["n"],
            "a": frozen["curve"]["a"],
            "dimension": frozen["factor_base"]["dimension"],
            "evidence_class": "public_synthetic_candidate_pending_independent_replay",
        }
    smoke = frozen["smoke"]
    return {
        "targets": smoke["targets"],
        "n": smoke["curve_n"],
        "a": smoke["curve_a"],
        "dimension": 4,
        "evidence_class": "operational_smoke",
    }


def environment() -> dict[str, str]:
    return {
        "CARGO_BUILD_JOBS": "1",
        "LANG": "C",
        "LC_ALL": "C",
        "MKL_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "PATH": "/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin",
        "RAYON_NUM_THREADS": "1",
        "TMPDIR": "/tmp",
        "TZ": "UTC",
        "VECLIB_MAXIMUM_THREADS": "1",
    }


def source_binding() -> dict[str, Any]:
    paths = [REPO / relative for relative in SOURCE_RELATIVE_PATHS]
    if custody.file_identity(REPO / "Cargo.lock")["sha256"] != custody.file_identity(LOCK)["sha256"]:
        raise Stage23Error("workspace Cargo.lock differs from the frozen lock")
    state = custody.git_state()
    return {
        "schema": "koblitz_unknown_scalar_source_binding.v1",
        "git": state,
        "direct_sources": {
            str(path.relative_to(REPO)): custody.file_identity(path, str(path.relative_to(REPO)))
            for path in paths
        },
        "dependency_lock": custody.file_identity(LOCK, "frozen Cargo.lock"),
        "python": custody.tool_identity(Path(sys.executable), "Python"),
        "cargo": custody.tool_identity(Path(shutil.which("cargo") or ""), "Cargo"),
        "rustc": custody.tool_identity(Path(shutil.which("rustc") or ""), "rustc"),
    }


def host_binding() -> dict[str, Any]:
    return {
        "schema": "koblitz_unknown_scalar_host_binding.v1",
        "platform": platform.platform(),
        "uname": list(platform.uname()),
        "logical_cpus": os.cpu_count(),
        "scientific_environment": environment(),
    }


def expected_inner_command(args: argparse.Namespace) -> list[str]:
    command = [
        str(Path(sys.executable).resolve()),
        str(RUNNER_SOURCE),
        "run",
        "--profile",
        args.profile,
        "--output",
        str(args.output.resolve()),
        "--meter",
        str(args.meter.resolve()),
    ]
    if args.allow_dirty:
        command.append("--allow-dirty")
    return command


def require_frozen_meter(path: Path) -> None:
    if path.resolve() != METER.resolve():
        raise Stage23Error("Stage-23 requires the bound repository process meter")


def plan(args: argparse.Namespace) -> dict[str, Any]:
    frozen = protocol(args.profile)
    require_frozen_meter(args.meter)
    output = custody.safe_new_directory(args.output, "Stage-23 planned output")
    source = source_binding()
    if args.profile == "production" and (source["git"]["dirty"] or args.allow_dirty):
        raise Stage23Error("production requires a clean checkout without --allow-dirty")
    inner = expected_inner_command(args)
    outer_root = output.with_name(output.name + "-outer")
    outer = custody.meter_command(
        meter=args.meter,
        cwd=REPO,
        timeout=float(frozen["execution"]["whole_driver_watchdog_seconds"]),
        stdout=outer_root / "driver.stdout",
        stderr=outer_root / "driver.stderr",
        metrics=outer_root / "driver.metrics.json",
        command=inner,
    )
    return {
        "schema": "koblitz_unknown_scalar_execution_plan.v1",
        "profile": args.profile,
        "output": str(output),
        "protocol_sha256": custody.file_identity(PROTOCOL)["sha256"],
        "source_revision": source["git"],
        "expected_targets": profile_values(args.profile, frozen)["targets"],
        "expected_child_processes": 4 + 2 * profile_values(args.profile, frozen)["targets"],
        "inner_command": inner,
        "outer_command": outer,
        "outer_root": str(outer_root),
    }


def run_metered(
    *, root: Path, name: str, command: list[str], timeout: float, meter: Path,
    env: dict[str, str], inputs: list[dict[str, Any]], expect_json: bool,
) -> dict[str, Any]:
    task = root / "tasks" / name
    task.mkdir(parents=True)
    stdout = task / "stdout"
    stderr = task / "stderr"
    metrics_path = task / "metrics.json"
    intent = {
        "schema": "koblitz_unknown_scalar_process_intent.v1",
        "name": name,
        "command": command,
        "watchdog_seconds": timeout,
        "environment": env,
        "inputs": inputs,
        "meter": custody.file_identity(meter, "process meter"),
        "meter_exclusive_create": True,
        "stdin": "devnull",
        "close_fds": True,
    }
    custody.write_json_new(task / "intent.json", intent)
    wrapper = custody.meter_command(
        meter=meter, cwd=REPO, timeout=timeout, stdout=stdout, stderr=stderr,
        metrics=metrics_path, command=command,
    )
    subprocess.run(wrapper, cwd=REPO, env=env, stdin=subprocess.DEVNULL, close_fds=True, check=False)
    process = custody.process_metrics(metrics_path, f"{name} metrics")
    if process["command"] != command or process["watchdog_seconds"] != timeout:
        raise Stage23Error(f"{name} metrics differ from its command or watchdog")
    process_ok = bool(
        process["returncode"] == 0
        and not process["timed_out"]
        and not process["orphan_group_terminated"]
    )
    parsed = None
    parse_error = None
    if expect_json and process_ok:
        try:
            parsed = load_json(stdout, f"{name} stdout")
            custody.write_json_new(task / "result.json", parsed)
        except Exception as error:  # retained as an incomplete task
            parse_error = str(error)
    receipt = {
        "schema": "koblitz_unknown_scalar_process_receipt.v1",
        "name": name,
        "process": process,
        "intent": custody.file_identity(task / "intent.json"),
        "stdout": custody.file_identity(stdout),
        "stderr": custody.file_identity(stderr),
        "metrics": custody.file_identity(metrics_path),
        "inputs": inputs,
        "result": custody.file_identity(task / "result.json") if parsed is not None else None,
        "parse_error": parse_error,
        "terminal_status": (
            "complete" if process_ok and (parsed is not None or not expect_json) else "incomplete"
        ),
    }
    custody.write_json_new(task / "receipt.json", receipt)
    return {"receipt": receipt, "result": parsed}


def task_resource(task: dict[str, Any]) -> dict[str, Any]:
    process = task["receipt"]["process"]
    return {
        "returncode": process["returncode"],
        "timed_out": process["timed_out"],
        "orphan_group_terminated": process["orphan_group_terminated"],
        **process["metrics"],
    }


def aggregate_resources(tasks: list[dict[str, Any]]) -> dict[str, Any]:
    resources = [task_resource(task) for task in tasks]
    return {
        "processes": len(resources),
        "summed_user_seconds": sum(row["user_seconds"] for row in resources),
        "summed_system_seconds": sum(row["system_seconds"] for row in resources),
        "total_core_seconds": sum(row["total_core_seconds"] for row in resources),
        "summed_process_wall_seconds": sum(row["wall_seconds"] for row in resources),
        "peak_process_rss_bytes": max((row["peak_rss_bytes"] for row in resources), default=0),
        "single_core_elapsed_seconds": None,
        "single_core_seconds_legacy_alias": "total_core_seconds",
        "aggregate_parallel_rss_bytes": None,
    }


def require_exact_keys(value: Any, expected: set[str], context: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != expected:
        raise Stage23Error(f"{context} uses an unexpected schema")
    return value


def validate_bound_file(
    identity: Any, expected_path: Path, run_root: Path, context: str
) -> dict[str, Any]:
    value = custody.validate_file_identity(identity, context)
    if Path(value["path"]).resolve() != expected_path.resolve():
        raise Stage23Error(f"{context} points at the wrong path")
    try:
        expected_path.resolve().relative_to(run_root)
    except ValueError as error:
        raise Stage23Error(f"{context} escapes the sealed run root") from error
    return value


def process_completed_cleanly(process: dict[str, Any]) -> bool:
    return bool(
        process["returncode"] == 0
        and not process["timed_out"]
        and not process["orphan_group_terminated"]
    )


def reconstruct_task(
    *,
    run_root: Path,
    name: str,
    command: list[str],
    watchdog: float,
    expected_environment: dict[str, str],
    inputs: list[dict[str, Any]],
    meter_identity: dict[str, Any],
    expect_json: bool,
) -> dict[str, Any]:
    task_root = run_root / "tasks" / name
    if not task_root.is_dir() or task_root.is_symlink():
        raise Stage23Error(f"missing task directory {name}")
    receipt_path = task_root / "receipt.json"
    receipt = require_exact_keys(
        load_json(receipt_path, f"{name} receipt"),
        {
            "schema",
            "name",
            "process",
            "intent",
            "stdout",
            "stderr",
            "metrics",
            "inputs",
            "result",
            "parse_error",
            "terminal_status",
        },
        f"{name} receipt",
    )
    if receipt["schema"] != "koblitz_unknown_scalar_process_receipt.v1" or receipt["name"] != name:
        raise Stage23Error(f"{name} receipt identity changed")

    intent_path = task_root / "intent.json"
    stdout_path = task_root / "stdout"
    stderr_path = task_root / "stderr"
    metrics_path = task_root / "metrics.json"
    validate_bound_file(receipt["intent"], intent_path, run_root, f"{name} intent")
    validate_bound_file(receipt["stdout"], stdout_path, run_root, f"{name} stdout")
    validate_bound_file(receipt["stderr"], stderr_path, run_root, f"{name} stderr")
    validate_bound_file(receipt["metrics"], metrics_path, run_root, f"{name} metrics")

    intent = require_exact_keys(
        load_json(intent_path, f"{name} intent"),
        {
            "schema",
            "name",
            "command",
            "watchdog_seconds",
            "environment",
            "inputs",
            "meter",
            "meter_exclusive_create",
            "stdin",
            "close_fds",
        },
        f"{name} intent",
    )
    expected_intent = {
        "schema": "koblitz_unknown_scalar_process_intent.v1",
        "name": name,
        "command": command,
        "watchdog_seconds": watchdog,
        "environment": expected_environment,
        "inputs": inputs,
        "meter": meter_identity,
        "meter_exclusive_create": True,
        "stdin": "devnull",
        "close_fds": True,
    }
    if intent != expected_intent:
        raise Stage23Error(f"{name} intent differs from the reconstructed task")
    if receipt["inputs"] != inputs:
        raise Stage23Error(f"{name} receipt inputs differ from its intent")

    raw_process = custody.process_metrics(metrics_path, f"{name} raw metrics")
    if receipt["process"] != raw_process:
        raise Stage23Error(f"{name} receipt process differs from raw metrics")
    if raw_process["command"] != command or raw_process["watchdog_seconds"] != watchdog:
        raise Stage23Error(f"{name} raw metrics differ from the reconstructed command or watchdog")
    clean = process_completed_cleanly(raw_process)

    parsed: dict[str, Any] | None = None
    expected_parse_error: str | None = None
    if expect_json and clean:
        try:
            parsed = load_json(stdout_path, f"{name} stdout")
        except Exception as error:
            expected_parse_error = str(error)

    result_path = task_root / "result.json"
    if parsed is not None:
        validate_bound_file(receipt["result"], result_path, run_root, f"{name} result")
        if load_json(result_path, f"{name} result") != parsed:
            raise Stage23Error(f"{name} result differs from parsed stdout")
    elif receipt["result"] is not None or result_path.exists() or result_path.is_symlink():
        raise Stage23Error(f"{name} retained a result for an incomplete process")

    expected_terminal = "complete" if clean and (parsed is not None or not expect_json) else "incomplete"
    if receipt["parse_error"] != expected_parse_error or receipt["terminal_status"] != expected_terminal:
        raise Stage23Error(f"{name} parse or terminal status is not derivable")
    expected_files = {"intent.json", "stdout", "stderr", "metrics.json", "receipt.json"}
    if parsed is not None:
        expected_files.add("result.json")
    actual_files = {path.name for path in task_root.iterdir()}
    if actual_files != expected_files or any(not path.is_file() or path.is_symlink() for path in task_root.iterdir()):
        raise Stage23Error(f"{name} task file inventory changed")
    return {"receipt": receipt, "result": parsed}


def validate_discovery_result(
    result: dict[str, Any], curve_a: int, values: dict[str, Any], frozen: dict[str, Any]
) -> None:
    if (
        result.get("schema") != "koblitz_public_factor_base_discovery.v1"
        or result.get("n") != values["n"]
        or result.get("a") != curve_a
        or result.get("m") != 2
        or result.get("requested_dimension") != values["dimension"]
    ):
        raise Stage23Error(f"curve-a={curve_a} discovery parameters changed")
    forbidden = result.get("forbidden_inputs")
    if (
        not isinstance(forbidden, dict)
        or set(forbidden)
        != {
            "target_constructed",
            "target_subgroup_enumerated",
            "discrete_log_labels_constructed",
            "relation_yield_used",
            "solver_timing_used",
        }
        or any(value is not False for value in forbidden.values())
    ):
        raise Stage23Error(f"curve-a={curve_a} discovery used a forbidden input")
    candidates = result.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise Stage23Error(f"curve-a={curve_a} discovery emitted no candidates")
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise Stage23Error(f"curve-a={curve_a} discovery candidate is not an object")
        if (
            not isinstance(candidate.get("divisor_indices"), list)
            or not candidate["divisor_indices"]
            or not all(type(index) is int and index >= 0 for index in candidate["divisor_indices"])
            or candidate.get("dimension") != values["dimension"]
        ):
            raise Stage23Error(f"curve-a={curve_a} discovery candidate is malformed")
        for field in (
            "divisor_polynomial",
            "abscissae",
            "rational_points",
            "signed_frobenius_orbits_before_projection",
            "projected_signed_frobenius_orbits",
        ):
            if type(candidate.get(field)) is not int or candidate[field] < 0:
                raise Stage23Error(f"curve-a={curve_a} discovery field {field} is invalid")
    if result.get("selected") not in candidates:
        raise Stage23Error(f"curve-a={curve_a} discovery selection is not a candidate")
    if curve_a == frozen["curve"]["a"]:
        selected_base = [
            candidate
            for candidate in candidates
            if candidate.get("divisor_indices") == frozen["factor_base"]["divisor_indices"]
        ]
        if not selected_base:
            raise Stage23Error("public discovery omitted the frozen K0 divisor base")
        if values["n"] == frozen["curve"]["n"]:
            candidate = selected_base[0]
            expected = frozen["factor_base"]
            if (
                candidate.get("divisor_polynomial") != expected["divisor_polynomial"]
                or candidate.get("rational_points") != expected["rational_points"]
                or candidate.get("signed_frobenius_orbits_before_projection")
                != expected["signed_frobenius_orbits_before_projection"]
                or candidate.get("projected_signed_frobenius_orbits")
                != expected["projected_signed_frobenius_columns"]
            ):
                raise Stage23Error("public discovery changed the frozen production factor base")


def validate_run_bindings(
    run_root: Path, summary: dict[str, Any], seal: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    inputs_root = run_root / "inputs"
    if not inputs_root.is_dir() or inputs_root.is_symlink():
        raise Stage23Error("run inputs directory is missing")
    if {path.name for path in inputs_root.iterdir()} != {"protocol.json", "source.json", "host.json"}:
        raise Stage23Error("run inputs inventory changed")
    protocol_path = inputs_root / "protocol.json"
    source_path = inputs_root / "source.json"
    host_path = inputs_root / "host.json"
    archived_protocol = validate_protocol_value(
        load_json(protocol_path, "archived protocol"), summary["profile"]
    )
    protocol_identity = custody.file_identity(protocol_path, "archived protocol")
    if (
        protocol_identity["sha256"] != summary.get("protocol_sha256")
        or protocol_identity["sha256"] != seal.get("protocol_sha256")
        or protocol_identity["sha256"] != custody.file_identity(PROTOCOL)["sha256"]
    ):
        raise Stage23Error("archived protocol identity changed")

    source = require_exact_keys(
        load_json(source_path, "source binding"),
        {"schema", "git", "direct_sources", "dependency_lock", "python", "cargo", "rustc"},
        "source binding",
    )
    if source["schema"] != "koblitz_unknown_scalar_source_binding.v1":
        raise Stage23Error("wrong source binding schema")
    if source["git"] != summary.get("source_revision") or source["git"].get("commit") != seal.get("source_commit"):
        raise Stage23Error("source revision binding changed")
    if source["git"] != custody.git_state():
        raise Stage23Error("checkout state differs from the bound source state")
    if set(source["direct_sources"]) != set(SOURCE_RELATIVE_PATHS):
        raise Stage23Error("direct source inventory changed")
    for relative in SOURCE_RELATIVE_PATHS:
        expected_path = (REPO / relative).resolve()
        identity = custody.validate_file_identity(source["direct_sources"][relative], relative)
        if Path(identity["path"]).resolve() != expected_path:
            raise Stage23Error(f"source binding path changed for {relative}")
    if source["dependency_lock"] != source["direct_sources"][str(LOCK.relative_to(REPO))]:
        raise Stage23Error("frozen dependency-lock binding changed")
    custody.validate_tool_identity(source["python"], "bound Python")
    custody.validate_tool_identity(source["cargo"], "bound Cargo")
    custody.validate_tool_identity(source["rustc"], "bound rustc")

    host = require_exact_keys(
        load_json(host_path, "host binding"),
        {"schema", "platform", "uname", "logical_cpus", "scientific_environment"},
        "host binding",
    )
    if host["schema"] != "koblitz_unknown_scalar_host_binding.v1":
        raise Stage23Error("wrong host binding schema")
    if host != summary.get("host") or host["scientific_environment"] != environment():
        raise Stage23Error("host or scientific environment binding changed")
    meter_identity = source["direct_sources"][str(METER.relative_to(REPO))]
    return archived_protocol, source, host, {
        "protocol": protocol_identity,
        "source": custody.file_identity(source_path, "source binding"),
        "host": custody.file_identity(host_path, "host binding"),
        "meter": meter_identity,
    }


def validate_binaries(
    run_root: Path, claimed: Any, build_complete: bool
) -> dict[str, dict[str, Any]]:
    binaries_root = run_root / "binaries"
    if not binaries_root.is_dir() or binaries_root.is_symlink():
        raise Stage23Error("run binaries directory is missing")
    names = {"koblitz_public_factor_base_discovery", "koblitz_unknown_scalar_panel"}
    if not build_complete:
        if claimed != {} or any(binaries_root.iterdir()):
            raise Stage23Error("an incomplete build retained unbound executables")
        return {}
    if not isinstance(claimed, dict) or set(claimed) != names:
        raise Stage23Error("binary identity inventory changed")
    identities: dict[str, dict[str, Any]] = {}
    for name in names:
        expected_path = (binaries_root / name).resolve()
        identity = custody.validate_tool_identity(claimed[name], name)
        if Path(identity["path"]).resolve() != expected_path:
            raise Stage23Error(f"{name} binary path changed")
        built = custody.file_identity(
            run_root / "build-target" / "release" / "examples" / name,
            f"built {name}",
        )
        if (built["bytes"], built["sha256"]) != (identity["bytes"], identity["sha256"]):
            raise Stage23Error(f"copied {name} differs from the fresh build output")
        identities[name] = identity
    if {path.name for path in binaries_root.iterdir()} != names:
        raise Stage23Error("binary directory inventory changed")
    return identities


def reconstruct_row(
    index: int,
    target: dict[str, Any],
    ic_task: dict[str, Any],
    rho_task: dict[str, Any],
    profile: str,
) -> dict[str, Any]:
    ic_validation_error = None
    rho_validation_error = None
    try:
        ic_complete = bool(
            ic_task["result"] and validate_ic(ic_task["result"], target, profile)
        )
    except Stage23Error as error:
        ic_complete = False
        ic_validation_error = str(error)
    try:
        rho_complete = bool(
            rho_task["result"] and validate_rho(rho_task["result"], target, profile)
        )
    except Stage23Error as error:
        rho_complete = False
        rho_validation_error = str(error)
    recovered_match = bool(
        ic_complete
        and rho_complete
        and ic_task["result"]["report"]["recovered_scalar"]
        == rho_task["result"]["report"]["recovered_scalar"]
    )
    return {
        "index": index,
        "target_id": target["target_id"],
        "target": target["point"],
        "ic": {
            "complete": ic_complete,
            "status": ic_task["result"].get("status") if ic_task["result"] else "process_incomplete",
            "resources": task_resource(ic_task),
            "result": ic_task["receipt"]["result"],
            "validation_error": ic_validation_error,
        },
        "rho": {
            "complete": rho_complete,
            "status": rho_task["result"].get("status") if rho_task["result"] else "process_incomplete",
            "resources": task_resource(rho_task),
            "result": rho_task["receipt"]["result"],
            "validation_error": rho_validation_error,
        },
        "recovered_scalars_match": recovered_match,
        "row_complete": ic_complete and rho_complete and recovered_match,
        "incomplete_is_inconclusive": not (ic_complete and rho_complete and recovered_match),
    }


def compute_ratios(tasks: list[dict[str, Any]], rows: list[dict[str, Any]]) -> dict[str, Any]:
    all_complete = bool(rows) and all(row["row_complete"] for row in rows)
    ic_core = sum(row["ic"]["resources"]["total_core_seconds"] for row in rows)
    rho_core = sum(row["rho"]["resources"]["total_core_seconds"] for row in rows)
    discovery_core = sum(
        task_resource(task)["total_core_seconds"]
        for task in tasks
        if "discovery-" in task["receipt"]["name"]
    )
    return {
        "online_ic_over_rho_core": ic_core / rho_core if rho_core else None,
        "setup_charged_ic_over_rho_core": (ic_core + discovery_core) / rho_core if rho_core else None,
        "verdict": "finite_complete_panel" if all_complete else "incomplete_inconclusive",
    }


def recursively_forbid_scalar_labels(value: Any) -> None:
    forbidden = {"secret", "target_scalar", "known_scalar", "planted_scalar"}
    if isinstance(value, dict):
        overlap = forbidden.intersection(value)
        if overlap:
            raise Stage23Error(f"forbidden target scalar label keys: {sorted(overlap)}")
        for child in value.values():
            recursively_forbid_scalar_labels(child)
    elif isinstance(value, list):
        for child in value:
            recursively_forbid_scalar_labels(child)


def validate_targets(result: dict[str, Any], profile: str, expected: int) -> list[dict[str, Any]]:
    if result.get("schema") != "koblitz_unknown_scalar_target_panel.v1" or result.get("status") != "complete":
        raise Stage23Error("target panel is not complete")
    if result.get("target_scalar_constructed_or_recorded") is not False:
        raise Stage23Error("target panel constructed a scalar label")
    if result.get("factor_base_log_labels_constructed_or_recorded") is not False:
        raise Stage23Error("target panel constructed factor-base log labels")
    recursively_forbid_scalar_labels(result)
    identity = result.get("identity", {})
    if custody.json_blake3(identity) != result.get("identity_blake3"):
        raise Stage23Error("target panel identity hash is invalid")
    if identity.get("profile") != profile or len(identity.get("targets", [])) != expected:
        raise Stage23Error("target panel profile or count changed")
    targets = identity["targets"]
    ids = [target.get("target_id") for target in targets]
    points = [(target.get("point", {}).get("x"), target.get("point", {}).get("y")) for target in targets]
    if len(set(ids)) != expected or len(set(points)) != expected:
        raise Stage23Error("target panel contains duplicate identities or points")
    return targets


def validate_ic(
    result: dict[str, Any], target: dict[str, Any], profile: str | None = None
) -> bool:
    if result.get("schema") != "koblitz_unknown_scalar_ic_result.v1":
        raise Stage23Error("wrong IC result schema")
    if result.get("status") not in {"complete_verified", "incomplete"}:
        raise Stage23Error("unknown IC terminal status")
    if profile is not None and result.get("profile") != profile:
        raise Stage23Error("IC profile mismatch")
    if result.get("target_id") != target["target_id"] or result.get("target") != target["point"]:
        raise Stage23Error("IC target mismatch")
    if result.get("seed") != target.get("ic_seed"):
        raise Stage23Error("IC seed differs from the frozen target-panel seed")
    if result.get("target_scalar_constructed_or_supplied") is not False:
        raise Stage23Error("IC process received a scalar label")
    if result.get("factor_base_logs_constructed_or_supplied") is not False:
        raise Stage23Error("IC process used factor-base log labels")
    report = result.get("report", {})
    attempts = result.get("attempt_records", [])
    matrix = result.get("relation_matrix", [])
    if report.get("trials") != len(attempts) or report.get("relations") != len(matrix):
        raise Stage23Error("IC attempt or matrix inventory is incomplete")
    if sum(report.get("outcomes", {}).values()) != len(attempts):
        raise Stage23Error("IC outcome partition does not cover every trial")
    if [row.get("trial") for row in attempts] != list(range(1, len(attempts) + 1)):
        raise Stage23Error("IC attempt sequence is not exact and contiguous")
    outcome_counts = {name: 0 for name in report["outcomes"]}
    for row in attempts:
        disposition = row.get("disposition")
        if disposition not in outcome_counts:
            raise Stage23Error("IC attempt uses an unknown disposition")
        outcome_counts[disposition] += 1
    if outcome_counts != report["outcomes"]:
        raise Stage23Error("IC outcome counts differ from attempt records")
    if outcome_counts.get("relation_found") != len(matrix):
        raise Stage23Error("IC relation count differs from relation-found attempts")
    if sum(row.get("conflicts", 0) for row in attempts) != report.get("sat_conflicts"):
        raise Stage23Error("IC conflict total differs from attempt records")
    if sum(row.get("solver_calls", 0) for row in attempts) != report.get("sat_calls"):
        raise Stage23Error("IC solver-call total differs from attempt records")
    if sum(row.get("models", 0) for row in attempts) != report.get("sat_models"):
        raise Stage23Error("IC model total differs from attempt records")
    if custody.json_blake3(attempts) != result.get("attempt_records_blake3"):
        raise Stage23Error("IC attempt-record hash is invalid")
    if custody.json_blake3(matrix) != result.get("relation_matrix_blake3"):
        raise Stage23Error("IC relation-matrix hash is invalid")
    progress = result.get("progress", [])
    if custody.json_blake3(progress) != result.get("progress_blake3"):
        raise Stage23Error("IC progress hash is invalid")
    attempt_events = [row for row in progress if row.get("event") == "relation_attempt_finished"]
    if len(attempt_events) != len(attempts):
        raise Stage23Error("IC progress omits relation attempts")
    if [row.get("trial") for row in attempt_events] != list(range(1, len(attempts) + 1)):
        raise Stage23Error("IC progress attempt sequence is not exact")
    rank_history = report.get("rank_history", [])
    rank_events = [row for row in progress if row.get("event") == "matrix_rank"]
    if (
        not rank_history
        or report.get("rank_checks") != len(rank_history)
        or report.get("linear_solve_attempts") != len(rank_history)
        or len(rank_events) != len(rank_history)
        or report.get("matrix_rows") != len(matrix)
        or type(report.get("matrix_columns")) is not int
        or type(report.get("terminal_matrix_rank")) is not int
        or not 0
        <= report["terminal_matrix_rank"]
        <= min(report["matrix_rows"], report["matrix_columns"])
    ):
        raise Stage23Error("IC rank history or progress is incomplete")
    terminal = rank_history[-1]
    terminal_event = rank_events[-1]
    for field in ("rows", "columns", "rank"):
        report_field = "terminal_matrix_rank" if field == "rank" else f"matrix_{field}"
        if terminal.get(field) != report.get(report_field) or terminal_event.get(field) != report.get(report_field):
            raise Stage23Error("IC terminal rank evidence is inconsistent")
    if result.get("status") == "complete_verified":
        if not report.get("recovered_scalar_point_verified") or report.get("recovered_scalar") is None:
            raise Stage23Error("completed IC row lacks point verification")
        if report.get("direct_relation") or report.get("sat_invalid_models") != 0:
            raise Stage23Error("completed IC row bypassed the matrix or admitted an invalid model")
        if (
            report.get("matrix_rows") != len(matrix)
            or report.get("terminal_matrix_rank") != report.get("matrix_columns")
            or not report.get("rank_history")
            or report["rank_history"][-1].get("candidate_verified") is not True
        ):
            raise Stage23Error("completed IC row lacks matrix/rank evidence")
        return True
    if report.get("recovered_scalar_point_verified") or report.get("recovered_scalar") is not None:
        raise Stage23Error("incomplete IC row retained a completed scalar claim")
    return False


def validate_rho(
    result: dict[str, Any], target: dict[str, Any], profile: str | None = None
) -> bool:
    if result.get("schema") != "koblitz_unknown_scalar_rho_result.v1":
        raise Stage23Error("wrong rho result schema")
    if result.get("status") not in {"complete_verified", "incomplete"}:
        raise Stage23Error("unknown rho terminal status")
    actual_profile = result.get("profile")
    if profile is not None and actual_profile != profile:
        raise Stage23Error("rho profile mismatch")
    if actual_profile not in {"production", "smoke"}:
        raise Stage23Error("rho profile is unknown")
    if result.get("target_id") != target["target_id"] or result.get("target") != target["point"]:
        raise Stage23Error("rho target mismatch")
    if result.get("seed") != target.get("rho_seed"):
        raise Stage23Error("rho seed differs from the frozen target-panel seed")
    if result.get("target_scalar_constructed_or_supplied") is not False:
        raise Stage23Error("rho process received a scalar label")
    report = result.get("report", {})
    charges = result.get("charges", {})
    if report.get("jump_table_rebuilds") != report.get("restarts_attempted"):
        raise Stage23Error("rho restart and jump-table counts differ")
    if charges.get("setup_scalar_multiplications") != 34 * report.get("jump_table_rebuilds", -1):
        raise Stage23Error("rho setup scalar-multiplication ledger does not balance")
    if charges.get("setup_group_additions") != 17 * report.get("jump_table_rebuilds", -1):
        raise Stage23Error("rho setup addition ledger does not balance")
    if charges.get("coefficient_draws") != charges.get("setup_scalar_multiplications"):
        raise Stage23Error("rho coefficient-draw ledger does not balance")
    if charges.get("walk_group_additions") != charges.get("partition_hashes"):
        raise Stage23Error("rho walk-addition ledger does not balance")
    if charges.get("walk_group_additions") != 3 * report.get("iterations", -1):
        raise Stage23Error("rho Floyd-walk addition ledger does not balance")
    if charges.get("canonicalizations") != (
        report.get("restarts_attempted", -1) + charges.get("walk_group_additions", -1)
    ):
        raise Stage23Error("rho canonicalization ledger does not balance")
    curve_n = 23 if actual_profile == "production" else 7
    if (
        charges.get("frobenius_maps") != charges.get("negations_examined")
        or charges.get("frobenius_maps", 0) > charges.get("canonicalizations", 0) * curve_n
        or charges.get("frobenius_maps", 0) % curve_n != 0
    ):
        raise Stage23Error("rho automorphism ledger does not balance")
    if charges.get("collisions", 0) < charges.get("failed_collisions", 0):
        raise Stage23Error("rho collision ledger does not balance")
    progress = result.get("progress", [])
    if custody.json_blake3(progress) != result.get("progress_blake3"):
        raise Stage23Error("rho progress hash is invalid")
    if sum(row.get("event") == "rho_restart_started" for row in progress) != report.get("restarts_attempted"):
        raise Stage23Error("rho progress omits restarts")
    if sum(row.get("event") == "rho_jump_table_ready" for row in progress) != report.get("jump_table_rebuilds"):
        raise Stage23Error("rho progress omits jump tables")
    if sum(row.get("event") == "rho_collision" for row in progress) != charges.get("collisions"):
        raise Stage23Error("rho progress omits collisions")
    if not progress or progress[-1].get("event") != "rho_finished":
        raise Stage23Error("rho progress lacks its terminal event")
    timing = result.get("timing_ns", {})
    timing_fields = (
        "target_and_subgroup_validation",
        "rho_setup",
        "rho_walk",
        "candidate_verification",
        "end_to_end",
    )
    if any(type(timing.get(field)) is not int or timing[field] < 0 for field in timing_fields):
        raise Stage23Error("rho timing ledger is incomplete")
    if sum(timing[field] for field in timing_fields[:-1]) > timing["end_to_end"]:
        raise Stage23Error("rho timing stages overlap or exceed end-to-end time")
    if result.get("status") == "complete_verified":
        if not report.get("recovered_scalar_point_verified") or report.get("recovered_scalar") is None:
            raise Stage23Error("completed rho row lacks point verification")
        if report.get("exhausted") is not False:
            raise Stage23Error("completed rho row is marked exhausted")
        return True
    if (
        report.get("recovered_scalar_point_verified")
        or report.get("recovered_scalar") is not None
        or report.get("exhausted") is not True
    ):
        raise Stage23Error("incomplete rho row retained a completed scalar claim")
    return False


def finish_run(
    *, output: Path, profile: str, frozen: dict[str, Any], source: dict[str, Any],
    host: dict[str, Any], tasks: list[dict[str, Any]], rows: list[dict[str, Any]],
    targets: dict[str, Any] | None, binaries: dict[str, Any], expected_inner: list[str],
    critical_failure: str | None,
) -> dict[str, Any]:
    resources = aggregate_resources(tasks)
    all_complete = bool(rows) and len(rows) == profile_values(profile, frozen)["targets"] and all(row["row_complete"] for row in rows)
    ratios = compute_ratios(tasks, rows)
    summary = {
        "schema": RUN_SUMMARY_SCHEMA,
        "status": "complete_verified_panel" if all_complete else "incomplete_panel",
        "profile": profile,
        "evidence_class": profile_values(profile, frozen)["evidence_class"],
        "protocol_sha256": custody.file_identity(PROTOCOL)["sha256"],
        "source_revision": source["git"],
        "host": host,
        "critical_failure": critical_failure,
        "target_panel": targets,
        "binaries": binaries,
        "rows": rows,
        "completed_rows": sum(row["row_complete"] for row in rows),
        "expected_rows": profile_values(profile, frozen)["targets"],
        "resources": resources,
        "ratios": ratios,
        "full_cost_gate_passed": False,
        "independent_external_reproduction_satisfied": False,
        "koblitz_index_calculus_sota": False,
        "outer_driver_accounting": "required_before_admission",
        "claim_boundary": frozen["claim_boundary"],
    }
    custody.write_json_new(output / "run-summary.json", summary)
    frozen_inventory = custody.inventory(output, {"run-seal.json"})
    seal_payload = {
        "schema": RUN_SEAL_SCHEMA,
        "status": "outputs_frozen",
        "profile": profile,
        "protocol_sha256": summary["protocol_sha256"],
        "source_commit": source["git"]["commit"],
        "summary": custody.file_identity(output / "run-summary.json"),
        "inventory": frozen_inventory,
        "inventory_sha256": custody.canonical_sha256(frozen_inventory),
        "expected_inner_command": expected_inner,
        "panel_complete": all_complete,
        "scientific_measurement_admitted": False,
    }
    seal = dict(seal_payload)
    seal["seal_payload_sha256"] = custody.canonical_sha256(seal_payload)
    custody.write_json_new(output / "run-seal.json", seal)
    return seal


def run(args: argparse.Namespace) -> dict[str, Any]:
    frozen = protocol(args.profile)
    require_frozen_meter(args.meter)
    values = profile_values(args.profile, frozen)
    output = custody.safe_new_directory(args.output, "Stage-23 run output")
    output.mkdir(parents=True)
    (output / "tasks").mkdir()
    (output / "binaries").mkdir()
    (output / "inputs").mkdir()
    source = source_binding()
    if args.profile == "production" and (source["git"]["dirty"] or args.allow_dirty):
        raise Stage23Error("production requires a clean checkout without --allow-dirty")
    host = host_binding()
    shutil.copyfile(PROTOCOL, output / "inputs/protocol.json")
    custody.write_json_new(output / "inputs/source.json", source)
    custody.write_json_new(output / "inputs/host.json", host)
    env = environment()
    tasks: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    binaries: dict[str, Any] = {}
    targets_result: dict[str, Any] | None = None
    critical_failure = None
    target_dir = output / "build-target"
    cargo = source["cargo"]["path"]
    build = run_metered(
        root=output,
        name="00-build",
        command=[cargo, "build", "--release", "--locked", "--jobs", "1", "--target-dir", str(target_dir), "--example", "koblitz_public_factor_base_discovery", "--example", "koblitz_unknown_scalar_panel"],
        timeout=float(frozen["execution"]["build_watchdog_seconds"]),
        meter=args.meter,
        env=env,
        inputs=[custody.file_identity(output / "inputs/source.json"), custody.file_identity(output / "inputs/protocol.json")],
        expect_json=False,
    )
    tasks.append(build)
    if build["receipt"]["terminal_status"] != "complete":
        critical_failure = "build"
    else:
        for name in ("koblitz_public_factor_base_discovery", "koblitz_unknown_scalar_panel"):
            source_binary = target_dir / "release/examples" / name
            destination = output / "binaries" / name
            shutil.copyfile(source_binary, destination)
            destination.chmod(0o755)
            binaries[name] = custody.tool_identity(destination, name)

    if critical_failure is None:
        for curve_a in (0, 1):
            task = run_metered(
                root=output,
                name=f"0{curve_a + 1}-discovery-a{curve_a}",
                command=[binaries["koblitz_public_factor_base_discovery"]["path"], str(values["n"]), str(curve_a), "2", str(values["dimension"])],
                timeout=float(frozen["execution"]["discovery_watchdog_seconds"]),
                meter=args.meter,
                env=env,
                inputs=[custody.file_identity(output / "inputs/protocol.json")],
                expect_json=True,
            )
            tasks.append(task)
            if task["receipt"]["terminal_status"] != "complete":
                critical_failure = f"discovery-a{curve_a}"
                break
            try:
                validate_discovery_result(task["result"], curve_a, values, frozen)
            except Stage23Error:
                critical_failure = f"discovery-a{curve_a}-invalid"
                break

    if critical_failure is None:
        target_task = run_metered(
            root=output,
            name="03-targets",
            command=[binaries["koblitz_unknown_scalar_panel"]["path"], "targets", args.profile],
            timeout=float(frozen["execution"]["target_generation_watchdog_seconds"]),
            meter=args.meter,
            env=env,
            inputs=[custody.file_identity(output / "inputs/protocol.json")],
            expect_json=True,
        )
        tasks.append(target_task)
        if target_task["receipt"]["terminal_status"] != "complete":
            critical_failure = "targets"
        else:
            targets_result = target_task["result"]
            try:
                validate_targets(targets_result, args.profile, values["targets"])
            except Stage23Error:
                critical_failure = "targets-invalid"

    if critical_failure is None and targets_result is not None:
        targets = targets_result["identity"]["targets"]
        for index, target in enumerate(targets, 1):
            common = [target["target_id"], target["point"]["x"], target["point"]["y"]]
            ic_task = run_metered(
                root=output,
                name=f"{2 * index + 2:02d}-row-{index:02d}-ic",
                command=[binaries["koblitz_unknown_scalar_panel"]["path"], "ic", args.profile, *common, str(target["ic_seed"])],
                timeout=float(frozen["execution"]["ic_watchdog_seconds"]),
                meter=args.meter,
                env=env,
                inputs=[custody.file_identity(output / "tasks/03-targets/result.json")],
                expect_json=True,
            )
            tasks.append(ic_task)
            rho_task = run_metered(
                root=output,
                name=f"{2 * index + 3:02d}-row-{index:02d}-rho",
                command=[binaries["koblitz_unknown_scalar_panel"]["path"], "rho", args.profile, *common, str(target["rho_seed"])],
                timeout=float(frozen["execution"]["rho_watchdog_seconds"]),
                meter=args.meter,
                env=env,
                inputs=[custody.file_identity(output / "tasks/03-targets/result.json")],
                expect_json=True,
            )
            tasks.append(rho_task)
            rows.append(reconstruct_row(index, target, ic_task, rho_task, args.profile))

    return finish_run(
        output=output,
        profile=args.profile,
        frozen=frozen,
        source=source,
        host=host,
        tasks=tasks,
        rows=rows,
        targets=targets_result,
        binaries=binaries,
        expected_inner=expected_inner_command(args),
        critical_failure=critical_failure,
    )


def reconstruct_run(
    run_root: Path, summary: dict[str, Any], seal: dict[str, Any]
) -> dict[str, Any]:
    require_exact_keys(
        summary,
        {
            "schema",
            "status",
            "profile",
            "evidence_class",
            "protocol_sha256",
            "source_revision",
            "host",
            "critical_failure",
            "target_panel",
            "binaries",
            "rows",
            "completed_rows",
            "expected_rows",
            "resources",
            "ratios",
            "full_cost_gate_passed",
            "independent_external_reproduction_satisfied",
            "koblitz_index_calculus_sota",
            "outer_driver_accounting",
            "claim_boundary",
        },
        "run summary",
    )
    if summary["schema"] != RUN_SUMMARY_SCHEMA or summary["profile"] not in {"production", "smoke"}:
        raise Stage23Error("wrong run summary schema or profile")
    frozen, source, _host, bindings = validate_run_bindings(run_root, summary, seal)
    values = profile_values(summary["profile"], frozen)
    if summary["expected_rows"] != values["targets"]:
        raise Stage23Error("summary expected-row count changed")

    full_names = ["00-build", "01-discovery-a0", "02-discovery-a1", "03-targets"]
    for index in range(1, values["targets"] + 1):
        full_names.extend(
            [f"{2 * index + 2:02d}-row-{index:02d}-ic", f"{2 * index + 3:02d}-row-{index:02d}-rho"]
        )
    tasks_root = run_root / "tasks"
    if not tasks_root.is_dir() or tasks_root.is_symlink():
        raise Stage23Error("run tasks directory is missing")
    children = list(tasks_root.iterdir())
    if any(not path.is_dir() or path.is_symlink() for path in children):
        raise Stage23Error("tasks root contains a non-directory entry")
    actual_names = sorted(path.name for path in children)
    if not actual_names or actual_names != full_names[: len(actual_names)]:
        raise Stage23Error("task directories are not an exact execution prefix")

    expected_environment = environment()
    tasks: list[dict[str, Any]] = []
    task_by_name: dict[str, dict[str, Any]] = {}

    def add_task(
        name: str,
        command: list[str],
        watchdog: float,
        inputs: list[dict[str, Any]],
        expect_json: bool,
    ) -> dict[str, Any]:
        if name not in actual_names:
            raise Stage23Error(f"sealed run omitted required task {name}")
        task = reconstruct_task(
            run_root=run_root,
            name=name,
            command=command,
            watchdog=watchdog,
            expected_environment=expected_environment,
            inputs=inputs,
            meter_identity=bindings["meter"],
            expect_json=expect_json,
        )
        tasks.append(task)
        task_by_name[name] = task
        return task

    build = add_task(
        "00-build",
        [
            source["cargo"]["path"],
            "build",
            "--release",
            "--locked",
            "--jobs",
            "1",
            "--target-dir",
            str((run_root / "build-target").resolve()),
            "--example",
            "koblitz_public_factor_base_discovery",
            "--example",
            "koblitz_unknown_scalar_panel",
        ],
        float(frozen["execution"]["build_watchdog_seconds"]),
        [bindings["source"], bindings["protocol"]],
        False,
    )
    build_complete = build["receipt"]["terminal_status"] == "complete"
    binaries = validate_binaries(run_root, summary["binaries"], build_complete)
    critical_failure: str | None = None
    targets_result: dict[str, Any] | None = None
    targets: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []

    if not build_complete:
        critical_failure = "build"
    else:
        for curve_a in (0, 1):
            name = f"0{curve_a + 1}-discovery-a{curve_a}"
            discovery = add_task(
                name,
                [
                    binaries["koblitz_public_factor_base_discovery"]["path"],
                    str(values["n"]),
                    str(curve_a),
                    "2",
                    str(values["dimension"]),
                ],
                float(frozen["execution"]["discovery_watchdog_seconds"]),
                [bindings["protocol"]],
                True,
            )
            if discovery["receipt"]["terminal_status"] != "complete":
                critical_failure = f"discovery-a{curve_a}"
                break
            try:
                validate_discovery_result(discovery["result"], curve_a, values, frozen)
            except Stage23Error:
                critical_failure = f"discovery-a{curve_a}-invalid"
                break

    if critical_failure is None:
        target_task = add_task(
            "03-targets",
            [binaries["koblitz_unknown_scalar_panel"]["path"], "targets", summary["profile"]],
            float(frozen["execution"]["target_generation_watchdog_seconds"]),
            [bindings["protocol"]],
            True,
        )
        if target_task["receipt"]["terminal_status"] != "complete":
            critical_failure = "targets"
        else:
            targets_result = target_task["result"]
            try:
                targets = validate_targets(targets_result, summary["profile"], values["targets"])
            except Stage23Error:
                critical_failure = "targets-invalid"

    if critical_failure is None:
        target_identity = task_by_name["03-targets"]["receipt"]["result"]
        for index, target in enumerate(targets, 1):
            common = [target["target_id"], target["point"]["x"], target["point"]["y"]]
            ic_name = f"{2 * index + 2:02d}-row-{index:02d}-ic"
            rho_name = f"{2 * index + 3:02d}-row-{index:02d}-rho"
            ic_task = add_task(
                ic_name,
                [
                    binaries["koblitz_unknown_scalar_panel"]["path"],
                    "ic",
                    summary["profile"],
                    *common,
                    str(target["ic_seed"]),
                ],
                float(frozen["execution"]["ic_watchdog_seconds"]),
                [target_identity],
                True,
            )
            rho_task = add_task(
                rho_name,
                [
                    binaries["koblitz_unknown_scalar_panel"]["path"],
                    "rho",
                    summary["profile"],
                    *common,
                    str(target["rho_seed"]),
                ],
                float(frozen["execution"]["rho_watchdog_seconds"]),
                [target_identity],
                True,
            )
            rows.append(reconstruct_row(index, target, ic_task, rho_task, summary["profile"]))

    if actual_names != [task["receipt"]["name"] for task in tasks]:
        raise Stage23Error("task graph continued after a critical failure or omitted a required task")
    resources = aggregate_resources(tasks)
    ratios = compute_ratios(tasks, rows)
    all_complete = (
        len(rows) == values["targets"] and all(row["row_complete"] for row in rows)
    )
    expected_status = "complete_verified_panel" if all_complete else "incomplete_panel"
    if summary["critical_failure"] != critical_failure:
        raise Stage23Error("summary critical-failure state is not derivable")
    if summary["target_panel"] != targets_result:
        raise Stage23Error("summary target panel differs from the target task")
    if summary["rows"] != rows:
        raise Stage23Error("summary rows differ from reconstructed task results")
    if summary["resources"] != resources or summary["ratios"] != ratios:
        raise Stage23Error("summary resources or ratios are not derivable")
    if (
        summary["status"] != expected_status
        or summary["completed_rows"] != sum(row["row_complete"] for row in rows)
        or summary["evidence_class"] != values["evidence_class"]
        or summary["claim_boundary"] != frozen["claim_boundary"]
        or summary["outer_driver_accounting"] != "required_before_admission"
        or summary["full_cost_gate_passed"] is not False
        or summary["independent_external_reproduction_satisfied"] is not False
        or summary["koblitz_index_calculus_sota"] is not False
    ):
        raise Stage23Error("summary panel state or claim boundary is not derivable")
    if (
        seal.get("profile") != summary["profile"]
        or seal.get("panel_complete") is not all_complete
        or seal.get("scientific_measurement_admitted") is not False
    ):
        raise Stage23Error("run seal panel state is not derivable")

    runner_path = source["direct_sources"][str(RUNNER_SOURCE.relative_to(REPO))]["path"]
    base_inner = [
        source["python"]["path"],
        runner_path,
        "run",
        "--profile",
        summary["profile"],
        "--output",
        str(run_root),
        "--meter",
        bindings["meter"]["path"],
    ]
    allowed_inner = [base_inner]
    if summary["profile"] == "smoke":
        allowed_inner.append([*base_inner, "--allow-dirty"])
    if seal.get("expected_inner_command") not in allowed_inner:
        raise Stage23Error("sealed inner command is not reconstructible")
    if summary["profile"] == "production" and summary["source_revision"]["dirty"]:
        raise Stage23Error("production run used a dirty source state")
    return {
        "frozen": frozen,
        "source": source,
        "tasks": tasks,
        "rows": rows,
        "resources": resources,
        "all_complete": all_complete,
    }


def verify(args: argparse.Namespace) -> dict[str, Any]:
    run_root = args.run_root.resolve()
    seal = require_exact_keys(
        load_json(run_root / "run-seal.json", "run seal"),
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
        "run seal",
    )
    if seal.get("schema") != RUN_SEAL_SCHEMA or seal.get("status") != "outputs_frozen":
        raise Stage23Error("run seal is not terminal")
    payload = dict(seal)
    claimed_payload = payload.pop("seal_payload_sha256", None)
    if claimed_payload != custody.canonical_sha256(payload):
        raise Stage23Error("run seal payload hash is invalid")
    actual_inventory = custody.inventory(run_root, {"run-seal.json"})
    if actual_inventory != seal["inventory"] or custody.canonical_sha256(actual_inventory) != seal["inventory_sha256"]:
        raise Stage23Error("run inventory changed after sealing")
    validate_bound_file(seal["summary"], run_root / "run-summary.json", run_root, "run summary")
    summary = load_json(run_root / "run-summary.json", "run summary")
    reconstructed = reconstruct_run(run_root, summary, seal)
    outer = custody.process_metrics(args.outer_metrics, "outer metrics")
    if outer["command"] != seal["expected_inner_command"]:
        raise Stage23Error("outer receipt does not enclose the frozen inner command")
    if outer["watchdog_seconds"] != float(reconstructed["frozen"]["execution"]["whole_driver_watchdog_seconds"]):
        raise Stage23Error("outer receipt watchdog differs from the frozen protocol")
    if not process_completed_cleanly(outer):
        raise Stage23Error("outer driver did not complete cleanly")
    child_resources = reconstructed["resources"]
    if outer["metrics"]["total_core_seconds"] + 1e-9 < child_resources["total_core_seconds"]:
        raise Stage23Error("outer CPU does not enclose child CPU")
    if outer["metrics"]["wall_seconds"] + 1e-9 < child_resources["summed_process_wall_seconds"]:
        raise Stage23Error("outer wall does not enclose sequential child wall")
    if outer["metrics"]["peak_rss_bytes"] < child_resources["peak_process_rss_bytes"]:
        raise Stage23Error("outer RSS does not enclose child RSS")
    output = custody.safe_new_directory(args.output, "Stage-23 verification output")
    output.mkdir(parents=True)
    production = summary["profile"] == "production"
    verification = {
        "schema": VERIFICATION_SCHEMA,
        "status": "run_custody_and_derivable_relationships_verified",
        "production": production,
        "run_root": str(run_root),
        "run_seal": custody.file_identity(run_root / "run-seal.json"),
        "run_inventory_sha256": seal["inventory_sha256"],
        "outer_metrics": custody.file_identity(args.outer_metrics),
        "outer_driver": outer,
        "panel_status": summary["status"],
        "completed_rows": summary["completed_rows"],
        "expected_rows": summary["expected_rows"],
        "resources": {"charged_children": child_resources, "whole_driver": outer["metrics"]},
        "scientific_measurement_admitted": False,
        "measurement_admission_status": "pending_independent_mathematical_payload_replay" if production else "operational_smoke_only",
        "pending_independent_payload_replay": [
            "domain-separated target generation and subgroup membership",
            "public factor-base discovery and materialization",
            "every relation-attempt target, decomposition witness, and relation row",
            "relation-matrix rank and modular linear algebra",
            "IC and rho recovered-scalar point verification",
            "signed-Frobenius rho operation ledger",
        ],
        "external_portable_verification_satisfied": False,
        "external_portable_verification_blocker": "run identities retain absolute checkout paths; no archive-relative source and executable replay exists",
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "claim_boundary": summary["claim_boundary"],
    }
    custody.write_json_new(output / "verification.json", verification)
    inventory = custody.inventory(output, {"verification-seal.json"})
    seal_payload = {
        "schema": VERIFICATION_SEAL_SCHEMA,
        "status": "verification_frozen",
        "verification": custody.file_identity(output / "verification.json"),
        "inventory": inventory,
        "inventory_sha256": custody.canonical_sha256(inventory),
    }
    verification_seal = dict(seal_payload)
    verification_seal["seal_payload_sha256"] = custody.canonical_sha256(seal_payload)
    custody.write_json_new(output / "verification-seal.json", verification_seal)
    return verification_seal


def self_test() -> dict[str, Any]:
    frozen = protocol("production")
    assert profile_values("production", frozen)["targets"] == 5
    assert profile_values("smoke", frozen)["targets"] == 2
    try:
        recursively_forbid_scalar_labels({"nested": [{"target_scalar": 1}]})
    except Stage23Error:
        pass
    else:
        raise AssertionError("forbidden scalar label was accepted")
    return {"self_test": "pass", "checks": 3, "protocol_sha256": custody.file_identity(PROTOCOL)["sha256"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ("plan", "run"):
        child = subparsers.add_parser(name)
        child.add_argument("--profile", choices=("production", "smoke"), default="production")
        child.add_argument("--output", type=Path, required=True)
        child.add_argument("--meter", type=Path, default=METER)
        child.add_argument("--allow-dirty", action="store_true")
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--run-root", type=Path, required=True)
    verify_parser.add_argument("--outer-metrics", type=Path, required=True)
    verify_parser.add_argument("--output", type=Path, required=True)
    subparsers.add_parser("self-test")
    args = parser.parse_args()
    try:
        if args.command == "plan":
            result = plan(args)
        elif args.command == "run":
            result = run(args)
        elif args.command == "verify":
            result = verify(args)
        else:
            result = self_test()
        print(json.dumps(result, indent=2, sort_keys=True))
    except (Stage23Error, custody.Stage21Error, OSError, ValueError, KeyError) as error:
        parser.exit(1, f"stage23: {error}\n")


if __name__ == "__main__":
    main()
