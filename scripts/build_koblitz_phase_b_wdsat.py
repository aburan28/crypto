#!/usr/bin/env python3
"""Build, meter, seal, and validate the frozen Phase-B WDSat binary."""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
from typing import Any

import run_koblitz_blind_pdp_phase_b as phase_b


BUILD_RECEIPT_SCHEMA = "koblitz_pdp_phase_b_wdsat_build_receipt.v2"
BUILD_SEAL_SCHEMA = "koblitz_pdp_phase_b_wdsat_build_seal.v2"
SOURCE_INVENTORY_SCHEMA = "koblitz_pdp_phase_b_wdsat_source_inventory.v2"
RECIPE_SCHEMA = "koblitz_pdp_phase_b_wdsat_build_recipe.v1"
EXPECTED_ROLES = ["source-copy", "config-install", "clean", "build", "binary-copy"]
CLAIM_BOUNDARY = (
    "Clean source-copy, config-install, clean, build, and binary-copy costs only; "
    "exact blind-target export sizing is charged and checked separately before solving"
)
FROZEN_RECIPE = {
    "schema": RECIPE_SCHEMA,
    "one_process_at_a_time": True,
    "stages": [
        {"role": "source-copy", "tool": "cp", "cwd": "capsule", "arguments": ["-R", "source/src", "capsule/build/src"]},
        {"role": "config-install", "tool": "cp", "cwd": "capsule", "arguments": ["frozen-config", "capsule/build/src/config.h"]},
        {"role": "clean", "tool": "make", "cwd": "capsule/build", "arguments": ["-C", "src", "clean"]},
        {"role": "build", "tool": "make", "cwd": "capsule/build", "arguments": ["-C", "src"]},
        {"role": "binary-copy", "tool": "cp", "cwd": "capsule", "arguments": ["capsule/build/wdsat_solver", "capsule/wdsat_solver"]},
    ],
}

RECEIPT_FIELDS = {
    "schema", "status", "protocol_payload_sha256", "source_commit", "source_dirty",
    "config", "config_sha256", "binary_identity", "binary_sha256", "limits",
    "build_context", "child_environment", "frozen_recipe",
    "configured_source_inventory_path", "build_driver_identity",
    "process_meter_identity", "toolchain", "build_processes", "resources",
    "claim_boundary", "receipt_payload_sha256",
}
SEAL_FIELDS = {
    "schema", "status", "protocol_payload_sha256", "source_commit", "config_sha256",
    "binary_path", "binary_identity", "binary_sha256", "receipt_path", "receipt_sha256",
    "implementation_state", "inventory", "inventory_sha256", "seal_payload_sha256",
}
PROCESS_FIELDS = {
    "command", "returncode", "watchdog_seconds", "timed_out", "orphan_group_terminated",
    "metrics", "role", "intent_path", "stdout_path", "stderr_path", "metrics_path",
    "stdout_bytes", "stdout_sha256", "stderr_bytes", "stderr_sha256", "inputs",
}
INTENT_FIELDS = {
    "schema", "role", "command", "cwd", "watchdog_seconds", "stdin", "close_fds",
    "environment", "inputs", "executable_sha256",
}
METRICS_FIELDS = {
    "command", "returncode", "watchdog_seconds", "timed_out",
    "orphan_group_terminated", "metrics",
}
RESOURCE_FIELDS = {
    "wall_seconds", "user_seconds", "system_seconds", "total_core_seconds",
    "single_core_seconds", "peak_rss_bytes", "meter",
}


def git_state(source: Path) -> tuple[str, list[str]]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=source, text=True,
        capture_output=True, check=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"], cwd=source,
        text=True, capture_output=True, check=True,
    ).stdout.splitlines()
    phase_b.require_hex40(commit, "WDSat source commit")
    return commit, status


def expected_configured_source_inventory(source: Path, config_bytes: bytes) -> list[dict[str, Any]]:
    completed = subprocess.run(
        ["git", "ls-files", "-z", "--", "src"],
        cwd=source,
        capture_output=True,
        check=True,
    )
    paths = [item.decode() for item in completed.stdout.split(b"\0") if item]
    inventory = []
    for path_text in paths:
        relative = Path(path_text).relative_to("src")
        data = config_bytes if relative.as_posix() == "config.h" else phase_b.regular_file_bytes(
            source / path_text, "tracked WDSat source input"
        )
        inventory.append(
            {"path": relative.as_posix(), "bytes": len(data), "sha256": phase_b.sha256_bytes(data)}
        )
    if not any(item["path"] == "config.h" for item in inventory):
        inventory.append(
            {"path": "config.h", "bytes": len(config_bytes), "sha256": phase_b.sha256_bytes(config_bytes)}
        )
    inventory.sort(key=lambda item: item["path"])
    return inventory


def current_implementation_state() -> dict[str, Any]:
    return phase_b.git_state()


def validate_implementation_state(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != {"commit", "dirty", "porcelain"}:
        raise phase_b.PhaseBError(f"{context} uses an invalid implementation-state schema")
    phase_b.require_hex40(value.get("commit"), f"{context} commit")
    if not isinstance(value.get("dirty"), bool):
        raise phase_b.PhaseBError(f"{context} dirty flag is invalid")
    porcelain = value.get("porcelain")
    if not isinstance(porcelain, list) or not all(isinstance(line, str) for line in porcelain):
        raise phase_b.PhaseBError(f"{context} porcelain status is invalid")
    if value["dirty"] != bool(porcelain):
        raise phase_b.PhaseBError(f"{context} dirty flag differs from its porcelain status")
    return value


def parse_config(config: str) -> dict[str, int]:
    macros: dict[str, int] = {}
    for line in config.splitlines():
        match = re.fullmatch(r"\s*#define\s+(\S+)\s+(\d+)\s*", line)
        if match is not None:
            macros[match.group(1)] = int(match.group(2))
    names = {
        "max_anf_id": "__MAX_ANF_ID__", "max_degree": "__MAX_DEGREE__",
        "max_id": "__MAX_ID__", "max_buffer_size": "__MAX_BUFFER_SIZE__",
        "max_eq": "__MAX_EQ__", "max_eq_size": "__MAX_EQ_SIZE__",
        "max_xeq": "__MAX_XEQ__", "max_xeq_size": "__MAX_XEQ_SIZE__",
    }
    try:
        return {field: macros[macro] for field, macro in names.items()}
    except KeyError as error:
        raise phase_b.PhaseBError(f"frozen WDSat config lacks {error.args[0]}") from error


def host_identity(path: Path, label: str, *, executable: bool) -> dict[str, Any]:
    try:
        invocation = path.absolute()
        resolved = invocation.resolve(strict=True)
        metadata = resolved.stat()
    except OSError as error:
        raise phase_b.PhaseBError(f"cannot identify {label} {path}: {error}") from error
    if not stat.S_ISREG(metadata.st_mode) or (executable and not os.access(invocation, os.X_OK)):
        raise phase_b.PhaseBError(f"{label} is not an executable regular file: {path}")
    try:
        data = resolved.read_bytes()
    except OSError as error:
        raise phase_b.PhaseBError(f"cannot read {label} {resolved}: {error}") from error
    return {
        "path": str(invocation), "resolved_path": str(resolved),
        "bytes": len(data), "sha256": phase_b.sha256_bytes(data),
    }


def portable_binary_identity(path: Path) -> dict[str, Any]:
    identity = phase_b.executable_identity(path, "final WDSat binary")
    return {"path": "wdsat_solver", "bytes": identity["bytes"], "sha256": identity["sha256"]}


def same_content(left: Any, right: Any) -> bool:
    return isinstance(left, dict) and isinstance(right, dict) and all(
        left.get(field) == right.get(field) for field in ("bytes", "sha256")
    )


def validate_recorded_identity(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != {"path", "resolved_path", "bytes", "sha256"}:
        raise phase_b.PhaseBError(f"{context} uses an invalid identity schema")
    if not all(isinstance(value[field], str) and Path(value[field]).is_absolute() for field in ("path", "resolved_path")):
        raise phase_b.PhaseBError(f"{context} paths must be absolute")
    if isinstance(value.get("bytes"), bool) or not isinstance(value.get("bytes"), int) or value["bytes"] <= 0:
        raise phase_b.PhaseBError(f"{context} byte count is invalid")
    phase_b.require_hex64(value.get("sha256"), f"{context} hash")
    return value


def resource_summary(processes: list[dict[str, Any]]) -> dict[str, Any]:
    metrics = [process["metrics"] for process in processes]
    total_core_seconds = round(math.fsum(item["total_core_seconds"] for item in metrics), 12)
    return {
        "total_core_seconds": total_core_seconds,
        "single_core_seconds": total_core_seconds,
        "single_core_seconds_alias_of": "total_core_seconds",
        "single_core_elapsed_seconds": None,
        "summed_process_wall_seconds": round(math.fsum(item["wall_seconds"] for item in metrics), 12),
        "largest_child_peak_rss_bytes": max((item["peak_rss_bytes"] for item in metrics), default=0),
        "aggregate_process_tree_peak_rss_bytes": None,
    }


def successful_process(record: dict[str, Any], role: str) -> dict[str, Any]:
    if record["timed_out"] or record["returncode"] != 0 or record["orphan_group_terminated"]:
        raise phase_b.PhaseBError(
            f"WDSat {role} failed: returncode={record['returncode']} "
            f"timed_out={record['timed_out']} orphan_group_terminated={record['orphan_group_terminated']}"
        )
    return phase_b.compact_process(record)


def self_hashed(value: dict[str, Any], field: str) -> dict[str, Any]:
    result = deepcopy(value)
    if field in result:
        raise phase_b.PhaseBError(f"refusing to overwrite self-hash field {field}")
    result[field] = phase_b.canonical_sha256(result)
    return result


def validate_self_hash(value: dict[str, Any], field: str, context: str) -> None:
    payload = dict(value)
    actual = payload.pop(field, None)
    phase_b.require_hex64(actual, f"{context} self-hash")
    if phase_b.canonical_sha256(payload) != actual:
        raise phase_b.PhaseBError(f"{context} payload self-hash differs")


def validate_resources(metrics: Any, context: str) -> dict[str, Any]:
    if not isinstance(metrics, dict) or set(metrics) != RESOURCE_FIELDS:
        raise phase_b.PhaseBError(f"{context} lacks the exact resource schema")
    for field in RESOURCE_FIELDS - {"meter", "peak_rss_bytes"}:
        value = metrics[field]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise phase_b.PhaseBError(f"{context} has invalid {field}")
    if isinstance(metrics["peak_rss_bytes"], bool) or not isinstance(metrics["peak_rss_bytes"], int) or metrics["peak_rss_bytes"] < 0:
        raise phase_b.PhaseBError(f"{context} has invalid peak RSS")
    if (
        metrics["meter"] != "fresh-process getrusage(RUSAGE_CHILDREN)"
        or not math.isclose(metrics["total_core_seconds"], metrics["user_seconds"] + metrics["system_seconds"], rel_tol=0, abs_tol=1e-9)
        or metrics["single_core_seconds"] != metrics["total_core_seconds"]
    ):
        raise phase_b.PhaseBError(f"{context} has inconsistent resource accounting")
    return metrics


def local_artifact_name(value: Any, expected: str, context: str) -> str:
    if value != expected or Path(value).is_absolute() or len(Path(value).parts) != 1:
        raise phase_b.PhaseBError(f"{context} must be the local artifact {expected}")
    return value


def validate_source_inventory(path: Path, receipt: dict[str, Any]) -> dict[str, Any]:
    value, _ = phase_b.read_json(path, "configured WDSat source inventory")
    fields = {
        "schema", "source_commit", "config_sha256", "inventory",
        "inventory_sha256", "record_payload_sha256",
    }
    if not isinstance(value, dict) or set(value) != fields or value.get("schema") != SOURCE_INVENTORY_SCHEMA:
        raise phase_b.PhaseBError("configured WDSat source inventory uses an unexpected schema")
    validate_self_hash(value, "record_payload_sha256", "configured WDSat source inventory")
    inventory = value.get("inventory")
    if not isinstance(inventory, list) or phase_b.canonical_sha256(inventory) != value.get("inventory_sha256"):
        raise phase_b.PhaseBError("configured WDSat source inventory hash differs")
    if value.get("source_commit") != receipt["source_commit"] or value.get("config_sha256") != receipt["config_sha256"]:
        raise phase_b.PhaseBError("configured WDSat source inventory binding differs")
    source_root = path.parent / "build" / "src"
    for item in inventory:
        if not isinstance(item, dict) or set(item) != {"path", "bytes", "sha256"}:
            raise phase_b.PhaseBError("configured WDSat source inventory contains an invalid record")
        relative = Path(item["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise phase_b.PhaseBError("configured WDSat source inventory contains an unsafe path")
        data = phase_b.regular_file_bytes(source_root / relative, "configured WDSat source input")
        if len(data) != item["bytes"] or phase_b.sha256_bytes(data) != item["sha256"]:
            raise phase_b.PhaseBError("configured WDSat source input changed during the build")
    by_path = {item.get("path"): item for item in inventory if isinstance(item, dict)}
    config = by_path.get("config.h")
    if not isinstance(config, dict) or config.get("sha256") != receipt["config_sha256"]:
        raise phase_b.PhaseBError("retained WDSat config differs from the frozen config")
    if "makefile" not in by_path:
        raise phase_b.PhaseBError("retained WDSat source lacks src/makefile")
    return value


def expected_recipe_processes(
    receipt: dict[str, Any], source_inventory: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    context = receipt["build_context"]
    source = Path(context["source_checkout"])
    config = Path(context["frozen_config_path"])
    capsule = Path(context["original_capsule_root"])
    build_root = capsule / "build"
    cp = receipt["toolchain"]["cp"]["resolved_path"]
    make = receipt["toolchain"]["make"]["resolved_path"]
    configured_inputs = [
        {
            "path": str(build_root / "src" / item["path"]),
            "bytes": item["bytes"],
            "sha256": item["sha256"],
        }
        for item in source_inventory
    ]
    config_input = next(item for item in configured_inputs if Path(item["path"]).name == "config.h")
    return [
        {
            "role": "source-copy",
            "command": [cp, "-R", str(source / "src"), str(build_root / "src")],
            "cwd": str(capsule),
            "tool_sha256": receipt["toolchain"]["cp"]["sha256"],
            "inputs": [],
        },
        {
            "role": "config-install",
            "command": [cp, str(config), str(build_root / "src/config.h")],
            "cwd": str(capsule),
            "tool_sha256": receipt["toolchain"]["cp"]["sha256"],
            "inputs": [],
        },
        {
            "role": "clean",
            "command": [make, "-C", "src", "clean"],
            "cwd": str(build_root),
            "tool_sha256": receipt["toolchain"]["make"]["sha256"],
            "inputs": [dict(config_input)],
        },
        {
            "role": "build",
            "command": [make, "-C", "src"],
            "cwd": str(build_root),
            "tool_sha256": receipt["toolchain"]["make"]["sha256"],
            "inputs": deepcopy(configured_inputs),
        },
        {
            "role": "binary-copy",
            "command": [cp, str(build_root / "wdsat_solver"), str(capsule / "wdsat_solver")],
            "cwd": str(capsule),
            "tool_sha256": receipt["toolchain"]["cp"]["sha256"],
            "inputs": [{
                "path": str(build_root / "wdsat_solver"),
                "bytes": receipt["binary_identity"]["bytes"],
                "sha256": receipt["binary_sha256"],
            }],
        },
    ]


def validate_processes(
    capsule: Path,
    receipt: dict[str, Any],
    protocol: dict[str, Any],
    source_inventory: list[dict[str, Any]],
) -> None:
    processes = receipt.get("build_processes")
    if not isinstance(processes, list) or [item.get("role") if isinstance(item, dict) else None for item in processes] != EXPECTED_ROLES:
        raise phase_b.PhaseBError("WDSat receipt lacks the required separately metered stages")
    if receipt.get("frozen_recipe") != FROZEN_RECIPE:
        raise phase_b.PhaseBError("WDSat receipt differs from the frozen WDSat build recipe")
    expected_rows = expected_recipe_processes(receipt, source_inventory)
    environment = receipt["child_environment"]
    timeout = protocol["execution"]["per_process_watchdog_seconds"]
    for process, expected in zip(processes, expected_rows, strict=True):
        role = expected["role"]
        if not isinstance(process, dict) or set(process) != PROCESS_FIELDS:
            raise phase_b.PhaseBError(f"WDSat {role} process record uses an unexpected schema")
        expected_names = {
            "intent_path": f"{role}.intent.json", "stdout_path": f"{role}.stdout",
            "stderr_path": f"{role}.stderr", "metrics_path": f"{role}.metrics.json",
        }
        for field, name in expected_names.items():
            local_artifact_name(process.get(field), name, f"WDSat {role} {field}")
        if (
            process.get("command") != expected["command"]
            or process.get("inputs") != expected["inputs"]
            or process.get("watchdog_seconds") != timeout
            or process.get("returncode") != 0
            or process.get("timed_out") is not False
            or process.get("orphan_group_terminated") is not False
        ):
            raise phase_b.PhaseBError(f"WDSat {role} differs from the frozen WDSat build recipe")
        validate_resources(process.get("metrics"), f"WDSat {role}")
        intent, _ = phase_b.read_json(capsule / process["intent_path"], f"WDSat {role} intent")
        expected_intent = {
            "schema": "koblitz_pdp_phase_b_child_intent.v1",
            "role": role,
            "command": expected["command"],
            "cwd": expected["cwd"],
            "watchdog_seconds": timeout,
            "stdin": "devnull",
            "close_fds": True,
            "environment": environment,
            "inputs": expected["inputs"],
            "executable_sha256": expected["tool_sha256"],
        }
        if not isinstance(intent, dict) or set(intent) != INTENT_FIELDS or intent != expected_intent:
            raise phase_b.PhaseBError(f"WDSat {role} differs from its exact raw intent")
        metrics, _ = phase_b.read_json(capsule / process["metrics_path"], f"WDSat {role} metrics")
        raw_projection = {
            key: process[key]
            for key in ("command", "returncode", "watchdog_seconds", "timed_out", "orphan_group_terminated", "metrics")
        }
        if not isinstance(metrics, dict) or set(metrics) != METRICS_FIELDS or metrics != raw_projection:
            raise phase_b.PhaseBError(f"WDSat {role} process differs from its exact raw metrics")
        validate_resources(metrics["metrics"], f"WDSat {role} raw metrics")
        for stream in ("stdout", "stderr"):
            data = phase_b.regular_file_bytes(capsule / process[f"{stream}_path"], f"WDSat {role} {stream}")
            if len(data) != process[f"{stream}_bytes"] or phase_b.sha256_bytes(data) != process[f"{stream}_sha256"]:
                raise phase_b.PhaseBError(f"WDSat {role} {stream} differs from its process record")


def validate_capsule(
    seal_path: Path,
    protocol: dict[str, Any],
    expected_binary_identity: dict[str, Any],
    expected_implementation_state: dict[str, Any],
    require_clean_implementation: bool = True,
) -> dict[str, Any]:
    """Validate one WDSat capsule and return a relocation-safe normalized view."""
    phase_b.validate_protocol(protocol)
    if not isinstance(require_clean_implementation, bool):
        raise phase_b.PhaseBError("clean-implementation policy must be boolean")
    seal_path = seal_path.resolve(strict=True)
    if seal_path.name != "build-seal.json":
        raise phase_b.PhaseBError("WDSat capsule seal must be named build-seal.json")
    capsule = seal_path.parent
    seal, seal_bytes = phase_b.read_json(seal_path, "WDSat build seal")
    if not isinstance(seal, dict) or set(seal) != SEAL_FIELDS:
        raise phase_b.PhaseBError("WDSat build seal uses an unexpected field set")
    validate_self_hash(seal, "seal_payload_sha256", "WDSat build seal")
    if seal.get("schema") != BUILD_SEAL_SCHEMA or seal.get("status") != "build_frozen":
        raise phase_b.PhaseBError("WDSat build seal uses an unexpected schema or status")
    expected_state = validate_implementation_state(expected_implementation_state, "expected implementation state")
    sealed_state = validate_implementation_state(seal.get("implementation_state"), "sealed implementation state")
    if sealed_state != expected_state:
        raise phase_b.PhaseBError("WDSat build seal differs from the expected implementation state")
    if require_clean_implementation and (sealed_state["dirty"] or sealed_state["porcelain"]):
        raise phase_b.PhaseBError("production WDSat custody requires a clean implementation state")
    local_artifact_name(seal.get("binary_path"), "wdsat_solver", "WDSat binary path")
    local_artifact_name(seal.get("receipt_path"), "receipt.json", "WDSat receipt path")
    inventory = seal.get("inventory")
    if not isinstance(inventory, list) or phase_b.canonical_sha256(inventory) != seal.get("inventory_sha256"):
        raise phase_b.PhaseBError("WDSat build seal inventory hash differs")
    if phase_b.all_regular_inventory(capsule, {"build-seal.json"}) != inventory:
        raise phase_b.PhaseBError("WDSat capsule inventory differs from its build seal")

    receipt_path = capsule / seal["receipt_path"]
    receipt, receipt_bytes = phase_b.read_json(receipt_path, "WDSat build receipt")
    if phase_b.sha256_bytes(receipt_bytes) != seal.get("receipt_sha256"):
        raise phase_b.PhaseBError("WDSat build receipt hash differs from its seal")
    if not isinstance(receipt, dict) or set(receipt) != RECEIPT_FIELDS:
        raise phase_b.PhaseBError("WDSat build receipt uses an unexpected field set")
    validate_self_hash(receipt, "receipt_payload_sha256", "WDSat build receipt")
    if receipt.get("schema") != BUILD_RECEIPT_SCHEMA or receipt.get("status") != "completed":
        raise phase_b.PhaseBError("WDSat build receipt uses an unexpected schema or status")
    protocol_hash = phase_b.canonical_sha256(protocol)
    if receipt.get("protocol_payload_sha256") != protocol_hash or seal.get("protocol_payload_sha256") != protocol_hash:
        raise phase_b.PhaseBError("WDSat capsule differs from the requested Phase-B protocol")
    policy = protocol["wdsat_build"]
    if (
        receipt.get("source_commit") != policy["source_commit"]
        or receipt.get("source_dirty") is not False
        or seal.get("source_commit") != policy["source_commit"]
    ):
        raise phase_b.PhaseBError("WDSat capsule does not bind the clean frozen source commit")
    config = receipt.get("config")
    if not isinstance(config, str) or phase_b.sha256_bytes(config.encode()) != receipt.get("config_sha256"):
        raise phase_b.PhaseBError("WDSat receipt config text differs from its hash")
    if (
        receipt["config_sha256"] != policy["frozen_config_sha256"]
        or seal.get("config_sha256") != policy["frozen_config_sha256"]
        or parse_config(config) != policy["limits"]
        or receipt.get("limits") != policy["limits"]
    ):
        raise phase_b.PhaseBError("WDSat capsule differs from the frozen config or limits")
    binary = receipt.get("binary_identity")
    if not isinstance(binary, dict) or set(binary) != {"path", "bytes", "sha256"} or binary.get("path") != "wdsat_solver":
        raise phase_b.PhaseBError("WDSat receipt binary identity uses an unexpected schema")
    phase_b.require_hex64(binary.get("sha256"), "WDSat receipt binary hash")
    if (
        binary.get("sha256") != receipt.get("binary_sha256")
        or binary != seal.get("binary_identity")
        or binary.get("sha256") != seal.get("binary_sha256")
    ):
        raise phase_b.PhaseBError("WDSat binary identity differs between receipt and seal")
    actual_binary = portable_binary_identity(capsule / "wdsat_solver")
    if actual_binary != binary or not same_content(binary, expected_binary_identity):
        raise phase_b.PhaseBError("WDSat capsule does not bind the requested executable")

    context = receipt.get("build_context")
    if (
        not isinstance(context, dict)
        or set(context) != {"source_checkout", "frozen_config_path", "original_capsule_root"}
        or not all(isinstance(path, str) and Path(path).is_absolute() for path in context.values())
    ):
        raise phase_b.PhaseBError("WDSat receipt uses an invalid build context")
    if receipt.get("claim_boundary") != CLAIM_BOUNDARY:
        raise phase_b.PhaseBError("WDSat build claim boundary changed")
    driver = validate_recorded_identity(receipt.get("build_driver_identity"), "WDSat build driver")
    meter = validate_recorded_identity(receipt.get("process_meter_identity"), "WDSat process meter")
    if not same_content(driver, host_identity(Path(__file__), "trusted WDSat build driver", executable=False)):
        raise phase_b.PhaseBError("WDSat build driver differs from the trusted current implementation")
    if not same_content(meter, host_identity(phase_b.DEFAULT_METER, "trusted process meter", executable=False)):
        raise phase_b.PhaseBError("WDSat process meter differs from the trusted current implementation")
    toolchain = receipt.get("toolchain")
    if not isinstance(toolchain, dict) or set(toolchain) != {"cp", "make", "cc"}:
        raise phase_b.PhaseBError("WDSat build toolchain uses an unexpected schema")
    for name, identity in toolchain.items():
        validate_recorded_identity(identity, f"WDSat {name} identity")
        if host_identity(Path(identity["path"]), name, executable=True) != identity:
            raise phase_b.PhaseBError(f"WDSat {name} differs from the recorded build tool")
    expected_environment = phase_b.safe_child_environment()
    expected_environment["CC"] = toolchain["cc"]["resolved_path"]
    if receipt.get("child_environment") != expected_environment:
        raise phase_b.PhaseBError("WDSat build environment differs from the frozen environment")
    local_artifact_name(
        receipt.get("configured_source_inventory_path"),
        "configured-source-inventory.json",
        "configured source inventory path",
    )
    source_record = validate_source_inventory(
        capsule / receipt["configured_source_inventory_path"], receipt
    )
    validate_processes(capsule, receipt, protocol, source_record["inventory"])
    if receipt.get("resources") != resource_summary(receipt["build_processes"]):
        raise phase_b.PhaseBError("WDSat aggregate build resources differ from the process records")
    markers = phase_b.forbidden_material(protocol)
    phase_b.assert_no_forbidden_structured(receipt, markers, "WDSat build receipt")
    phase_b.assert_no_forbidden_text(json.dumps(receipt, sort_keys=True), markers, "WDSat build receipt")

    return {
        "schema": BUILD_RECEIPT_SCHEMA,
        "seal_sha256": phase_b.sha256_bytes(seal_bytes),
        "receipt_bytes": len(receipt_bytes),
        "receipt_sha256": phase_b.sha256_bytes(receipt_bytes),
        "protocol_payload_sha256": protocol_hash,
        "source_commit": receipt["source_commit"],
        "config_sha256": receipt["config_sha256"],
        "config": receipt["config"],
        "binary_path": "wdsat_solver",
        "binary_identity": dict(binary),
        "binary_sha256": binary["sha256"],
        "limits": deepcopy(receipt["limits"]),
        "build_processes": deepcopy(receipt["build_processes"]),
        "resources": deepcopy(receipt["resources"]),
        "claim_boundary": receipt["claim_boundary"],
        "implementation_state": deepcopy(sealed_state),
        "inventory_sha256": seal["inventory_sha256"],
    }


def copy_capsule(
    seal_path: Path,
    destination: Path,
    protocol: dict[str, Any],
    expected_binary_identity: dict[str, Any],
    expected_implementation_state: dict[str, Any],
    require_clean_implementation: bool = True,
) -> dict[str, Any]:
    """Validate before copying, then validate both source and relocated copy."""
    source_seal = seal_path.resolve(strict=True)
    source_root = source_seal.parent
    destination = destination.resolve()
    if destination.exists() or destination.is_symlink():
        raise phase_b.PhaseBError(f"WDSat capsule destination must be new: {destination}")
    try:
        destination.relative_to(source_root)
    except ValueError:
        pass
    else:
        raise phase_b.PhaseBError("WDSat capsule destination must not be inside its source")
    before = validate_capsule(
        source_seal, protocol, expected_binary_identity, expected_implementation_state,
        require_clean_implementation,
    )
    shutil.copytree(source_root, destination, copy_function=shutil.copy2, symlinks=False)
    copied = validate_capsule(
        destination / source_seal.name, protocol, expected_binary_identity,
        expected_implementation_state, require_clean_implementation,
    )
    after = validate_capsule(
        source_seal, protocol, expected_binary_identity, expected_implementation_state,
        require_clean_implementation,
    )
    if before != copied or before != after:
        raise phase_b.PhaseBError("WDSat capsule changed while making its relocated copy")
    return copied


def build(
    protocol_path: Path,
    source: Path,
    config_path: Path,
    output: Path,
    meter: Path,
) -> dict[str, Any]:
    protocol_path = protocol_path.resolve()
    source = source.resolve()
    config_path = config_path.resolve()
    output = output.resolve()
    meter = meter.resolve()
    protocol, _ = phase_b.read_json(protocol_path, "Phase-B protocol")
    phase_b.validate_protocol(protocol)
    policy = protocol["wdsat_build"]
    if output.exists() or output.is_symlink():
        raise phase_b.PhaseBError(f"WDSat build output must be new: {output}")
    source_metadata = source.lstat()
    if stat.S_ISLNK(source_metadata.st_mode) or not stat.S_ISDIR(source_metadata.st_mode):
        raise phase_b.PhaseBError("WDSat source must be a real directory")
    if not (source / "src/makefile").is_file():
        raise phase_b.PhaseBError("WDSat source lacks src/makefile")
    commit, status = git_state(source)
    if commit != policy["source_commit"] or status:
        raise phase_b.PhaseBError(
            f"WDSat source must be clean commit {policy['source_commit']}; got {commit} with {status}"
        )
    config_bytes = phase_b.regular_file_bytes(config_path, "frozen WDSat config")
    if phase_b.sha256_bytes(config_bytes) != policy["frozen_config_sha256"]:
        raise phase_b.PhaseBError("WDSat config bytes differ from the frozen protocol")
    try:
        config_text = config_bytes.decode()
    except UnicodeDecodeError as error:
        raise phase_b.PhaseBError("WDSat config is not UTF-8 text") from error
    if parse_config(config_text) != policy["limits"]:
        raise phase_b.PhaseBError("WDSat config macros differ from the frozen protocol limits")
    implementation_before = validate_implementation_state(
        current_implementation_state(), "current implementation state"
    )
    cp_path = Path(shutil.which("cp") or "/bin/cp").resolve()
    make_path = Path(shutil.which("make") or "/usr/bin/make").resolve()
    cc_path = Path(shutil.which("cc") or "/usr/bin/cc").resolve()
    toolchain = {
        "cp": host_identity(cp_path, "cp", executable=True),
        "make": host_identity(make_path, "make", executable=True),
        "cc": host_identity(cc_path, "cc", executable=True),
    }
    driver_identity = host_identity(Path(__file__), "WDSat build driver", executable=False)
    meter_identity = host_identity(meter, "WDSat process meter", executable=False)
    trusted_meter = host_identity(phase_b.DEFAULT_METER, "trusted process meter", executable=False)
    if not same_content(meter_identity, trusted_meter):
        raise phase_b.PhaseBError("WDSat build meter differs from the trusted Phase-B process meter")

    output.mkdir(parents=False)
    build_root = output / "build"
    build_root.mkdir()
    environment = phase_b.safe_child_environment()
    environment["CC"] = str(cc_path)
    markers = phase_b.forbidden_material(protocol)
    timeout = protocol["execution"]["per_process_watchdog_seconds"]
    processes: list[dict[str, Any]] = []

    def measured(role: str, command: list[str], cwd: Path, inputs: list[Path], tool: str) -> dict[str, Any]:
        record = phase_b.run_metered(
            role=role,
            command=command,
            cwd=cwd,
            task_root=output,
            input_paths=inputs,
            timeout=timeout,
            meter=meter,
            environment=environment,
            markers=markers,
            # System build tools may intentionally be hard-linked. Their exact
            # content identities are checked before and after the build rather
            # than passed through the single-link artifact guard.
            expected_executable_sha256=toolchain[tool]["sha256"],
            allow_hardlinked_executable=True,
        )
        processes.append(successful_process(record, role))
        return record

    measured(
        "source-copy",
        [str(cp_path), "-R", str(source / "src"), str(build_root / "src")],
        output, [], "cp",
    )
    measured(
        "config-install",
        [str(cp_path), str(config_path), str(build_root / "src/config.h")],
        output, [], "cp",
    )
    if phase_b.regular_file_bytes(build_root / "src/config.h", "installed WDSat config") != config_bytes:
        raise phase_b.PhaseBError("installed WDSat config differs from its frozen source")
    measured(
        "clean", [str(make_path), "-C", "src", "clean"], build_root,
        [build_root / "src/config.h"], "make",
    )
    source_inventory = phase_b.all_regular_inventory(build_root / "src")
    expected_inventory = expected_configured_source_inventory(source, config_bytes)
    if source_inventory != expected_inventory:
        raise phase_b.PhaseBError(
            "cleaned WDSat build source differs from the tracked commit plus frozen config"
        )
    source_record = self_hashed(
        {
            "schema": SOURCE_INVENTORY_SCHEMA,
            "source_commit": commit,
            "config_sha256": policy["frozen_config_sha256"],
            "inventory": source_inventory,
            "inventory_sha256": phase_b.canonical_sha256(source_inventory),
        },
        "record_payload_sha256",
    )
    phase_b.write_json_new(output / "configured-source-inventory.json", source_record)
    source_inputs = [build_root / "src" / item["path"] for item in source_inventory]
    measured(
        "build", [str(make_path), "-C", "src"], build_root,
        source_inputs, "make",
    )
    built_binary = build_root / "wdsat_solver"
    phase_b.regular_file_bytes(built_binary, "built WDSat binary")
    final_binary = output / "wdsat_solver"
    measured(
        "binary-copy", [str(cp_path), str(built_binary), str(final_binary)], output,
        [built_binary], "cp",
    )
    os.chmod(final_binary, 0o755)
    binary_identity = portable_binary_identity(final_binary)
    if phase_b.sha256_file(built_binary, "built WDSat binary") != binary_identity["sha256"]:
        raise phase_b.PhaseBError("final WDSat binary differs from the compiled binary")
    if phase_b.regular_file_bytes(build_root / "src/config.h", "retained WDSat config") != config_bytes:
        raise phase_b.PhaseBError("WDSat build changed the frozen config")
    validate_source_inventory(output / "configured-source-inventory.json", {
        "source_commit": commit,
        "config_sha256": policy["frozen_config_sha256"],
    })
    if git_state(source) != (commit, status):
        raise phase_b.PhaseBError("WDSat source changed during the measured build")
    if current_implementation_state() != implementation_before:
        raise phase_b.PhaseBError("Phase-B implementation state changed during the WDSat build")
    for name, identity in toolchain.items():
        if host_identity(Path(identity["path"]), name, executable=True) != identity:
            raise phase_b.PhaseBError(f"WDSat build tool changed during execution: {name}")

    receipt = self_hashed(
        {
            "schema": BUILD_RECEIPT_SCHEMA,
            "status": "completed",
            "protocol_payload_sha256": phase_b.canonical_sha256(protocol),
            "source_commit": commit,
            "source_dirty": False,
            "config": config_text,
            "config_sha256": policy["frozen_config_sha256"],
            "binary_identity": binary_identity,
            "binary_sha256": binary_identity["sha256"],
            "limits": policy["limits"],
            "build_context": {
                "source_checkout": str(source),
                "frozen_config_path": str(config_path),
                "original_capsule_root": str(output),
            },
            "child_environment": environment,
            "frozen_recipe": FROZEN_RECIPE,
            "configured_source_inventory_path": "configured-source-inventory.json",
            "build_driver_identity": driver_identity,
            "process_meter_identity": meter_identity,
            "toolchain": toolchain,
            "build_processes": processes,
            "resources": resource_summary(processes),
            "claim_boundary": CLAIM_BOUNDARY,
        },
        "receipt_payload_sha256",
    )
    phase_b.write_json_new(output / "receipt.json", receipt)
    inventory = phase_b.all_regular_inventory(output, {"build-seal.json"})
    build_seal = self_hashed(
        {
            "schema": BUILD_SEAL_SCHEMA,
            "status": "build_frozen",
            "protocol_payload_sha256": phase_b.canonical_sha256(protocol),
            "source_commit": commit,
            "config_sha256": policy["frozen_config_sha256"],
            "binary_path": "wdsat_solver",
            "binary_identity": binary_identity,
            "binary_sha256": binary_identity["sha256"],
            "receipt_path": "receipt.json",
            "receipt_sha256": phase_b.sha256_file(output / "receipt.json", "WDSat build receipt"),
            "implementation_state": implementation_before,
            "inventory": inventory,
            "inventory_sha256": phase_b.canonical_sha256(inventory),
        },
        "seal_payload_sha256",
    )
    phase_b.write_json_new(output / "build-seal.json", build_seal)
    validate_capsule(
        output / "build-seal.json", protocol,
        phase_b.executable_identity(final_binary, "final WDSat binary"),
        implementation_before, require_clean_implementation=False,
    )
    return build_seal


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=phase_b.DEFAULT_PROTOCOL)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--meter", type=Path, default=phase_b.DEFAULT_METER)
    args = parser.parse_args()
    try:
        protocol, _ = phase_b.read_json(args.protocol.resolve(), "Phase-B protocol")
        phase_b.validate_protocol(protocol)
        config = args.config.resolve() if args.config is not None else (
            phase_b.REPO / protocol["wdsat_build"]["frozen_config_path"]
        ).resolve()
        result = build(
            args.protocol.resolve(), args.source.resolve(), config,
            args.output.resolve(), args.meter.resolve(),
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, phase_b.PhaseBError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
