#!/usr/bin/env python3
"""Stage and run the blinded Koblitz PDP Phase-B SAT panel.

The solver-facing command deliberately has no ground-truth argument.  It reads
an authenticated, sanitized root containing only ``blind/bundle.json`` and a
class-blind binding.  A separate program opens the Phase-A ground truth only
after this program has sealed every solver output.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import statistics
from typing import Any, Iterable

import run_koblitz_pdp_matrix as matrix


REPO = Path(__file__).resolve().parents[1]
DEFAULT_PROTOCOL = (
    REPO
    / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
    / "stage-20-balanced-pdp-phase-b-protocol.json"
)
DEFAULT_METER = Path(__file__).with_name("process_meter.py")

PROTOCOL_SCHEMA = "koblitz_balanced_pdp_phase_b_protocol.v1"
BINDING_SCHEMA = "koblitz_pdp_phase_b_solver_input_binding.v1"
BLIND_BUNDLE_SCHEMA = "koblitz_pdp_blind_bundle.v1"
BLIND_INSTANCE_SCHEMA = "koblitz_pdp_blind_instance.v1"
RUN_PLAN_SCHEMA = "koblitz_pdp_phase_b_execution_plan.v1"
TASK_SCHEMA = "koblitz_pdp_phase_b_task_result.v1"
RUN_SUMMARY_SCHEMA = "koblitz_pdp_phase_b_run_summary.v1"
RUN_SEAL_SCHEMA = "koblitz_pdp_phase_b_run_seal.v1"
PRODUCTION_EVIDENCE_CLASS = "internal_measurement_with_bound_tool_builds_pending_outer_receipt"

BLIND_TOP_LEVEL_FIELDS = {"schema", "scope", "instance_count", "instances"}
BLIND_INSTANCE_FIELDS = {
    "schema",
    "blind_instance_id",
    "source_nonce",
    "source_system_id",
    "cell_id",
    "n",
    "ell",
    "m",
    "basis",
    "curve_a",
    "factor_index",
    "target",
}
TARGET_FIELDS = {"x", "y"}
BINDING_FIELDS = {
    "schema",
    "scope",
    "source_phase_a_seal_sha256",
    "source_phase_a_protocol_sha256",
    "phase_a_source_revision",
    "blind_bundle_path",
    "blind_bundle_bytes",
    "blind_bundle_sha256",
    "blind_bundle_schema",
    "blind_instance_count",
    "selected_cell_ids",
    "binding_payload_sha256",
}
DEFAULT_FORBIDDEN_MATERIAL = {
    "target_class",
    "witness_indices",
    "unconditional_draw_index",
    "selection_priority_hex",
    "planted_points",
    "sealed-oracle",
    "oracle-ledger",
    "oracle",
    "decomposable",
    "nondecomposable",
}
BACKENDS = ("native-xor", "wdsat", "cryptominisat")
HEX64 = re.compile(r"^[0-9a-f]{64}$")
HEX40 = re.compile(r"^[0-9a-f]{40}$")
BLIND_ID = re.compile(r"^b-[0-9a-f]{64}$")
SOURCE_SYSTEM_ID = re.compile(r"^x-[0-9a-f]{64}$")


class PhaseBError(RuntimeError):
    """A fail-closed Phase-B contract error."""


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise PhaseBError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def parse_json_bytes(data: bytes, context: str) -> Any:
    try:
        return json.loads(data, object_pairs_hook=reject_duplicate_pairs)
    except PhaseBError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PhaseBError(f"invalid JSON in {context}: {error}") from error


def regular_file_bytes(path: Path, context: str) -> bytes:
    try:
        metadata = path.lstat()
    except OSError as error:
        raise PhaseBError(f"cannot stat {context} {path}: {error}") from error
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        raise PhaseBError(f"{context} must be a regular non-symlink file: {path}")
    if metadata.st_nlink != 1:
        raise PhaseBError(f"{context} must not be hard-linked: {path}")
    try:
        return path.read_bytes()
    except OSError as error:
        raise PhaseBError(f"cannot read {context} {path}: {error}") from error


def read_json(path: Path, context: str) -> tuple[Any, bytes]:
    data = regular_file_bytes(path, context)
    return parse_json_bytes(data, context), data


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path, context: str = "file") -> str:
    return sha256_bytes(regular_file_bytes(path, context))


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def canonical_sha256(value: Any) -> str:
    return sha256_bytes(canonical_bytes(value))


def pretty_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def write_new(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags, 0o644)
    except OSError as error:
        raise PhaseBError(f"refusing to overwrite {path}: {error}") from error
    try:
        with os.fdopen(descriptor, "wb", closefd=True) as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
    except Exception:
        try:
            path.unlink()
        except OSError:
            pass
        raise


def write_json_new(path: Path, value: Any) -> None:
    write_new(path, pretty_bytes(value))


def require_hex64(value: Any, context: str) -> str:
    if not isinstance(value, str) or HEX64.fullmatch(value) is None:
        raise PhaseBError(f"{context} must be 64 lowercase hexadecimal characters")
    return value


def require_hex40(value: Any, context: str) -> str:
    if not isinstance(value, str) or HEX40.fullmatch(value) is None:
        raise PhaseBError(f"{context} must be 40 lowercase hexadecimal characters")
    return value


def forbidden_material(protocol: dict[str, Any]) -> set[str]:
    configured = protocol.get("source_contract", {}).get("forbidden_solver_material", [])
    if not isinstance(configured, list) or not all(isinstance(item, str) for item in configured):
        raise PhaseBError("source_contract.forbidden_solver_material must be a string list")
    result = {item.lower() for item in configured}
    if not DEFAULT_FORBIDDEN_MATERIAL.issubset(result):
        raise PhaseBError("protocol removed a mandatory forbidden-material marker")
    return result


def walk_keys(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for key, child in value.items():
            yield key
            yield from walk_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk_keys(child)


def assert_no_forbidden_structured(value: Any, markers: set[str], context: str) -> None:
    for key in walk_keys(value):
        if key.lower() in markers:
            raise PhaseBError(f"{context} contains forbidden field {key!r}")


def assert_no_forbidden_text(text: str, markers: set[str], context: str) -> None:
    lowered = text.lower()
    found = sorted(marker for marker in markers if marker in lowered)
    if found:
        raise PhaseBError(f"{context} contains forbidden material: {', '.join(found)}")


def validate_protocol(protocol: Any) -> dict[str, Any]:
    if not isinstance(protocol, dict) or protocol.get("schema") != PROTOCOL_SCHEMA:
        raise PhaseBError("unsupported Phase-B protocol schema")
    if protocol.get("status") != "frozen_before_solver_execution":
        raise PhaseBError("Phase-B protocol is not frozen before execution")
    phase_a = protocol.get("phase_a_binding")
    if not isinstance(phase_a, dict):
        raise PhaseBError("phase_a_binding must be an object")
    require_hex40(phase_a.get("source_revision"), "phase_a_binding.source_revision")
    for field in ("protocol_sha256", "seal_sha256", "blind_bundle_sha256"):
        require_hex64(phase_a.get(field), f"phase_a_binding.{field}")
    if phase_a.get("blind_bundle_schema") != BLIND_BUNDLE_SCHEMA:
        raise PhaseBError("unexpected bound blind bundle schema")
    if phase_a.get("blind_instance_schema") != BLIND_INSTANCE_SCHEMA:
        raise PhaseBError("unexpected bound blind instance schema")
    if not isinstance(phase_a.get("instance_count"), int) or phase_a["instance_count"] <= 0:
        raise PhaseBError("bound instance count must be positive")
    root_policy = protocol.get("solver_root")
    if not isinstance(root_policy, dict):
        raise PhaseBError("solver_root policy must be an object")
    if sorted(root_policy.get("exact_files", [])) != ["binding.json", "blind/bundle.json"]:
        raise PhaseBError("solver root exact file inventory changed")
    if root_policy.get("symlinks_allowed") is not False or root_policy.get(
        "unexpected_files_allowed"
    ) is not False:
        raise PhaseBError("solver root must reject symlinks and unexpected files")
    if root_policy.get("binding_schema") != BINDING_SCHEMA:
        raise PhaseBError("solver root binding schema changed")
    cells = protocol.get("cells")
    if not isinstance(cells, list) or not cells:
        raise PhaseBError("protocol must contain at least one cell")
    ids: set[str] = set()
    expected_instances = 0
    required_cell_fields = {
        "id",
        "n",
        "ell",
        "m",
        "basis",
        "curve_a",
        "factor_index",
        "expected_instances",
    }
    for cell in cells:
        if not isinstance(cell, dict) or set(cell) != required_cell_fields:
            raise PhaseBError("every protocol cell must use the exact Phase-B cell schema")
        if not isinstance(cell["id"], str) or cell["id"] in ids:
            raise PhaseBError("protocol cell ids must be unique strings")
        ids.add(cell["id"])
        if (
            not isinstance(cell["n"], int)
            or not 3 <= cell["n"] <= 63
            or not isinstance(cell["ell"], int)
            or not 1 <= cell["ell"] <= cell["n"]
            or cell["m"] != 3
            or cell["basis"] not in {"standard", "ggmp"}
            or cell["curve_a"] not in {0, 1}
            or not isinstance(cell["factor_index"], int)
            or cell["factor_index"] < 0
            or not isinstance(cell["expected_instances"], int)
            or cell["expected_instances"] <= 0
        ):
            raise PhaseBError(f"invalid parameters for protocol cell {cell['id']}")
        expected_instances += cell["expected_instances"]
    if expected_instances != phase_a["instance_count"]:
        raise PhaseBError("cell quotas do not equal the bound instance count")
    execution = protocol.get("execution")
    if not isinstance(execution, dict):
        raise PhaseBError("execution policy must be an object")
    if execution.get("backend_order") != list(BACKENDS):
        raise PhaseBError("backend order must be native-xor, WDSat, CryptoMiniSat")
    if execution.get("expected_backend_runs") != expected_instances * len(BACKENDS):
        raise PhaseBError("expected backend run count is inconsistent")
    if execution.get("one_process_at_a_time") is not True:
        raise PhaseBError("Phase B must run one process at a time")
    if execution.get("per_process_watchdog_seconds") != 120:
        raise PhaseBError("Phase B requires the frozen 120-second watchdog")
    if execution.get("single_thread_requested") is not True:
        raise PhaseBError("Phase B must request one solver thread")
    if execution.get("unknown_or_timeout_is_unsat") is not False:
        raise PhaseBError("Unknown and timeout must remain inconclusive")
    if not isinstance(execution.get("native_conflict_budget"), int) or execution[
        "native_conflict_budget"
    ] <= 0:
        raise PhaseBError("native conflict budget must be positive")
    build = protocol.get("wdsat_build")
    build_limit_names = {
        "max_anf_id",
        "max_degree",
        "max_id",
        "max_buffer_size",
        "max_eq",
        "max_eq_size",
        "max_xeq",
        "max_xeq_size",
    }
    if not isinstance(build, dict):
        raise PhaseBError("protocol lacks its frozen WDSat build policy")
    require_hex40(build.get("source_commit"), "wdsat_build.source_commit")
    require_hex64(build.get("frozen_config_sha256"), "wdsat_build.frozen_config_sha256")
    if not isinstance(build.get("frozen_config_path"), str):
        raise PhaseBError("wdsat_build.frozen_config_path must be a string")
    if (
        not isinstance(build.get("limits"), dict)
        or set(build["limits"]) != build_limit_names
        or not all(isinstance(value, int) and value > 0 for value in build["limits"].values())
    ):
        raise PhaseBError("protocol has invalid frozen WDSat limits")
    tool_builds = protocol.get("tool_build_accounting")
    if (
        not isinstance(tool_builds, dict)
        or tool_builds.get("full_cost_gate_passed") is not False
        or tool_builds.get("rust_release_build", {}).get("status") != "not_run"
        or tool_builds.get("rust_release_build", {}).get("required_for_full_cost_gate")
        is not True
        or tool_builds.get("cryptominisat_build", {}).get("status") != "not_run"
        or tool_builds.get("cryptominisat_build", {}).get("required_for_full_cost_gate")
        is not True
    ):
        raise PhaseBError(
            "protocol must keep uncharged Rust and CryptoMiniSat builds outside the full-cost gate"
        )
    require_hex40(
        tool_builds["cryptominisat_build"].get("source_commit"),
        "tool_build_accounting.cryptominisat_build.source_commit",
    )
    forbidden_material(protocol)
    return protocol


def cell_map(protocol: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {cell["id"]: cell for cell in protocol["cells"]}


def validate_blind_bundle(bundle: Any, protocol: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(bundle, dict) or set(bundle) != BLIND_TOP_LEVEL_FIELDS:
        raise PhaseBError("blind bundle uses an unexpected top-level schema")
    if bundle.get("schema") != BLIND_BUNDLE_SCHEMA:
        raise PhaseBError("blind bundle schema changed")
    instances = bundle.get("instances")
    if not isinstance(instances, list) or bundle.get("instance_count") != len(instances):
        raise PhaseBError("blind bundle instance count is inconsistent")
    if len(instances) != protocol["phase_a_binding"]["instance_count"]:
        raise PhaseBError("blind bundle count differs from its frozen binding")
    markers = forbidden_material(protocol)
    assert_no_forbidden_structured(bundle, markers, "blind bundle")
    assert_no_forbidden_text(json.dumps(bundle, sort_keys=True), markers, "blind bundle")
    cells = cell_map(protocol)
    counts: Counter[str] = Counter()
    ids: set[str] = set()
    targets: set[tuple[str, int, int]] = set()
    for index, instance in enumerate(instances):
        context = f"blind instance {index}"
        if not isinstance(instance, dict) or set(instance) != BLIND_INSTANCE_FIELDS:
            raise PhaseBError(f"{context} uses an unexpected field set")
        if instance.get("schema") != BLIND_INSTANCE_SCHEMA:
            raise PhaseBError(f"{context} has an unexpected schema")
        blind_id = instance.get("blind_instance_id")
        source_system = instance.get("source_system_id")
        if not isinstance(blind_id, str) or BLIND_ID.fullmatch(blind_id) is None:
            raise PhaseBError(f"{context} has an invalid blind id")
        if not isinstance(source_system, str) or SOURCE_SYSTEM_ID.fullmatch(source_system) is None:
            raise PhaseBError(f"{context} has an invalid source-system id")
        if blind_id in ids:
            raise PhaseBError(f"duplicate blind id {blind_id}")
        ids.add(blind_id)
        nonce = instance.get("source_nonce")
        if not isinstance(nonce, int) or not 0 <= nonce < 1 << 64:
            raise PhaseBError(f"{context} has an invalid source nonce")
        cell_id = instance.get("cell_id")
        if cell_id not in cells:
            raise PhaseBError(f"{context} refers to an unknown cell")
        cell = cells[cell_id]
        for field in ("n", "ell", "m", "basis", "curve_a", "factor_index"):
            if instance.get(field) != cell[field]:
                raise PhaseBError(f"{context} differs from frozen cell field {field}")
        target = instance.get("target")
        if not isinstance(target, dict) or set(target) != TARGET_FIELDS:
            raise PhaseBError(f"{context} target uses an unexpected schema")
        try:
            x = int(target["x"])
            y = int(target["y"])
        except (TypeError, ValueError) as error:
            raise PhaseBError(f"{context} target coordinates are not base-10 integers") from error
        if str(x) != target["x"] or str(y) != target["y"] or x < 0 or y < 0:
            raise PhaseBError(f"{context} target coordinates are not canonical unsigned decimals")
        if x >= 1 << cell["n"] or y >= 1 << cell["n"]:
            raise PhaseBError(f"{context} target coordinate exceeds the field width")
        target_key = (cell_id, x, y)
        if target_key in targets:
            raise PhaseBError(f"duplicate target in cell {cell_id}")
        targets.add(target_key)
        counts[cell_id] += 1
    for cell_id, cell in cells.items():
        if counts[cell_id] != cell["expected_instances"]:
            raise PhaseBError(f"{cell_id} has {counts[cell_id]} instances, expected {cell['expected_instances']}")
    return bundle


def binding_payload(binding: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in binding.items() if key != "binding_payload_sha256"}


def validate_binding(binding: Any, protocol: dict[str, Any], bundle_bytes: bytes) -> dict[str, Any]:
    if not isinstance(binding, dict) or set(binding) != BINDING_FIELDS:
        raise PhaseBError("solver input binding uses an unexpected field set")
    if binding.get("schema") != BINDING_SCHEMA:
        raise PhaseBError("solver input binding schema changed")
    markers = forbidden_material(protocol)
    assert_no_forbidden_structured(binding, markers, "solver input binding")
    assert_no_forbidden_text(json.dumps(binding, sort_keys=True), markers, "solver input binding")
    phase_a = protocol["phase_a_binding"]
    expected = {
        "source_phase_a_seal_sha256": phase_a["seal_sha256"],
        "source_phase_a_protocol_sha256": phase_a["protocol_sha256"],
        "phase_a_source_revision": phase_a["source_revision"],
        "blind_bundle_path": "blind/bundle.json",
        "blind_bundle_bytes": len(bundle_bytes),
        "blind_bundle_sha256": sha256_bytes(bundle_bytes),
        "blind_bundle_schema": phase_a["blind_bundle_schema"],
        "blind_instance_count": phase_a["instance_count"],
        "selected_cell_ids": sorted(cell_map(protocol)),
    }
    for field, value in expected.items():
        if binding.get(field) != value:
            raise PhaseBError(f"solver input binding differs at {field}")
    if binding["blind_bundle_sha256"] != phase_a["blind_bundle_sha256"]:
        raise PhaseBError("solver input binding does not match the frozen bundle hash")
    require_hex64(binding.get("binding_payload_sha256"), "binding_payload_sha256")
    if canonical_sha256(binding_payload(binding)) != binding["binding_payload_sha256"]:
        raise PhaseBError("solver input binding self-hash is invalid")
    return binding


def solver_root_inventory(root: Path) -> set[str]:
    try:
        root_metadata = root.lstat()
    except OSError as error:
        raise PhaseBError(f"cannot stat solver root {root}: {error}") from error
    if stat.S_ISLNK(root_metadata.st_mode) or not stat.S_ISDIR(root_metadata.st_mode):
        raise PhaseBError("solver root must be a real directory")
    files: set[str] = set()
    directories: set[str] = set()
    for current, names, filenames in os.walk(root, followlinks=False):
        current_path = Path(current)
        for name in names:
            path = current_path / name
            metadata = path.lstat()
            if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
                raise PhaseBError(f"solver root contains a non-directory or symlink: {path}")
            directories.add(path.relative_to(root).as_posix())
        for name in filenames:
            path = current_path / name
            regular_file_bytes(path, "solver-root file")
            files.add(path.relative_to(root).as_posix())
    if directories != {"blind"}:
        raise PhaseBError(f"solver root directory inventory changed: {sorted(directories)}")
    if files != {"binding.json", "blind/bundle.json"}:
        raise PhaseBError(f"solver root file inventory changed: {sorted(files)}")
    return files


def load_solver_root(root: Path, protocol: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    solver_root_inventory(root)
    bundle, bundle_bytes = read_json(root / "blind/bundle.json", "blind bundle")
    validate_blind_bundle(bundle, protocol)
    binding, _ = read_json(root / "binding.json", "solver input binding")
    validate_binding(binding, protocol, bundle_bytes)
    return binding, bundle


def validate_phase_a_seal(seal: Any, seal_bytes: bytes, protocol: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(seal, dict) or seal.get("schema") != "koblitz_pdp_phase_a_seal.v1":
        raise PhaseBError("source Phase-A seal schema changed")
    phase_a = protocol["phase_a_binding"]
    if sha256_bytes(seal_bytes) != phase_a["seal_sha256"]:
        raise PhaseBError("source Phase-A seal bytes do not match the frozen commitment")
    if seal.get("status") != "sealed_before_solver_execution" or seal.get(
        "execution_status"
    ) != "not_run":
        raise PhaseBError("source Phase-A seal is not in its pre-solver state")
    if seal.get("protocol_sha256") != phase_a["protocol_sha256"]:
        raise PhaseBError("source Phase-A protocol commitment changed")
    if seal.get("blind_bundle_sha256") != phase_a["blind_bundle_sha256"]:
        raise PhaseBError("source Phase-A blind commitment changed")
    implementation = seal.get("implementation", {})
    revision = implementation.get("source_revision", {}) if isinstance(implementation, dict) else {}
    if revision.get("commit") != phase_a["source_revision"] or revision.get("dirty") is not False:
        raise PhaseBError("source Phase-A implementation provenance is not the frozen clean revision")
    if sorted(seal.get("selected_cell_ids", [])) != sorted(cell_map(protocol)):
        raise PhaseBError("source Phase-A seal cell selection changed")
    return seal


def stage_solver_root(protocol_path: Path, seal_path: Path, bundle_path: Path, output: Path) -> dict[str, Any]:
    protocol, _ = read_json(protocol_path, "Phase-B protocol")
    validate_protocol(protocol)
    if output.exists() or output.is_symlink():
        raise PhaseBError(f"solver-root output must be new: {output}")
    seal, seal_bytes = read_json(seal_path, "source Phase-A seal")
    validate_phase_a_seal(seal, seal_bytes, protocol)
    bundle, bundle_bytes = read_json(bundle_path, "source blind bundle")
    validate_blind_bundle(bundle, protocol)
    phase_a = protocol["phase_a_binding"]
    if sha256_bytes(bundle_bytes) != phase_a["blind_bundle_sha256"]:
        raise PhaseBError("source blind bundle bytes do not match the frozen commitment")
    if seal.get("blind_bundle_sha256") != sha256_bytes(bundle_bytes):
        raise PhaseBError("source Phase-A seal does not authenticate the supplied blind bundle")
    payload = {
        "schema": BINDING_SCHEMA,
        "scope": "Sanitized class-blind input binding for the separately metered Phase-B solver run",
        "source_phase_a_seal_sha256": phase_a["seal_sha256"],
        "source_phase_a_protocol_sha256": phase_a["protocol_sha256"],
        "phase_a_source_revision": phase_a["source_revision"],
        "blind_bundle_path": "blind/bundle.json",
        "blind_bundle_bytes": len(bundle_bytes),
        "blind_bundle_sha256": sha256_bytes(bundle_bytes),
        "blind_bundle_schema": phase_a["blind_bundle_schema"],
        "blind_instance_count": len(bundle["instances"]),
        "selected_cell_ids": sorted(cell_map(protocol)),
    }
    binding = dict(payload)
    binding["binding_payload_sha256"] = canonical_sha256(payload)
    assert_no_forbidden_structured(binding, forbidden_material(protocol), "staged binding")
    assert_no_forbidden_text(
        json.dumps(binding, sort_keys=True), forbidden_material(protocol), "staged binding"
    )
    output.mkdir(parents=False)
    (output / "blind").mkdir()
    write_new(output / "blind/bundle.json", bundle_bytes)
    write_json_new(output / "binding.json", binding)
    load_solver_root(output, protocol)
    return binding


def git_state() -> dict[str, Any]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPO, text=True, capture_output=True, check=True
    ).stdout.strip()
    status_lines = subprocess.run(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.splitlines()
    return {"commit": commit, "dirty": bool(status_lines), "porcelain": status_lines}


def require_phase_a_ancestor(commit: str) -> None:
    completed = subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPO,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise PhaseBError("the solver implementation does not descend from the bound Phase-A revision")


def require_git_state(expected: dict[str, Any], context: str) -> None:
    current = git_state()
    if current != expected:
        raise PhaseBError(f"repository state changed during {context}")


def executable_identity(path: Path, label: str, *, executable: bool = True) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    data = regular_file_bytes(resolved, label)
    if executable and not os.access(resolved, os.X_OK):
        raise PhaseBError(f"{label} is not executable: {resolved}")
    return {"path": str(resolved), "bytes": len(data), "sha256": sha256_bytes(data)}


def validate_wdsat_build_receipt(
    value: Any,
    receipt_bytes: bytes,
    wdsat_identity: dict[str, Any],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    required_fields = {
        "schema",
        "status",
        "source_commit",
        "source_dirty",
        "config",
        "config_sha256",
        "binary_sha256",
        "limits",
        "build_processes",
        "claim_boundary",
    }
    if not isinstance(value, dict) or set(value) != required_fields:
        raise PhaseBError("WDSat build receipt uses an unexpected field set")
    if value.get("schema") != "koblitz_pdp_phase_b_wdsat_build_receipt.v1":
        raise PhaseBError("WDSat build receipt schema changed")
    if value.get("status") != "completed" or value.get("source_dirty") is not False:
        raise PhaseBError("WDSat build receipt is not a clean completed build")
    require_hex40(value.get("source_commit"), "WDSat source commit")
    frozen_build = protocol["wdsat_build"]
    if value["source_commit"] != frozen_build["source_commit"]:
        raise PhaseBError("WDSat build receipt uses the wrong source commit")
    config = value.get("config")
    if not isinstance(config, str) or sha256_bytes(config.encode()) != require_hex64(
        value.get("config_sha256"), "WDSat config hash"
    ):
        raise PhaseBError("WDSat build receipt config text does not match its hash")
    if value["config_sha256"] != frozen_build["frozen_config_sha256"]:
        raise PhaseBError("WDSat build receipt uses the wrong frozen config")
    if require_hex64(value.get("binary_sha256"), "WDSat binary hash") != wdsat_identity[
        "sha256"
    ]:
        raise PhaseBError("WDSat build receipt does not bind the requested executable")
    limit_names = {
        "max_anf_id",
        "max_degree",
        "max_id",
        "max_buffer_size",
        "max_eq",
        "max_eq_size",
        "max_xeq",
        "max_xeq_size",
    }
    limits = value.get("limits")
    if (
        not isinstance(limits, dict)
        or set(limits) != limit_names
        or not all(isinstance(item, int) and item > 0 for item in limits.values())
    ):
        raise PhaseBError("WDSat build receipt has invalid compile-time limits")
    if limits != frozen_build["limits"]:
        raise PhaseBError("WDSat build receipt limits differ from the frozen protocol")
    macro_names = {
        "max_anf_id": "__MAX_ANF_ID__",
        "max_degree": "__MAX_DEGREE__",
        "max_id": "__MAX_ID__",
        "max_buffer_size": "__MAX_BUFFER_SIZE__",
        "max_eq": "__MAX_EQ__",
        "max_eq_size": "__MAX_EQ_SIZE__",
        "max_xeq": "__MAX_XEQ__",
        "max_xeq_size": "__MAX_XEQ_SIZE__",
    }
    parsed_macros: dict[str, int] = {}
    for line in config.splitlines():
        match = re.fullmatch(r"\s*#define\s+(\S+)\s+(\d+)\s*", line)
        if match is not None:
            parsed_macros[match.group(1)] = int(match.group(2))
    for field, macro in macro_names.items():
        if parsed_macros.get(macro) != limits[field]:
            raise PhaseBError(f"WDSat build receipt limit {field} differs from {macro}")
    processes = value.get("build_processes")
    expected_roles = ["source-copy", "config-install", "clean", "build", "binary-copy"]
    if (
        not isinstance(processes, list)
        or [item.get("role") if isinstance(item, dict) else None for item in processes]
        != expected_roles
    ):
        raise PhaseBError("WDSat build receipt lacks separately charged build processes")
    required_metrics = {
        "wall_seconds",
        "user_seconds",
        "system_seconds",
        "total_core_seconds",
        "single_core_seconds",
        "peak_rss_bytes",
        "meter",
    }
    for process in processes:
        if (
            not isinstance(process, dict)
            or set(process)
            != {
                "role",
                "command",
                "returncode",
                "timed_out",
                "orphan_group_terminated",
                "metrics",
                "stdout_sha256",
                "stderr_sha256",
            }
            or process["returncode"] != 0
            or process["timed_out"] is not False
            or process["orphan_group_terminated"] is not False
            or not isinstance(process["role"], str)
            or not isinstance(process["command"], list)
            or not all(isinstance(argument, str) for argument in process["command"])
            or require_hex64(process["stdout_sha256"], "WDSat build stdout hash")
            != process["stdout_sha256"]
            or require_hex64(process["stderr_sha256"], "WDSat build stderr hash")
            != process["stderr_sha256"]
            or not isinstance(process["metrics"], dict)
            or set(process["metrics"]) != required_metrics
        ):
            raise PhaseBError("WDSat build receipt contains an invalid process record")
        metrics = process["metrics"]
        numeric = required_metrics - {"meter", "peak_rss_bytes"}
        if (
            any(
                isinstance(metrics[name], bool)
                or not isinstance(metrics[name], (int, float))
                or not math.isfinite(metrics[name])
                or metrics[name] < 0
                for name in numeric
            )
            or not isinstance(metrics["peak_rss_bytes"], int)
            or isinstance(metrics["peak_rss_bytes"], bool)
            or metrics["peak_rss_bytes"] < 0
            or metrics["meter"] != "fresh-process getrusage(RUSAGE_CHILDREN)"
            or not math.isclose(
                metrics["total_core_seconds"],
                metrics["user_seconds"] + metrics["system_seconds"],
                rel_tol=0,
                abs_tol=1e-9,
            )
            or metrics["single_core_seconds"] != metrics["total_core_seconds"]
        ):
            raise PhaseBError("WDSat build receipt contains invalid process resources")
    markers = forbidden_material(protocol)
    assert_no_forbidden_structured(value, markers, "WDSat build receipt")
    assert_no_forbidden_text(json.dumps(value, sort_keys=True), markers, "WDSat build receipt")
    return {
        "receipt_bytes": len(receipt_bytes),
        "receipt_sha256": sha256_bytes(receipt_bytes),
        "source_commit": value["source_commit"],
        "config_sha256": value["config_sha256"],
        "config": config,
        "binary_sha256": value["binary_sha256"],
        "limits": value["limits"],
        "build_processes": value["build_processes"],
        "claim_boundary": value["claim_boundary"],
    }


def wdsat_requirements(instance_root: Path, manifest: dict[str, Any]) -> dict[str, int]:
    cms = manifest["exports"]["cryptominisat_xor_dimacs"]
    n_source = int(manifest["source_variables"])
    max_degree = int(manifest["source_max_degree"])
    max_id = int(cms["variables"])
    max_eq = int(cms["cnf_clauses"])
    max_xeq = int(cms["xor_rows"])
    max_terms = int(manifest["source_max_monomials_per_equation"])
    source_xor_atoms = matrix.wdsat_xor_atom_count(instance_root / "instance.anf")
    return {
        "max_anf_id": n_source + 1,
        "max_degree": max_degree + 1,
        "max_id": max_id,
        "max_buffer_size": max(5000, max_terms * 8, max_id * 16, source_xor_atoms + 1),
        "max_eq": max(64, max_eq + 16),
        "max_eq_size": max_degree + 2,
        "max_xeq": max(8, max_xeq + 2),
        "max_xeq_size": max(max_id + 1, max_terms + 2),
    }


def require_wdsat_capacity(
    receipt: dict[str, Any], requirements_by_task: list[dict[str, Any]]
) -> dict[str, Any]:
    maxima = {
        name: max(item["requirements"][name] for item in requirements_by_task)
        for name in receipt["limits"]
    }
    insufficient = {
        name: {"required": maxima[name], "available": receipt["limits"][name]}
        for name in maxima
        if maxima[name] > receipt["limits"][name]
    }
    if insufficient:
        raise PhaseBError(f"WDSat binary is undersized for frozen exports: {insufficient}")
    by_cell: dict[str, dict[str, int]] = {}
    for item in requirements_by_task:
        cell = item["cell_id"]
        if cell not in by_cell:
            by_cell[cell] = dict(item["requirements"])
        else:
            for name, required in item["requirements"].items():
                by_cell[cell][name] = max(by_cell[cell][name], required)
    return {"panel_maxima": maxima, "per_cell_maxima": by_cell, "capacity_verified": True}


def safe_child_environment() -> dict[str, str]:
    return {
        "PATH": "/usr/bin:/bin:/usr/sbin:/sbin",
        "LANG": "C",
        "LC_ALL": "C",
        "TZ": "UTC",
        "TMPDIR": "/tmp",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "RAYON_NUM_THREADS": "1",
    }


def path_within(path: Path, root: Path) -> bool:
    try:
        path.resolve(strict=True).relative_to(root.resolve(strict=True))
        return True
    except (OSError, ValueError):
        return False


def input_receipts(paths: list[Path], task_root: Path, markers: set[str]) -> list[dict[str, Any]]:
    receipts = []
    for path in paths:
        if not path_within(path, task_root):
            raise PhaseBError(f"child input escapes its task root: {path}")
        data = regular_file_bytes(path, "child input")
        try:
            text = data.decode()
        except UnicodeDecodeError:
            text = ""
        assert_no_forbidden_text(text, markers, f"child input {path.name}")
        if path.suffix == ".json":
            value = parse_json_bytes(data, f"child input {path.name}")
            assert_no_forbidden_structured(value, markers, f"child input {path.name}")
        receipts.append(
            {
                "path": str(path.resolve()),
                "bytes": len(data),
                "sha256": sha256_bytes(data),
            }
        )
    return receipts


def assert_launch_boundary(
    command: list[str], environment: dict[str, str], inputs: list[dict[str, Any]], markers: set[str]
) -> None:
    assert_no_forbidden_text("\0".join(command), markers, "child command")
    assert_no_forbidden_text(
        "\0".join(f"{key}={value}" for key, value in sorted(environment.items())),
        markers,
        "child environment",
    )
    assert_no_forbidden_structured(inputs, markers, "child input receipt")


def run_metered(
    *,
    role: str,
    command: list[str],
    cwd: Path,
    task_root: Path,
    input_paths: list[Path],
    timeout: int,
    meter: Path,
    environment: dict[str, str],
    markers: set[str],
    expected_executable_sha256: str | None = None,
    allow_hardlinked_executable: bool = False,
) -> dict[str, Any]:
    def launch_executable_sha256() -> str:
        if not allow_hardlinked_executable:
            return sha256_file(Path(command[0]), f"{role} executable")
        resolved = Path(command[0]).resolve(strict=True)
        metadata = resolved.stat()
        if not stat.S_ISREG(metadata.st_mode) or not os.access(resolved, os.X_OK):
            raise PhaseBError(f"{role} executable is not a regular executable: {resolved}")
        return sha256_bytes(resolved.read_bytes())

    if not isinstance(allow_hardlinked_executable, bool):
        raise PhaseBError(f"{role} hard-link policy must be boolean")
    if expected_executable_sha256 is not None:
        require_hex64(expected_executable_sha256, f"{role} executable identity")
        if launch_executable_sha256() != expected_executable_sha256:
            raise PhaseBError(f"{role} executable changed before launch")
    intent_path = task_root / f"{role}.intent.json"
    stdout_path = task_root / f"{role}.stdout"
    stderr_path = task_root / f"{role}.stderr"
    metrics_path = task_root / f"{role}.metrics.json"
    for path in (intent_path, stdout_path, stderr_path, metrics_path):
        if path.exists() or path.is_symlink():
            raise PhaseBError(f"child artifact path already exists: {path}")
    inputs = input_receipts(input_paths, task_root, markers)
    assert_launch_boundary(command, environment, inputs, markers)
    intent = {
        "schema": "koblitz_pdp_phase_b_child_intent.v1",
        "role": role,
        "command": command,
        "cwd": str(cwd.resolve()),
        "watchdog_seconds": timeout,
        "stdin": "devnull",
        "close_fds": True,
        "environment": environment,
        "inputs": inputs,
        "executable_sha256": expected_executable_sha256,
    }
    write_json_new(intent_path, intent)
    meter_command = [
        str(Path(sys.executable).resolve()),
        str(meter.resolve()),
        "--cwd",
        str(cwd.resolve()),
        "--timeout",
        str(timeout),
        "--stdout",
        str(stdout_path.resolve()),
        "--stderr",
        str(stderr_path.resolve()),
        "--metrics",
        str(metrics_path.resolve()),
        "--",
        *command,
    ]
    completed = subprocess.run(
        meter_command,
        cwd=REPO,
        env=environment,
        stdin=subprocess.DEVNULL,
        close_fds=True,
        check=False,
    )
    if completed.returncode != 0 or not metrics_path.is_file():
        raise PhaseBError(f"process meter failed for {role} with return code {completed.returncode}")
    metrics, _ = read_json(metrics_path, f"{role} process metrics")
    if not isinstance(metrics, dict) or metrics.get("command") != command:
        raise PhaseBError(f"{role} process receipt does not bind its child command")
    if metrics.get("watchdog_seconds") != timeout:
        raise PhaseBError(f"{role} process receipt changed its watchdog")
    if (
        not isinstance(metrics.get("returncode"), int)
        or not isinstance(metrics.get("timed_out"), bool)
        or not isinstance(metrics.get("orphan_group_terminated"), bool)
    ):
        raise PhaseBError(f"{role} process receipt has invalid terminal fields")
    if metrics["orphan_group_terminated"]:
        raise PhaseBError(f"{role} left descendant processes after its leader exited")
    resources = metrics.get("metrics")
    required_resources = {
        "wall_seconds",
        "user_seconds",
        "system_seconds",
        "total_core_seconds",
        "single_core_seconds",
        "peak_rss_bytes",
        "meter",
    }
    if not isinstance(resources, dict) or set(resources) != required_resources:
        raise PhaseBError(f"{role} process receipt lacks required resources")
    numeric_names = {
        "wall_seconds",
        "user_seconds",
        "system_seconds",
        "total_core_seconds",
        "single_core_seconds",
    }
    if any(
        isinstance(resources[name], bool)
        or not isinstance(resources[name], (int, float))
        or not math.isfinite(resources[name])
        or resources[name] < 0
        for name in numeric_names
    ):
        raise PhaseBError(f"{role} process receipt has invalid numeric resources")
    if (
        not isinstance(resources["peak_rss_bytes"], int)
        or isinstance(resources["peak_rss_bytes"], bool)
        or resources["peak_rss_bytes"] < 0
        or resources["meter"] != "fresh-process getrusage(RUSAGE_CHILDREN)"
        or not math.isclose(
            resources["total_core_seconds"],
            resources["user_seconds"] + resources["system_seconds"],
            rel_tol=0,
            abs_tol=1e-9,
        )
        or resources["single_core_seconds"] != resources["total_core_seconds"]
    ):
        raise PhaseBError(f"{role} process receipt has inconsistent resource accounting")
    stdout = regular_file_bytes(stdout_path, f"{role} stdout")
    stderr = regular_file_bytes(stderr_path, f"{role} stderr")
    if expected_executable_sha256 is not None and launch_executable_sha256() != expected_executable_sha256:
        raise PhaseBError(f"{role} executable changed during launch")
    record = dict(metrics)
    record.update(
        {
            "role": role,
            "intent_path": intent_path.name,
            "stdout_path": stdout_path.name,
            "stderr_path": stderr_path.name,
            "metrics_path": metrics_path.name,
            "stdout_bytes": len(stdout),
            "stdout_sha256": sha256_bytes(stdout),
            "stderr_bytes": len(stderr),
            "stderr_sha256": sha256_bytes(stderr),
            "inputs": inputs,
            "stdout_text": stdout.decode(errors="replace"),
            "stderr_text": stderr.decode(errors="replace"),
        }
    )
    return record


def compact_process(record: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in record.items()
        if key not in {"stdout_text", "stderr_text"}
    }


def validate_export_manifest(
    manifest: Any,
    instance: dict[str, Any],
    instance_root: Path,
    protocol: dict[str, Any],
) -> dict[str, Any]:
    require_exact_instance_inventory(instance_root)
    if not isinstance(manifest, dict):
        raise PhaseBError("exporter manifest must be an object")
    markers = forbidden_material(protocol)
    assert_no_forbidden_structured(manifest, markers, "exporter manifest")
    assert_no_forbidden_text(json.dumps(manifest, sort_keys=True), markers, "exporter manifest")
    exact = {
        "n": instance["n"],
        "ell": instance["ell"],
        "m": instance["m"],
        "seed": instance["source_nonce"],
        "blind_instance_id": instance["blind_instance_id"],
        "target_mode": "explicit_affine",
        "curve_a": instance["curve_a"],
        "target": instance["target"],
    }
    for field, value in exact.items():
        if manifest.get(field) != value:
            raise PhaseBError(f"exporter manifest differs from blind instance at {field}")
    if manifest.get("kind") != "binary_koblitz_pdp_cross_solver_instance":
        raise PhaseBError("exporter manifest kind changed")
    predicate = manifest.get("factor_base_predicate")
    if not isinstance(predicate, dict):
        raise PhaseBError("exporter manifest lacks its algebraic factor-base predicate")
    expected_predicate = "polynomial_subspace" if instance["basis"] == "standard" else "ggmp_linearised_kernel"
    if predicate.get("kind") != expected_predicate:
        raise PhaseBError("exporter manifest used the wrong factor-base construction")
    if (
        predicate.get("enumerates_target_subgroup") is not False
        or predicate.get("uses_discrete_log_labels") is not False
    ):
        raise PhaseBError("exporter manifest violates the public algebraic factor-base contract")
    if instance["basis"] == "ggmp" and predicate.get("factor_index") != instance["factor_index"]:
        raise PhaseBError("exporter manifest used the wrong GGMP factor index")
    basis_bits = manifest.get("factor_base_basis_bitmasks")
    if not isinstance(basis_bits, list) or len(basis_bits) != instance["ell"]:
        raise PhaseBError("exporter manifest used the wrong factor-base dimension")
    if instance["basis"] == "standard" and basis_bits != [
        str(1 << index) for index in range(instance["ell"])
    ]:
        raise PhaseBError("standard factor-base basis is not the frozen polynomial subspace")
    if manifest.get("native_sat", {}).get("status") != "not_run_in_export_process":
        raise PhaseBError("explicit exporter ran the native solver")
    if manifest.get("direct_meet_in_the_middle", {}).get("status") != "not_run_in_export_process":
        raise PhaseBError("explicit exporter ran direct MITM")
    source = manifest.get("source_instance")
    identity = source.get("identity") if isinstance(source, dict) else None
    if not isinstance(source, dict) or source.get("schema") != "koblitz_pdp_source_instance.v2":
        raise PhaseBError("explicit exporter did not produce source identity v2")
    require_hex64(source.get("id_blake3"), "source_instance.id_blake3")
    if not isinstance(identity, dict) or identity.get("schema") != "koblitz_pdp_source_identity.v2":
        raise PhaseBError("explicit exporter did not produce source identity payload v2")
    for field in (
        "n",
        "ell",
        "m",
        "seed",
        "blind_instance_id",
        "target_mode",
        "curve_a",
        "target",
        "representation",
        "source_variables",
        "source_equations",
        "exports",
    ):
        if identity.get(field) != manifest.get(field):
            raise PhaseBError(f"source identity v2 differs from its manifest at {field}")
    try:
        snapshot = matrix.source_artifact_snapshot(instance_root, manifest)
    except (OSError, ValueError) as error:
        raise PhaseBError(f"exported source artifact verification failed: {error}") from error
    manifest_bytes = regular_file_bytes(instance_root / "manifest.json", "explicit-target manifest")
    return {
        "manifest": {
            "path": "manifest.json",
            "bytes": len(manifest_bytes),
            "sha256": sha256_bytes(manifest_bytes),
        },
        "exports": snapshot,
    }


def frozen_source_snapshot(instance_root: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    require_exact_instance_inventory(instance_root)
    manifest_bytes = regular_file_bytes(instance_root / "manifest.json", "explicit-target manifest")
    return {
        "manifest": {
            "path": "manifest.json",
            "bytes": len(manifest_bytes),
            "sha256": sha256_bytes(manifest_bytes),
        },
        "exports": matrix.source_artifact_snapshot(instance_root, manifest),
    }


def require_source_unchanged(
    instance_root: Path, manifest: dict[str, Any], before: dict[str, Any], context: str
) -> None:
    try:
        after = frozen_source_snapshot(instance_root, manifest)
    except (OSError, ValueError) as error:
        raise PhaseBError(f"source custody failed after {context}: {error}") from error
    if after != before:
        raise PhaseBError(f"source artifacts changed during {context}")


def require_exact_instance_inventory(instance_root: Path) -> None:
    try:
        metadata = instance_root.lstat()
    except OSError as error:
        raise PhaseBError(f"cannot stat exported instance root: {error}") from error
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise PhaseBError("exported instance root must be a real directory")
    expected = {"manifest.json", "instance.anf", "instance.xor.cnf", "instance.magma"}
    actual: set[str] = set()
    for child in instance_root.iterdir():
        data = regular_file_bytes(child, "exported instance artifact")
        if not data:
            # Empty source files are most likely a truncated producer artifact.
            raise PhaseBError(f"exported instance artifact is empty: {child.name}")
        actual.add(child.name)
    if actual != expected:
        raise PhaseBError(f"exported instance inventory changed: {sorted(actual)}")


def manifest_source_paths(instance_root: Path, manifest: dict[str, Any]) -> list[Path]:
    exports = manifest.get("exports", {})
    return [
        instance_root / exports[name]["path"]
        for name in ("wdsat_anf", "cryptominisat_xor_dimacs", "magma_boolean_f4")
    ]


def process_record_for_matrix(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "stdout": record["stdout_text"],
        "stderr": record["stderr_text"],
        "returncode": record["returncode"],
        "timed_out": record["timed_out"],
        "metrics": record["metrics"],
        "command": record["command"],
    }


def source_verification_result(
    record: dict[str, Any], manifest: dict[str, Any]
) -> dict[str, Any]:
    try:
        report = parse_json_bytes(record["stdout_text"].encode(), "source verification stdout")
    except PhaseBError:
        report = None
    expected_id = manifest["source_instance"]["id_blake3"]
    artifacts = report.get("source_artifacts") if isinstance(report, dict) else None
    valid_artifacts = (
        isinstance(artifacts, dict)
        and set(artifacts)
        == {"wdsat_anf", "cryptominisat_xor_dimacs", "magma_boolean_f4"}
        and all(isinstance(value, dict) and value.get("valid") is True for value in artifacts.values())
    )
    accepted = (
        not record["timed_out"]
        and record["returncode"] == 0
        and isinstance(report, dict)
        and report.get("schema") == "koblitz_pdp_source_verification.v1"
        and report.get("status") == "verified"
        and report.get("source_instance_id") == expected_id
        and report.get("source_instance_verified") is True
        and report.get("regenerated_source_exact") is True
        and valid_artifacts
    )
    if not accepted:
        raise PhaseBError("isolated source verifier rejected the explicit-target source")
    return {
        "status": "verified",
        "source_instance_id": expected_id,
        "source_instance_verified": True,
        "regenerated_source_exact": True,
        "source_artifacts": artifacts,
        "process": compact_process(record),
        "report": report,
    }


def native_result(record: dict[str, Any], manifest: dict[str, Any]) -> dict[str, Any]:
    row = matrix.isolated_backend_status(process_record_for_matrix(record), "native-sat", manifest)
    row["solver"] = "native-xor"
    row["process"] = compact_process(record)
    return row


def validate_external_model(
    *,
    solver: str,
    assignment: list[bool],
    manifest: dict[str, Any],
    manifest_path: Path,
    backend: Path,
    task_root: Path,
    instance_root: Path,
    timeout: int,
    meter: Path,
    environment: dict[str, str],
    markers: set[str],
    backend_sha256: str,
) -> dict[str, Any]:
    assignment_path = task_root / f"{solver}.source-model.json"
    assignment_bytes = (json.dumps(assignment, separators=(",", ":")) + "\n").encode()
    write_new(assignment_path, assignment_bytes)
    record = run_metered(
        role=f"{solver}.point-validation",
        command=[
            str(backend.resolve()),
            "validate-model",
            str(manifest_path.resolve()),
            str(assignment_path.resolve()),
        ],
        cwd=instance_root,
        task_root=task_root,
        input_paths=[manifest_path, *manifest_source_paths(instance_root, manifest), assignment_path],
        timeout=timeout,
        meter=meter,
        environment=environment,
        markers=markers,
        expected_executable_sha256=backend_sha256,
    )
    status = matrix.assignment_validation_status(
        process_record_for_matrix(record), manifest, assignment_path
    )
    return {
        "status": status["status"],
        "assignment_path": assignment_path.name,
        "assignment_bytes": len(assignment_bytes),
        "assignment_sha256": sha256_bytes(assignment_bytes),
        "process": compact_process(record),
    }


def external_result(
    *,
    solver: str,
    record: dict[str, Any],
    manifest: dict[str, Any],
    manifest_path: Path,
    backend: Path,
    task_root: Path,
    instance_root: Path,
    timeout: int,
    meter: Path,
    environment: dict[str, str],
    markers: set[str],
    backend_sha256: str,
) -> dict[str, Any]:
    matrix_record = process_record_for_matrix(record)
    row = matrix.solver_status(matrix_record, solver, instance_root, manifest)
    lines = [line.strip() for line in record["stdout_text"].splitlines() if line.strip()]
    if not record["timed_out"] and solver == "wdsat":
        model = matrix.parse_wdsat_model(record["stdout_text"], manifest["source_variables"])
        model_lines = [
            line
            for line in lines
            if len(line) >= manifest["source_variables"] and set(line) <= {"0", "1"}
        ]
        exact_unsat = [line for line in lines if line == "UNSAT"]
        exact_unknown = [line for line in lines if line in {"UNKNOWN", "s UNKNOWN"}]
        if model is not None:
            if (
                len(model_lines) != 1
                or exact_unsat
                or exact_unknown
                or record["returncode"] != 0
            ):
                row["status"] = "solver_error"
        else:
            if len(exact_unsat) == 1 and not exact_unknown and record["returncode"] == 0:
                row["status"] = "unsat"
            elif len(exact_unknown) == 1 and not exact_unsat and record["returncode"] == 0:
                row["status"] = "unknown_inconclusive"
            else:
                row["status"] = "solver_error"
    elif not record["timed_out"] and solver == "cryptominisat":
        sat_lines = [line for line in lines if line == "s SATISFIABLE"]
        unsat_lines = [line for line in lines if line == "s UNSATISFIABLE"]
        unknown_lines = [line for line in lines if line in {"UNKNOWN", "s UNKNOWN"}]
        terminals = len(sat_lines) + len(unsat_lines) + len(unknown_lines)
        if terminals != 1:
            row["status"] = "solver_error"
        elif unknown_lines:
            row["status"] = (
                "unknown_inconclusive" if record["returncode"] == 0 else "solver_error"
            )
        elif sat_lines and record["returncode"] != 10:
            row["status"] = "sat_invalid_terminal_status"
        elif unsat_lines and record["returncode"] != 20:
            row["status"] = "unsat_invalid_terminal_status"
    validation = None
    if row["status"] == "sat_source_model_unverified_point_witness":
        assignment = matrix.parsed_source_assignment(matrix_record, solver, manifest)
        if assignment is None:
            row["status"] = "sat_invalid_model"
        else:
            validation = validate_external_model(
                solver=solver,
                assignment=assignment,
                manifest=manifest,
                manifest_path=manifest_path,
                backend=backend,
                task_root=task_root,
                instance_root=instance_root,
                timeout=timeout,
                meter=meter,
                environment=environment,
                markers=markers,
                backend_sha256=backend_sha256,
            )
            row["source_witness_valid"] = validation["status"] == "valid_point_witness"
            row["status"] = (
                "sat"
                if validation["status"] == "valid_point_witness"
                else "sat_nonlifting_model_inconclusive"
                if validation["status"] == "nonlifting_source_model"
                else "sat_point_validation_timeout_inconclusive"
                if validation["status"] == "timeout_inconclusive"
                else "sat_point_validation_error"
            )
    row["process"] = compact_process(record)
    row["point_witness_validation"] = validation
    return row


def task_directory(output: Path, index: int, blind_id: str) -> Path:
    return output / "tasks" / f"{index:06d}-{blind_id}"


def export_one_task(
    *,
    index: int,
    instance: dict[str, Any],
    output: Path,
    protocol: dict[str, Any],
    tools: dict[str, Path],
    identities: dict[str, dict[str, Any]],
    environment: dict[str, str],
) -> dict[str, Any]:
    task_root = task_directory(output, index, instance["blind_instance_id"])
    task_root.mkdir(parents=True)
    markers = forbidden_material(protocol)
    timeout = protocol["execution"]["per_process_watchdog_seconds"]
    conflict_budget = protocol["execution"]["native_conflict_budget"]
    instance_root = task_root / "instance"
    exporter_command = [
        str(tools["exporter"].resolve()),
        str(instance["n"]),
        str(instance["ell"]),
        instance["basis"],
        str(instance["source_nonce"]),
        str(conflict_budget),
        str(instance_root.resolve()),
        str(instance["curve_a"]),
        str(instance["factor_index"]),
        "--target-x",
        instance["target"]["x"],
        "--target-y",
        instance["target"]["y"],
        "--blind-instance-id",
        instance["blind_instance_id"],
        "--export-only",
    ]
    export_record = run_metered(
        role="export",
        command=exporter_command,
        cwd=REPO,
        task_root=task_root,
        input_paths=[],
        timeout=timeout,
        meter=tools["meter"],
        environment=environment,
        markers=markers,
        expected_executable_sha256=identities["exporter"]["sha256"],
    )
    if export_record["timed_out"] or export_record["returncode"] != 0:
        raise PhaseBError(
            f"explicit exporter failed for {instance['blind_instance_id']}: "
            f"returncode={export_record['returncode']} timed_out={export_record['timed_out']}"
        )
    manifest_path = instance_root / "manifest.json"
    manifest, _ = read_json(manifest_path, "explicit-target manifest")
    source_before = validate_export_manifest(manifest, instance, instance_root, protocol)
    source_verification_record = run_metered(
        role="source-verification",
        command=[
            str(tools["backend"].resolve()),
            "verify-source",
            str(manifest_path.resolve()),
        ],
        cwd=instance_root,
        task_root=task_root,
        input_paths=[manifest_path, *manifest_source_paths(instance_root, manifest)],
        timeout=timeout,
        meter=tools["meter"],
        environment=environment,
        markers=markers,
        expected_executable_sha256=identities["backend"]["sha256"],
    )
    require_source_unchanged(
        instance_root, manifest, source_before, "isolated source verification"
    )
    source_verification = source_verification_result(source_verification_record, manifest)
    requirements = wdsat_requirements(instance_root, manifest)
    export_result: dict[str, Any] = {
        "schema": TASK_SCHEMA,
        "ordinal": index,
        "blind_instance_id": instance["blind_instance_id"],
        "source_system_id": instance["source_system_id"],
        "cell_id": instance["cell_id"],
        "target": instance["target"],
        "source_instance_id": manifest["source_instance"]["id_blake3"],
        "export": compact_process(export_record),
        "source_artifacts_before": source_before,
        "source_verification": source_verification,
        "wdsat_requirements": requirements,
    }
    write_json_new(task_root / "export-result.json", export_result)
    return {
        "index": index,
        "instance": instance,
        "task_root": task_root,
        "instance_root": instance_root,
        "manifest_path": manifest_path,
        "manifest": manifest,
        "source_before": source_before,
        "export_result": export_result,
        "requirements": requirements,
    }


def solve_one_task(
    *,
    prepared: dict[str, Any],
    protocol: dict[str, Any],
    tools: dict[str, Path],
    identities: dict[str, dict[str, Any]],
    environment: dict[str, str],
) -> dict[str, Any]:
    index = prepared["index"]
    instance = prepared["instance"]
    task_root = prepared["task_root"]
    instance_root = prepared["instance_root"]
    manifest_path = prepared["manifest_path"]
    manifest = prepared["manifest"]
    source_before = prepared["source_before"]
    markers = forbidden_material(protocol)
    timeout = protocol["execution"]["per_process_watchdog_seconds"]
    conflict_budget = protocol["execution"]["native_conflict_budget"]
    task_result: dict[str, Any] = dict(prepared["export_result"])
    task_result["backends"] = []
    native_record = run_metered(
        role="native-xor",
        command=[
            str(tools["backend"].resolve()),
            "native-sat",
            str(manifest_path.resolve()),
            str(conflict_budget),
        ],
        cwd=instance_root,
        task_root=task_root,
        input_paths=[manifest_path, *manifest_source_paths(instance_root, manifest)],
        timeout=timeout,
        meter=tools["meter"],
        environment=environment,
        markers=markers,
        expected_executable_sha256=identities["backend"]["sha256"],
    )
    require_source_unchanged(instance_root, manifest, source_before, "native-XOR execution")
    task_result["backends"].append(native_result(native_record, manifest))
    branch_variables = ",".join(str(value) for value in range(1, 3 * instance["ell"] + 1))
    wdsat_record = run_metered(
        role="wdsat",
        command=[
            str(tools["wdsat"].resolve()),
            "-i",
            str((instance_root / "instance.anf").resolve()),
            "-g",
            branch_variables,
        ],
        cwd=instance_root,
        task_root=task_root,
        input_paths=[instance_root / "instance.anf"],
        timeout=timeout,
        meter=tools["meter"],
        environment=environment,
        markers=markers,
        expected_executable_sha256=identities["wdsat"]["sha256"],
    )
    require_source_unchanged(instance_root, manifest, source_before, "WDSat execution")
    task_result["backends"].append(
        external_result(
            solver="wdsat",
            record=wdsat_record,
            manifest=manifest,
            manifest_path=manifest_path,
            backend=tools["backend"],
            task_root=task_root,
            instance_root=instance_root,
            timeout=timeout,
            meter=tools["meter"],
            environment=environment,
            markers=markers,
            backend_sha256=identities["backend"]["sha256"],
        )
    )
    require_source_unchanged(
        instance_root, manifest, source_before, "WDSat point-witness validation"
    )
    cms_record = run_metered(
        role="cryptominisat",
        command=[
            str(tools["cryptominisat"].resolve()),
            "--verb",
            "1",
            "--threads",
            "1",
            str((instance_root / "instance.xor.cnf").resolve()),
        ],
        cwd=instance_root,
        task_root=task_root,
        input_paths=[instance_root / "instance.xor.cnf"],
        timeout=timeout,
        meter=tools["meter"],
        environment=environment,
        markers=markers,
        expected_executable_sha256=identities["cryptominisat"]["sha256"],
    )
    require_source_unchanged(instance_root, manifest, source_before, "CryptoMiniSat execution")
    task_result["backends"].append(
        external_result(
            solver="cryptominisat",
            record=cms_record,
            manifest=manifest,
            manifest_path=manifest_path,
            backend=tools["backend"],
            task_root=task_root,
            instance_root=instance_root,
            timeout=timeout,
            meter=tools["meter"],
            environment=environment,
            markers=markers,
            backend_sha256=identities["backend"]["sha256"],
        )
    )
    require_source_unchanged(
        instance_root, manifest, source_before, "CryptoMiniSat point-witness validation"
    )
    source_after = frozen_source_snapshot(instance_root, manifest)
    if source_after != source_before:
        raise PhaseBError(f"a solver modified source artifacts for {instance['blind_instance_id']}")
    if [row["solver"] for row in task_result["backends"]] != list(BACKENDS):
        raise PhaseBError("backend result order changed")
    task_result["source_artifacts_after"] = source_after
    task_result["source_artifacts_unchanged"] = True
    task_result["status"] = "terminal_outputs_recorded"
    write_json_new(task_root / "task-result.json", task_result)
    return task_result


def all_regular_inventory(root: Path, excluded: set[str] | None = None) -> list[dict[str, Any]]:
    excluded = excluded or set()
    inventory: list[dict[str, Any]] = []
    for current, names, filenames in os.walk(root, followlinks=False):
        current_path = Path(current)
        for name in names:
            path = current_path / name
            metadata = path.lstat()
            if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
                raise PhaseBError(f"result tree contains an unsafe directory entry: {path}")
        for name in filenames:
            path = current_path / name
            relative = path.relative_to(root).as_posix()
            if relative in excluded:
                continue
            data = regular_file_bytes(path, "result artifact")
            inventory.append({"path": relative, "bytes": len(data), "sha256": sha256_bytes(data)})
    inventory.sort(key=lambda item: item["path"])
    return inventory


def task_processes(task: dict[str, Any]) -> list[dict[str, Any]]:
    processes = [task["export"], task["source_verification"]["process"]]
    for row in task["backends"]:
        processes.append(row["process"])
        validation = row.get("point_witness_validation")
        if isinstance(validation, dict):
            processes.append(validation["process"])
    return processes


def process_metrics_from_tasks(
    output: Path, task_results: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    expected_paths: set[Path] = set()
    records = []
    for task in task_results:
        root = task_directory(output, task["ordinal"], task["blind_instance_id"])
        for process in task_processes(task):
            relative = Path(process["metrics_path"])
            if relative.is_absolute() or len(relative.parts) != 1:
                raise PhaseBError("task process metrics path is not one local filename")
            path = root / relative
            if path in expected_paths:
                raise PhaseBError("task records reuse one process metrics file")
            expected_paths.add(path)
            record, _ = read_json(path, "process metrics")
            if (
                record.get("command") != process.get("command")
                or record.get("returncode") != process.get("returncode")
                or record.get("timed_out") != process.get("timed_out")
                or record.get("metrics") != process.get("metrics")
            ):
                raise PhaseBError("task process summary differs from its metrics file")
            records.append(record)
    actual_paths = set(output.glob("tasks/**/*.metrics.json"))
    if actual_paths != expected_paths:
        raise PhaseBError("process metric files differ from the exact task-record inventory")
    return records


def summarize_run(
    output: Path,
    protocol: dict[str, Any],
    selected_count: int,
    full_count: int,
    task_results: list[dict[str, Any]],
    wdsat_build: dict[str, Any],
    tool_builds: dict[str, Any] | None = None,
) -> dict[str, Any]:
    tool_builds = tool_builds or {}
    statuses: dict[str, Counter[str]] = {backend: Counter() for backend in BACKENDS}
    conflict_values: dict[str, list[int]] = {backend: [] for backend in BACKENDS}
    for task in task_results:
        for row in task["backends"]:
            statuses[row["solver"]][row["status"]] += 1
            if row.get("conflicts") is not None:
                conflict_values[row["solver"]].append(row["conflicts"])
    process_records = process_metrics_from_tasks(output, task_results)
    resources = [record["metrics"] for record in process_records]
    backend_outcomes = sum(sum(counter.values()) for counter in statuses.values())
    build_metrics = [process["metrics"] for process in wdsat_build["build_processes"]]
    return {
        "schema": RUN_SUMMARY_SCHEMA,
        "claim_boundary": "Balanced public toy PDP solver evidence only; no end-to-end index calculus or SOTA claim",
        "selected_instances": selected_count,
        "frozen_instances": full_count,
        "task_records": len(task_results),
        "backend_outcomes": backend_outcomes,
        "expected_selected_backend_outcomes": selected_count * len(BACKENDS),
        "full_frozen_backend_outcomes": protocol["execution"]["expected_backend_runs"],
        "status_counts": {backend: dict(sorted(counter.items())) for backend, counter in statuses.items()},
        "conflicts": {
            backend: {
                "reported": len(values),
                "values": sorted(values),
                "values_sha256": canonical_sha256(sorted(values)),
                "sum": sum(values),
                "min": min(values) if values else None,
                "median": statistics.median(values) if values else None,
                "max": max(values) if values else None,
            }
            for backend, values in conflict_values.items()
        },
        "process_receipts": len(process_records),
        "charged_process_resources": {
            "total_core_seconds": round(math.fsum(item["total_core_seconds"] for item in resources), 12),
            "single_core_seconds": round(math.fsum(item["total_core_seconds"] for item in resources), 12),
            "single_core_seconds_alias_of": "total_core_seconds",
            "single_core_elapsed_seconds": None,
            "summed_process_wall_seconds": round(math.fsum(item["wall_seconds"] for item in resources), 12),
            "peak_rss_bytes": max((item["peak_rss_bytes"] for item in resources), default=0),
            "peak_rss_scope": "maximum fresh-process high-water mark, not aggregate parallel memory",
        },
        "separately_charged_wdsat_build": {
            "receipt_sha256": wdsat_build["receipt_sha256"],
            "total_core_seconds": round(math.fsum(item["total_core_seconds"] for item in build_metrics), 12),
            "single_core_seconds": round(math.fsum(item["total_core_seconds"] for item in build_metrics), 12),
            "single_core_seconds_alias_of": "total_core_seconds",
            "single_core_elapsed_seconds": None,
            "summed_process_wall_seconds": round(math.fsum(item["wall_seconds"] for item in build_metrics), 12),
            "peak_rss_bytes": max((item["peak_rss_bytes"] for item in build_metrics), default=0),
            "peak_rss_scope": "maximum build-child high-water mark, not aggregate parallel memory",
        },
        "separately_charged_tool_builds": {
            name: {"receipt_sha256": value["receipt_sha256"], **value["receipt"]["resources"]}
            for name, value in sorted(tool_builds.items())
        },
        "rust_and_cryptominisat_build_receipts_bound": set(tool_builds) == {"rust", "cryptominisat"},
        "selection_complete": len(task_results) == selected_count
        and backend_outcomes == selected_count * len(BACKENDS),
        "full_panel_complete": selected_count == full_count
        and len(task_results) == full_count
        and backend_outcomes == protocol["execution"]["expected_backend_runs"],
        "truth_scoring_status": "withheld_until_separate_post_run_step",
        "native_xor_validation_accounting": "Each native backend process receipt includes source regeneration, solving, source-model checking, and exact lifted point-witness checking",
        "full_cost_gate_passed": False,
        "full_cost_blockers": tool_build_cost_blockers(tool_builds),
    }


def tool_build_cost_blockers(tool_builds: dict[str, Any]) -> list[str]:
    missing = [f"fresh metered {name} build receipt not yet bound" for name in ("rust", "cryptominisat") if name not in tool_builds]
    return missing + [
        "measured single-core elapsed time is unavailable; legacy single_core_seconds aliases total CPU",
        "outer build-driver overhead and simultaneous aggregate build memory require separate measurement",
        "Phase-B costs exclude Phase-A preparation and the remaining end-to-end index-calculus stages",
    ]


def run_panel(args: argparse.Namespace) -> dict[str, Any]:
    import build_koblitz_phase_b_tools as tool_builder
    import build_koblitz_phase_b_wdsat as wdsat_builder

    protocol, protocol_bytes = read_json(args.protocol, "Phase-B protocol")
    validate_protocol(protocol)
    require_phase_a_ancestor(protocol["phase_a_binding"]["source_revision"])
    binding, bundle = load_solver_root(args.solver_root, protocol)
    if args.output.exists() or args.output.is_symlink():
        raise PhaseBError(f"run output must be new: {args.output}")
    state = git_state()
    full_count = len(bundle["instances"])
    if args.max_instances is not None and args.max_instances <= 0:
        raise PhaseBError("--max-instances must be positive")
    selected = bundle["instances"][: args.max_instances] if args.max_instances else bundle["instances"]
    full_selection = len(selected) == full_count
    if state["dirty"] and not args.allow_dirty:
        raise PhaseBError("solver execution requires a clean checkout; dirty runs need --allow-dirty")
    tools = {
        "exporter": args.exporter.resolve(strict=True),
        "backend": args.backend.resolve(strict=True),
        "wdsat": args.wdsat.resolve(strict=True),
        "cryptominisat": args.cryptominisat.resolve(strict=True),
        "meter": args.meter.resolve(strict=True),
    }
    identities = {
        name: executable_identity(path, name, executable=name != "meter")
        for name, path in tools.items()
    }
    wdsat_build = wdsat_builder.validate_capsule(
        args.wdsat_build_seal.resolve(),
        protocol,
        identities["wdsat"],
        state,
        require_clean_implementation=full_selection and not args.allow_dirty,
    )
    tool_receipt_paths = {
        name: path for name, path in (
            ("rust", args.rust_build_receipt), ("cryptominisat", args.cryptominisat_build_receipt)
        ) if path is not None
    }
    if full_selection and not args.allow_dirty and set(tool_receipt_paths) != {"rust", "cryptominisat"}:
        raise PhaseBError("a full production selection requires both Rust and CryptoMiniSat build receipts")
    source_objects = tool_builder.rust_source_objects(REPO) if "rust" in tool_receipt_paths else None
    tool_builds = {
        name: tool_builder.validate_receipt(
            path.resolve(), name, identities,
            source_objects if name == "rust" else None,
            state,
            require_clean_implementation=full_selection and not args.allow_dirty,
        )
        for name, path in tool_receipt_paths.items()
    }
    environment = safe_child_environment()
    markers = forbidden_material(protocol)
    assert_launch_boundary([], environment, [], markers)
    evidence_class = (
        PRODUCTION_EVIDENCE_CLASS
        if full_selection and not state["dirty"] and not args.allow_dirty
        else "operational_smoke"
    )
    plan = {
        "schema": RUN_PLAN_SCHEMA,
        "created_at": now(),
        "protocol_path": str(args.protocol.resolve()),
        "protocol_sha256": sha256_bytes(protocol_bytes),
        "source_revision": state,
        "allow_dirty_requested": args.allow_dirty,
        "solver_input_binding": binding,
        "tool_identities": identities,
        "wdsat_build": wdsat_build,
        "wdsat_build_capsule_bound": True,
        "tool_build_accounting": protocol["tool_build_accounting"],
        "additional_tool_builds": tool_builds,
        "rust_build_source_objects": source_objects,
        "child_environment": environment,
        "stdin": "devnull",
        "close_fds": True,
        "watchdog_seconds": protocol["execution"]["per_process_watchdog_seconds"],
        "single_thread_requested": True,
        "instance_order": "authenticated_bundle_presentation_order",
        "selected_instance_count": len(selected),
        "full_instance_count": full_count,
        "selected_blind_instance_ids": [item["blind_instance_id"] for item in selected],
        "backend_order": list(BACKENDS),
        "planned_backend_outcomes": len(selected) * len(BACKENDS),
        "evidence_class": evidence_class,
        "truth_scoring_status": "withheld_until_separate_post_run_step",
        "outer_expected_command": [
            str(Path(sys.executable).resolve()),
            str(Path(__file__).resolve()),
            *sys.argv[1:],
        ],
    }
    assert_no_forbidden_structured(plan, markers, "execution plan")
    assert_no_forbidden_text(json.dumps(plan, sort_keys=True), markers, "execution plan")
    if args.plan:
        return plan
    args.output.mkdir(parents=True)
    (args.output / "tasks").mkdir()
    write_json_new(args.output / "execution-plan.json", plan)
    frozen_wdsat = wdsat_builder.copy_capsule(
        args.wdsat_build_seal.resolve(),
        args.output / "tool-builds" / "wdsat",
        protocol,
        identities["wdsat"],
        state,
        require_clean_implementation=full_selection and not args.allow_dirty,
    )
    if frozen_wdsat != wdsat_build:
        raise PhaseBError("WDSat build capsule changed while freezing the run")
    for name, path in tool_receipt_paths.items():
        destination = args.output / "tool-builds" / name
        tool_builder.copy_capsule(path.resolve(), destination)
        if tool_builder.validate_receipt(
            destination / "receipt.json", name, identities,
            source_objects if name == "rust" else None,
            state,
            require_clean_implementation=full_selection and not args.allow_dirty,
        ) != tool_builds[name]:
            raise PhaseBError("tool build receipt changed while freezing the run capsule")
    prepared_tasks = []
    for index, instance in enumerate(selected):
        prepared_tasks.append(
            export_one_task(
                index=index,
                instance=instance,
                output=args.output,
                protocol=protocol,
                tools=tools,
                identities=identities,
                environment=environment,
            )
        )
    requirements_by_task = [
        {
            "blind_instance_id": item["instance"]["blind_instance_id"],
            "cell_id": item["instance"]["cell_id"],
            "requirements": item["requirements"],
        }
        for item in prepared_tasks
    ]
    capacity = require_wdsat_capacity(wdsat_build, requirements_by_task)
    require_git_state(state, "the export-only sizing pass")
    export_inventory = {
        "schema": "koblitz_pdp_phase_b_frozen_export_inventory.v1",
        "status": "all_selected_explicit_targets_exported_before_solver_execution",
        "selected_instances": len(prepared_tasks),
        "source_instance_ids": [
            item["manifest"]["source_instance"]["id_blake3"] for item in prepared_tasks
        ],
        "wdsat_capacity": capacity,
        "source_artifacts": [
            {
                "blind_instance_id": item["instance"]["blind_instance_id"],
                "source_instance_id": item["manifest"]["source_instance"]["id_blake3"],
                "artifacts": item["source_before"],
            }
            for item in prepared_tasks
        ],
    }
    write_json_new(args.output / "export-inventory.json", export_inventory)
    task_results = [
        solve_one_task(
            prepared=item,
            protocol=protocol,
            tools=tools,
            identities=identities,
            environment=environment,
        )
        for item in prepared_tasks
    ]
    require_git_state(state, "solver execution")
    summary = summarize_run(
        args.output, protocol, len(selected), full_count, task_results, wdsat_build, tool_builds
    )
    write_json_new(args.output / "run-summary.json", summary)
    inventory = all_regular_inventory(args.output, {"run-seal.json"})
    seal_payload = {
        "schema": RUN_SEAL_SCHEMA,
        "status": "solver_outputs_frozen",
        "created_at": now(),
        "protocol_sha256": sha256_bytes(protocol_bytes),
        "source_phase_a_seal_sha256": binding["source_phase_a_seal_sha256"],
        "blind_bundle_sha256": binding["blind_bundle_sha256"],
        "selected_instance_count": len(selected),
        "full_instance_count": full_count,
        "backend_outcomes": summary["backend_outcomes"],
        "full_panel_complete": summary["full_panel_complete"],
        "inventory": inventory,
        "inventory_sha256": canonical_sha256(inventory),
        "truth_scoring_status": "withheld_until_separate_post_run_step",
    }
    seal = dict(seal_payload)
    seal["seal_payload_sha256"] = canonical_sha256(seal_payload)
    write_json_new(args.output / "run-seal.json", seal)
    return seal


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    stage = subparsers.add_parser("stage", help="create a sanitized solver-only input root")
    stage.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    stage.add_argument("--phase-a-seal", type=Path, required=True)
    stage.add_argument("--blind-bundle", type=Path, required=True)
    stage.add_argument("--output", type=Path, required=True)
    run = subparsers.add_parser("run", help="run and seal the blinded solver panel")
    run.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    run.add_argument("--solver-root", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--exporter", type=Path, required=True)
    run.add_argument("--backend", type=Path, required=True)
    run.add_argument("--wdsat", type=Path, required=True)
    run.add_argument("--wdsat-build-seal", type=Path, required=True)
    run.add_argument("--cryptominisat", type=Path, required=True)
    run.add_argument("--rust-build-receipt", type=Path)
    run.add_argument("--cryptominisat-build-receipt", type=Path)
    run.add_argument("--meter", type=Path, default=DEFAULT_METER)
    run.add_argument("--max-instances", type=int)
    run.add_argument("--allow-dirty", action="store_true")
    run.add_argument("--plan", action="store_true")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.command == "stage":
            result = stage_solver_root(
                args.protocol.resolve(),
                args.phase_a_seal.resolve(),
                args.blind_bundle.resolve(),
                args.output.resolve(),
            )
        else:
            args.protocol = args.protocol.resolve()
            args.solver_root = args.solver_root.resolve()
            args.output = args.output.resolve()
            result = run_panel(args)
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, PhaseBError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
