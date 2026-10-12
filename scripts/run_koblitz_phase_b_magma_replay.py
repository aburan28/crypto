#!/usr/bin/env python3
"""Prepare, run, seal, and score the licensed-host Phase-B Magma replay."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import stat
import subprocess
import sys
from typing import Any

import koblitz_phase_b_magma_core as legacy
import run_koblitz_blind_pdp_phase_b as phase_b
import verify_stage20_phase_b_terminal_evidence as terminal


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
DEFAULT_PROTOCOL = STAGE / "stage-22-phase-b-magma-replay-protocol.json"
DEFAULT_BUNDLE = STAGE / "stage-20-phase-b-terminal-evidence-successor-04-20260910"
DEFAULT_RUN_ROOT = Path("/Volumes/SSD990/koblitz-balanced-pdp-phase-b-run-successor-01-20260910")
DEFAULT_SCORE_SEAL = Path("/Volumes/SSD990/koblitz-balanced-pdp-phase-b-score-successor-01-20260910/score-seal.json")
DEFAULT_METER = REPO / "scripts/process_meter.py"
DEFAULT_BACKEND = REPO / "target/release/examples/koblitz_pdp_backend"

PROTOCOL_SCHEMA = "koblitz_phase_b_magma_replay_protocol.v1"
PACKET_MANIFEST_SCHEMA = "koblitz_phase_b_magma_replay_packet.v1"
PACKET_SEAL_SCHEMA = "koblitz_phase_b_magma_replay_packet_seal.v1"
PREFLIGHT_SCHEMA = "koblitz_phase_b_magma_replay_preflight.v1"
TASK_SCHEMA = "koblitz_phase_b_magma_replay_task.v1"
SUMMARY_SCHEMA = "koblitz_phase_b_magma_replay_summary.v1"
RETURN_SEAL_SCHEMA = "koblitz_phase_b_magma_replay_return_seal.v1"
SCORE_SCHEMA = "koblitz_phase_b_magma_replay_score.v1"
SCORE_SEAL_SCHEMA = "koblitz_phase_b_magma_replay_score_seal.v1"
HEX64 = re.compile(r"^[0-9a-f]{64}$")
BLIND_ID = re.compile(r"^b-[0-9a-f]{64}$")
FORBIDDEN_INPUT_FIELDS = {
    "target_class", "witness_indices", "planted_points", "oracle",
    "oracle-ledger", "sealed-oracle", "unconditional_draw_index",
    "selection_priority_hex", "decomposable", "nondecomposable",
}


class ReplayError(RuntimeError):
    """A replay packet, execution, or score violates its frozen contract."""


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def regular_bytes(path: Path, context: str) -> bytes:
    try:
        metadata = path.lstat()
    except OSError as error:
        raise ReplayError(f"cannot stat {context} {path}: {error}") from error
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        raise ReplayError(f"{context} must be a regular non-symlink file: {path}")
    if metadata.st_nlink != 1:
        raise ReplayError(f"{context} must not be hard-linked: {path}")
    try:
        return path.read_bytes()
    except OSError as error:
        raise ReplayError(f"cannot read {context} {path}: {error}") from error


def sha256_file(path: Path, context: str = "file") -> str:
    return hashlib.sha256(regular_bytes(path, context)).hexdigest()


def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ReplayError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def read_json(path: Path, context: str) -> tuple[dict[str, Any], bytes]:
    data = regular_bytes(path, context)
    try:
        value = json.loads(data, object_pairs_hook=reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ReplayError(f"invalid JSON in {context}: {error}") from error
    if not isinstance(value, dict):
        raise ReplayError(f"{context} must contain a JSON object")
    return value, data


def pretty_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def write_new(path: Path, data: bytes, mode: int = 0o644) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags, mode)
    except OSError as error:
        raise ReplayError(f"refusing to overwrite {path}: {error}") from error
    with os.fdopen(descriptor, "wb", closefd=True) as output:
        output.write(data)
        output.flush()
        os.fsync(output.fileno())


def write_json_new(path: Path, value: Any) -> None:
    write_new(path, pretty_bytes(value))


def file_record(path: Path, root: Path, context: str) -> dict[str, Any]:
    data = regular_bytes(path, context)
    try:
        relative = path.resolve(strict=True).relative_to(root.resolve(strict=True)).as_posix()
    except (OSError, ValueError) as error:
        raise ReplayError(f"{context} escapes its root: {path}") from error
    return {"path": relative, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def inventory(root: Path, excluded: set[str] | None = None) -> list[dict[str, Any]]:
    excluded = excluded or set()
    records: list[dict[str, Any]] = []
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root).as_posix()
        metadata = path.lstat()
        if stat.S_ISLNK(metadata.st_mode):
            raise ReplayError(f"artifact tree contains symlink {relative}")
        if stat.S_ISDIR(metadata.st_mode):
            continue
        if relative in excluded:
            continue
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
            raise ReplayError(f"artifact tree contains inadmissible file {relative}")
        data = path.read_bytes()
        records.append({"path": relative, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
    records.sort(key=lambda record: record["path"])
    return records


def with_self_hash(value: dict[str, Any], field: str) -> dict[str, Any]:
    result = dict(value)
    result[field] = canonical_sha256(result)
    return result


def validate_self_hash(value: dict[str, Any], field: str, context: str) -> None:
    payload = {key: child for key, child in value.items() if key != field}
    if value.get(field) != canonical_sha256(payload):
        raise ReplayError(f"{context} self-hash is invalid")


def require_hash(value: Any, context: str) -> str:
    if not isinstance(value, str) or HEX64.fullmatch(value) is None:
        raise ReplayError(f"{context} must be a lowercase SHA-256 digest")
    return value


def validate_protocol(protocol: dict[str, Any]) -> dict[str, Any]:
    if protocol.get("schema") != PROTOCOL_SCHEMA or protocol.get("status") != "frozen_before_external_execution":
        raise ReplayError("unsupported or unfrozen Magma replay protocol")
    source = protocol.get("source")
    packet = protocol.get("packet")
    execution = protocol.get("execution")
    preflight = protocol.get("preflight")
    return_policy = protocol.get("return")
    if not all(isinstance(value, dict) for value in (source, packet, execution, preflight, return_policy)):
        raise ReplayError("Magma replay protocol is incomplete")
    for field in (
        "terminal_bundle_seal_sha256", "blind_bundle_sha256",
        "phase_b_protocol_sha256", "terminal_run_seal_sha256",
        "terminal_run_inventory_sha256", "terminal_score_seal_sha256",
        "terminal_score_sha256",
    ):
        require_hash(source.get(field), f"source.{field}")
    if packet.get("instance_count") != 160 or packet.get("order") != "exact blind-bundle order":
        raise ReplayError("production packet must contain all 160 instances in blind order")
    if packet.get("copy_source_exports_byte_for_byte") is not True or packet.get("ground_truth_allowed") is not False or packet.get("class_oracle_or_known_witness_labels_allowed") is not False:
        raise ReplayError("production solver packet widened its information boundary")
    if (
        execution.get("parallel_workers") != 1
        or execution.get("magma_threads") != 1
        or execution.get("gpu_enabled") is not False
        or execution.get("watchdog_seconds_per_process") != 120
        or execution.get("attempts_per_phase") != 1
        or execution.get("retry_allowed") is not False
        or execution.get("unknown_or_timeout_is_unsat") is not False
        or execution.get("unknown_or_timeout_is_inconclusive") is not True
        or "Al := \"Direct\"" not in execution.get("f4_algorithm", "")
        or "Dense := false" not in execution.get("f4_algorithm", "")
        or "Nthreads := 1" not in execution.get("f4_algorithm", "")
    ):
        raise ReplayError("production execution policy changed")
    if not all(preflight.get(field) is True for field in (
        "bind_magma_executable_hash_and_version", "bind_minisat_executable_hash_and_version",
        "bind_backend_executable_hash_and_version", "bind_process_meter_hash",
        "bind_runner_hash_and_repository_state", "bind_host_identity",
        "require_license_access_statement",
    )):
        raise ReplayError("licensed-host preflight policy changed")
    require_hash(preflight.get("expected_backend_sha256"), "preflight.expected_backend_sha256")
    require_hash(preflight.get("expected_process_meter_sha256"), "preflight.expected_process_meter_sha256")
    if not isinstance(preflight.get("expected_backend_bytes"), int) or preflight["expected_backend_bytes"] <= 0:
        raise ReplayError("preflight expected backend byte count is invalid")
    if return_policy.get("write_once") is not True or return_policy.get("score_only_after_return_seal") is not True:
        raise ReplayError("external return custody policy changed")
    return protocol


def walk_keys(value: Any):
    if isinstance(value, dict):
        for key, child in value.items():
            yield key
            yield from walk_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk_keys(child)


def reject_ground_truth(value: Any, context: str) -> None:
    found = {key.lower() for key in walk_keys(value)} & FORBIDDEN_INPUT_FIELDS
    if found:
        raise ReplayError(f"{context} exposes forbidden solver fields: {sorted(found)}")
    serialized = json.dumps(value, sort_keys=True).lower()
    leaked = {marker for marker in FORBIDDEN_INPUT_FIELDS if marker in serialized}
    if leaked:
        raise ReplayError(f"{context} exposes forbidden solver values: {sorted(leaked)}")


def validate_terminal_bundle(protocol: dict[str, Any], bundle: Path) -> dict[str, Any]:
    bundle = bundle.resolve(strict=True)
    source = protocol["source"]
    if sha256_file(bundle / "bundle-seal.json", "terminal bundle seal") != source["terminal_bundle_seal_sha256"]:
        raise ReplayError("terminal successor-04 bundle seal differs from the replay protocol")
    try:
        verified = terminal.verify(bundle)
    except terminal.EvidenceError as error:
        raise ReplayError(f"terminal successor-04 bundle failed verification: {error}") from error
    paths = {
        "blind_bundle": bundle / "inputs/blind-bundle.json",
        "phase_b_protocol": bundle / "inputs/phase-b-protocol.json",
        "run_seal": bundle / "run/run-seal.json",
        "score_seal": bundle / "score/score-seal.json",
    }
    expected = {
        "blind_bundle": source["blind_bundle_sha256"],
        "phase_b_protocol": source["phase_b_protocol_sha256"],
        "run_seal": source["terminal_run_seal_sha256"],
        "score_seal": source["terminal_score_seal_sha256"],
    }
    for name, path in paths.items():
        if sha256_file(path, name) != expected[name]:
            raise ReplayError(f"terminal bundle {name} differs from the replay protocol")
    blind, _ = read_json(paths["blind_bundle"], "blind bundle")
    phase_b_protocol, _ = phase_b.read_json(paths["phase_b_protocol"], "archived Phase-B protocol")
    try:
        phase_b.validate_protocol(phase_b_protocol)
        phase_b.validate_blind_bundle(blind, phase_b_protocol)
    except phase_b.PhaseBError as error:
        raise ReplayError(f"archived blind source failed semantic validation: {error}") from error
    run_seal, _ = read_json(paths["run_seal"], "terminal run seal")
    score_seal, _ = read_json(paths["score_seal"], "terminal score seal")
    reject_ground_truth(blind, "blind bundle")
    if blind.get("instance_count") != 160 or len(blind.get("instances", [])) != 160:
        raise ReplayError("blind bundle does not contain exactly 160 instances")
    if run_seal.get("inventory_sha256") != source["terminal_run_inventory_sha256"]:
        raise ReplayError("terminal run inventory differs from the replay protocol")
    magma_rows = [row for row in run_seal.get("inventory", []) if row.get("path", "").endswith("/instance/instance.magma")]
    if len(magma_rows) != 160:
        raise ReplayError("terminal run seal does not bind exactly 160 Magma exports")
    if score_seal.get("score_sha256") != source["terminal_score_sha256"]:
        raise ReplayError("terminal score seal differs from the frozen score identity")
    return {
        "verified": verified,
        "bundle": bundle,
        "blind": blind,
        "run_seal": run_seal,
        "score_seal": score_seal,
        "phase_b_protocol": phase_b_protocol,
        "paths": paths,
    }


def dry_run(protocol_path: Path, bundle: Path) -> dict[str, Any]:
    protocol, protocol_bytes = read_json(protocol_path.resolve(), "Magma replay protocol")
    validate_protocol(protocol)
    source = validate_terminal_bundle(protocol, bundle)
    return {
        "schema": "koblitz_phase_b_magma_replay_dry_run.v1",
        "status": "ready_for_source_export_copy",
        "protocol_sha256": hashlib.sha256(protocol_bytes).hexdigest(),
        "bundle_seal_sha256": protocol["source"]["terminal_bundle_seal_sha256"],
        "blind_instance_count": len(source["blind"]["instances"]),
        "sealed_magma_export_count": 160,
        "run_inventory_sha256": source["run_seal"]["inventory_sha256"],
        "licensed_magma_executed": False,
        "claim_boundary": protocol["claim_boundary"],
    }


def validate_source_run(run_root: Path, source: dict[str, Any], protocol: dict[str, Any]) -> dict[str, dict[str, Any]]:
    run_root = run_root.resolve(strict=True)
    if run_root.is_symlink() or not run_root.is_dir():
        raise ReplayError("terminal Phase-B run root must be a real directory")
    run_seal_path = run_root / "run-seal.json"
    if sha256_file(run_seal_path, "terminal run seal") != protocol["source"]["terminal_run_seal_sha256"]:
        raise ReplayError("live terminal run seal differs from successor-04")
    live_seal, _ = read_json(run_seal_path, "live terminal run seal")
    if live_seal != source["run_seal"]:
        raise ReplayError("live terminal run seal bytes do not match successor-04")
    observed = inventory(run_root, {"run-seal.json"})
    if observed != live_seal.get("inventory") or canonical_sha256(observed) != live_seal.get("inventory_sha256"):
        raise ReplayError("live terminal run inventory does not match its seal")
    return {row["path"]: row for row in observed}


def validate_score_seal(score_seal_path: Path, source: dict[str, Any], protocol: dict[str, Any]) -> None:
    if sha256_file(score_seal_path.resolve(strict=True), "terminal score seal") != protocol["source"]["terminal_score_seal_sha256"]:
        raise ReplayError("terminal score seal differs from the replay protocol")
    value, raw = read_json(score_seal_path.resolve(strict=True), "terminal score seal")
    if raw != regular_bytes(source["paths"]["score_seal"], "bundled score seal"):
        raise ReplayError("terminal score seal differs from successor-04")
    if value.get("score_sha256") != protocol["source"]["terminal_score_sha256"]:
        raise ReplayError("terminal score identity differs from the replay protocol")


def source_task_records(
    source: dict[str, Any], run_root: Path, sealed_inventory: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for ordinal, blind in enumerate(source["blind"]["instances"]):
        blind_id = blind.get("blind_instance_id")
        if not isinstance(blind_id, str) or BLIND_ID.fullmatch(blind_id) is None:
            raise ReplayError(f"blind instance {ordinal} has an invalid identifier")
        task_name = f"{ordinal:06d}-{blind_id}"
        source_prefix = f"tasks/{task_name}/instance"
        magma_relative = f"{source_prefix}/instance.magma"
        manifest_relative = f"{source_prefix}/manifest.json"
        magma_path = run_root / magma_relative
        manifest_path = run_root / manifest_relative
        for relative, path, label in (
            (magma_relative, magma_path, "Magma export"),
            (manifest_relative, manifest_path, "export manifest"),
        ):
            record = sealed_inventory.get(relative)
            if record is None or file_record(path, run_root, label) != record:
                raise ReplayError(f"{task_name}: {label} differs from the terminal run seal")
        manifest, _ = read_json(manifest_path, f"{task_name} manifest")
        reject_ground_truth(manifest, f"{task_name} manifest")
        try:
            phase_b.validate_export_manifest(
                manifest, blind, manifest_path.parent, source["phase_b_protocol"]
            )
        except phase_b.PhaseBError as error:
            raise ReplayError(f"{task_name}: source export failed semantic validation: {error}") from error
        source_identity = manifest.get("source_instance", {}).get("identity", {})
        expected = {
            "blind_instance_id": blind_id,
            "seed": blind["source_nonce"],
            "n": blind["n"],
            "ell": blind["ell"],
            "m": blind["m"],
            "curve_a": blind["curve_a"],
            "target": blind["target"],
        }
        observed = {
            "blind_instance_id": manifest.get("blind_instance_id"),
            "seed": manifest.get("seed"),
            "n": manifest.get("n"),
            "ell": manifest.get("ell"),
            "m": manifest.get("m"),
            "curve_a": manifest.get("curve_a"),
            "target": manifest.get("target"),
        }
        if observed != expected or source_identity.get("blind_instance_id") != blind_id:
            raise ReplayError(f"{task_name}: export manifest differs from the blind instance")
        try:
            legacy.split_magma_source(magma_path)
        except legacy.ExternalMagmaError as error:
            raise ReplayError(f"{task_name}: unsafe Magma source: {error}") from error
        magma_text = regular_bytes(magma_path, f"{task_name} Magma source").decode()
        forbidden_text = {item for item in FORBIDDEN_INPUT_FIELDS if item in magma_text.lower()}
        if forbidden_text:
            raise ReplayError(f"{task_name}: Magma source exposes forbidden labels {sorted(forbidden_text)}")
        source_variables = manifest.get("source_variables")
        source_blake3 = manifest.get("source_instance", {}).get("id_blake3")
        if not isinstance(source_variables, int) or source_variables <= 0 or not isinstance(source_blake3, str) or HEX64.fullmatch(source_blake3) is None:
            raise ReplayError(f"{task_name}: source identity is incomplete")
        tasks.append({
            "ordinal": ordinal,
            "id": task_name,
            "blind_instance_id": blind_id,
            "source_system_id": blind["source_system_id"],
            "source_nonce": blind["source_nonce"],
            "source_instance_blake3": source_blake3,
            "source_variables": source_variables,
            "cell_id": blind["cell_id"],
            "n": blind["n"],
            "ell": blind["ell"],
            "m": blind["m"],
            "basis": blind["basis"],
            "curve_a": blind["curve_a"],
            "factor_index": blind["factor_index"],
            "target": blind["target"],
            "source_files": {
                "magma_boolean_f4": dict(sealed_inventory[magma_relative]),
                "manifest": dict(sealed_inventory[manifest_relative]),
            },
        })
    if len(tasks) != 160 or len({task["blind_instance_id"] for task in tasks}) != 160:
        raise ReplayError("source task selection is not the exact 160-instance blind order")
    return tasks


def build_packet(
    protocol_path: Path,
    bundle: Path,
    run_root: Path,
    score_seal_path: Path,
    output: Path,
) -> dict[str, Any]:
    protocol_path = protocol_path.resolve(strict=True)
    protocol, protocol_bytes = read_json(protocol_path, "Magma replay protocol")
    validate_protocol(protocol)
    source = validate_terminal_bundle(protocol, bundle)
    run_root = run_root.resolve(strict=True)
    sealed_inventory = validate_source_run(run_root, source, protocol)
    validate_score_seal(score_seal_path, source, protocol)
    tasks = source_task_records(source, run_root, sealed_inventory)
    output = output.resolve()
    for protected, label in ((run_root, "terminal run"), (source["bundle"], "terminal evidence bundle")):
        try:
            output.relative_to(protected)
        except ValueError:
            pass
        else:
            raise ReplayError(f"replay packet output must not reside inside the {label}")
    if output.exists() or output.is_symlink():
        raise ReplayError(f"replay packet output must be new: {output}")
    output.mkdir(parents=True)
    write_new(output / "protocol.json", protocol_bytes)
    packet_tasks: list[dict[str, Any]] = []
    for task in tasks:
        destination = output / "instances" / task["id"]
        destination.mkdir(parents=True)
        files = {}
        for name, filename in (("magma_boolean_f4", "instance.magma"), ("manifest", "manifest.json")):
            source_path = run_root / task["source_files"][name]["path"]
            destination_path = destination / filename
            write_new(destination_path, regular_bytes(source_path, f"{task['id']} {name}"))
            files[name] = file_record(destination_path, output, f"packet {task['id']} {name}")
            if files[name]["bytes"] != task["source_files"][name]["bytes"] or files[name]["sha256"] != task["source_files"][name]["sha256"]:
                raise ReplayError(f"{task['id']}: copied {name} differs from its source")
        packet_task = {key: value for key, value in task.items() if key != "source_files"}
        packet_task["files"] = files
        packet_tasks.append(packet_task)
    manifest = with_self_hash({
        "schema": PACKET_MANIFEST_SCHEMA,
        "status": "solver_inputs_frozen",
        "created_at": now(),
        "protocol_path": "protocol.json",
        "protocol_sha256": hashlib.sha256(protocol_bytes).hexdigest(),
        "source_terminal_bundle_seal_sha256": protocol["source"]["terminal_bundle_seal_sha256"],
        "source_terminal_run_seal_sha256": protocol["source"]["terminal_run_seal_sha256"],
        "source_terminal_run_inventory_sha256": protocol["source"]["terminal_run_inventory_sha256"],
        "source_terminal_score_seal_sha256": protocol["source"]["terminal_score_seal_sha256"],
        "order": "exact blind-bundle order",
        "instance_count": 160,
        "tasks": packet_tasks,
        "ground_truth_included": False,
        "execution_policy": protocol["execution"],
        "claim_boundary": protocol["claim_boundary"],
    }, "manifest_payload_sha256")
    write_json_new(output / "packet-manifest.json", manifest)
    sealed_files = inventory(output, {"packet-seal.json"})
    seal = with_self_hash({
        "schema": PACKET_SEAL_SCHEMA,
        "status": "external_solver_packet_frozen",
        "manifest_path": "packet-manifest.json",
        "manifest_sha256": sha256_file(output / "packet-manifest.json", "packet manifest"),
        "instance_count": 160,
        "magma_input_count": 160,
        "inventory": sealed_files,
        "inventory_sha256": canonical_sha256(sealed_files),
        "claim_boundary": protocol["claim_boundary"],
    }, "seal_payload_sha256")
    write_json_new(output / "packet-seal.json", seal)
    validate_packet(output, protocol)
    return seal


def resolve_packet_file(packet: Path, record: dict[str, Any], context: str) -> Path:
    relative = Path(str(record.get("path", "")))
    if relative.is_absolute() or ".." in relative.parts or not relative.parts:
        raise ReplayError(f"unsafe path for {context}")
    path = packet / relative
    if file_record(path, packet, context) != record:
        raise ReplayError(f"{context} differs from its packet identity")
    return path


def absolute_path_has_suffix(value: Any, relative: str) -> bool:
    if not isinstance(value, str):
        return False
    path = Path(value)
    suffix = Path(relative)
    return path.is_absolute() and len(path.parts) >= len(suffix.parts) and path.parts[-len(suffix.parts):] == suffix.parts


def validate_packet(packet: Path, expected_protocol: dict[str, Any] | None = None) -> dict[str, Any]:
    packet = packet.resolve(strict=True)
    if packet.is_symlink() or not packet.is_dir():
        raise ReplayError("replay packet must be a real directory")
    seal, seal_bytes = read_json(packet / "packet-seal.json", "packet seal")
    if seal.get("schema") != PACKET_SEAL_SCHEMA or seal.get("status") != "external_solver_packet_frozen":
        raise ReplayError("replay packet seal is not terminal")
    validate_self_hash(seal, "seal_payload_sha256", "packet seal")
    actual_inventory = inventory(packet, {"packet-seal.json"})
    if actual_inventory != seal.get("inventory") or canonical_sha256(actual_inventory) != seal.get("inventory_sha256"):
        raise ReplayError("replay packet inventory changed after sealing")
    manifest_path = packet / "packet-manifest.json"
    manifest, manifest_bytes = read_json(manifest_path, "packet manifest")
    if sha256_file(manifest_path, "packet manifest") != seal.get("manifest_sha256"):
        raise ReplayError("packet manifest differs from its seal")
    if manifest.get("schema") != PACKET_MANIFEST_SCHEMA or manifest.get("status") != "solver_inputs_frozen":
        raise ReplayError("packet manifest has an unexpected schema or status")
    validate_self_hash(manifest, "manifest_payload_sha256", "packet manifest")
    protocol, protocol_bytes = read_json(packet / "protocol.json", "packet protocol")
    validate_protocol(protocol)
    if expected_protocol is not None and protocol != expected_protocol:
        raise ReplayError("packet protocol differs from the requested protocol")
    if hashlib.sha256(protocol_bytes).hexdigest() != manifest.get("protocol_sha256"):
        raise ReplayError("packet protocol differs from its manifest")
    tasks = manifest.get("tasks")
    if manifest.get("instance_count") != 160 or seal.get("instance_count") != 160 or not isinstance(tasks, list) or len(tasks) != 160:
        raise ReplayError("replay packet does not contain exactly 160 tasks")
    if manifest.get("ground_truth_included") is not False:
        raise ReplayError("replay packet exposes ground truth")
    for ordinal, task in enumerate(tasks):
        if not isinstance(task, dict) or task.get("ordinal") != ordinal or task.get("id") != f"{ordinal:06d}-{task.get('blind_instance_id')}":
            raise ReplayError("replay packet task order changed")
        reject_ground_truth(task, f"packet task {ordinal}")
        magma = resolve_packet_file(packet, task.get("files", {}).get("magma_boolean_f4", {}), f"task {ordinal} Magma input")
        manifest_file = resolve_packet_file(packet, task.get("files", {}).get("manifest", {}), f"task {ordinal} manifest")
        decoded, _ = read_json(manifest_file, f"task {ordinal} manifest")
        reject_ground_truth(decoded, f"task {ordinal} manifest")
        if decoded.get("blind_instance_id") != task.get("blind_instance_id") or decoded.get("source_variables") != task.get("source_variables"):
            raise ReplayError(f"task {ordinal} manifest identity changed")
        try:
            legacy.split_magma_source(magma)
        except legacy.ExternalMagmaError as error:
            raise ReplayError(f"task {ordinal} Magma template changed: {error}") from error
    return {
        "packet_seal_sha256": hashlib.sha256(seal_bytes).hexdigest(),
        "packet_inventory_sha256": seal["inventory_sha256"],
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "protocol": protocol,
        "manifest": manifest,
        "packet": packet,
    }


def absolute_identity(path: Path, context: str) -> dict[str, Any]:
    data = regular_bytes(path.resolve(strict=True), context)
    return {"path": str(path.resolve(strict=True)), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def repository_state() -> dict[str, Any]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPO, text=True,
        capture_output=True, check=True,
    ).stdout.strip()
    status = subprocess.run(
        [
            "git", "status", "--porcelain=v1", "--",
            "scripts/run_koblitz_phase_b_magma_replay.py",
            "scripts/koblitz_phase_b_magma_core.py",
            "scripts/run_koblitz_pdp_matrix.py",
            "scripts/process_meter.py",
            "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-22-phase-b-magma-replay-protocol.json",
        ],
        cwd=REPO, text=True, capture_output=True, check=True,
    ).stdout.splitlines()
    if re.fullmatch(r"[0-9a-f]{40}", commit) is None:
        raise ReplayError("cannot bind the replay source revision")
    return {"commit": commit, "dirty": bool(status), "porcelain": status}


def tool_identity(value: str | Path, options: list[list[str]], context: str) -> dict[str, Any]:
    try:
        identity = legacy.binary_identity(value, options)
    except legacy.ExternalMagmaError as error:
        raise ReplayError(f"{context} preflight failed: {error}") from error
    if identity.get("version_returncode") != 0 or identity.get("version") == "unreported":
        raise ReplayError(f"{context} did not return a successful version probe")
    return identity


def preflight(
    packet: Path,
    magma: str | Path,
    minisat: str | Path,
    backend: str | Path,
    meter: Path,
    license_access_statement: str,
    *,
    synthetic_test_mode: bool = False,
) -> dict[str, Any]:
    validated = validate_packet(packet)
    statement = license_access_statement.strip()
    if len(statement) < 12 or any(marker in statement.lower() for marker in ("todo", "placeholder", "unavailable")):
        raise ReplayError("licensed-host preflight requires a concrete license-access statement")
    magma_identity = tool_identity(magma, [["--version"]], "Magma")
    minisat_identity = tool_identity(minisat, [["--version"], ["-h"]], "minisat")
    backend_identity = tool_identity(backend, [["--version"]], "point validator")
    if Path(minisat_identity["path"]).name != "minisat":
        raise ReplayError("the SAT(G) helper must resolve to an executable named minisat")
    recognized = legacy.recognized_magma_version(magma_identity["version"])
    magma_identity["recognized_version_banner"] = recognized
    if not synthetic_test_mode and not recognized:
        raise ReplayError("Magma version output is not a recognized licensed-host V2 banner")
    preflight_policy = validated["protocol"]["preflight"]
    if not synthetic_test_mode and (
        backend_identity["sha256"] != preflight_policy["expected_backend_sha256"]
        or backend_identity["bytes"] != preflight_policy["expected_backend_bytes"]
    ):
        raise ReplayError("point validator differs from the frozen Phase-B backend executable")
    meter = meter.resolve(strict=True)
    if sha256_file(meter, "process meter") != preflight_policy["expected_process_meter_sha256"]:
        raise ReplayError("process meter differs from the frozen implementation")
    source_state = repository_state()
    if not synthetic_test_mode and source_state["dirty"]:
        raise ReplayError("licensed production replay requires a clean runner source state")
    result = with_self_hash({
        "schema": PREFLIGHT_SCHEMA,
        "status": "synthetic_test_ready" if synthetic_test_mode else "licensed_host_ready",
        "created_at": now(),
        "packet_seal_sha256": validated["packet_seal_sha256"],
        "packet_inventory_sha256": validated["packet_inventory_sha256"],
        "tools": {
            "magma": magma_identity,
            "minisat": minisat_identity,
            "backend": backend_identity,
            "meter": absolute_identity(meter, "process meter"),
        },
        "host": legacy.host_identity(),
        "source": {
            "repository": str(REPO),
            "state": source_state,
            "runner": absolute_identity(Path(__file__), "Magma replay runner"),
            "magma_core": absolute_identity(REPO / "scripts/koblitz_phase_b_magma_core.py", "Magma execution core"),
            "matrix_contract": absolute_identity(REPO / "scripts/run_koblitz_pdp_matrix.py", "matrix contract"),
        },
        "license_access_statement": statement,
        "synthetic_test_mode": synthetic_test_mode,
        "execution_policy": validated["protocol"]["execution"],
    }, "preflight_payload_sha256")
    return result


def validate_preflight_record(
    value: dict[str, Any], packet_info: dict[str, Any]
) -> None:
    validate_self_hash(value, "preflight_payload_sha256", "licensed-host preflight")
    synthetic = value.get("synthetic_test_mode")
    if type(synthetic) is not bool:
        raise ReplayError("licensed-host preflight has an invalid synthetic-mode flag")
    expected_status = "synthetic_test_ready" if synthetic else "licensed_host_ready"
    if value.get("schema") != PREFLIGHT_SCHEMA or value.get("status") != expected_status:
        raise ReplayError("licensed-host preflight has an invalid schema or status")
    if (
        value.get("packet_seal_sha256") != packet_info["packet_seal_sha256"]
        or value.get("packet_inventory_sha256") != packet_info["packet_inventory_sha256"]
        or value.get("execution_policy") != packet_info["protocol"]["execution"]
    ):
        raise ReplayError("licensed-host preflight is bound to a different packet or policy")
    statement = value.get("license_access_statement")
    if not isinstance(statement, str) or len(statement.strip()) < 12 or any(
        marker in statement.lower() for marker in ("todo", "placeholder", "unavailable")
    ):
        raise ReplayError("licensed-host preflight lacks a concrete license-access statement")
    tools = value.get("tools")
    if not isinstance(tools, dict) or set(tools) != {"magma", "minisat", "backend", "meter"}:
        raise ReplayError("licensed-host preflight tool inventory changed")
    for name in ("magma", "minisat", "backend"):
        identity = tools[name]
        if (
            not isinstance(identity, dict)
            or not isinstance(identity.get("path"), str)
            or not Path(identity["path"]).is_absolute()
            or not isinstance(identity.get("bytes"), int)
            or identity["bytes"] <= 0
            or not isinstance(identity.get("version"), str)
            or not identity["version"].strip()
            or identity.get("version_returncode") != 0
        ):
            raise ReplayError(f"licensed-host preflight {name} identity is incomplete")
        require_hash(identity.get("sha256"), f"licensed-host {name} identity")
    if Path(tools["minisat"]["path"]).name != "minisat":
        raise ReplayError("licensed-host preflight did not bind an executable named minisat")
    recognized = tools["magma"].get("recognized_version_banner")
    recomputed_recognition = legacy.recognized_magma_version(tools["magma"]["version"])
    if (
        type(recognized) is not bool
        or recognized is not recomputed_recognition
        or (not synthetic and recognized is not True)
    ):
        raise ReplayError("licensed-host preflight lacks a recognized Magma version banner")
    meter = tools["meter"]
    if (
        not isinstance(meter, dict)
        or require_hash(meter.get("sha256"), "licensed-host meter identity")
        != packet_info["protocol"]["preflight"]["expected_process_meter_sha256"]
        or not isinstance(meter.get("bytes"), int)
        or meter["bytes"] <= 0
    ):
        raise ReplayError("licensed-host preflight process-meter identity changed")
    if not synthetic and (
        tools["backend"]["sha256"]
        != packet_info["protocol"]["preflight"]["expected_backend_sha256"]
        or tools["backend"]["bytes"]
        != packet_info["protocol"]["preflight"]["expected_backend_bytes"]
    ):
        raise ReplayError("licensed-host point validator differs from the frozen backend")
    source = value.get("source")
    if not isinstance(source, dict) or set(source) != {"repository", "state", "runner", "magma_core", "matrix_contract"}:
        raise ReplayError("licensed-host preflight source inventory changed")
    state = source["state"]
    if (
        not isinstance(state, dict)
        or set(state) != {"commit", "dirty", "porcelain"}
        or re.fullmatch(r"[0-9a-f]{40}", str(state.get("commit"))) is None
        or type(state.get("dirty")) is not bool
        or not isinstance(state.get("porcelain"), list)
        or state["dirty"] != bool(state["porcelain"])
        or (not synthetic and state["dirty"])
    ):
        raise ReplayError("licensed-host preflight repository state is invalid")
    trusted_sources = {
        "runner": Path(__file__),
        "magma_core": REPO / "scripts/koblitz_phase_b_magma_core.py",
        "matrix_contract": REPO / "scripts/run_koblitz_pdp_matrix.py",
    }
    for name, path in trusted_sources.items():
        recorded = source[name]
        trusted = absolute_identity(path, f"trusted {name}")
        if not isinstance(recorded, dict) or any(recorded.get(field) != trusted[field] for field in ("bytes", "sha256")):
            raise ReplayError(f"licensed-host preflight {name} source differs from the trusted implementation")
    host = value.get("host")
    required_host = {"node", "platform", "system", "release", "machine", "cpu_model", "logical_cpus", "python", "python_executable"}
    if not isinstance(host, dict) or set(host) != required_host or not all(
        isinstance(host[field], str) and host[field] for field in required_host - {"logical_cpus"}
    ) or not isinstance(host["logical_cpus"], int) or host["logical_cpus"] <= 0:
        raise ReplayError("licensed-host preflight host identity is incomplete")


def receipt_reference(output: Path, receipt_path: Path) -> dict[str, Any]:
    record = file_record(receipt_path, output, "process receipt")
    receipt, _ = read_json(receipt_path, "process receipt")
    record["process"] = receipt["process"]
    return record


def write_attempt_start(
    task_root: Path,
    role: str,
    task: dict[str, Any],
    command: list[str],
    tool: dict[str, Any],
    watchdog: int,
) -> None:
    record = with_self_hash({
        "schema": "koblitz_phase_b_magma_replay_attempt_start.v1",
        "role": role,
        "task_ordinal": task["ordinal"],
        "blind_instance_id": task["blind_instance_id"],
        "source_magma_sha256": task["files"]["magma_boolean_f4"]["sha256"],
        "command": command,
        "tool_sha256": tool["sha256"],
        "watchdog_seconds": watchdog,
        "attempt_number": 1,
        "retry_allowed": False,
        "created_at": now(),
    }, "start_payload_sha256")
    write_json_new(task_root / f"{role}-start.json", record)


def task_input_records(packet: Path, task: dict[str, Any]) -> dict[str, Any]:
    records = {}
    for name, record in task["files"].items():
        path = resolve_packet_file(packet, record, f"{task['id']} {name}")
        records[name] = file_record(path, packet, f"{task['id']} {name}")
    return records


def identity_marker_source(task: dict[str, Any]) -> str:
    blind = task["blind_instance_id"].removeprefix("b-")
    source_hash = task["files"]["magma_boolean_f4"]["sha256"]
    return (
        f'printf "KOBLITZ_REPLAY_TASK_ORDINAL={task["ordinal"]}\\n";\n'
        f'printf "KOBLITZ_REPLAY_BLIND_ID_HI={blind[:32]}\\n";\n'
        f'printf "KOBLITZ_REPLAY_BLIND_ID_LO={blind[32:]}\\n";\n'
        f'printf "KOBLITZ_REPLAY_SOURCE_SHA256_HI={source_hash[:32]}\\n";\n'
        f'printf "KOBLITZ_REPLAY_SOURCE_SHA256_LO={source_hash[32:]}\\n";\n'
    )


def add_identity_markers(script: str, task: dict[str, Any]) -> str:
    if not script.endswith("quit;\n"):
        raise ReplayError(f"{task['id']}: generated Magma script lacks terminal quit")
    return script.removesuffix("quit;\n") + identity_marker_source(task) + "quit;\n"


def validate_identity_markers(output: str, task: dict[str, Any]) -> None:
    observed: dict[str, str] = {}
    for line in output.splitlines():
        if not line.startswith("KOBLITZ_REPLAY_") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in observed:
            raise ReplayError(f"{task['id']}: duplicate Magma identity marker {key}")
        observed[key] = value.strip()
    blind = task["blind_instance_id"].removeprefix("b-")
    source_hash = task["files"]["magma_boolean_f4"]["sha256"]
    expected = {
        "KOBLITZ_REPLAY_TASK_ORDINAL": str(task["ordinal"]),
        "KOBLITZ_REPLAY_BLIND_ID_HI": blind[:32],
        "KOBLITZ_REPLAY_BLIND_ID_LO": blind[32:],
        "KOBLITZ_REPLAY_SOURCE_SHA256_HI": source_hash[:32],
        "KOBLITZ_REPLAY_SOURCE_SHA256_LO": source_hash[32:],
    }
    if observed != expected:
        raise ReplayError(f"{task['id']}: Magma output identity markers differ")


def execute_one(
    packet_info: dict[str, Any],
    task: dict[str, Any],
    output: Path,
    preflight_record: dict[str, Any],
) -> dict[str, Any]:
    packet = packet_info["packet"]
    task_root = output / "tasks" / task["id"]
    task_root.mkdir(parents=True)
    before = task_input_records(packet, task)
    magma_path = resolve_packet_file(packet, task["files"]["magma_boolean_f4"], "Magma input")
    manifest_path = resolve_packet_file(packet, task["files"]["manifest"], "export manifest")
    manifest, _ = read_json(manifest_path, "export manifest")
    tools = preflight_record["tools"]
    timeout = packet_info["protocol"]["execution"]["watchdog_seconds_per_process"]
    meter = Path(tools["meter"]["path"])
    environment = os.environ.copy()
    environment.update({
        "PATH": str(Path(tools["minisat"]["path"]).parent) + os.pathsep + environment.get("PATH", ""),
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
    })
    primary_script = task_root / "primary-f4.magma"
    raw_magma = regular_bytes(magma_path, "sealed Magma input").decode()
    write_new(primary_script, add_identity_markers(raw_magma, task).encode())
    f4_command = [tools["magma"]["path"], "-t", "1", "-b", str(primary_script)]
    f4_attempt = task_root / "primary-f4"
    write_attempt_start(task_root, "primary-f4", task, f4_command, tools["magma"], timeout)
    try:
        f4_receipt = legacy.run_metered(meter, f4_command, task_root, timeout, f4_attempt, environment)
        f4_status, f4_terminal = legacy.classify_f4(f4_receipt)
        if not f4_receipt["process"]["timed_out"]:
            validate_identity_markers(legacy.receipt_stdout(f4_receipt), task)
    except legacy.ExternalMagmaError as error:
        raise ReplayError(f"{task['id']}: primary F4 execution failed: {error}") from error
    result: dict[str, Any] = {
        "schema": TASK_SCHEMA,
        "ordinal": task["ordinal"],
        "id": task["id"],
        "blind_instance_id": task["blind_instance_id"],
        "cell_id": task["cell_id"],
        "inputs_before": before,
        "primary_f4": {
            "status": f4_status,
            "terminal": f4_terminal,
            "script": file_record(primary_script, output, "identity-bound primary F4 script"),
            "source_magma_sha256": task["files"]["magma_boolean_f4"]["sha256"],
            "receipt": receipt_reference(output, f4_attempt / "receipt.json"),
        },
        "f4_plus_sat_g": None,
        "point_validation": None,
        "contradictions": [],
        "conflicts": None,
        "conflicts_semantics": "Magma does not expose SAT-style conflict counts for this path",
        "attempts": {"primary_f4": 1, "f4_plus_sat_g": 0, "point_validation": 0},
    }
    if f4_status == "unsat":
        result["final_status"] = "unsat"
    elif f4_status != "sat_basis_certificate_unverified_model":
        result["final_status"] = f4_status if f4_status.endswith("inconclusive") else f"{f4_status}_inconclusive"
    else:
        witness_path = task_root / "f4-plus-sat-g.magma"
        try:
            witness_text = add_identity_markers(
                legacy.render_witness_script(magma_path, [], int(task["source_variables"])), task
            )
        except legacy.ExternalMagmaError as error:
            raise ReplayError(f"{task['id']}: cannot render SAT(G) script: {error}") from error
        write_new(witness_path, witness_text.encode())
        model_command = [tools["magma"]["path"], "-t", "1", "-b", str(witness_path)]
        model_attempt = task_root / "f4-plus-sat-g"
        write_attempt_start(task_root, "f4-plus-sat-g", task, model_command, tools["magma"], timeout)
        model_receipt = legacy.run_metered(meter, model_command, task_root, timeout, model_attempt, environment)
        repeated_status, repeated_terminal = legacy.classify_f4(model_receipt)
        if not model_receipt["process"]["timed_out"]:
            validate_identity_markers(legacy.receipt_stdout(model_receipt), task)
        model_terminal = legacy.parse_model_terminal(
            legacy.receipt_stdout(model_receipt), int(task["source_variables"]), 0
        )
        result["attempts"]["f4_plus_sat_g"] = 1
        result["f4_plus_sat_g"] = {
            "f4_status": repeated_status,
            "f4_terminal": repeated_terminal,
            "model_terminal": model_terminal,
            "script": file_record(witness_path, output, "F4 plus SAT(G) script"),
            "receipt": receipt_reference(output, model_attempt / "receipt.json"),
        }
        if model_receipt["process"]["timed_out"]:
            result["final_status"] = "f4_plus_sat_g_timeout_inconclusive"
        elif f4_terminal is None or repeated_terminal is None or not legacy.same_f4_terminal(f4_terminal, repeated_terminal):
            result["contradictions"].append("F4 plus SAT(G) did not reproduce the primary F4 terminal")
            result["final_status"] = "f4_repetition_contradiction"
        elif repeated_status != "sat_basis_certificate_unverified_model":
            result["final_status"] = "f4_plus_sat_g_solver_error_inconclusive"
        elif model_terminal is None:
            result["final_status"] = "sat_g_terminal_missing_inconclusive"
        elif model_terminal["status"] == "unsat":
            result["contradictions"].append("proper F4 basis but SAT(G) returned UNSAT")
            result["final_status"] = "sat_g_contradiction"
        else:
            assignment = model_terminal["assignment"]
            assignment_path = task_root / "assignment.json"
            write_new(assignment_path, (json.dumps(assignment, separators=(",", ":")) + "\n").encode())
            validation_command = [
                tools["backend"]["path"], "validate-model", str(manifest_path), str(assignment_path)
            ]
            validation_attempt = task_root / "point-validation"
            write_attempt_start(
                task_root, "point-validation", task, validation_command,
                tools["backend"], timeout,
            )
            validation_receipt = legacy.run_metered(
                meter, validation_command, task_root, timeout, validation_attempt, environment
            )
            validation_status, report = legacy.validation_status(validation_receipt, manifest, assignment_path)
            result["attempts"]["point_validation"] = 1
            result["point_validation"] = {
                "status": validation_status,
                "report": report,
                "assignment": file_record(assignment_path, output, "model assignment"),
                "receipt": receipt_reference(output, validation_attempt / "receipt.json"),
            }
            if validation_status == "valid_point_witness":
                result["final_status"] = "sat"
            elif validation_status == "nonlifting_source_model":
                result["final_status"] = "nonlifting_model_inconclusive"
            elif validation_status == "timeout_inconclusive":
                result["final_status"] = "point_validation_timeout_inconclusive"
            else:
                result["final_status"] = f"point_validation_{validation_status}_inconclusive"
    after = task_input_records(packet, task)
    result["inputs_after"] = after
    result["inputs_unchanged"] = before == after
    if not result["inputs_unchanged"]:
        result["contradictions"].append("sealed solver input changed during execution")
        result["final_status"] = "input_custody_failure"
    result["finished_at"] = now()
    write_json_new(task_root / "task-result.json", result)
    return result


def resource_totals(rows: list[dict[str, Any]]) -> dict[str, Any]:
    metrics = [row["process"]["metrics"] for row in rows]
    return {
        "processes": len(metrics),
        "total_core_seconds": sum(float(item["total_core_seconds"]) for item in metrics),
        "single_core_seconds": sum(float(item["single_core_seconds"]) for item in metrics),
        "single_core_seconds_alias_of": "total_core_seconds",
        "single_core_elapsed_seconds": None,
        "summed_process_wall_seconds": sum(float(item["wall_seconds"]) for item in metrics),
        "largest_process_peak_rss_bytes": max((int(item["peak_rss_bytes"]) for item in metrics), default=0),
        "aggregate_parallel_peak_rss_bytes": None,
    }


def summarize(results: list[dict[str, Any]], preflight_record: dict[str, Any], expected: int) -> dict[str, Any]:
    status_counts = Counter(str(result.get("final_status", "missing")) for result in results)
    phase_rows: dict[str, list[dict[str, Any]]] = {
        "primary_f4": [], "f4_plus_sat_g": [], "point_validation": [],
    }
    for result in results:
        for field in phase_rows:
            value = result.get(field)
            if isinstance(value, dict) and isinstance(value.get("receipt"), dict):
                phase_rows[field].append(value["receipt"])
    all_rows = [row for values in phase_rows.values() for row in values]
    per_cell: dict[str, Any] = {}
    for cell_id in sorted({str(result.get("cell_id")) for result in results}):
        cell_results = [result for result in results if result.get("cell_id") == cell_id]
        cell_phases = {
            name: [
                result[name]["receipt"]
                for result in cell_results
                if isinstance(result.get(name), dict)
                and isinstance(result[name].get("receipt"), dict)
            ]
            for name in phase_rows
        }
        per_cell[cell_id] = {
            "tasks": len(cell_results),
            "statuses": dict(sorted(Counter(str(result.get("final_status")) for result in cell_results).items())),
            "resources": {
                **{name: resource_totals(rows) for name, rows in cell_phases.items()},
                "all_metered_processes": resource_totals(
                    [row for rows in cell_phases.values() for row in rows]
                ),
            },
            "conflicts": None,
        }
    synthetic = preflight_record.get("synthetic_test_mode") is True
    full = len(results) == expected == 160
    return {
        "schema": SUMMARY_SCHEMA,
        "status": "all_inputs_terminal" if full else "partial_test_selection_terminal",
        "expected_tasks": 160,
        "selected_tasks": expected,
        "terminal_task_records": len(results),
        "status_counts": dict(sorted(status_counts.items())),
        "one_attempt_no_retry": all(
            all(isinstance(count, int) and 0 <= count <= 1 for count in result.get("attempts", {}).values())
            for result in results
        ),
        "resources": {
            **{name: resource_totals(rows) for name, rows in phase_rows.items()},
            "all_metered_processes": resource_totals(all_rows),
        },
        "per_cell": per_cell,
        "conflicts": None,
        "conflicts_semantics": "unavailable from the Magma F4 and SAT(G) interfaces",
        "unknown_and_timeout_are_inconclusive": True,
        "licensed_magma_execution_completed": full and not synthetic,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def run(
    packet: Path,
    output: Path,
    magma: str | Path,
    minisat: str | Path,
    backend: str | Path,
    meter: Path,
    license_access_statement: str,
    *,
    synthetic_test_mode: bool = False,
    max_tasks: int | None = None,
) -> dict[str, Any]:
    packet_info = validate_packet(packet)
    preflight_record = preflight(
        packet, magma, minisat, backend, meter, license_access_statement,
        synthetic_test_mode=synthetic_test_mode,
    )
    tasks = packet_info["manifest"]["tasks"]
    if max_tasks is not None:
        if not synthetic_test_mode or not isinstance(max_tasks, int) or not 1 <= max_tasks <= len(tasks):
            raise ReplayError("a partial task selection is permitted only in synthetic test mode")
        tasks = tasks[:max_tasks]
    output = output.resolve()
    try:
        output.relative_to(packet_info["packet"])
    except ValueError:
        pass
    else:
        raise ReplayError("external return must not reside inside the sealed solver packet")
    if output.exists() or output.is_symlink():
        raise ReplayError(f"external return path must be new; retries and resume are forbidden: {output}")
    output.mkdir(parents=True)
    write_json_new(output / "preflight.json", preflight_record)
    binding = with_self_hash({
        "schema": "koblitz_phase_b_magma_replay_packet_binding.v1",
        "packet_seal_sha256": packet_info["packet_seal_sha256"],
        "packet_inventory_sha256": packet_info["packet_inventory_sha256"],
        "packet_manifest_sha256": packet_info["manifest_sha256"],
        "selected_task_ids": [task["id"] for task in tasks],
    }, "binding_payload_sha256")
    write_json_new(output / "packet-binding.json", binding)
    results: list[dict[str, Any]] = []
    for task in tasks:
        try:
            result = execute_one(packet_info, task, output, preflight_record)
        except (ReplayError, legacy.ExternalMagmaError, OSError) as error:
            task_root = output / "tasks" / task["id"]
            task_root.mkdir(parents=True, exist_ok=True)
            attempts = {
                "primary_f4": int((task_root / "primary-f4-start.json").exists()),
                "f4_plus_sat_g": int((task_root / "f4-plus-sat-g-start.json").exists()),
                "point_validation": int((task_root / "point-validation-start.json").exists()),
            }
            retained_phases: dict[str, Any] = {}
            for field, directory in (
                ("primary_f4", "primary-f4"),
                ("f4_plus_sat_g", "f4-plus-sat-g"),
                ("point_validation", "point-validation"),
            ):
                receipt_path = task_root / directory / "receipt.json"
                retained_phases[field] = (
                    {"status": "retained_after_execution_error", "receipt": receipt_reference(output, receipt_path)}
                    if receipt_path.is_file()
                    else None
                )
            rechecked_inputs = task_input_records(packet_info["packet"], task)
            result = {
                "schema": TASK_SCHEMA,
                "ordinal": task["ordinal"],
                "id": task["id"],
                "blind_instance_id": task["blind_instance_id"],
                "cell_id": task["cell_id"],
                "final_status": "execution_error_inconclusive",
                "attempts": attempts,
                **retained_phases,
                "contradictions": [],
                "conflicts": None,
                "conflicts_semantics": "Magma does not expose SAT-style conflict counts for this path",
                "error": str(error),
                "inputs_rechecked_after_error": rechecked_inputs,
                "inputs_unchanged": True,
                "finished_at": now(),
            }
            if not (task_root / "task-result.json").exists():
                write_json_new(task_root / "task-result.json", result)
        results.append(result)
    summary = summarize(results, preflight_record, len(tasks))
    summary["finished_at"] = now()
    write_json_new(output / "summary.json", summary)
    sealed_files = inventory(output, {"return-seal.json"})
    return_claim = dict(packet_info["protocol"]["claim_boundary"])
    return_claim["licensed_magma_execution_completed"] = summary[
        "licensed_magma_execution_completed"
    ]
    seal = with_self_hash({
        "schema": RETURN_SEAL_SCHEMA,
        "status": "external_results_frozen",
        "created_at": now(),
        "packet_seal_sha256": packet_info["packet_seal_sha256"],
        "packet_inventory_sha256": packet_info["packet_inventory_sha256"],
        "preflight_sha256": sha256_file(output / "preflight.json", "preflight"),
        "summary_sha256": sha256_file(output / "summary.json", "summary"),
        "selected_tasks": len(tasks),
        "terminal_task_records": len(results),
        "one_attempt_no_retry": summary["one_attempt_no_retry"],
        "inventory": sealed_files,
        "inventory_sha256": canonical_sha256(sealed_files),
        "claim_boundary": return_claim,
    }, "seal_payload_sha256")
    write_json_new(output / "return-seal.json", seal)
    validate_return(packet, output)
    return seal


def local_receipt(return_root: Path, reference: dict[str, Any], context: str) -> dict[str, Any]:
    relative = Path(str(reference.get("path", "")))
    if relative.is_absolute() or ".." in relative.parts or not relative.parts:
        raise ReplayError(f"unsafe process receipt path for {context}")
    path = return_root / relative
    observed = file_record(path, return_root, context)
    if any(observed[field] != reference.get(field) for field in ("path", "bytes", "sha256")):
        raise ReplayError(f"{context} receipt identity changed")
    receipt, _ = read_json(path, context)
    process = receipt.get("process")
    if process != reference.get("process"):
        raise ReplayError(f"{context} embedded process record changed")
    try:
        legacy.validate_meter_record(process, process.get("command"), 120)
    except legacy.ExternalMagmaError as error:
        raise ReplayError(f"{context} process receipt is invalid: {error}") from error
    launcher = receipt.get("meter_launcher_command")
    command = process["command"]
    if (
        receipt.get("meter_launcher_returncode") != 0
        or not isinstance(launcher, list)
        or len(launcher) <= len(command)
        or launcher[-len(command):] != command
    ):
        raise ReplayError(f"{context} meter launcher differs from its child command")
    attempt = path.parent
    metrics, _ = read_json(attempt / "metrics.json", f"{context} raw metrics")
    if metrics != process:
        raise ReplayError(f"{context} raw metrics differ from the receipt")
    for stream in ("stdout", "stderr"):
        raw = regular_bytes(attempt / stream, f"{context} {stream}")
        record = receipt.get(stream)
        if not isinstance(record, dict) or len(raw) != record.get("bytes") or hashlib.sha256(raw).hexdigest() != record.get("sha256"):
            raise ReplayError(f"{context} raw {stream} differs from the receipt")
    return receipt


def receipt_stream(return_root: Path, reference: dict[str, Any], stream: str, context: str) -> str:
    relative = Path(reference["path"])
    data = regular_bytes(return_root / relative.parent / stream, f"{context} {stream}")
    try:
        return data.decode()
    except UnicodeDecodeError as error:
        raise ReplayError(f"{context} {stream} is not UTF-8") from error


def validate_attempt_start(
    return_root: Path,
    result: dict[str, Any],
    task: dict[str, Any],
    role: str,
    tool: dict[str, Any],
    receipt: dict[str, Any] | None,
) -> None:
    path = return_root / "tasks" / task["id"] / f"{role}-start.json"
    count = result["attempts"][role.replace("-", "_")]
    if count == 0:
        if path.exists() or path.is_symlink():
            raise ReplayError(f"task {task['id']} has an uncounted {role} attempt start")
        return
    start, _ = read_json(path, f"task {task['id']} {role} attempt start")
    validate_self_hash(start, "start_payload_sha256", f"task {task['id']} {role} attempt start")
    if (
        start.get("schema") != "koblitz_phase_b_magma_replay_attempt_start.v1"
        or start.get("role") != role
        or start.get("task_ordinal") != task["ordinal"]
        or start.get("blind_instance_id") != task["blind_instance_id"]
        or start.get("source_magma_sha256") != task["files"]["magma_boolean_f4"]["sha256"]
        or start.get("tool_sha256") != tool["sha256"]
        or start.get("watchdog_seconds") != 120
        or start.get("attempt_number") != 1
        or start.get("retry_allowed") is not False
        or (receipt is not None and start.get("command") != receipt["process"]["command"])
    ):
        raise ReplayError(f"task {task['id']} {role} attempt start changed")


def classify_local_f4(receipt: dict[str, Any], stdout: str, stderr: str) -> tuple[str, dict[str, Any] | None]:
    process = receipt["process"]
    if process["timed_out"]:
        return "timeout_inconclusive", None
    terminal_record = legacy.MATRIX.parse_magma_terminal(stdout)
    if process["returncode"] != 0 or "Runtime error" in stdout + stderr or "User error" in stdout + stderr:
        return "solver_error", terminal_record
    if terminal_record is None:
        return "terminal_certificate_missing", None
    if terminal_record["terminal_status"] == "unsat":
        return "unsat", terminal_record
    return "sat_basis_certificate_unverified_model", terminal_record


def derive_final_status(
    result: dict[str, Any],
    primary_status: str | None,
    primary_terminal: dict[str, Any] | None,
    model_status: str | None,
    model_terminal: dict[str, Any] | None,
    parsed_model: dict[str, Any] | None,
    validation_status: str | None,
) -> str:
    if "error" in result:
        return "execution_error_inconclusive"
    if primary_status is None:
        raise ReplayError(f"task {result.get('id')} lacks a primary F4 receipt")
    if primary_status == "unsat":
        return "unsat"
    if primary_status == "timeout_inconclusive":
        return "timeout_inconclusive"
    if primary_status != "sat_basis_certificate_unverified_model":
        return primary_status if primary_status.endswith("inconclusive") else f"{primary_status}_inconclusive"
    if model_status is None:
        raise ReplayError(f"task {result.get('id')} lacks F4 plus SAT(G) after a proper basis")
    if model_status == "timeout_inconclusive":
        return "f4_plus_sat_g_timeout_inconclusive"
    if primary_terminal is None or model_terminal is None or not legacy.same_f4_terminal(primary_terminal, model_terminal):
        return "f4_repetition_contradiction"
    if model_status != "sat_basis_certificate_unverified_model":
        return "f4_plus_sat_g_solver_error_inconclusive"
    if parsed_model is None:
        return "sat_g_terminal_missing_inconclusive"
    if parsed_model["status"] == "unsat":
        return "sat_g_contradiction"
    if validation_status is None:
        raise ReplayError(f"task {result.get('id')} lacks point validation for its SAT(G) model")
    if validation_status == "valid_point_witness":
        return "sat"
    if validation_status == "nonlifting_source_model":
        return "nonlifting_model_inconclusive"
    if validation_status == "timeout_inconclusive":
        return "point_validation_timeout_inconclusive"
    return f"point_validation_{validation_status}_inconclusive"


def validate_return(packet: Path, return_root: Path) -> dict[str, Any]:
    packet_info = validate_packet(packet)
    return_root = return_root.resolve(strict=True)
    seal, seal_bytes = read_json(return_root / "return-seal.json", "external return seal")
    if seal.get("schema") != RETURN_SEAL_SCHEMA or seal.get("status") != "external_results_frozen":
        raise ReplayError("external return is not sealed")
    validate_self_hash(seal, "seal_payload_sha256", "external return seal")
    actual = inventory(return_root, {"return-seal.json"})
    if actual != seal.get("inventory") or canonical_sha256(actual) != seal.get("inventory_sha256"):
        raise ReplayError("external return inventory changed after sealing")
    if seal.get("packet_seal_sha256") != packet_info["packet_seal_sha256"] or seal.get("packet_inventory_sha256") != packet_info["packet_inventory_sha256"]:
        raise ReplayError("external return is bound to a different solver packet")
    preflight_record, _ = read_json(return_root / "preflight.json", "licensed-host preflight")
    validate_preflight_record(preflight_record, packet_info)
    if sha256_file(return_root / "preflight.json", "preflight") != seal.get("preflight_sha256"):
        raise ReplayError("preflight differs from the return seal")
    binding, _ = read_json(return_root / "packet-binding.json", "packet binding")
    validate_self_hash(binding, "binding_payload_sha256", "packet binding")
    if (
        preflight_record.get("schema") != PREFLIGHT_SCHEMA
        or preflight_record.get("packet_seal_sha256") != packet_info["packet_seal_sha256"]
        or preflight_record.get("packet_inventory_sha256") != packet_info["packet_inventory_sha256"]
        or binding.get("packet_seal_sha256") != packet_info["packet_seal_sha256"]
        or binding.get("packet_inventory_sha256") != packet_info["packet_inventory_sha256"]
        or binding.get("packet_manifest_sha256") != packet_info["manifest_sha256"]
    ):
        raise ReplayError("preflight or packet binding refers to a different solver packet")
    task_ids = binding.get("selected_task_ids")
    if not isinstance(task_ids, list) or len(task_ids) != seal.get("selected_tasks"):
        raise ReplayError("external return task selection is inconsistent")
    task_map = {task["id"]: task for task in packet_info["manifest"]["tasks"]}
    full_order = [task["id"] for task in packet_info["manifest"]["tasks"]]
    if len(task_ids) != len(set(task_ids)) or task_ids != full_order[: len(task_ids)]:
        raise ReplayError("external return task selection is not a unique blind-order prefix")
    if preflight_record.get("synthetic_test_mode") is not True and task_ids != full_order:
        raise ReplayError("licensed production return must contain all 160 tasks in blind order")
    results = []
    for task_id in task_ids:
        if task_id not in task_map:
            raise ReplayError(f"external return contains unknown task {task_id}")
        result, _ = read_json(return_root / "tasks" / task_id / "task-result.json", f"task {task_id}")
        if (
            result.get("schema") != TASK_SCHEMA
            or result.get("id") != task_id
            or result.get("ordinal") != task_map[task_id]["ordinal"]
            or result.get("blind_instance_id") != task_map[task_id]["blind_instance_id"]
            or result.get("cell_id") != task_map[task_id]["cell_id"]
        ):
            raise ReplayError(f"task {task_id} identity changed")
        attempts = result.get("attempts")
        if not isinstance(attempts, dict) or set(attempts) != {"primary_f4", "f4_plus_sat_g", "point_validation"} or any(type(count) is not int or not 0 <= count <= 1 for count in attempts.values()):
            raise ReplayError(f"task {task_id} violates one-attempt/no-retry policy")
        tools = preflight_record.get("tools", {})
        phase_receipts: dict[str, dict[str, Any] | None] = {}
        observed_primary_status: str | None = None
        observed_primary_terminal: dict[str, Any] | None = None
        observed_model_status: str | None = None
        observed_model_terminal: dict[str, Any] | None = None
        observed_model: dict[str, Any] | None = None
        observed_validation_status: str | None = None
        for field, tool in (("primary_f4", "magma"), ("f4_plus_sat_g", "magma"), ("point_validation", "backend")):
            phase = result.get(field)
            if not isinstance(phase, dict) or not isinstance(phase.get("receipt"), dict):
                phase_receipts[field] = None
                continue
            receipt = local_receipt(return_root, phase["receipt"], f"{task_id} {field}")
            phase_receipts[field] = receipt
            command = receipt["process"]["command"]
            if not isinstance(command, list) or not command or command[0] != tools[tool]["path"]:
                raise ReplayError(f"task {task_id} {field} used a different executable")
            if tool == "magma" and (len(command) != 5 or command[1:4] != ["-t", "1", "-b"]):
                raise ReplayError(f"task {task_id} {field} changed Magma thread or batch flags")
            if tool == "backend" and (len(command) != 4 or command[1] != "validate-model"):
                raise ReplayError(f"task {task_id} changed point-validation command")
            if tool == "magma" and not receipt["process"]["timed_out"]:
                validate_identity_markers(
                    receipt_stream(return_root, phase["receipt"], "stdout", f"{task_id} {field}"),
                    task_map[task_id],
                )
        for role, field, tool in (
            ("primary-f4", "primary_f4", "magma"),
            ("f4-plus-sat-g", "f4_plus_sat_g", "magma"),
            ("point-validation", "point_validation", "backend"),
        ):
            validate_attempt_start(
                return_root, result, task_map[task_id], role, tools[tool],
                phase_receipts[field],
            )
        primary = result.get("primary_f4")
        if isinstance(primary, dict) and primary.get("status") != "retained_after_execution_error":
            primary_script = return_root / Path(primary["script"]["path"])
            source_magma = resolve_packet_file(
                packet_info["packet"], task_map[task_id]["files"]["magma_boolean_f4"],
                f"{task_id} Magma input",
            )
            expected_primary = add_identity_markers(
                regular_bytes(source_magma, f"{task_id} Magma input").decode(), task_map[task_id]
            ).encode()
            if regular_bytes(primary_script, f"{task_id} primary script") != expected_primary or file_record(primary_script, return_root, f"{task_id} primary script") != primary["script"]:
                raise ReplayError(f"task {task_id} identity-bound primary script changed")
            command = phase_receipts["primary_f4"]["process"]["command"]
            if command[:4] != [tools["magma"]["path"], "-t", "1", "-b"] or len(command) != 5 or not absolute_path_has_suffix(command[4], primary["script"]["path"]):
                raise ReplayError(f"task {task_id} primary F4 command used a different script")
            stdout = receipt_stream(return_root, primary["receipt"], "stdout", f"{task_id} primary F4")
            stderr = receipt_stream(return_root, primary["receipt"], "stderr", f"{task_id} primary F4")
            status, terminal_record = classify_local_f4(phase_receipts["primary_f4"], stdout, stderr)
            observed_primary_status = status
            observed_primary_terminal = terminal_record
            if primary.get("status") != status or primary.get("terminal") != terminal_record:
                raise ReplayError(f"task {task_id} primary F4 interpretation changed")
        model_phase = result.get("f4_plus_sat_g")
        if isinstance(model_phase, dict) and model_phase.get("status") != "retained_after_execution_error":
            model_script = return_root / Path(model_phase["script"]["path"])
            expected_model = add_identity_markers(
                legacy.render_witness_script(
                    source_magma, [], int(task_map[task_id]["source_variables"])
                ),
                task_map[task_id],
            ).encode()
            if regular_bytes(model_script, f"{task_id} F4 plus SAT(G) script") != expected_model or file_record(model_script, return_root, f"{task_id} F4 plus SAT(G) script") != model_phase["script"]:
                raise ReplayError(f"task {task_id} F4 plus SAT(G) script changed")
            command = phase_receipts["f4_plus_sat_g"]["process"]["command"]
            if command[:4] != [tools["magma"]["path"], "-t", "1", "-b"] or len(command) != 5 or not absolute_path_has_suffix(command[4], model_phase["script"]["path"]):
                raise ReplayError(f"task {task_id} F4 plus SAT(G) command used a different script")
            stdout = receipt_stream(return_root, model_phase["receipt"], "stdout", f"{task_id} F4 plus SAT(G)")
            stderr = receipt_stream(return_root, model_phase["receipt"], "stderr", f"{task_id} F4 plus SAT(G)")
            status, terminal_record = classify_local_f4(phase_receipts["f4_plus_sat_g"], stdout, stderr)
            parsed_model = legacy.parse_model_terminal(stdout, int(task_map[task_id]["source_variables"]), 0)
            observed_model_status = status
            observed_model_terminal = terminal_record
            observed_model = parsed_model
            if (
                model_phase.get("f4_status") != status
                or model_phase.get("f4_terminal") != terminal_record
                or model_phase.get("model_terminal") != parsed_model
            ):
                raise ReplayError(f"task {task_id} F4 plus SAT(G) interpretation changed")
        validation_phase = result.get("point_validation")
        if isinstance(validation_phase, dict) and validation_phase.get("status") != "retained_after_execution_error":
            assignment_record = validation_phase.get("assignment")
            if not isinstance(assignment_record, dict):
                raise ReplayError(f"task {task_id} point validation lacks its assignment")
            assignment_path = return_root / Path(assignment_record["path"])
            if file_record(assignment_path, return_root, f"{task_id} assignment") != assignment_record:
                raise ReplayError(f"task {task_id} assignment identity changed")
            try:
                assignment_value = json.loads(regular_bytes(assignment_path, f"{task_id} assignment"))
            except json.JSONDecodeError as error:
                raise ReplayError(f"task {task_id} assignment is invalid JSON") from error
            if not isinstance(observed_model, dict) or assignment_value != observed_model.get("assignment"):
                raise ReplayError(f"task {task_id} validator assignment differs from SAT(G)")
            stdout = receipt_stream(return_root, validation_phase["receipt"], "stdout", f"{task_id} validation")
            stderr = receipt_stream(return_root, validation_phase["receipt"], "stderr", f"{task_id} validation")
            manifest_path = resolve_packet_file(
                packet_info["packet"], task_map[task_id]["files"]["manifest"],
                f"{task_id} manifest",
            )
            manifest, _ = read_json(manifest_path, f"{task_id} manifest")
            process = phase_receipts["point_validation"]["process"]
            command = process["command"]
            if (
                len(command) != 4
                or command[:2] != [tools["backend"]["path"], "validate-model"]
                or not absolute_path_has_suffix(
                    command[2], task_map[task_id]["files"]["manifest"]["path"]
                )
                or not absolute_path_has_suffix(command[3], assignment_record["path"])
            ):
                raise ReplayError(f"task {task_id} point validator used different inputs")
            interpreted = legacy.MATRIX.assignment_validation_status({
                "stdout": stdout, "stderr": stderr,
                "returncode": process["returncode"], "timed_out": process["timed_out"],
                "metrics": process["metrics"], "command": process["command"],
            }, manifest, assignment_path)
            observed_validation_status = interpreted["status"]
            if validation_phase.get("status") != interpreted["status"] or validation_phase.get("report") != interpreted.get("report"):
                raise ReplayError(f"task {task_id} point-validation interpretation changed")
        derived_status = derive_final_status(
            result,
            observed_primary_status,
            observed_primary_terminal,
            observed_model_status,
            observed_model_terminal,
            observed_model,
            observed_validation_status,
        )
        if result.get("final_status") != derived_status:
            raise ReplayError(
                f"task {task_id} final status differs from its verified process evidence"
            )
        if "error" not in result:
            model_attempted = int(observed_primary_status == "sat_basis_certificate_unverified_model")
            validation_attempted = int(
                model_attempted == 1
                and observed_model_status == "sat_basis_certificate_unverified_model"
                and observed_primary_terminal is not None
                and observed_model_terminal is not None
                and legacy.same_f4_terminal(observed_primary_terminal, observed_model_terminal)
                and isinstance(observed_model, dict)
                and observed_model.get("status") == "sat"
            )
            expected_attempts = {
                "primary_f4": 1,
                "f4_plus_sat_g": model_attempted,
                "point_validation": validation_attempted,
            }
            if attempts != expected_attempts:
                raise ReplayError(f"task {task_id} phase sequence differs from its verified evidence")
        elif attempts["f4_plus_sat_g"] > attempts["primary_f4"] or attempts["point_validation"] > attempts["f4_plus_sat_g"]:
            raise ReplayError(f"task {task_id} failed-attempt sequence is impossible")
        expected_contradictions: list[str] = []
        if derived_status == "f4_repetition_contradiction":
            expected_contradictions.append("F4 plus SAT(G) did not reproduce the primary F4 terminal")
        if derived_status == "sat_g_contradiction":
            expected_contradictions.append("proper F4 basis but SAT(G) returned UNSAT")
        if result.get("contradictions") != expected_contradictions:
            raise ReplayError(f"task {task_id} contradiction record differs from verified evidence")
        results.append(result)
    summary, _ = read_json(return_root / "summary.json", "external return summary")
    if sha256_file(return_root / "summary.json", "summary") != seal.get("summary_sha256") or summary != summarize(results, preflight_record, len(task_ids)) | {"finished_at": summary.get("finished_at")}:
        raise ReplayError("external return summary does not reconstruct from sealed tasks")
    if seal.get("terminal_task_records") != len(results) or seal.get("one_attempt_no_retry") is not True:
        raise ReplayError("external return terminal count or attempt policy changed")
    expected_claim = dict(packet_info["protocol"]["claim_boundary"])
    expected_claim["licensed_magma_execution_completed"] = summary[
        "licensed_magma_execution_completed"
    ]
    if seal.get("claim_boundary") != expected_claim:
        raise ReplayError("external return widened or suppressed its claim boundary")
    return {
        "return_seal_sha256": hashlib.sha256(seal_bytes).hexdigest(),
        "return_inventory_sha256": seal["inventory_sha256"],
        "packet": packet_info,
        "preflight": preflight_record,
        "summary": summary,
        "results": results,
    }


def score(
    packet: Path,
    return_root: Path,
    truth_score_path: Path,
    truth_score_seal_path: Path,
    output: Path,
) -> dict[str, Any]:
    validated = validate_return(packet, return_root)
    protocol = validated["packet"]["protocol"]
    if sha256_file(truth_score_seal_path.resolve(strict=True), "truth score seal") != protocol["source"]["terminal_score_seal_sha256"]:
        raise ReplayError("post-seal scorer received a different truth score seal")
    truth_seal, _ = read_json(truth_score_seal_path.resolve(strict=True), "truth score seal")
    truth, truth_bytes = read_json(truth_score_path.resolve(strict=True), "truth score")
    if hashlib.sha256(truth_bytes).hexdigest() != truth_seal.get("score_sha256") or truth_seal.get("score_sha256") != protocol["source"]["terminal_score_sha256"]:
        raise ReplayError("post-seal truth score differs from its frozen seal")
    classes: dict[str, str] = {}
    for row in truth.get("rows", []):
        blind_id = row.get("blind_instance_id")
        target_class = row.get("target_class")
        if target_class not in {"decomposable", "nondecomposable"}:
            raise ReplayError("truth score contains an invalid target class")
        if blind_id in classes and classes[blind_id] != target_class:
            raise ReplayError("truth score disagrees across backends")
        classes[blind_id] = target_class
    if len(classes) != 160:
        raise ReplayError("truth score does not classify exactly 160 blind instances")
    rows = []
    counts: Counter[str] = Counter()
    for result in validated["results"]:
        blind_id = result["blind_instance_id"]
        target_class = classes.get(blind_id)
        if target_class is None:
            raise ReplayError(f"truth score lacks {blind_id}")
        solver_status = result.get("final_status")
        if solver_status == "sat":
            classification = "true_positive" if target_class == "decomposable" else "false_positive"
        elif solver_status == "unsat":
            classification = "true_negative" if target_class == "nondecomposable" else "false_negative"
        else:
            classification = "inconclusive"
        counts[classification] += 1
        rows.append({
            "ordinal": result["ordinal"],
            "blind_instance_id": blind_id,
            "cell_id": result["cell_id"],
            "solver_status": solver_status,
            "target_class": target_class,
            "classification": classification,
        })
    output = output.resolve()
    for protected, label in ((validated["packet"]["packet"], "solver packet"), (Path(return_root).resolve(), "external return")):
        try:
            output.relative_to(protected)
        except ValueError:
            pass
        else:
            raise ReplayError(f"post-seal score must not reside inside the {label}")
    if output.exists() or output.is_symlink():
        raise ReplayError(f"post-seal score output must be new: {output}")
    output.mkdir(parents=True)
    per_cell_classifications: dict[str, dict[str, int]] = {}
    for row in rows:
        cell_counts = per_cell_classifications.setdefault(row["cell_id"], {})
        label = row["classification"]
        cell_counts[label] = cell_counts.get(label, 0) + 1
    scored = with_self_hash({
        "schema": SCORE_SCHEMA,
        "status": "post_return_truth_scored",
        "created_at": now(),
        "packet_seal_sha256": validated["packet"]["packet_seal_sha256"],
        "return_seal_sha256": validated["return_seal_sha256"],
        "truth_score_sha256": protocol["source"]["terminal_score_sha256"],
        "selected_instances": len(rows),
        "full_panel_scored": len(rows) == 160,
        "classifications": dict(sorted(counts.items())),
        "per_cell_classifications": {
            cell: dict(sorted(values.items()))
            for cell, values in sorted(per_cell_classifications.items())
        },
        "rows": rows,
        "licensed_magma_execution_completed": validated["summary"]["licensed_magma_execution_completed"],
        "independent_external_reproduction_satisfied": False,
        "novelty_review_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }, "score_payload_sha256")
    write_json_new(output / "score.json", scored)
    score_inventory = inventory(output, {"score-seal.json"})
    score_claim = dict(protocol["claim_boundary"])
    score_claim["licensed_magma_execution_completed"] = scored[
        "licensed_magma_execution_completed"
    ]
    seal = with_self_hash({
        "schema": SCORE_SEAL_SCHEMA,
        "status": "post_return_score_frozen",
        "score_path": "score.json",
        "score_sha256": sha256_file(output / "score.json", "Magma replay score"),
        "return_seal_sha256": validated["return_seal_sha256"],
        "inventory": score_inventory,
        "inventory_sha256": canonical_sha256(score_inventory),
        "claim_boundary": score_claim,
    }, "seal_payload_sha256")
    write_json_new(output / "score-seal.json", seal)
    return seal


def license_statement(path: Path) -> str:
    try:
        return regular_bytes(path.resolve(strict=True), "license-access statement").decode().strip()
    except UnicodeDecodeError as error:
        raise ReplayError("license-access statement must be UTF-8 text") from error


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    dry = subparsers.add_parser("dry-run", help="verify committed sources without Magma")
    dry.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    dry.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    prepare = subparsers.add_parser("prepare", help="create the sealed 160-input solver packet")
    prepare.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    prepare.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    prepare.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    prepare.add_argument("--score-seal", type=Path)
    prepare.add_argument("--output", type=Path, required=True)
    verify_packet_parser = subparsers.add_parser("verify-packet")
    verify_packet_parser.add_argument("--packet", type=Path, required=True)
    for name in ("preflight", "run"):
        command = subparsers.add_parser(name)
        command.add_argument("--packet", type=Path, required=True)
        command.add_argument("--magma", required=True)
        command.add_argument("--minisat", required=True)
        command.add_argument("--backend", type=Path, default=DEFAULT_BACKEND)
        command.add_argument("--meter", type=Path, default=DEFAULT_METER)
        command.add_argument("--license-access-statement-file", type=Path, required=True)
        command.add_argument("--synthetic-test-mode", action="store_true", help=argparse.SUPPRESS)
        if name == "run":
            command.add_argument("--output", type=Path, required=True)
            command.add_argument("--max-tasks", type=int, help=argparse.SUPPRESS)
    verify_return_parser = subparsers.add_parser("verify-return")
    verify_return_parser.add_argument("--packet", type=Path, required=True)
    verify_return_parser.add_argument("--return-root", type=Path, required=True)
    score_parser = subparsers.add_parser("score")
    score_parser.add_argument("--packet", type=Path, required=True)
    score_parser.add_argument("--return-root", type=Path, required=True)
    score_parser.add_argument("--truth-score", type=Path, default=DEFAULT_BUNDLE / "score/score.json")
    score_parser.add_argument("--truth-score-seal", type=Path, default=DEFAULT_BUNDLE / "score/score-seal.json")
    score_parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.command == "dry-run":
            result = dry_run(args.protocol, args.bundle)
        elif args.command == "prepare":
            score_seal = args.score_seal or args.bundle / "score/score-seal.json"
            result = build_packet(args.protocol, args.bundle, args.run_root, score_seal, args.output)
        elif args.command == "verify-packet":
            validated = validate_packet(args.packet)
            result = {
                "status": "pass", "packet_seal_sha256": validated["packet_seal_sha256"],
                "packet_inventory_sha256": validated["packet_inventory_sha256"],
                "instance_count": validated["manifest"]["instance_count"],
            }
        elif args.command == "preflight":
            result = preflight(
                args.packet, args.magma, args.minisat, args.backend, args.meter,
                license_statement(args.license_access_statement_file),
                synthetic_test_mode=args.synthetic_test_mode,
            )
        elif args.command == "run":
            result = run(
                args.packet, args.output, args.magma, args.minisat, args.backend,
                args.meter, license_statement(args.license_access_statement_file),
                synthetic_test_mode=args.synthetic_test_mode, max_tasks=args.max_tasks,
            )
        elif args.command == "verify-return":
            validated = validate_return(args.packet, args.return_root)
            result = {
                "status": "pass", "return_seal_sha256": validated["return_seal_sha256"],
                "return_inventory_sha256": validated["return_inventory_sha256"],
                "terminal_task_records": len(validated["results"]),
            }
        else:
            result = score(
                args.packet, args.return_root, args.truth_score,
                args.truth_score_seal, args.output,
            )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (ReplayError, terminal.EvidenceError, legacy.ExternalMagmaError, OSError, subprocess.SubprocessError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
