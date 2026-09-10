#!/usr/bin/env python3
"""Verify Stage 19 source, Git, transport, crash, and no-retry custody fail closed."""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import stat
import subprocess
import sys
from typing import Any
import urllib.parse


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_ARTIFACT = HERE / "stage-19-magma-calculator-panel-amendment-02-20260910"
SUMMARY_SCHEMA = "koblitz_magma_calculator_stage19_summary.v2"


class VerificationError(RuntimeError):
    """The retained panel is inconsistent, ambiguous, or outside its claim boundary."""


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise VerificationError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RENDERER = load_module("stage19_renderer_for_verifier", HERE / "render_stage19_magma_calculator_panel.py")
RUNNER = load_module("stage19_runner_for_verifier", HERE / "run_stage19_magma_calculator_panel.py")
CHILD = RUNNER.CHILD
STAGE15 = RENDERER.STAGE15


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise VerificationError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise VerificationError(f"expected JSON object in {path}")
    return value


def regular_record(path: Path, relative_to: Path | None = None) -> dict:
    try:
        info = path.lstat()
    except FileNotFoundError as error:
        raise VerificationError(f"missing retained file {path}") from error
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise VerificationError(f"retained path is not a single-link regular file: {path}")
    data = path.read_bytes()
    return {
        "path": str(path.relative_to(relative_to)) if relative_to is not None else str(path),
        "bytes": len(data),
        "sha256": sha256_bytes(data),
    }


def parse_time(value: Any, label: str) -> datetime:
    if not isinstance(value, str) or not value.endswith("Z"):
        raise VerificationError(f"{label} is not a UTC timestamp")
    try:
        return datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise VerificationError(f"{label} is not an ISO timestamp") from error


def git(*args: str) -> bytes:
    result = subprocess.run(["git", *args], cwd=REPO, capture_output=True, check=False)
    if result.returncode != 0:
        raise VerificationError(f"git {' '.join(args)} failed: {result.stderr.decode(errors='replace').strip()}")
    return result.stdout


def git_tree_entry(revision: str, relative: str) -> tuple[str, str]:
    line = git("ls-tree", revision, "--", relative).decode().strip()
    if not line or "\t" not in line:
        raise VerificationError(f"missing Git tree entry {revision}:{relative}")
    metadata, observed = line.split("\t", 1)
    fields = metadata.split()
    if len(fields) != 3 or fields[1] != "blob" or observed != relative:
        raise VerificationError(f"unexpected Git tree entry {revision}:{relative}")
    return fields[0], fields[2]


def verify_prepared(
    artifact: Path, allow_external_test_artifact: bool = False
) -> tuple[dict, dict[str, bytes]]:
    if not allow_external_test_artifact and artifact.resolve() != DEFAULT_ARTIFACT.resolve():
        raise VerificationError("Stage 19 verifier requires the repository artifact")
    expected_plan, rendered_inputs = RENDERER.build_plan()
    plan_path = artifact / "plan.json"
    plan = read_json(plan_path)
    if plan != expected_plan or regular_record(plan_path)["sha256"] != sha256_bytes(canonical_bytes(expected_plan)):
        raise VerificationError("plan.json differs from deterministic source rendering")
    input_dir = artifact / "inputs"
    if input_dir.is_symlink() or not input_dir.is_dir():
        raise VerificationError("named-input directory is missing or a symlink")
    expected_names = {Path(task["named_input"]["path"]).name for task in plan["tasks"]}
    actual = list(input_dir.iterdir())
    if {path.name for path in actual} != expected_names:
        raise VerificationError("named-input inventory is missing, duplicated, or unexpected")
    for task in plan["tasks"]:
        path = artifact / task["named_input"]["path"]
        expected = rendered_inputs[path.name]
        observed = regular_record(path)
        if observed["sha256"] != sha256_bytes(expected) or path.read_bytes() != expected:
            raise VerificationError(f"named input changed for {task['id']}")
    allowed = {
        "plan.json", "inputs", "prepared-summary.json", "execution-manifest.json",
        "run.json", "attempts", "summary.json",
    }
    unexpected = {path.name for path in artifact.iterdir()} - allowed
    if unexpected:
        raise VerificationError(f"unexpected top-level artifact entries: {sorted(unexpected)}")
    regular_record(artifact / "prepared-summary.json")
    return plan, rendered_inputs


def verify_execution_binding(binding: dict) -> dict:
    if not isinstance(binding, dict) or binding.get("schema") != "koblitz_magma_calculator_stage19_execution_binding.v1":
        raise VerificationError("run lacks the execution binding")
    body = dict(binding)
    claimed_hash = body.pop("binding_sha256", None)
    if claimed_hash != sha256_bytes(canonical_bytes(body)):
        raise VerificationError("execution-binding hash changed")
    first = binding.get("first_preexecution_commit")
    execution = binding.get("execution_commit")
    if not all(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{40}", value) for value in (first, execution)):
        raise VerificationError("execution binding contains an invalid commit")
    if git("rev-parse", f"{first}^{{tree}}").decode().strip() != binding.get("first_preexecution_tree"):
        raise VerificationError("first pre-execution tree changed")
    if git("rev-parse", f"{execution}^{{tree}}").decode().strip() != binding.get("execution_tree"):
        raise VerificationError("execution tree changed")
    manifest_record = binding.get("execution_manifest")
    if not isinstance(manifest_record, dict):
        raise VerificationError("execution manifest record is missing")
    manifest_path = REPO / str(manifest_record.get("path", ""))
    observed_manifest = regular_record(manifest_path, REPO)
    if any(observed_manifest.get(key) != manifest_record.get(key) for key in ("path", "bytes", "sha256")):
        raise VerificationError("current execution-manifest bytes changed")
    manifest = read_json(manifest_path)
    if manifest.get("schema") != RUNNER.EXECUTION_MANIFEST_SCHEMA or manifest.get("first_preexecution_commit") != first:
        raise VerificationError("execution manifest identity changed")
    delta = [line for line in git("diff", "--name-status", first, execution).decode().splitlines() if line]
    if delta != [f"A\t{observed_manifest['path']}"] or binding.get("tree_delta") != [
        {"status": "A", "path": observed_manifest["path"]}
    ]:
        raise VerificationError("execution commit delta is not exactly the manifest addition")
    mode, blob = git_tree_entry(execution, observed_manifest["path"])
    if mode != manifest_record.get("git_mode") or blob != manifest_record.get("git_blob_oid"):
        raise VerificationError("execution-manifest Git blob changed")
    if git("show", f"{execution}:{observed_manifest['path']}") != manifest_path.read_bytes():
        raise VerificationError("current execution manifest differs from its execution commit")
    records = manifest.get("relevant_blobs")
    if not isinstance(records, list) or binding.get("relevant_blob_count") != len(records):
        raise VerificationError("execution relevant-blob count changed")
    if {record.get("path") for record in records if isinstance(record, dict)} != RUNNER.expected_execution_bound_paths():
        raise VerificationError("execution manifest required path set changed")
    ca_bundle = CHILD.ca_bundle_record()
    if (
        manifest.get("external_dependencies", {}).get("ca_bundle") != ca_bundle
        or binding.get("ca_bundle") != ca_bundle
    ):
        raise VerificationError("execution CA bundle identity changed")
    archived_probe = RUNNER.TLS_PROBE.archived_records()
    if (
        manifest.get("archived_tls_get_probe") != archived_probe
    ):
        raise VerificationError("execution archived TLS GET probe binding changed")
    try:
        fresh_probe = RUNNER.TLS_PROBE.verify_fresh(
            manifest.get("fresh_preexecution_tls_get_probe")
        )
    except RUNNER.TLS_PROBE.ProbeError as error:
        raise VerificationError(str(error)) from error
    if binding.get("fresh_preexecution_tls_get_probe") != fresh_probe:
        raise VerificationError("execution fresh TLS GET probe binding changed")
    for record in records:
        relative = record.get("path") if isinstance(record, dict) else None
        if not isinstance(relative, str):
            raise VerificationError("malformed relevant-blob record")
        base_mode, base_blob = git_tree_entry(first, relative)
        execution_mode, execution_blob = git_tree_entry(execution, relative)
        current = regular_record(REPO / relative)
        if (
            base_mode != record.get("git_mode") or execution_mode != record.get("git_mode")
            or base_blob != record.get("git_blob_oid") or execution_blob != record.get("git_blob_oid")
            or current["bytes"] != record.get("bytes") or current["sha256"] != record.get("sha256")
            or git("show", f"{execution}:{relative}") != (REPO / relative).read_bytes()
        ):
            raise VerificationError(f"execution-bound current blob changed: {relative}")
    return manifest


def expected_request(plan: dict, task: dict, input_bytes: bytes, binding: dict) -> dict:
    service = plan["service"]
    body = urllib.parse.urlencode({service["form_field"]: input_bytes.decode("ascii")}).encode("ascii")
    return {
        "endpoint": service["endpoint"], "method": "POST", "form_field": service["form_field"],
        "content_type": "application/x-www-form-urlencoded", "user_agent": service["user_agent"],
        "input": task["named_input"], "body_bytes": len(body), "body_sha256": sha256_bytes(body),
        "hard_parent_watchdog_seconds": service["client_timeout_seconds"],
        "response_byte_limit": service["max_response_bytes"],
        "ca_bundle": binding["ca_bundle"],
        "execution_commit": binding["execution_commit"], "execution_tree": binding["execution_tree"],
        "execution_binding_sha256": binding["binding_sha256"],
        "execution_manifest_sha256": binding["execution_manifest"]["sha256"],
    }


def classify_response(path: Path, status: int | None, final_url: str | None, plan: dict, task: dict) -> dict:
    endpoint = plan["service"]["endpoint"]
    if status != 200:
        return {"outcome": "http_error", "clean_terminal": False}
    if final_url != endpoint:
        return {"outcome": "final_url_mismatch", "clean_terminal": False}
    try:
        response = STAGE15.parse_calculator_xml(path)
    except STAGE15.VerificationError as error:
        return {"outcome": "invalid_calculator_xml", "clean_terminal": False, "diagnostic": str(error)}
    if response["service"]["warning"] is not None or response["service"]["alert"] is not None:
        return {"outcome": "service_diagnostic", "clean_terminal": False, "service": response["service"]}
    prefix = (
        f"KOBLITZ_MAGMA_TASK_ID={task['id']}\n"
        f"KOBLITZ_MAGMA_SOURCE_SHA256={task['source_instance_sha256']}\n"
    )
    if not response["output"].startswith(prefix):
        return {
            "outcome": "identity_or_terminal_mismatch", "clean_terminal": False,
            "diagnostic": "response task/source identity markers are absent or mismatched",
            "service": response["service"],
        }
    remainder = response["output"][len(prefix):]
    if "KOBLITZ_MAGMA_TASK_ID=" in remainder or "KOBLITZ_MAGMA_SOURCE_SHA256=" in remainder:
        return {
            "outcome": "identity_or_terminal_mismatch", "clean_terminal": False,
            "diagnostic": "response contains duplicate identity markers", "service": response["service"],
        }
    try:
        STAGE15.require_clean_f4_output(remainder, int(task["seed"]))
    except STAGE15.VerificationError as error:
        return {"outcome": "identity_or_terminal_mismatch", "clean_terminal": False, "diagnostic": str(error), "service": response["service"]}
    terminal = STAGE15.parse_magma_terminal(remainder)
    if terminal is None:
        return {"outcome": "identity_or_terminal_mismatch", "clean_terminal": False, "diagnostic": "strict seven-marker parser rejected response", "service": response["service"]}
    if terminal["terminal_status"] != "sat":
        return {"outcome": "planted_unsat_contradiction", "clean_terminal": False, "terminal": terminal, "service": response["service"]}
    return {
        "outcome": "clean_f4_terminal", "clean_terminal": True,
        "classification": "sat_basis_certificate_unverified_model", "response_identity_verified": True,
        "terminal": terminal, "service": response["service"],
    }


def expected_transport_command(artifact: Path, task: dict) -> list[str]:
    return [
        str(Path(sys.executable).resolve()), str((HERE / "post_stage19_magma_calculator_request.py").resolve()),
        "--attempt", str((artifact / "attempts" / RUNNER.task_stem(task)).resolve()),
        "--execute-child",
    ]


def verify_envelope(attempt: Path, plan: dict, task: dict, request: dict) -> tuple[dict, Path | None]:
    path = attempt / "transport-envelope.json"
    envelope = read_json(path)
    if envelope.get("schema") != CHILD.ENVELOPE_SCHEMA:
        raise VerificationError(f"{task['id']}: envelope schema changed")
    child_request = {
        "endpoint": request["endpoint"], "method": "POST", "form_field": request["form_field"],
        "user_agent": request["user_agent"], "input_bytes": task["named_input"]["bytes"],
        "input_sha256": task["named_input"]["sha256"], "body_bytes": request["body_bytes"],
        "body_sha256": request["body_sha256"], "socket_timeout_seconds": CHILD.SOCKET_TIMEOUT_SECONDS,
        "ca_bundle": request["ca_bundle"],
    }
    if envelope.get("request") != child_request or envelope.get("response_byte_limit") != plan["service"]["max_response_bytes"]:
        raise VerificationError(f"{task['id']}: envelope request or byte limit changed")
    http = envelope.get("http")
    if not isinstance(http, dict) or set(http) != {
        "status", "final_url", "selected_headers", "transport_error"
    }:
        raise VerificationError(f"{task['id']}: malformed HTTP envelope")
    status = http["status"]
    if status is not None and (not isinstance(status, int) or isinstance(status, bool)):
        raise VerificationError(f"{task['id']}: invalid HTTP status")
    headers = http["selected_headers"]
    if (
        not isinstance(headers, dict)
        or set(headers) - set(CHILD.SAFE_RESPONSE_HEADERS)
        or any(not isinstance(value, str) for value in headers.values())
    ):
        raise VerificationError(f"{task['id']}: unsafe HTTP header custody")
    started = parse_time(envelope.get("started_at"), "transport started_at")
    if parse_time(envelope.get("finished_at"), "transport finished_at") < started:
        raise VerificationError(f"{task['id']}: negative transport chronology")
    response = envelope.get("response")
    if response is None:
        return envelope, None
    response_path = attempt / str(response.get("path", ""))
    observed = regular_record(response_path, attempt)
    if any(observed.get(key) != response.get(key) for key in ("path", "bytes", "sha256")):
        raise VerificationError(f"{task['id']}: envelope response changed")
    if response.get("response_byte_limit") != plan["service"]["max_response_bytes"]:
        raise VerificationError(f"{task['id']}: response limit changed")
    expected_name = (
        "response.partial"
        if response.get("body_limit_exceeded") is True
        else "response.xml" if status == 200 else "response.body"
    )
    if (
        response.get("path") != expected_name
        or response.get("complete") is not (response.get("body_limit_exceeded") is False)
        or response.get("bytes", CHILD.MAX_RESPONSE_BYTES + 1) > CHILD.MAX_RESPONSE_BYTES
        or (
            response.get("body_limit_exceeded") is True
            and response.get("bytes") != CHILD.MAX_RESPONSE_BYTES
        )
    ):
        raise VerificationError(f"{task['id']}: response completeness/path contract changed")
    return envelope, response_path


def verify_normal_receipt(artifact: Path, plan: dict, task: dict, attempt: Path, start: dict, receipt: dict) -> dict:
    process = receipt.get("transport_process")
    if not isinstance(process, dict):
        raise VerificationError(f"{task['id']}: normal receipt lacks transport process")
    metrics = process.get("metrics")
    command = expected_transport_command(artifact, task)
    if not isinstance(metrics, dict) or metrics.get("command") != command or metrics.get("watchdog_seconds") != 75.0:
        raise VerificationError(f"{task['id']}: transport command/watchdog changed")
    if type(metrics.get("timed_out")) is not bool or type(metrics.get("orphan_group_terminated")) is not bool:
        raise VerificationError(f"{task['id']}: transport terminal flags are invalid")
    resources = metrics.get("metrics")
    if not isinstance(resources, dict) or resources.get("meter") != "fresh-process getrusage(RUSAGE_CHILDREN)":
        raise VerificationError(f"{task['id']}: transport resource meter changed")
    for key in ("wall_seconds", "user_seconds", "system_seconds", "total_core_seconds", "single_core_seconds"):
        value = resources.get(key)
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value < 0:
            raise VerificationError(f"{task['id']}: invalid transport metric {key}")
    if not math.isclose(resources["total_core_seconds"], resources["user_seconds"] + resources["system_seconds"], rel_tol=1e-9, abs_tol=1e-9):
        raise VerificationError(f"{task['id']}: transport core-seconds do not add up")
    expected_files = {
        "stdout": "transport.stdout", "stderr": "transport.stderr", "metrics_file": "transport-metrics.json"
    }
    for key, filename in expected_files.items():
        if process.get(key) != regular_record(attempt / filename, attempt):
            raise VerificationError(f"{task['id']}: transport {key} changed")
    retained = [
        regular_record(path, attempt) for path in sorted(attempt.iterdir())
        if path.name != "receipt.json"
    ]
    if receipt.get("retained_files") != retained:
        raise VerificationError(f"{task['id']}: retained transport inventory changed")
    authorization_paths = [
        path for path in (
            attempt / "launch-authorization.json",
            attempt / "launch-authorization.consumed.json",
        ) if path.is_file()
    ]
    if len(authorization_paths) != 1:
        raise VerificationError(f"{task['id']}: launch authorization is missing or duplicated")
    if receipt.get("launch_authorization") != regular_record(authorization_paths[0], attempt):
        raise VerificationError(f"{task['id']}: launch authorization receipt changed")
    expected_authorization = CHILD.expected_launch_authorization(
        attempt, plan, task, start, artifact=artifact
    )
    if authorization_paths[0].read_bytes() != CHILD.canonical_bytes(expected_authorization):
        raise VerificationError(f"{task['id']}: launch authorization content changed")
    envelope_path = attempt / "transport-envelope.json"
    envelope_record = regular_record(envelope_path, attempt) if envelope_path.is_file() else None
    if receipt.get("transport_envelope") != envelope_record:
        raise VerificationError(f"{task['id']}: retained envelope was not bound")
    response_paths = [
        path for path in (attempt / "response.xml", attempt / "response.body", attempt / "response.partial")
        if path.is_file()
    ]
    if len(response_paths) > 1:
        raise VerificationError(f"{task['id']}: multiple retained response bodies")
    response_record = regular_record(response_paths[0], attempt) if response_paths else None
    if receipt.get("response") != response_record:
        raise VerificationError(f"{task['id']}: retained response was not bound")
    classification: dict
    envelope_record = receipt.get("transport_envelope")
    response_record = receipt.get("response")
    if metrics.get("timed_out") is True:
        classification = {"outcome": "transport_timeout", "clean_terminal": False}
        expected_boundary = "transport_timeout_with_retained_files"
    elif metrics.get("orphan_group_terminated") is True:
        classification = {"outcome": "transport_orphan_terminated", "clean_terminal": False}
        expected_boundary = "transport_orphan_with_retained_files"
    elif metrics.get("returncode") != 0:
        classification = {"outcome": "transport_child_error", "clean_terminal": False}
        expected_boundary = "transport_child_error_with_retained_files"
    elif (attempt / "transport-envelope.json").is_file():
        expected_boundary = None
        envelope, response_path = verify_envelope(attempt, plan, task, start["request"])
        if envelope_record != regular_record(attempt / "transport-envelope.json", attempt):
            raise VerificationError(f"{task['id']}: envelope record changed")
        if response_path is not None and response_record != regular_record(response_path, attempt):
            raise VerificationError(f"{task['id']}: response record changed")
        if isinstance(envelope.get("response"), dict) and envelope["response"].get("body_limit_exceeded") is True:
            classification = {"outcome": "response_body_limit_exceeded", "clean_terminal": False}
        elif response_path is None:
            classification = {"outcome": "transport_error", "clean_terminal": False}
        else:
            classification = classify_response(response_path, envelope["http"]["status"], envelope["http"]["final_url"], plan, task)
    else:
        expected_boundary = None
        classification = {"outcome": "transport_envelope_missing", "clean_terminal": False}
    for key, value in classification.items():
        if receipt.get(key) != value:
            raise VerificationError(f"{task['id']}: normal receipt classification changed at {key}")
    if receipt.get("recovery_boundary") != expected_boundary:
        raise VerificationError(f"{task['id']}: transport recovery boundary changed")
    temporary_files = [path for path in attempt.iterdir() if path.name.startswith(".") and ".tmp-" in path.name]
    if classification["clean_terminal"] is True and temporary_files:
        raise VerificationError(f"{task['id']}: clean receipt retained an atomic temporary")
    if classification["clean_terminal"] is True and authorization_paths[0].name != "launch-authorization.consumed.json":
        raise VerificationError(f"{task['id']}: clean receipt did not consume launch authorization")
    return classification


def verify_recovery_receipt(task: dict, attempt: Path, receipt: dict) -> dict:
    expected_records = [regular_record(path, attempt) for path in sorted(attempt.iterdir()) if path.name != "receipt.json"]
    if receipt.get("recovery_files") != expected_records or receipt.get("retained_files") != expected_records:
        raise VerificationError(f"{task['id']}: recovery file inventory changed")
    names = {record["path"] for record in expected_records}
    responses = names & {"response.xml", "response.body", "response.partial"}
    temporary = any(name.startswith(".") and ".tmp-" in name for name in names)
    if temporary:
        outcome = "interrupted_during_atomic_write"
    elif responses and "transport-envelope.json" in names:
        outcome = "interrupted_with_retained_response"
    elif responses:
        outcome = "interrupted_with_unbound_response"
    elif "transport-envelope.json" in names:
        outcome = "interrupted_with_incomplete_transport_envelope"
    elif "transport-metrics.json" in names:
        outcome = "interrupted_after_transport_process"
    elif read_json(attempt / "attempt-start.json").get("recovered_empty_attempt_directory") is True:
        outcome = "interrupted_before_attempt_start"
    else:
        outcome = "interrupted_before_transport_receipt"
    if receipt.get("outcome") != outcome or receipt.get("recovery_boundary") != outcome:
        raise VerificationError(f"{task['id']}: recovery boundary changed")
    envelope_record = (
        regular_record(attempt / "transport-envelope.json", attempt)
        if "transport-envelope.json" in names else None
    )
    response_paths = [attempt / name for name in responses]
    response_record = regular_record(response_paths[0], attempt) if len(response_paths) == 1 else None
    if receipt.get("transport_envelope") != envelope_record or receipt.get("response") != response_record:
        raise VerificationError(f"{task['id']}: recovery omitted retained response/envelope custody")
    authorization_paths = [
        path for path in (
            attempt / "launch-authorization.json",
            attempt / "launch-authorization.consumed.json",
        ) if path.is_file()
    ]
    authorization_record = regular_record(authorization_paths[0], attempt) if len(authorization_paths) == 1 else None
    if receipt.get("launch_authorization") != authorization_record:
        raise VerificationError(f"{task['id']}: recovery launch authorization changed")
    return {"outcome": outcome, "clean_terminal": False}


def verify_attempt(artifact: Path, plan: dict, task: dict, attempt: Path, binding: dict) -> dict:
    if attempt.is_symlink() or not attempt.is_dir():
        raise VerificationError(f"{task['id']}: attempt is missing or a symlink")
    files = list(attempt.iterdir())
    if any(path.is_symlink() or (path.is_file() and path.stat().st_nlink != 1) for path in files):
        raise VerificationError(f"{task['id']}: attempt contains a symlink or hardlink")
    allowed_files = {
        "attempt-start.json", "request-body.bin", "launch-authorization.json",
        "launch-authorization.consumed.json", "transport.stdout", "transport.stderr",
        "transport-metrics.json", "transport-envelope.json", "response.xml", "response.body",
        "response.partial", "receipt.json",
    }
    unexpected = [path for path in files if path.name not in allowed_files]
    if any(not (path.name.startswith(".") and ".tmp-" in path.name) for path in unexpected):
        raise VerificationError(f"{task['id']}: unexpected attempt file")
    if not (attempt / "attempt-start.json").is_file() or not (attempt / "receipt.json").is_file():
        raise VerificationError(f"{task['id']}: incomplete attempt custody")
    start = read_json(attempt / "attempt-start.json")
    receipt = read_json(attempt / "receipt.json")
    input_bytes = (artifact / task["named_input"]["path"]).read_bytes()
    request = expected_request(plan, task, input_bytes, binding)
    base_start = {
        "schema": RUNNER.ATTEMPT_SCHEMA, "task_id": task["id"], "ordinal": task["ordinal"],
        "attempt_ordinal": 1, "started_at": start.get("started_at"), "request": request,
        "execution_commit": binding["execution_commit"], "execution_tree": binding["execution_tree"],
        "execution_binding_sha256": binding["binding_sha256"], "retry_permitted": False,
        "launch_nonce_sha256": start.get("launch_nonce_sha256"),
        "recovered_empty_attempt_directory": start.get("recovered_empty_attempt_directory"),
    }
    if start != base_start or type(start["recovered_empty_attempt_directory"]) is not bool:
        raise VerificationError(f"{task['id']}: attempt-start receipt changed")
    nonce_hash = start["launch_nonce_sha256"]
    if start["recovered_empty_attempt_directory"]:
        if nonce_hash is not None:
            raise VerificationError(f"{task['id']}: recovered empty attempt invented a nonce")
    elif not isinstance(nonce_hash, str) or re.fullmatch(r"[0-9a-f]{64}", nonce_hash) is None:
        raise VerificationError(f"{task['id']}: launch nonce hash is malformed")
    started = parse_time(start["started_at"], f"{task['id']} started_at")
    required = {
        "schema": RUNNER.RECEIPT_SCHEMA, "task_id": task["id"], "ordinal": task["ordinal"],
        "attempt_ordinal": 1, "started_at": start["started_at"], "request": request,
        "execution_commit": binding["execution_commit"], "execution_tree": binding["execution_tree"],
        "execution_binding_sha256": binding["binding_sha256"], "retry_permitted": False,
    }
    if any(receipt.get(key) != value for key, value in required.items()):
        raise VerificationError(f"{task['id']}: final receipt identity changed")
    if parse_time(receipt.get("finished_at"), f"{task['id']} finished_at") < started:
        raise VerificationError(f"{task['id']}: receipt chronology is negative")
    if (attempt / "request-body.bin").is_file():
        body = urllib.parse.urlencode({plan["service"]["form_field"]: input_bytes.decode("ascii")}).encode("ascii")
        if (attempt / "request-body.bin").read_bytes() != body:
            raise VerificationError(f"{task['id']}: request body changed")
    if "recovery_files" in receipt:
        classification = verify_recovery_receipt(task, attempt, receipt)
    else:
        classification = verify_normal_receipt(artifact, plan, task, attempt, start, receipt)
    if receipt.get("clean_terminal") is not classification["clean_terminal"]:
        raise VerificationError(f"{task['id']}: clean-terminal flag changed")
    if receipt.get("claim_admitted") is not (classification["clean_terminal"] is True):
        raise VerificationError(f"{task['id']}: claim-admission flag changed")
    return {
        "id": task["id"], "ordinal": task["ordinal"], "started_at": receipt["started_at"],
        "finished_at": receipt["finished_at"], "outcome": classification["outcome"],
        "clean_f4_terminal": classification["clean_terminal"],
        "classification": classification.get("classification"), "terminal": classification.get("terminal"),
        "service": classification.get("service"),
        "attempt_receipt_sha256": sha256_bytes((attempt / "receipt.json").read_bytes()),
    }


def summary_body(plan: dict, cases: list[dict], artifact_status: str, execution_binding: dict | None) -> dict:
    clean = sum(case["clean_f4_terminal"] is True for case in cases)
    outcomes: dict[str, int] = {}
    for case in cases:
        outcomes[case["outcome"]] = outcomes.get(case["outcome"], 0) + 1
    no_retry = bool(cases) and all(case["ordinal"] == index for index, case in enumerate(cases, 1))
    return {
        "schema": SUMMARY_SCHEMA, "artifact_status": artifact_status,
        "plan_sha256": sha256_bytes(canonical_bytes(plan)),
        "execution_binding_sha256": None if execution_binding is None else execution_binding["binding_sha256"],
        "expected_new_requests": 10, "attempted_new_requests": len(cases), "unattempted_new_requests": 10 - len(cases),
        "one_attempt_no_retry_verified": no_retry,
        "one_attempt_no_retry_status": "verified" if no_retry else "not_applicable",
        "cross_attempt_spacing_verified": len(cases) > 1,
        "clean_new_f4_terminals": clean, "operational_or_invalid_receipts": len(cases) - clean,
        "outcome_counts": dict(sorted(outcomes.items())), "retained_prior_f4_terminals": 5,
        "verified_public_service_f4_terminal_coverage": 5 + clean,
        "maximum_public_service_f4_terminal_coverage": 15,
        "five_n59_tasks_excluded_by_input_cap": True, "licensed_twenty_cell_matrix_executed": False,
        "magma_point_witnesses_validated": 0, "process_scoped_solver_resources_complete": False,
        "service_headers_admitted_as_solver_process_metrics": False, "performance_ranking_admitted": False,
        "full_solver_matrix_gate_passed": False, "cases": cases,
        "claim": (
            "This artifact binds a bounded official-calculator F4 panel. Clean receipts are source-equivalent "
            "basis terminals without Magma-linked point witnesses or solver-process CPU/RSS. It does not execute "
            "the licensed twenty-cell matrix and supports no performance, novelty, end-to-end, or SOTA claim."
        ),
        "claim_boundary": plan["claim_boundary"],
    }


def verify(artifact: Path, execution_binding_override: dict | None = None) -> dict:
    plan, _ = verify_prepared(
        artifact, allow_external_test_artifact=execution_binding_override is not None
    )
    attempts_root = artifact / "attempts"
    run_path = artifact / "run.json"
    if not attempts_root.exists():
        if run_path.exists() or (artifact / "summary.json").exists():
            raise VerificationError("runtime receipt exists without attempts")
        return summary_body(plan, [], "prepared_not_executed", None)
    if attempts_root.is_symlink() or not attempts_root.is_dir() or not run_path.is_file():
        raise VerificationError("attempts/run custody is malformed")
    run = read_json(run_path)
    binding = execution_binding_override or run.get("execution_binding")
    if execution_binding_override is None:
        verify_execution_binding(binding)
    expected_names = [RUNNER.task_stem(task) for task in plan["tasks"]]
    entries = list(attempts_root.iterdir())
    names = {path.name for path in entries}
    if names != set(expected_names[:len(entries)]) or len(entries) > 10:
        raise VerificationError("attempts are not the exact contiguous seed-major prefix")
    cases = [
        verify_attempt(artifact, plan, task, attempts_root / RUNNER.task_stem(task), binding)
        for task in plan["tasks"][:len(entries)]
    ]
    spacing = float(plan["execution_policy"]["inter_request_delay_seconds"])
    for index in range(1, len(cases)):
        if parse_time(cases[index]["started_at"], "attempt started_at") < parse_time(
            cases[index - 1]["finished_at"], "prior finished_at"
        ) + timedelta(seconds=spacing):
            raise VerificationError("attempts violate global sequential spacing")
    expected_index = {
        task["id"]: {
            "path": str(Path("attempts") / RUNNER.task_stem(task) / "receipt.json"),
            "sha256": cases[index]["attempt_receipt_sha256"],
            "outcome": cases[index]["outcome"],
            "clean_terminal": cases[index]["clean_f4_terminal"],
        }
        for index, task in enumerate(plan["tasks"][:len(cases)])
    }
    if (
        run.get("schema") != RUNNER.RUN_SCHEMA or run.get("plan_sha256") != sha256_bytes(canonical_bytes(plan))
        or run.get("protocol_sha256") != plan["protocol"]["sha256"]
        or run.get("task_order") != [task["id"] for task in plan["tasks"]]
        or run.get("execution_policy") != plan["execution_policy"] or run.get("execution_binding") != binding
        or run.get("receipts") != expected_index
    ):
        raise VerificationError("run identity, policy, binding, or receipt index changed")
    parse_time(run.get("started_at"), "run started_at")
    for index, stamp in enumerate(run.get("resumed_at", [])):
        parse_time(stamp, f"run resumed_at[{index}]")
    nonclean = [index for index, case in enumerate(cases) if case["clean_f4_terminal"] is not True]
    status = run.get("status")
    if status == "complete":
        if len(cases) != 10 or nonclean:
            raise VerificationError("complete requires all ten clean receipts")
        artifact_status = "ten_clean_requests_complete"
    elif status == "halted_after_nonclean_receipt":
        if nonclean != [len(cases) - 1]:
            raise VerificationError("halted requires the first and only non-clean receipt last")
        artifact_status = "halted_after_nonclean_receipt"
    elif status == "running":
        if nonclean:
            raise VerificationError("running state cannot pass a historical non-clean receipt")
        artifact_status = "interrupted_after_clean_prefix"
    else:
        raise VerificationError("unexpected run status")
    if status != "running":
        parse_time(run.get("finished_at"), "run finished_at")
    raw_receipts = [read_json(attempts_root / RUNNER.task_stem(task) / "receipt.json") for task in plan["tasks"][:len(cases)]]
    expected_runner_summary = RUNNER.summarize_receipts(plan, raw_receipts, status)
    if run.get("summary") != expected_runner_summary:
        raise VerificationError("run summary changed")
    if status != "running":
        if read_json(artifact / "summary.json") != expected_runner_summary:
            raise VerificationError("standalone run summary changed")
    return summary_body(plan, cases, artifact_status, binding)


def self_test() -> dict:
    plan, inputs = RENDERER.build_plan()
    task = plan["tasks"][0]
    data = inputs[Path(task["named_input"]["path"]).name]
    if task["id"].encode() not in data or task["source_instance_sha256"].encode() not in data:
        raise AssertionError("identity markers are absent from named input")
    if (
        plan["service"]["max_response_bytes"] != CHILD.MAX_RESPONSE_BYTES
        or CHILD.expected_ca_bundle_record()["sha256"] != CHILD.CA_BUNDLE_SHA256
    ):
        raise AssertionError("response byte-limit contract changed")
    return {"self_test": "pass", "expected_requests": 10, "identity_bound_inputs": 10, "response_byte_limit": CHILD.MAX_RESPONSE_BYTES}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--expected", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    try:
        summary = verify(args.artifact.resolve())
    except (VerificationError, RENDERER.RenderError, RUNNER.RunError, CHILD.ChildError) as error:
        parser.error(str(error))
    rendered = canonical_bytes(summary)
    if args.expected is not None and args.expected.read_bytes() != rendered:
        parser.error("recomputed Stage 19 summary differs from expected")
    if args.output is not None:
        RUNNER.atomic_write(args.output.resolve(), rendered)
    print(rendered.decode(), end="")


if __name__ == "__main__":
    main()
