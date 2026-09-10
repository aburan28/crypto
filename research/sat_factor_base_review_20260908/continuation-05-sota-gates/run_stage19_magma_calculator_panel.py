#!/usr/bin/env python3
"""Run the commit-bound Stage 19 calculator panel once, sequentially, without retry."""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import secrets
import stat
import subprocess
import sys
import time
from typing import Any
import urllib.parse


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RENDERER_PATH = HERE / "render_stage19_magma_calculator_panel.py"
CHILD_PATH = HERE / "post_stage19_magma_calculator_request.py"
DEFAULT_ARTIFACT = HERE / "stage-19-magma-calculator-panel-amendment-01-20260910"
METER = REPO / "scripts" / "process_meter.py"
RUN_SCHEMA = "koblitz_magma_calculator_stage19_run.v2"
ATTEMPT_SCHEMA = "koblitz_magma_calculator_stage19_attempt.v2"
RECEIPT_SCHEMA = "koblitz_magma_calculator_stage19_receipt.v2"
EXECUTION_MANIFEST_SCHEMA = "koblitz_magma_calculator_stage19_execution_manifest.v1"


class RunError(RuntimeError):
    """Execution state is not safe for another public-service request."""


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RunError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RENDERER = load_module("stage19_renderer_for_runner", RENDERER_PATH)
CHILD = load_module("stage19_child_for_runner", CHILD_PATH)
STAGE15 = RENDERER.STAGE15
atomic_write = CHILD.atomic_write


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def parse_time(value: Any, label: str) -> datetime:
    if not isinstance(value, str) or not value.endswith("Z"):
        raise RunError(f"{label} is not a UTC timestamp")
    try:
        return datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise RunError(f"{label} is not an ISO timestamp") from error


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise RunError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise RunError(f"expected JSON object in {path}")
    return value


def atomic_json(path: Path, value: Any) -> None:
    atomic_write(path, canonical_bytes(value))


def regular_record(path: Path, relative_to: Path | None = None) -> dict:
    try:
        info = path.lstat()
    except FileNotFoundError as error:
        raise RunError(f"missing required regular file: {path}") from error
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise RunError(f"not a single-link regular file: {path}")
    data = path.read_bytes()
    return {
        "path": str(path.relative_to(relative_to)) if relative_to is not None else str(path),
        "bytes": len(data),
        "sha256": sha256_bytes(data),
    }


def git(*args: str, check: bool = True) -> bytes:
    result = subprocess.run(["git", *args], cwd=REPO, capture_output=True, check=False)
    if check and result.returncode != 0:
        raise RunError(f"git {' '.join(args)} failed: {result.stderr.decode(errors='replace').strip()}")
    return result.stdout


def task_stem(task: dict) -> str:
    return Path(task["named_input"]["path"]).stem


def request_body(input_bytes: bytes, field: str) -> bytes:
    try:
        input_text = input_bytes.decode("ascii")
    except UnicodeDecodeError as error:
        raise RunError("named Magma input is not ASCII") from error
    return urllib.parse.urlencode({field: input_text}).encode("ascii")


def verify_prepared_artifact(
    artifact: Path, allow_external_test_artifact: bool = False
) -> tuple[dict, dict[str, bytes]]:
    if not allow_external_test_artifact and artifact.resolve() != DEFAULT_ARTIFACT.resolve():
        raise RunError("network execution is restricted to the repository Stage 19 artifact")
    expected, rendered_inputs = RENDERER.build_plan()
    plan_path = artifact / "plan.json"
    if plan_path.is_symlink() or not plan_path.is_file():
        raise RunError("prepared plan.json is missing or a symlink")
    plan = read_json(plan_path)
    if plan != expected or plan_path.read_bytes() != canonical_bytes(expected):
        raise RunError("prepared plan differs from deterministic source rendering")
    for task in plan["tasks"]:
        path = artifact / task["named_input"]["path"]
        expected_bytes = rendered_inputs[path.name]
        if regular_record(path)["sha256"] != sha256_bytes(expected_bytes) or path.read_bytes() != expected_bytes:
            raise RunError(f"prepared named input changed for {task['id']}")
    regular_record(artifact / "prepared-summary.json")
    return plan, rendered_inputs


def parse_git_tree_entry(revision: str, relative: str) -> tuple[str, str]:
    line = git("ls-tree", revision, "--", relative).decode().strip()
    if not line or "\t" not in line:
        raise RunError(f"missing Git tree entry {revision}:{relative}")
    metadata, observed = line.split("\t", 1)
    fields = metadata.split()
    if len(fields) != 3 or fields[1] != "blob" or observed != relative:
        raise RunError(f"unexpected Git tree entry {revision}:{relative}")
    return fields[0], fields[2]


def runtime_dirty_path_allowed(relative: str) -> bool:
    root = str(DEFAULT_ARTIFACT.relative_to(REPO))
    return relative in {f"{root}/run.json", f"{root}/summary.json"} or relative.startswith(f"{root}/attempts/")


def expected_execution_bound_paths() -> set[str]:
    fixed = {
        ".github/workflows/koblitz-sota-reproduction.yml",
        "scripts/process_meter.py",
        str((HERE / "STAGE19_RESULTS.md").relative_to(REPO)),
        str((HERE / "stage-19-amendment-01-zero-post-parent-check.json").relative_to(REPO)),
        str((HERE / "stage-19-amendment-01-summary.json").relative_to(REPO)),
        str((HERE / "verify_stage19_amendment01.py").relative_to(REPO)),
        str((HERE / "stage-19-magma-calculator-panel-protocol.json").relative_to(REPO)),
        str((HERE / "render_stage19_magma_calculator_panel.py").relative_to(REPO)),
        str((HERE / "post_stage19_magma_calculator_request.py").relative_to(REPO)),
        str((HERE / "run_stage19_magma_calculator_panel.py").relative_to(REPO)),
        str((HERE / "verify_stage19_magma_calculator_panel.py").relative_to(REPO)),
        str((HERE / "test_stage19_magma_calculator_panel.py").relative_to(REPO)),
        str((HERE / "prepare_stage19_magma_execution_manifest.py").relative_to(REPO)),
        str((DEFAULT_ARTIFACT / "plan.json").relative_to(REPO)),
        str((DEFAULT_ARTIFACT / "prepared-summary.json").relative_to(REPO)),
    }
    inputs = {
        str(path.relative_to(REPO)) for path in (DEFAULT_ARTIFACT / "inputs").glob("*.magma")
    }
    if len(inputs) != 10:
        raise RunError("execution binding requires exactly ten prepared inputs")
    original_artifact = HERE / "stage-19-magma-calculator-panel-20260909"
    original = {
        str(path.relative_to(REPO))
        for path in original_artifact.rglob("*")
        if path.is_file() and not path.is_symlink()
    }
    if len(original) != 22:
        raise RunError("execution binding requires the exact immutable 22-file original artifact")
    return fixed | inputs | original


def checkout_dirty_paths() -> list[str]:
    output = git("status", "--porcelain=v1", "--untracked-files=all").decode()
    paths = []
    for line in output.splitlines():
        if len(line) < 4:
            raise RunError("malformed git status output")
        path = line[3:]
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        if path.startswith('"'):
            raise RunError("quoted git status paths are unsupported in execution preflight")
        paths.append(path)
    return paths


def validate_execution_binding(artifact: Path, resume: bool) -> dict:
    if artifact.resolve() != DEFAULT_ARTIFACT.resolve():
        raise RunError("execution manifest applies only to the repository artifact")
    manifest_path = artifact / "execution-manifest.json"
    manifest_record = regular_record(manifest_path, REPO)
    manifest = read_json(manifest_path)
    if manifest.get("schema") != EXECUTION_MANIFEST_SCHEMA:
        raise RunError("unexpected execution-manifest schema")
    first = manifest.get("first_preexecution_commit")
    first_tree = manifest.get("first_preexecution_tree")
    if not isinstance(first, str) or not re.fullmatch(r"[0-9a-f]{40}", first):
        raise RunError("execution manifest lacks its first pre-execution commit")
    if git("rev-parse", f"{first}^{{tree}}").decode().strip() != first_tree:
        raise RunError("first pre-execution tree identity changed")
    head = git("rev-parse", "HEAD").decode().strip()
    head_tree = git("rev-parse", "HEAD^{tree}").decode().strip()
    ancestry = subprocess.run(["git", "merge-base", "--is-ancestor", first, head], cwd=REPO, check=False)
    if ancestry.returncode != 0:
        raise RunError("execution commit does not descend from the first pre-execution commit")
    delta_lines = [line for line in git("diff", "--name-status", first, head).decode().splitlines() if line]
    expected_delta = f"A\t{manifest_record['path']}"
    if delta_lines != [expected_delta]:
        raise RunError("second pre-execution commit delta is not exactly execution-manifest.json")
    manifest_mode, manifest_blob = parse_git_tree_entry(head, manifest_record["path"])
    if git("show", f"{head}:{manifest_record['path']}") != manifest_path.read_bytes():
        raise RunError("execution manifest is not the committed HEAD blob")
    records = manifest.get("relevant_blobs")
    if not isinstance(records, list):
        raise RunError("execution manifest relevant-blob inventory is incomplete")
    if {record.get("path") for record in records if isinstance(record, dict)} != expected_execution_bound_paths():
        raise RunError("execution manifest does not bind the exact required path set")
    seen = set()
    for record in records:
        if not isinstance(record, dict) or set(record) != {"path", "git_mode", "git_blob_oid", "bytes", "sha256"}:
            raise RunError("malformed execution-manifest blob record")
        relative = record["path"]
        if relative in seen:
            raise RunError("duplicate execution-manifest blob path")
        seen.add(relative)
        base_mode, base_blob = parse_git_tree_entry(first, relative)
        head_mode, head_blob = parse_git_tree_entry(head, relative)
        path = REPO / relative
        current = regular_record(path)
        if (
            base_mode != record["git_mode"] or head_mode != record["git_mode"]
            or base_blob != record["git_blob_oid"] or head_blob != record["git_blob_oid"]
            or current["bytes"] != record["bytes"] or current["sha256"] != record["sha256"]
            or git("show", f"{first}:{relative}") != path.read_bytes()
        ):
            raise RunError(f"execution-bound blob changed: {relative}")
    dirty = checkout_dirty_paths()
    disallowed = [path for path in dirty if not (resume and runtime_dirty_path_allowed(path))]
    if disallowed:
        raise RunError(f"checkout is dirty outside allowed runtime artifacts: {disallowed}")
    binding = {
        "schema": "koblitz_magma_calculator_stage19_execution_binding.v1",
        "first_preexecution_commit": first,
        "first_preexecution_tree": first_tree,
        "execution_commit": head,
        "execution_tree": head_tree,
        "tree_delta": [{"status": "A", "path": manifest_record["path"]}],
        "execution_manifest": {**manifest_record, "git_mode": manifest_mode, "git_blob_oid": manifest_blob},
        "relevant_blob_count": len(records),
        "relevant_path_list_sha256": manifest["relevant_path_list_sha256"],
        "checkout_clean_outside_runtime_artifacts": True,
    }
    binding["binding_sha256"] = sha256_bytes(canonical_bytes(binding))
    return binding


def expected_request(plan: dict, task: dict, input_bytes: bytes, binding: dict) -> dict:
    service = plan["service"]
    body = request_body(input_bytes, service["form_field"])
    return {
        "endpoint": service["endpoint"], "method": "POST", "form_field": service["form_field"],
        "content_type": "application/x-www-form-urlencoded", "user_agent": service["user_agent"],
        "input": task["named_input"], "body_bytes": len(body), "body_sha256": sha256_bytes(body),
        "hard_parent_watchdog_seconds": service["client_timeout_seconds"],
        "response_byte_limit": service["max_response_bytes"],
        "execution_commit": binding["execution_commit"], "execution_tree": binding["execution_tree"],
        "execution_binding_sha256": binding["binding_sha256"],
        "execution_manifest_sha256": binding["execution_manifest"]["sha256"],
    }


def parse_identity_bound_terminal(output: str, task: dict) -> tuple[dict | None, str | None]:
    prefix = (
        f"KOBLITZ_MAGMA_TASK_ID={task['id']}\n"
        f"KOBLITZ_MAGMA_SOURCE_SHA256={task['source_instance_sha256']}\n"
    )
    if not output.startswith(prefix):
        return None, "response task/source identity markers are absent or mismatched"
    remainder = output[len(prefix):]
    if "KOBLITZ_MAGMA_TASK_ID=" in remainder or "KOBLITZ_MAGMA_SOURCE_SHA256=" in remainder:
        return None, "response contains duplicate identity markers"
    try:
        STAGE15.require_clean_f4_output(remainder, int(task["seed"]))
    except STAGE15.VerificationError as error:
        return None, str(error)
    terminal = STAGE15.parse_magma_terminal(remainder)
    return (terminal, None) if terminal is not None else (None, "strict seven-marker F4 parser rejected the response")


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
    terminal, diagnostic = parse_identity_bound_terminal(response["output"], task)
    if terminal is None:
        return {"outcome": "identity_or_terminal_mismatch", "clean_terminal": False, "diagnostic": diagnostic, "service": response["service"]}
    if terminal["terminal_status"] != "sat":
        return {"outcome": "planted_unsat_contradiction", "clean_terminal": False, "terminal": terminal, "service": response["service"]}
    return {
        "outcome": "clean_f4_terminal", "clean_terminal": True,
        "classification": "sat_basis_certificate_unverified_model", "response_identity_verified": True,
        "terminal": terminal, "service": response["service"],
    }


def validate_transport_envelope(attempt: Path, task: dict, request: dict) -> tuple[dict, Path | None]:
    envelope_path = attempt / "transport-envelope.json"
    envelope = read_json(envelope_path)
    if envelope.get("schema") != CHILD.ENVELOPE_SCHEMA:
        raise RunError(f"{task['id']}: transport envelope schema changed")
    expected_child_request = {
        "endpoint": request["endpoint"], "method": "POST", "form_field": request["form_field"],
        "user_agent": request["user_agent"], "input_bytes": task["named_input"]["bytes"],
        "input_sha256": task["named_input"]["sha256"], "body_bytes": request["body_bytes"],
        "body_sha256": request["body_sha256"], "socket_timeout_seconds": CHILD.SOCKET_TIMEOUT_SECONDS,
    }
    if envelope.get("request") != expected_child_request or envelope.get("response_byte_limit") != CHILD.MAX_RESPONSE_BYTES:
        raise RunError(f"{task['id']}: transport request binding changed")
    http = envelope.get("http")
    if not isinstance(http, dict) or set(http) != {
        "status", "final_url", "selected_headers", "transport_error"
    }:
        raise RunError(f"{task['id']}: malformed transport HTTP envelope")
    status = http["status"]
    if status is not None and (not isinstance(status, int) or isinstance(status, bool)):
        raise RunError(f"{task['id']}: malformed HTTP status")
    headers = http["selected_headers"]
    if (
        not isinstance(headers, dict)
        or set(headers) - set(CHILD.SAFE_RESPONSE_HEADERS)
        or any(not isinstance(value, str) for value in headers.values())
    ):
        raise RunError(f"{task['id']}: unsafe retained HTTP headers")
    started = parse_time(envelope.get("started_at"), f"{task['id']} transport started_at")
    if parse_time(envelope.get("finished_at"), f"{task['id']} transport finished_at") < started:
        raise RunError(f"{task['id']}: negative transport chronology")
    response = envelope.get("response")
    if response is None:
        return envelope, None
    if not isinstance(response, dict) or set(response) != {"path", "bytes", "sha256", "complete", "body_limit_exceeded", "response_byte_limit"}:
        raise RunError(f"{task['id']}: malformed envelope response record")
    path = attempt / response["path"]
    observed = regular_record(path, attempt)
    if observed != {key: response[key] for key in ("path", "bytes", "sha256")} or response["response_byte_limit"] != CHILD.MAX_RESPONSE_BYTES:
        raise RunError(f"{task['id']}: retained response bytes or limit changed")
    expected_name = (
        "response.partial"
        if response["body_limit_exceeded"] is True
        else "response.xml" if status == 200 else "response.body"
    )
    if (
        response["path"] != expected_name
        or response["complete"] is not (response["body_limit_exceeded"] is False)
        or response["bytes"] > CHILD.MAX_RESPONSE_BYTES
        or (response["body_limit_exceeded"] is True and response["bytes"] != CHILD.MAX_RESPONSE_BYTES)
    ):
        raise RunError(f"{task['id']}: response completeness/path contract changed")
    return envelope, path


def transport_process_record(attempt: Path, expected_command: list[str], timeout: float) -> dict:
    metrics_path = attempt / "transport-metrics.json"
    metrics = read_json(metrics_path)
    if metrics.get("command") != expected_command or metrics.get("watchdog_seconds") != timeout:
        raise RunError("transport process command or watchdog changed")
    if type(metrics.get("timed_out")) is not bool or type(metrics.get("orphan_group_terminated")) is not bool:
        raise RunError("transport process lacks terminal booleans")
    if not isinstance(metrics.get("returncode"), int) or isinstance(metrics.get("returncode"), bool):
        raise RunError("transport process return code is invalid")
    resources = metrics.get("metrics")
    required = {
        "wall_seconds", "user_seconds", "system_seconds", "total_core_seconds",
        "single_core_seconds", "peak_rss_bytes", "meter",
    }
    if not isinstance(resources, dict) or set(resources) != required:
        raise RunError("transport process resource record is malformed")
    for key in ("wall_seconds", "user_seconds", "system_seconds", "total_core_seconds", "single_core_seconds"):
        value = resources[key]
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value < 0:
            raise RunError(f"invalid transport process metric {key}")
    if not isinstance(resources["peak_rss_bytes"], int) or isinstance(resources["peak_rss_bytes"], bool) or resources["peak_rss_bytes"] < 0:
        raise RunError("invalid transport peak RSS")
    if resources["meter"] != "fresh-process getrusage(RUSAGE_CHILDREN)":
        raise RunError("unexpected transport process meter")
    if not math.isclose(resources["total_core_seconds"], resources["user_seconds"] + resources["system_seconds"], rel_tol=1e-9, abs_tol=1e-9):
        raise RunError("transport total core-seconds do not add up")
    return {
        "metrics": metrics,
        "stdout": regular_record(attempt / "transport.stdout", attempt),
        "stderr": regular_record(attempt / "transport.stderr", attempt),
        "metrics_file": regular_record(metrics_path, attempt),
    }


def run_metered_transport(
    attempt: Path, plan: dict, task: dict, launch_nonce: str
) -> dict:
    command = [
        str(Path(sys.executable).resolve()), str(CHILD_PATH.resolve()),
        "--attempt", str(attempt.resolve()), "--execute-child",
    ]
    timeout = float(plan["service"]["client_timeout_seconds"])
    meter_command = [
        str(Path(sys.executable).resolve()), str(METER.resolve()), "--cwd", str(REPO.resolve()),
        "--timeout", str(timeout), "--stdout", str((attempt / "transport.stdout").resolve()),
        "--stderr", str((attempt / "transport.stderr").resolve()), "--metrics", str((attempt / "transport-metrics.json").resolve()),
        "--", *command,
    ]
    child_environment = os.environ.copy()
    child_environment[CHILD.LAUNCH_NONCE_ENV] = launch_nonce
    launcher = subprocess.run(
        meter_command, cwd=REPO, env=child_environment, check=False
    )
    if launcher.returncode != 0 or not (attempt / "transport-metrics.json").is_file():
        raise RunError(f"transport process meter failed with return code {launcher.returncode}")
    return transport_process_record(attempt, command, timeout)


def build_final_receipt(attempt: Path, plan: dict, task: dict, start: dict, process: dict) -> dict:
    metrics = process["metrics"]
    classification: dict = {"outcome": "transport_child_error", "clean_terminal": False}
    envelope_path = attempt / "transport-envelope.json"
    envelope_record = regular_record(envelope_path, attempt) if envelope_path.is_file() else None
    response_paths = [
        path for path in (attempt / "response.xml", attempt / "response.body", attempt / "response.partial")
        if path.is_file()
    ]
    if len(response_paths) > 1:
        raise RunError(f"{task['id']}: multiple retained response bodies")
    response_record = regular_record(response_paths[0], attempt) if response_paths else None
    authorization_paths = [
        path for path in (
            attempt / "launch-authorization.json",
            attempt / "launch-authorization.consumed.json",
        ) if path.is_file()
    ]
    if len(authorization_paths) != 1:
        raise RunError(f"{task['id']}: launch authorization custody is missing or duplicated")
    authorization_record = regular_record(authorization_paths[0], attempt)
    recovery_boundary = None
    if metrics["timed_out"]:
        classification = {"outcome": "transport_timeout", "clean_terminal": False}
        recovery_boundary = "transport_timeout_with_retained_files"
    elif metrics["orphan_group_terminated"]:
        classification = {"outcome": "transport_orphan_terminated", "clean_terminal": False}
        recovery_boundary = "transport_orphan_with_retained_files"
    elif metrics["returncode"] != 0:
        classification = {"outcome": "transport_child_error", "clean_terminal": False}
        recovery_boundary = "transport_child_error_with_retained_files"
    elif envelope_path.is_file():
        envelope, response_path = validate_transport_envelope(attempt, task, start["request"])
        if response_path is not None and response_record != regular_record(response_path, attempt):
            raise RunError(f"{task['id']}: envelope and retained response disagree")
        if isinstance(envelope.get("response"), dict) and envelope["response"].get("body_limit_exceeded") is True:
            classification = {"outcome": "response_body_limit_exceeded", "clean_terminal": False}
        elif response_path is None:
            classification = {"outcome": "transport_error", "clean_terminal": False}
        else:
            classification = classify_response(response_path, envelope["http"]["status"], envelope["http"]["final_url"], plan, task)
    else:
        classification = {"outcome": "transport_envelope_missing", "clean_terminal": False}
    receipt = {
        "schema": RECEIPT_SCHEMA, "task_id": task["id"], "ordinal": task["ordinal"], "attempt_ordinal": 1,
        "started_at": start["started_at"], "finished_at": now(), **classification,
        "request": start["request"], "execution_commit": start["execution_commit"],
        "execution_tree": start["execution_tree"], "execution_binding_sha256": start["execution_binding_sha256"],
        "launch_authorization": authorization_record,
        "transport_process": process, "transport_envelope": envelope_record, "response": response_record,
        "retained_files": recovery_file_inventory(attempt),
        "retry_permitted": False, "claim_admitted": classification["clean_terminal"] is True,
        "recovery_boundary": recovery_boundary,
    }
    atomic_json(attempt / "receipt.json", receipt)
    return receipt


def recovery_file_inventory(attempt: Path) -> list[dict]:
    records = []
    for path in sorted(attempt.iterdir()):
        if path.name == "receipt.json":
            continue
        if path.is_symlink() or not path.is_file():
            raise RunError(f"ambiguous recovery entry is not a regular file: {path}")
        records.append(regular_record(path, attempt))
    return records


def recover_interrupted(
    artifact: Path, attempt: Path, plan: dict, task: dict, binding: dict
) -> dict:
    start_path = attempt / "attempt-start.json"
    if start_path.is_file():
        start = read_json(start_path)
        started_at = start.get("started_at")
        request = start.get("request")
    else:
        started_at = now()
        input_bytes = (artifact / task["named_input"]["path"]).read_bytes()
        request = expected_request(plan, task, input_bytes, binding)
        start = {
            "schema": ATTEMPT_SCHEMA, "task_id": task["id"], "ordinal": task["ordinal"], "attempt_ordinal": 1,
            "started_at": started_at, "request": request, "execution_commit": binding["execution_commit"],
            "execution_tree": binding["execution_tree"], "execution_binding_sha256": binding["binding_sha256"],
            "launch_nonce_sha256": None,
            "retry_permitted": False, "recovered_empty_attempt_directory": True,
        }
        atomic_json(start_path, start)
    response_paths = [path for path in (attempt / "response.xml", attempt / "response.body", attempt / "response.partial") if path.exists()]
    temporary_paths = [path for path in attempt.iterdir() if path.name.startswith(".") and ".tmp-" in path.name]
    envelope_exists = (attempt / "transport-envelope.json").is_file()
    metrics_exists = (attempt / "transport-metrics.json").is_file()
    if temporary_paths:
        outcome = "interrupted_during_atomic_write"
    elif response_paths and envelope_exists:
        outcome = "interrupted_with_retained_response"
    elif response_paths:
        outcome = "interrupted_with_unbound_response"
    elif envelope_exists:
        outcome = "interrupted_with_incomplete_transport_envelope"
    elif metrics_exists:
        outcome = "interrupted_after_transport_process"
    elif start.get("recovered_empty_attempt_directory") is True:
        outcome = "interrupted_before_attempt_start"
    else:
        outcome = "interrupted_before_transport_receipt"
    authorization_paths = [
        path for path in (
            attempt / "launch-authorization.json",
            attempt / "launch-authorization.consumed.json",
        ) if path.is_file()
    ]
    authorization_record = (
        regular_record(authorization_paths[0], attempt) if len(authorization_paths) == 1 else None
    )
    envelope_record = regular_record(attempt / "transport-envelope.json", attempt) if envelope_exists else None
    response_record = regular_record(response_paths[0], attempt) if len(response_paths) == 1 else None
    retained = recovery_file_inventory(attempt)
    receipt = {
        "schema": RECEIPT_SCHEMA, "task_id": task["id"], "ordinal": task["ordinal"], "attempt_ordinal": 1,
        "started_at": started_at, "finished_at": now(), "outcome": outcome, "clean_terminal": False,
        "request": request, "execution_commit": binding["execution_commit"], "execution_tree": binding["execution_tree"],
        "execution_binding_sha256": binding["binding_sha256"],
        "launch_authorization": authorization_record, "transport_process": None,
        "transport_envelope": envelope_record, "response": response_record,
        "retained_files": retained, "recovery_files": retained,
        "retry_permitted": False, "claim_admitted": False, "recovery_boundary": outcome,
    }
    atomic_json(attempt / "receipt.json", receipt)
    return receipt


def audit_existing_receipt(
    artifact: Path, attempt: Path, plan: dict, task: dict, binding: dict
) -> dict:
    receipt = read_json(attempt / "receipt.json")
    allowed_files = {
        "attempt-start.json", "request-body.bin", "launch-authorization.json",
        "launch-authorization.consumed.json", "transport.stdout", "transport.stderr",
        "transport-metrics.json", "transport-envelope.json", "response.xml", "response.body",
        "response.partial", "receipt.json",
    }
    unexpected_files = [path for path in attempt.iterdir() if path.name not in allowed_files]
    if any(not (path.name.startswith(".") and ".tmp-" in path.name) for path in unexpected_files):
        raise RunError(f"{task['id']}: attempt contains an unexpected file")
    start = read_json(attempt / "attempt-start.json")
    input_bytes = (artifact / task["named_input"]["path"]).read_bytes()
    request = expected_request(plan, task, input_bytes, binding)
    required = {
        "schema": RECEIPT_SCHEMA, "task_id": task["id"], "ordinal": task["ordinal"], "attempt_ordinal": 1,
        "request": request, "execution_commit": binding["execution_commit"], "execution_tree": binding["execution_tree"],
        "execution_binding_sha256": binding["binding_sha256"], "retry_permitted": False,
    }
    if any(receipt.get(key) != value for key, value in required.items()):
        raise RunError(f"{task['id']}: existing receipt identity changed")
    if start.get("request") != request or start.get("task_id") != task["id"]:
        raise RunError(f"{task['id']}: attempt start identity changed")
    nonce_hash = start.get("launch_nonce_sha256")
    if start.get("recovered_empty_attempt_directory") is True:
        if nonce_hash is not None:
            raise RunError(f"{task['id']}: recovered empty attempt invented a nonce")
    elif not isinstance(nonce_hash, str) or re.fullmatch(r"[0-9a-f]{64}", nonce_hash) is None:
        raise RunError(f"{task['id']}: launch nonce hash is malformed")
    parse_time(receipt.get("started_at"), f"{task['id']} started_at")
    parse_time(receipt.get("finished_at"), f"{task['id']} finished_at")
    if receipt.get("claim_admitted") is not (receipt.get("clean_terminal") is True):
        raise RunError(f"{task['id']}: claim admission disagrees with terminal class")
    retained = recovery_file_inventory(attempt)
    if receipt.get("retained_files") != retained:
        raise RunError(f"{task['id']}: retained file inventory changed")
    authorization_paths = [
        path for path in (
            attempt / "launch-authorization.json",
            attempt / "launch-authorization.consumed.json",
        ) if path.is_file()
    ]
    observed_authorization = (
        regular_record(authorization_paths[0], attempt) if len(authorization_paths) == 1 else None
    )
    if receipt.get("launch_authorization") != observed_authorization:
        raise RunError(f"{task['id']}: launch authorization receipt changed")
    if len(authorization_paths) == 1:
        expected_authorization = CHILD.expected_launch_authorization(
            attempt, plan, task, start, artifact=artifact
        )
        if authorization_paths[0].read_bytes() != CHILD.canonical_bytes(expected_authorization):
            raise RunError(f"{task['id']}: launch authorization content changed")
    if "recovery_files" in receipt:
        if receipt.get("clean_terminal") is not False or receipt.get("claim_admitted") is not False:
            raise RunError(f"{task['id']}: recovery receipt claims a terminal")
        expected_recovery = retained
        if receipt.get("recovery_files") != expected_recovery:
            raise RunError(f"{task['id']}: recovery file inventory changed")
    else:
        process = receipt.get("transport_process")
        if not isinstance(process, dict) or not isinstance(process.get("metrics"), dict):
            raise RunError(f"{task['id']}: normal receipt lacks transport metrics")
        expected_command = [
            str(Path(sys.executable).resolve()), str(CHILD_PATH.resolve()),
            "--attempt", str(attempt.resolve()), "--execute-child",
        ]
        observed_process = transport_process_record(
            attempt, expected_command, float(plan["service"]["client_timeout_seconds"])
        )
        if process != observed_process:
            raise RunError(f"{task['id']}: transport process receipt changed")
        metrics = process["metrics"]
        envelope_record = (
            regular_record(attempt / "transport-envelope.json", attempt)
            if (attempt / "transport-envelope.json").is_file() else None
        )
        if receipt.get("transport_envelope") != envelope_record:
            raise RunError(f"{task['id']}: retained envelope was not bound")
        response_paths = [
            path for path in (
                attempt / "response.xml", attempt / "response.body", attempt / "response.partial"
            ) if path.is_file()
        ]
        if len(response_paths) > 1:
            raise RunError(f"{task['id']}: multiple retained response bodies")
        response_record = regular_record(response_paths[0], attempt) if response_paths else None
        if receipt.get("response") != response_record:
            raise RunError(f"{task['id']}: retained response was not bound")
        if metrics["timed_out"]:
            classification = {"outcome": "transport_timeout", "clean_terminal": False}
            expected_boundary = "transport_timeout_with_retained_files"
        elif metrics["orphan_group_terminated"]:
            classification = {"outcome": "transport_orphan_terminated", "clean_terminal": False}
            expected_boundary = "transport_orphan_with_retained_files"
        elif metrics["returncode"] != 0:
            classification = {"outcome": "transport_child_error", "clean_terminal": False}
            expected_boundary = "transport_child_error_with_retained_files"
        elif (attempt / "transport-envelope.json").is_file():
            expected_boundary = None
            envelope, response_path = validate_transport_envelope(
                attempt, task, request
            )
            if receipt.get("transport_envelope") != regular_record(
                attempt / "transport-envelope.json", attempt
            ):
                raise RunError(f"{task['id']}: envelope receipt changed")
            if isinstance(envelope.get("response"), dict) and envelope["response"].get(
                "body_limit_exceeded"
            ) is True:
                classification = {
                    "outcome": "response_body_limit_exceeded",
                    "clean_terminal": False,
                }
            elif response_path is None:
                classification = {"outcome": "transport_error", "clean_terminal": False}
            else:
                if receipt.get("response") != regular_record(response_path, attempt):
                    raise RunError(f"{task['id']}: response receipt changed")
                classification = classify_response(
                    response_path,
                    envelope["http"]["status"],
                    envelope["http"]["final_url"],
                    plan,
                    task,
                )
        else:
            expected_boundary = None
            classification = {"outcome": "transport_envelope_missing", "clean_terminal": False}
        if receipt.get("recovery_boundary") != expected_boundary:
            raise RunError(f"{task['id']}: transport recovery boundary changed")
        for key, value in classification.items():
            if receipt.get(key) != value:
                raise RunError(f"{task['id']}: retained classification changed at {key}")
    regular_record(attempt / "attempt-start.json")
    regular_record(attempt / "receipt.json")
    return receipt


def create_attempt(
    artifact: Path, plan: dict, task: dict, binding: dict
) -> tuple[Path, dict, str]:
    attempts_root = artifact / "attempts"
    attempts_root.mkdir(parents=True, exist_ok=True)
    if attempts_root.is_symlink() or not attempts_root.is_dir():
        raise RunError("attempts root is not a real directory")
    attempt = attempts_root / task_stem(task)
    attempt.mkdir()
    input_path = artifact / task["named_input"]["path"]
    input_bytes = input_path.read_bytes()
    request = expected_request(plan, task, input_bytes, binding)
    launch_nonce = secrets.token_hex(32)
    start = {
        "schema": ATTEMPT_SCHEMA, "task_id": task["id"], "ordinal": task["ordinal"], "attempt_ordinal": 1,
        "started_at": now(), "request": request, "execution_commit": binding["execution_commit"],
        "execution_tree": binding["execution_tree"], "execution_binding_sha256": binding["binding_sha256"],
        "launch_nonce_sha256": sha256_bytes(launch_nonce.encode()),
        "retry_permitted": False, "recovered_empty_attempt_directory": False,
    }
    atomic_json(attempt / "attempt-start.json", start)
    atomic_write(attempt / "request-body.bin", request_body(input_bytes, plan["service"]["form_field"]))
    authorization = CHILD.expected_launch_authorization(
        attempt, plan, task, start, artifact=artifact
    )
    atomic_json(attempt / "launch-authorization.json", authorization)
    return attempt, start, launch_nonce


def wait_for_global_spacing(previous_finished_at: str | None, seconds: float) -> None:
    if previous_finished_at is None:
        return
    deadline = parse_time(previous_finished_at, "previous receipt finished_at") + timedelta(seconds=seconds)
    remaining = (deadline - datetime.now(timezone.utc)).total_seconds()
    if remaining > 0:
        monotonic_deadline = time.monotonic() + remaining
        while time.monotonic() < monotonic_deadline:
            time.sleep(min(monotonic_deadline - time.monotonic(), 0.25))


def summarize_receipts(plan: dict, receipts: list[dict], status: str) -> dict:
    clean = sum(receipt.get("clean_terminal") is True for receipt in receipts)
    outcomes: dict[str, int] = {}
    for receipt in receipts:
        outcome = str(receipt.get("outcome", "missing"))
        outcomes[outcome] = outcomes.get(outcome, 0) + 1
    no_retry = bool(receipts) and all(receipt.get("attempt_ordinal") == 1 for receipt in receipts)
    return {
        "status": status, "expected_tasks": 10, "attempted_tasks": len(receipts), "unattempted_tasks": 10 - len(receipts),
        "clean_f4_terminals": clean, "operational_or_invalid_receipts": len(receipts) - clean,
        "outcome_counts": dict(sorted(outcomes.items())), "all_tasks_attempted": len(receipts) == 10,
        "one_attempt_no_retry_verified": no_retry, "one_attempt_no_retry_status": "verified" if no_retry else "not_applicable",
        "public_service_f4_terminal_coverage_after_run": 5 + clean, "maximum_public_service_f4_terminal_coverage": 15,
        "complete_requires_all_ten_clean": len(receipts) == 10 and clean == 10,
        "halted_state_terminal": status == "halted_after_nonclean_receipt",
        "licensed_twenty_cell_matrix_executed": False, "magma_point_witnesses_validated": 0,
        "process_scoped_solver_resources_complete": False, "performance_ranking_admitted": False,
        "claim_boundary": plan["claim_boundary"],
    }


def attempt_prefix(artifact: Path, plan: dict) -> list[tuple[dict, Path]]:
    root = artifact / "attempts"
    if not root.exists():
        return []
    if root.is_symlink() or not root.is_dir():
        raise RunError("attempts root is not a real directory")
    entries = list(root.iterdir())
    if any(path.is_symlink() or not path.is_dir() for path in entries):
        raise RunError("attempts root contains a non-directory or symlink")
    expected = [task_stem(task) for task in plan["tasks"]]
    names = {path.name for path in entries}
    if names != set(expected[:len(entries)]):
        raise RunError("existing attempts are not the exact contiguous seed-major prefix")
    return [(task, root / task_stem(task)) for task in plan["tasks"][:len(entries)]]


def receipt_index_entry(artifact: Path, attempt: Path, receipt: dict) -> dict:
    return {
        "path": str((attempt / "receipt.json").relative_to(artifact)),
        "sha256": sha256_bytes((attempt / "receipt.json").read_bytes()),
        "outcome": receipt["outcome"], "clean_terminal": receipt["clean_terminal"],
    }


def run(artifact: Path, resume: bool, execution_binding_override: dict | None = None) -> dict:
    plan, _ = verify_prepared_artifact(
        artifact, allow_external_test_artifact=execution_binding_override is not None
    )
    binding = execution_binding_override or validate_execution_binding(artifact, resume)
    run_path = artifact / "run.json"
    if run_path.exists() and not resume:
        raise RunError("run.json exists; explicit resume cannot retry or pass a terminal state")
    if not run_path.exists() and resume:
        raise RunError("cannot resume before a run has started")
    plan_hash = sha256_bytes(canonical_bytes(plan))
    if run_path.exists():
        state = read_json(run_path)
        if state.get("status") in {"complete", "halted_after_nonclean_receipt"}:
            raise RunError(f"run state {state['status']} is terminal")
        if (
            state.get("schema") != RUN_SCHEMA or state.get("plan_sha256") != plan_hash
            or state.get("task_order") != [task["id"] for task in plan["tasks"]]
            or state.get("execution_binding") != binding or not isinstance(state.get("receipts"), dict)
        ):
            raise RunError("saved run identity differs from the frozen execution binding")
        state.setdefault("resumed_at", []).append(now())
    else:
        if (artifact / "attempts").exists() or (artifact / "summary.json").exists():
            raise RunError("runtime artifacts exist without a run receipt")
        state = {
            "schema": RUN_SCHEMA, "plan_sha256": plan_hash, "protocol_sha256": plan["protocol"]["sha256"],
            "task_order": [task["id"] for task in plan["tasks"]], "execution_policy": plan["execution_policy"],
            "execution_binding": binding, "started_at": now(), "status": "running", "receipts": {},
        }
        atomic_json(run_path, state)

    receipts = []
    prefix = attempt_prefix(artifact, plan)
    indexed_ids = set(state["receipts"])
    prefix_ids = {task["id"] for task, _ in prefix}
    if indexed_ids - prefix_ids:
        raise RunError("run index names a task without an attempt directory")
    for task, attempt in prefix:
        receipt = (
            audit_existing_receipt(artifact, attempt, plan, task, binding)
            if (attempt / "receipt.json").is_file()
            else recover_interrupted(artifact, attempt, plan, task, binding)
        )
        entry = receipt_index_entry(artifact, attempt, receipt)
        indexed = state["receipts"].get(task["id"])
        if indexed is not None and indexed != entry:
            raise RunError(f"{task['id']}: run index differs from retained receipt")
        state["receipts"][task["id"]] = entry
        receipts.append(receipt)
    for index, receipt in enumerate(receipts):
        if index and parse_time(receipt["started_at"], "attempt started_at") < (
            parse_time(receipts[index - 1]["finished_at"], "prior finished_at")
            + timedelta(seconds=float(plan["execution_policy"]["inter_request_delay_seconds"]))
        ):
            raise RunError("retained attempts violate global sequential spacing")
    nonclean = [index for index, receipt in enumerate(receipts) if receipt["clean_terminal"] is not True]
    if nonclean:
        if nonclean != [len(receipts) - 1]:
            raise RunError("historical non-clean receipt is not the first and terminal non-clean result")
        state["status"] = "halted_after_nonclean_receipt"
        state["finished_at"] = now()
        state["summary"] = summarize_receipts(plan, receipts, state["status"])
        atomic_json(run_path, state)
        atomic_json(artifact / "summary.json", state["summary"])
        return state

    next_index = len(receipts)
    previous_finished = receipts[-1]["finished_at"] if receipts else None
    while next_index < len(plan["tasks"]):
        wait_for_global_spacing(previous_finished, float(plan["execution_policy"]["inter_request_delay_seconds"]))
        task = plan["tasks"][next_index]
        attempt, start, launch_nonce = create_attempt(artifact, plan, task, binding)
        process = run_metered_transport(attempt, plan, task, launch_nonce)
        receipt = build_final_receipt(attempt, plan, task, start, process)
        receipts.append(receipt)
        state["receipts"][task["id"]] = receipt_index_entry(artifact, attempt, receipt)
        previous_finished = receipt["finished_at"]
        state["summary"] = summarize_receipts(plan, receipts, "running")
        atomic_json(run_path, state)
        if receipt["clean_terminal"] is not True:
            state["status"] = "halted_after_nonclean_receipt"
            state["finished_at"] = now()
            state["summary"] = summarize_receipts(plan, receipts, state["status"])
            atomic_json(run_path, state)
            atomic_json(artifact / "summary.json", state["summary"])
            return state
        next_index += 1
    if len(receipts) != 10 or any(receipt["clean_terminal"] is not True for receipt in receipts):
        raise RunError("complete state requires ten clean receipts")
    state["status"] = "complete"
    state["finished_at"] = now()
    state["summary"] = summarize_receipts(plan, receipts, "complete")
    atomic_json(run_path, state)
    atomic_json(artifact / "summary.json", state["summary"])
    return state


def self_test() -> dict:
    plan, inputs = RENDERER.build_plan()
    first = plan["tasks"][0]
    data = inputs[Path(first["named_input"]["path"]).name]
    body = request_body(data, plan["service"]["form_field"])
    if urllib.parse.parse_qs(body.decode("ascii"), strict_parsing=True) != {"input": [data.decode("ascii")]}:
        raise AssertionError("form body does not decode to the exact named input")
    if plan["service"]["client_timeout_seconds"] != 75 or CHILD.MAX_RESPONSE_BYTES != 1_048_576:
        raise AssertionError("transport bounds changed")
    return {
        "self_test": "pass", "attempts_per_task": 1, "retries_per_task": 0,
        "hard_parent_watchdog_seconds": 75, "response_byte_limit": CHILD.MAX_RESPONSE_BYTES,
        "first_request_body_sha256": sha256_bytes(body),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    if not args.execute:
        parser.error("network submission requires the explicit --execute flag")
    try:
        result = run(args.artifact.resolve(), args.resume)
    except (RunError, RENDERER.RenderError, CHILD.ChildError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
