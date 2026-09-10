#!/usr/bin/env python3
"""One bounded Magma Calculator POST, intended to run under process_meter.py."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import stat
import subprocess
import tempfile
from typing import Any
import urllib.error
import urllib.parse
import urllib.request


ENDPOINT = "https://magma.maths.usyd.edu.au/xml/calculator.xml"
FORM_FIELD = "input"
USER_AGENT = "aburan28-crypto-koblitz-stage19/1"
SOCKET_TIMEOUT_SECONDS = 70.0
MAX_RESPONSE_BYTES = 1_048_576
ENVELOPE_SCHEMA = "koblitz_magma_calculator_stage19_transport_envelope.v1"
SAFE_RESPONSE_HEADERS = ("content-type", "content-length", "date", "server")
LAUNCH_NONCE_ENV = "KOBLITZ_STAGE19_LAUNCH_NONCE"
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_ARTIFACT = HERE / "stage-19-magma-calculator-panel-amendment-01-20260910"
PLAN = DEFAULT_ARTIFACT / "plan.json"
EXECUTION_MANIFEST = DEFAULT_ARTIFACT / "execution-manifest.json"
METER = REPO / "scripts" / "process_meter.py"
ATTEMPT_START_SCHEMA = "koblitz_magma_calculator_stage19_attempt.v2"
LAUNCH_AUTH_SCHEMA = "koblitz_magma_calculator_stage19_launch_authorization.v1"


class ChildError(RuntimeError):
    """The one-request child could not preserve its bounded transport result."""


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ChildError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ChildError(f"expected JSON object in {path}")
    return value


def regular_record(path: Path, relative_to: Path | None = None) -> dict:
    ensure_regular_single_link(path, "launch-bound file")
    data = path.read_bytes()
    return {
        "path": str(path.relative_to(relative_to)) if relative_to is not None else str(path),
        "bytes": len(data),
        "sha256": sha256_bytes(data),
    }


def git(*args: str) -> bytes:
    result = subprocess.run(["git", *args], cwd=REPO, capture_output=True, check=False)
    if result.returncode != 0:
        raise ChildError(f"git {' '.join(args)} failed")
    return result.stdout


def ensure_regular_single_link(path: Path, label: str) -> None:
    try:
        info = path.lstat()
    except FileNotFoundError as error:
        raise ChildError(f"{label} is missing: {path}") from error
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ChildError(f"{label} is not a single-link regular file: {path}")


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    parent_info = path.parent.lstat()
    if not stat.S_ISDIR(parent_info.st_mode) or path.parent.is_symlink():
        raise ChildError(f"atomic-write parent is not a real directory: {path.parent}")
    if path.exists() or path.is_symlink():
        ensure_regular_single_link(path, "atomic-write destination")
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.tmp-"
    )
    temporary = Path(temporary_name)
    try:
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise ChildError("atomic temporary is not a single-link regular file")
        with os.fdopen(descriptor, "wb", closefd=True) as handle:
            descriptor = -1
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        ensure_regular_single_link(temporary, "atomic temporary")
        os.replace(temporary, path)
        directory_fd = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def expected_launch_authorization(
    attempt: Path,
    plan: dict,
    task: dict,
    start: dict,
    artifact: Path = DEFAULT_ARTIFACT,
) -> dict:
    input_path = artifact / task["named_input"]["path"]
    plan_path = artifact / "plan.json"
    body_path = attempt / "request-body.bin"
    return {
        "schema": LAUNCH_AUTH_SCHEMA,
        "task_id": task["id"],
        "ordinal": task["ordinal"],
        "attempt_ordinal": 1,
        "launch_nonce_sha256": start["launch_nonce_sha256"],
        "attempt_path": str(attempt.relative_to(artifact)),
        "attempt_start": regular_record(attempt / "attempt-start.json", attempt),
        "plan": regular_record(plan_path, artifact),
        "input": regular_record(input_path, artifact),
        "request_body": regular_record(body_path, attempt),
        "request": start["request"],
        "execution": {
            "commit": start["execution_commit"],
            "tree": start["execution_tree"],
            "binding_sha256": start["execution_binding_sha256"],
            "manifest_sha256": start["request"]["execution_manifest_sha256"],
        },
        "meter": {
            **regular_record(METER, REPO),
            "hard_watchdog_seconds": 75,
            "orphan_process_group_must_be_terminated": True,
        },
        "child": {
            **regular_record(Path(__file__).resolve(), REPO),
            "argv": [
                str(Path(__file__).resolve()),
                "--attempt",
                str(attempt.resolve()),
                "--execute-child",
            ],
        },
        "single_use": True,
    }


def validate_and_consume_authorization(attempt: Path) -> tuple[dict, dict, bytes]:
    try:
        attempt = attempt.resolve(strict=True)
        attempts_root = (DEFAULT_ARTIFACT / "attempts").resolve(strict=True)
    except FileNotFoundError as error:
        raise ChildError("default Stage 19 attempt ledger is absent") from error
    try:
        attempt.relative_to(attempts_root)
    except ValueError as error:
        raise ChildError("attempt is outside the default committed Stage 19 artifact") from error
    if attempt.parent != attempts_root or attempt.is_symlink() or not attempt.is_dir():
        raise ChildError("attempt is not one direct real directory in the Stage 19 ledger")
    authorization_path = attempt / "launch-authorization.json"
    consumed_path = attempt / "launch-authorization.consumed.json"
    if consumed_path.exists() or consumed_path.is_symlink():
        raise ChildError("launch authorization was already consumed")
    ensure_regular_single_link(authorization_path, "launch authorization")
    start = read_json(attempt / "attempt-start.json")
    if start.get("schema") != ATTEMPT_START_SCHEMA or start.get("attempt_ordinal") != 1:
        raise ChildError("attempt-start ledger is not the exact one-attempt schema")
    nonce = os.environ.pop(LAUNCH_NONCE_ENV, None)
    if (
        not isinstance(nonce, str)
        or not re.fullmatch(r"[0-9a-f]{64}", nonce)
        or sha256_bytes(nonce.encode()) != start.get("launch_nonce_sha256")
    ):
        raise ChildError("inherited launch nonce does not match the parent-bound attempt ledger")
    plan = read_json(PLAN)
    tasks = plan.get("tasks")
    if not isinstance(tasks, list):
        raise ChildError("prepared plan lacks tasks")
    matches = [
        task for task in tasks
        if task.get("id") == start.get("task_id") and task.get("ordinal") == start.get("ordinal")
    ]
    if len(matches) != 1:
        raise ChildError("attempt-start task is absent or duplicated in the prepared plan")
    task = matches[0]
    if attempt.name != Path(task["named_input"]["path"]).stem:
        raise ChildError("attempt directory does not match the frozen task ordinal and input")
    body_path = attempt / "request-body.bin"
    expected = expected_launch_authorization(attempt, plan, task, start)
    observed = read_json(authorization_path)
    if observed != expected or authorization_path.read_bytes() != canonical_bytes(expected):
        raise ChildError("launch authorization differs from the exact parent ledger")
    input_path = DEFAULT_ARTIFACT / task["named_input"]["path"]
    input_bytes = input_path.read_bytes()
    body = urllib.parse.urlencode({FORM_FIELD: input_bytes.decode("ascii")}).encode("ascii")
    if body_path.read_bytes() != body or sha256_bytes(body) != start["request"].get("body_sha256"):
        raise ChildError("request body differs from the parent-bound form encoding")
    execution = start.get("execution_commit")
    if not isinstance(execution, str) or not re.fullmatch(r"[0-9a-f]{40}", execution):
        raise ChildError("attempt lacks an exact execution commit")
    if git("rev-parse", "HEAD").decode().strip() != execution:
        raise ChildError("child checkout is not at the authorized execution commit")
    manifest = read_json(EXECUTION_MANIFEST)
    if sha256_bytes(EXECUTION_MANIFEST.read_bytes()) != start["request"].get("execution_manifest_sha256"):
        raise ChildError("execution manifest hash differs from attempt-start")
    if git("show", f"{execution}:{EXECUTION_MANIFEST.relative_to(REPO)}") != EXECUTION_MANIFEST.read_bytes():
        raise ChildError("execution manifest is not the committed execution blob")
    relevant = {
        record.get("path"): record
        for record in manifest.get("relevant_blobs", [])
        if isinstance(record, dict)
    }
    for path in (PLAN, input_path, METER, Path(__file__).resolve()):
        relative = str(path.relative_to(REPO))
        record = relevant.get(relative)
        current = regular_record(path, REPO)
        if (
            not isinstance(record, dict)
            or current["bytes"] != record.get("bytes")
            or current["sha256"] != record.get("sha256")
            or git("show", f"{execution}:{relative}") != path.read_bytes()
        ):
            raise ChildError(f"launch-bound committed blob changed: {relative}")
    allowed = {
        "attempt-start.json", "request-body.bin", "launch-authorization.json",
        "transport.stdout", "transport.stderr",
    }
    unexpected = {path.name for path in attempt.iterdir()} - allowed
    if unexpected:
        raise ChildError(f"attempt contains pre-existing transport state: {sorted(unexpected)}")
    os.replace(authorization_path, consumed_path)
    directory_fd = os.open(attempt, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    ensure_regular_single_link(consumed_path, "consumed launch authorization")
    if consumed_path.read_bytes() != canonical_bytes(expected):
        raise ChildError("consumed launch authorization changed")
    return plan, task, input_bytes


def selected_headers(headers) -> dict:  # noqa: ANN001
    return {
        name: headers.get(name)
        for name in SAFE_RESPONSE_HEADERS
        if headers.get(name) is not None
    }


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: ANN001
        return None


def bounded_read(handle) -> tuple[bytes, bool]:  # noqa: ANN001
    data = handle.read(MAX_RESPONSE_BYTES + 1)
    if len(data) <= MAX_RESPONSE_BYTES:
        return data, False
    return data[:MAX_RESPONSE_BYTES], True


def perform_request(input_bytes: bytes) -> dict:
    try:
        input_text = input_bytes.decode("ascii")
    except UnicodeDecodeError as error:
        raise ChildError("Magma input is not ASCII") from error
    body = urllib.parse.urlencode({FORM_FIELD: input_text}).encode("ascii")
    request = urllib.request.Request(
        ENDPOINT,
        data=body,
        method="POST",
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": USER_AGENT,
        },
    )
    opener = urllib.request.build_opener(NoRedirect())
    try:
        with opener.open(request, timeout=SOCKET_TIMEOUT_SECONDS) as response:
            payload, exceeded = bounded_read(response)
            return {
                "http_status": int(response.status),
                "final_url": response.geturl(),
                "headers": selected_headers(response.headers),
                "body": payload,
                "body_limit_exceeded": exceeded,
                "transport_error": None,
                "request_body_bytes": len(body),
                "request_body_sha256": sha256_bytes(body),
            }
    except urllib.error.HTTPError as error:
        payload, exceeded = bounded_read(error)
        return {
            "http_status": int(error.code),
            "final_url": error.geturl(),
            "headers": selected_headers(error.headers),
            "body": payload,
            "body_limit_exceeded": exceeded,
            "transport_error": f"HTTPError: {error.reason}",
            "request_body_bytes": len(body),
            "request_body_sha256": sha256_bytes(body),
        }
    except (urllib.error.URLError, TimeoutError, socket.timeout, OSError) as error:
        return {
            "http_status": None,
            "final_url": None,
            "headers": {},
            "body": None,
            "body_limit_exceeded": False,
            "transport_error": f"{type(error).__name__}: {error}",
            "request_body_bytes": len(body),
            "request_body_sha256": sha256_bytes(body),
        }


def persist_result(output: Path, input_bytes: bytes, result: dict, started_at: str) -> dict:
    body = result["body"]
    response_record = None
    if body is not None:
        if result["body_limit_exceeded"]:
            response_name = "response.partial"
        elif result["http_status"] == 200:
            response_name = "response.xml"
        else:
            response_name = "response.body"
        response_path = output / response_name
        atomic_write(response_path, body)
        response_record = {
            "path": response_name,
            "bytes": len(body),
            "sha256": sha256_bytes(body),
            "complete": not result["body_limit_exceeded"],
            "body_limit_exceeded": result["body_limit_exceeded"],
            "response_byte_limit": MAX_RESPONSE_BYTES,
        }
    envelope = {
        "schema": ENVELOPE_SCHEMA,
        "started_at": started_at,
        "finished_at": now(),
        "request": {
            "endpoint": ENDPOINT,
            "method": "POST",
            "form_field": FORM_FIELD,
            "user_agent": USER_AGENT,
            "input_bytes": len(input_bytes),
            "input_sha256": sha256_bytes(input_bytes),
            "body_bytes": result["request_body_bytes"],
            "body_sha256": result["request_body_sha256"],
            "socket_timeout_seconds": SOCKET_TIMEOUT_SECONDS,
        },
        "http": {
            "status": result["http_status"],
            "final_url": result["final_url"],
            "selected_headers": result["headers"],
            "transport_error": result["transport_error"],
        },
        "response": response_record,
        "response_byte_limit": MAX_RESPONSE_BYTES,
    }
    atomic_write(output / "transport-envelope.json", canonical_bytes(envelope))
    return envelope


def execute(attempt: Path) -> dict:
    _, _, input_bytes = validate_and_consume_authorization(attempt)
    started_at = now()
    result = perform_request(input_bytes)
    return persist_result(attempt, input_bytes, result, started_at)


def self_test() -> dict:
    sample = b"1 + 1;\n"
    body = urllib.parse.urlencode({FORM_FIELD: sample.decode("ascii")}).encode("ascii")
    if urllib.parse.parse_qs(body.decode("ascii"), strict_parsing=True) != {
        FORM_FIELD: [sample.decode("ascii")]
    }:
        raise AssertionError("form encoding does not round-trip")
    return {
        "self_test": "pass",
        "response_byte_limit": MAX_RESPONSE_BYTES,
        "socket_timeout_seconds": SOCKET_TIMEOUT_SECONDS,
        "sample_body_sha256": sha256_bytes(body),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempt", type=Path)
    parser.add_argument("--execute-child", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    if not args.execute_child:
        parser.error("the metered parent must pass --execute-child")
    if args.attempt is None:
        parser.error("execution requires --attempt")
    try:
        envelope = execute(args.attempt)
    except ChildError as error:
        parser.error(str(error))
    print(json.dumps({"transport_envelope_sha256": sha256_bytes(canonical_bytes(envelope))}))


if __name__ == "__main__":
    main()
