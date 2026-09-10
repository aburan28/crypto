#!/usr/bin/env python3
"""Verify Stage 19 Amendment 01 and the immutable zero-POST failure ledger."""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import stat
import subprocess
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AMENDMENT = HERE / "stage-19-amendment-01-zero-post-parent-check.json"
ORIGINAL = HERE / "stage-19-magma-calculator-panel-20260909"
CORRECTED = HERE / "stage-19-magma-calculator-panel-amendment-01-20260910"
SUMMARY = HERE / "stage-19-amendment-01-summary.json"
SCHEMA = "koblitz_magma_calculator_stage19_amendment_verification.v1"


class VerificationError(RuntimeError):
    """Amendment 01 or the original zero-POST custody is inconsistent."""


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise VerificationError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


STAGE19 = load_module("stage19_corrected_verifier_for_amendment", HERE / "verify_stage19_magma_calculator_panel.py")


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


def file_record(path: Path, relative_to: Path) -> dict:
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise VerificationError(f"not a single-link regular file: {path}")
    data = path.read_bytes()
    return {
        "path": str(path.relative_to(relative_to)),
        "bytes": len(data),
        "sha256": sha256_bytes(data),
    }


def git_show(revision: str, relative: str) -> bytes:
    result = subprocess.run(
        ["git", "show", f"{revision}:{relative}"],
        cwd=REPO,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise VerificationError(f"cannot read committed source {revision}:{relative}")
    return result.stdout


def call_lines(function: ast.FunctionDef) -> dict[str, list[int]]:
    rows: dict[str, list[int]] = {}
    for node in ast.walk(function):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Name):
            name = node.func.id
        elif isinstance(node.func, ast.Attribute):
            parts = [node.func.attr]
            owner = node.func.value
            while isinstance(owner, ast.Attribute):
                parts.append(owner.attr)
                owner = owner.value
            if isinstance(owner, ast.Name):
                parts.append(owner.id)
            name = ".".join(reversed(parts))
        else:
            continue
        rows.setdefault(name, []).append(node.lineno)
    return rows


def verify_committed_call_order(amendment: dict) -> dict:
    proof = amendment["zero_post_proof"]
    revision = amendment["original_execution"]["execution_commit"]
    relative = proof["committed_child_path"]
    source = git_show(revision, relative)
    if sha256_bytes(source) != proof["committed_child_sha256"]:
        raise VerificationError("committed failing child source hash changed")
    tree = ast.parse(source.decode())
    functions = {
        node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)
    }
    validate = functions.get("validate_and_consume_authorization")
    execute = functions.get("execute")
    if validate is None or execute is None:
        raise VerificationError("committed child lacks the required functions")
    validation_calls = call_lines(validate)
    execution_calls = call_lines(execute)
    parent_line = validation_calls.get("validate_meter_parent", [None])[0]
    consume_line = validation_calls.get("os.replace", [None])[0]
    validate_line = execution_calls.get("validate_and_consume_authorization", [None])[0]
    post_line = execution_calls.get("perform_request", [None])[0]
    if not all(isinstance(value, int) for value in (parent_line, consume_line, validate_line, post_line)):
        raise VerificationError("cannot locate the load-bearing committed call order")
    if not parent_line < consume_line or not validate_line < post_line:
        raise VerificationError("committed child no longer proves failure before consume and POST")
    if parent_line != proof["traceback_call_line"]:
        raise VerificationError("traceback call line differs from committed AST")
    return {
        "execution_commit": revision,
        "committed_child_sha256": sha256_bytes(source),
        "validate_meter_parent_line": parent_line,
        "authorization_consume_line": consume_line,
        "authorization_validation_line": validate_line,
        "perform_request_line": post_line,
        "failure_precedes_authorization_consumption": True,
        "failure_precedes_perform_request": True,
    }


def verify_original(amendment: dict) -> dict:
    expected_files = amendment.get("immutable_files")
    if not isinstance(expected_files, list) or len(expected_files) != 10:
        raise VerificationError("amendment must bind ten immutable original files")
    observed_files = [file_record(ORIGINAL / row["path"], ORIGINAL) for row in expected_files]
    if observed_files != expected_files:
        raise VerificationError("original failure artifact bytes changed")
    attempt = ORIGINAL / "attempts" / "01-seed-2026091301-n31-l5-m3-ggmp-a0-f0"
    expected_attempt_names = {
        "attempt-start.json", "launch-authorization.json", "receipt.json",
        "request-body.bin", "transport-metrics.json", "transport.stderr", "transport.stdout",
    }
    if {path.name for path in attempt.iterdir()} != expected_attempt_names:
        raise VerificationError("original attempt inventory changed")
    forbidden = {
        "launch-authorization.consumed.json", "transport-envelope.json",
        "response.xml", "response.body", "response.partial",
    }
    if any((attempt / name).exists() or (attempt / name).is_symlink() for name in forbidden):
        raise VerificationError("original attempt contains evidence of authorization consume or response")
    metrics = read_json(attempt / "transport-metrics.json")
    charged = amendment["charged_transport_process"]
    resources = metrics["metrics"]
    checks = {
        "returncode": metrics["returncode"],
        "timed_out": metrics["timed_out"],
        "orphan_group_terminated": metrics["orphan_group_terminated"],
        "wall_seconds_exact": resources["wall_seconds"],
        "user_seconds": resources["user_seconds"],
        "system_seconds": resources["system_seconds"],
        "total_core_seconds": resources["total_core_seconds"],
        "peak_rss_bytes": resources["peak_rss_bytes"],
    }
    if checks != {key: charged[key] for key in checks}:
        raise VerificationError("charged original transport metrics changed")
    if not math.isclose(charged["wall_seconds_reported"], round(charged["wall_seconds_exact"], 6), abs_tol=1e-12):
        raise VerificationError("reported wall-time rounding changed")
    stderr = (attempt / "transport.stderr").read_text()
    proof = amendment["zero_post_proof"]
    for marker in (
        proof["failure"],
        "validate_meter_parent(attempt)",
        "PermissionError",
    ):
        if marker not in stderr:
            raise VerificationError(f"original stderr lacks zero-POST marker {marker!r}")
    receipt = read_json(attempt / "receipt.json")
    if (
        receipt.get("outcome") != "transport_child_error"
        or receipt.get("claim_admitted") is not False
        or receipt.get("clean_terminal") is not False
        or receipt.get("transport_envelope") is not None
        or receipt.get("response") is not None
    ):
        raise VerificationError("original receipt boundary changed")
    retained = [
        file_record(path, attempt)
        for path in sorted(attempt.iterdir())
        if path.name != "receipt.json"
    ]
    if receipt.get("retained_files") != retained:
        raise VerificationError("original receipt retained-file custody changed")
    prepared_names = {
        "plan.json",
        "prepared-summary.json",
        *(f"inputs/{path.name}" for path in (CORRECTED / "inputs").glob("*.magma")),
    }
    if len(prepared_names) != 12:
        raise VerificationError("corrected prepared artifact does not contain ten inputs")
    for relative in sorted(prepared_names):
        if (ORIGINAL / relative).read_bytes() != (CORRECTED / relative).read_bytes():
            raise VerificationError(f"original tracked prepared file changed: {relative}")
    run = read_json(ORIGINAL / "run.json")
    summary = read_json(ORIGINAL / "summary.json")
    if (
        run.get("status") != "halted_after_nonclean_receipt"
        or run.get("summary") != summary
        or run.get("receipts", {}).get(amendment["original_execution"]["task_id"], {}).get("sha256")
        != sha256_bytes((attempt / "receipt.json").read_bytes())
    ):
        raise VerificationError("original terminal run index changed")
    return {
        "immutable_files_verified": len(observed_files),
        "tracked_prepared_files_byte_identical_to_correction": len(prepared_names),
        "attempts": 1,
        "post_requests_started": 0,
        "authorization_present_unconsumed": True,
        "response_or_envelope_files": 0,
        "charged_transport_process": charged,
        "terminal_status": run["status"],
    }


def summarize() -> dict:
    amendment = read_json(AMENDMENT)
    if amendment.get("schema") != "koblitz_magma_calculator_stage19_amendment.v1" or amendment.get("amendment") != 1:
        raise VerificationError("unexpected Stage 19 amendment")
    source_order = verify_committed_call_order(amendment)
    original = verify_original(amendment)
    corrected = STAGE19.verify(CORRECTED)
    if corrected.get("artifact_status") != "prepared_not_executed":
        raise VerificationError("corrected artifact is not prepared-only")
    return {
        "schema": SCHEMA,
        "amendment_sha256": sha256_bytes(AMENDMENT.read_bytes()),
        "original_zero_post_failure": original,
        "committed_source_order": source_order,
        "corrected_artifact": {
            "path": str(CORRECTED.relative_to(REPO)),
            "plan_sha256": corrected["plan_sha256"],
            "status": corrected["artifact_status"],
            "attempted_new_requests": corrected["attempted_new_requests"],
        },
        "disposition": "verified_zero_post_operational_failure_with_fresh_corrected_execution_path",
        "claim_boundary": amendment["claim_boundary"],
    }


def self_test() -> dict:
    amendment = read_json(AMENDMENT)
    order = verify_committed_call_order(amendment)
    if not order["failure_precedes_perform_request"]:
        raise AssertionError("zero-POST source-order proof failed")
    return {"self_test": "pass", "amendment": 1, "zero_post_source_order": True}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    try:
        summary = summarize()
    except (VerificationError, STAGE19.VerificationError) as error:
        parser.error(str(error))
    rendered = canonical_bytes(summary)
    if args.expected is not None and args.expected.read_bytes() != rendered:
        parser.error("recomputed Amendment 01 summary differs from expected")
    if args.output is not None:
        STAGE19.RUNNER.atomic_write(args.output.resolve(), rendered)
    print(rendered.decode(), end="")


if __name__ == "__main__":
    main()
