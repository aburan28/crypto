#!/usr/bin/env python3
"""Verify the immutable Amendment 02 ledger and its exact wrapped-marker adjudication."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import stat
import subprocess
from typing import Any
import urllib.parse


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AMENDMENT = HERE / "stage-19-amendment-03-output-fold.json"
SUMMARY = HERE / "stage-19-amendment-03-summary.json"
PREDECESSOR = HERE / "stage-19-magma-calculator-panel-amendment-02-20260910"
PREDECESSOR_COMMIT = "f5dba0fc84a8d6c22c0098e8a8c3506ece55c4f1"
ATTEMPT_STEM = "01-seed-2026091301-n31-l5-m3-ggmp-a0-f0"
TASK_ID = "seed-2026091301/n31-l5-m3-ggmp-a0-f0"
SCHEMA = "koblitz_magma_calculator_stage19_amendment03_verification.v1"
FOLD_WIDTH = 79


class VerificationError(RuntimeError):
    """The predecessor bytes or exact presentation-fold adjudication changed."""


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise VerificationError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


STAGE15 = load_module(
    "stage15_for_stage19_amendment03", HERE / "verify_stage15_magma_calculator.py"
)
LEGACY_RENDERER = load_module(
    "legacy_renderer_for_stage19_amendment03",
    HERE / "render_stage19_magma_calculator_panel.py",
)


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


def regular_record(path: Path, relative_to: Path) -> dict:
    try:
        info = path.lstat()
    except FileNotFoundError as error:
        raise VerificationError(f"missing predecessor file {path}") from error
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise VerificationError(f"predecessor path is not a single-link regular file: {path}")
    data = path.read_bytes()
    return {
        "path": str(path.relative_to(relative_to)),
        "bytes": len(data),
        "sha256": sha256_bytes(data),
    }


def artifact_inventory(root: Path) -> list[dict]:
    if root.is_symlink() or not root.is_dir():
        raise VerificationError("Amendment 02 artifact is absent or a symlink")
    entries = list(root.rglob("*"))
    if any(path.is_symlink() for path in entries):
        raise VerificationError("Amendment 02 artifact contains a symlink")
    return [
        regular_record(path, root)
        for path in sorted(entries)
        if path.is_file()
    ]


def git(*args: str) -> bytes:
    result = subprocess.run(["git", *args], cwd=REPO, capture_output=True, check=False)
    if result.returncode != 0:
        raise VerificationError(
            f"git {' '.join(args)} failed: {result.stderr.decode(errors='replace').strip()}"
        )
    return result.stdout


def verify_committed_inventory(inventory: list[dict]) -> None:
    prefix = str(PREDECESSOR.relative_to(REPO)) + "/"
    tree_lines = [
        line for line in git("ls-tree", "-r", PREDECESSOR_COMMIT, "--", str(PREDECESSOR.relative_to(REPO))).decode().splitlines()
        if line
    ]
    if len(tree_lines) != len(inventory):
        raise VerificationError("Amendment 02 commit inventory count changed")
    by_path = {prefix + row["path"]: row for row in inventory}
    if len(by_path) != len(inventory):
        raise VerificationError("Amendment 02 inventory contains a duplicate path")
    for line in tree_lines:
        try:
            metadata, relative = line.split("\t", 1)
            mode, kind, _oid = metadata.split()
        except ValueError as error:
            raise VerificationError("malformed committed Amendment 02 tree entry") from error
        if mode != "100644" or kind != "blob" or relative not in by_path:
            raise VerificationError("committed Amendment 02 tree differs from current inventory")
        if git("show", f"{PREDECESSOR_COMMIT}:{relative}") != (
            REPO / relative
        ).read_bytes():
            raise VerificationError(f"Amendment 02 bytes differ from {PREDECESSOR_COMMIT}: {relative}")


def exact_folded_identity_lines(task_id: str, source_sha256: str) -> tuple[str, str, str]:
    if len(source_sha256) != 64 or set(source_sha256) - set("0123456789abcdef"):
        raise VerificationError("source SHA-256 is malformed")
    task = f"KOBLITZ_MAGMA_TASK_ID={task_id}"
    source = f"KOBLITZ_MAGMA_SOURCE_SHA256={source_sha256}"
    if len(source) <= FOLD_WIDTH:
        raise VerificationError("source marker no longer exercises the frozen fold")
    return task, source[:FOLD_WIDTH] + "\\", source[FOLD_WIDTH:]


def recover_exact_folded_output(output: str, task_id: str, source_sha256: str) -> str:
    """Accept only the frozen one-backslash fold at the exact expected source marker."""
    lines = output.splitlines()
    expected = exact_folded_identity_lines(task_id, source_sha256)
    if len(lines) != 11 or tuple(lines[:3]) != expected:
        raise VerificationError("response is not the exact expected one-backslash source-marker fold")
    if any("\\" in line for line in lines[3:]):
        raise VerificationError("response contains an unapproved continuation outside the source marker")
    unfolded_source = lines[1][:-1] + lines[2]
    exact_source = f"KOBLITZ_MAGMA_SOURCE_SHA256={source_sha256}"
    if unfolded_source != exact_source:
        raise VerificationError("fold reconstruction differs from the exact expected source marker")
    normalized_lines = [lines[0], unfolded_source, *lines[3:]]
    return "\n".join(normalized_lines) + "\n"


def verify_predecessor(amendment: dict) -> dict:
    inventory = artifact_inventory(PREDECESSOR)
    observed_inventory = {
        "files": len(inventory),
        "inventory_sha256": sha256_bytes(canonical_bytes(inventory)),
        "commit": PREDECESSOR_COMMIT,
    }
    if amendment.get("immutable_amendment02") != observed_inventory:
        raise VerificationError("Amendment 02 inventory binding changed")
    verify_committed_inventory(inventory)

    plan = read_json(PREDECESSOR / "plan.json")
    expected_plan, expected_inputs = LEGACY_RENDERER.build_plan()
    if plan != expected_plan or (PREDECESSOR / "plan.json").read_bytes() != canonical_bytes(expected_plan):
        raise VerificationError("Amendment 02 plan no longer regenerates from the Stage 13 sources")
    tasks = plan.get("tasks")
    if not isinstance(tasks, list) or len(tasks) != 10 or tasks[0].get("id") != TASK_ID:
        raise VerificationError("Amendment 02 plan first-task identity changed")
    task = tasks[0]
    attempt = PREDECESSOR / "attempts" / ATTEMPT_STEM
    start = read_json(attempt / "attempt-start.json")
    receipt = read_json(attempt / "receipt.json")
    envelope = read_json(attempt / "transport-envelope.json")
    run = read_json(PREDECESSOR / "run.json")
    summary = read_json(PREDECESSOR / "summary.json")
    if (
        start.get("task_id") != TASK_ID
        or start.get("ordinal") != 1
        or receipt.get("outcome") != "identity_or_terminal_mismatch"
        or receipt.get("clean_terminal") is not False
        or receipt.get("claim_admitted") is not False
        or receipt.get("retry_permitted") is not False
        or run.get("status") != "halted_after_nonclean_receipt"
        or run.get("summary") != summary
        or summary.get("attempted_tasks") != 1
        or summary.get("unattempted_tasks") != 9
        or summary.get("clean_f4_terminals") != 0
    ):
        raise VerificationError("Amendment 02 terminal non-admitted ledger semantics changed")

    input_path = PREDECESSOR / str(task["named_input"]["path"])
    input_bytes = input_path.read_bytes()
    if (
        regular_record(input_path, PREDECESSOR) != task["named_input"]
        or expected_inputs.get(input_path.name) != input_bytes
    ):
        raise VerificationError("Amendment 02 first named input changed")
    expected_body = urllib.parse.urlencode(
        {plan["service"]["form_field"]: input_bytes.decode("ascii")}
    ).encode("ascii")
    body = (attempt / "request-body.bin").read_bytes()
    if body != expected_body or sha256_bytes(body) != start["request"]["body_sha256"]:
        raise VerificationError("Amendment 02 retained POST body differs from its bound input")
    http = envelope.get("http", {})
    retained_response = envelope.get("response", {})
    response_path = attempt / "response.xml"
    if (
        http.get("status") != 200
        or http.get("final_url") != plan["service"]["endpoint"]
        or http.get("transport_error") is not None
        or retained_response.get("complete") is not True
        or retained_response.get("body_limit_exceeded") is not False
        or retained_response.get("sha256") != sha256_bytes(response_path.read_bytes())
    ):
        raise VerificationError("Amendment 02 successful HTTP response custody changed")

    response = STAGE15.parse_calculator_xml(response_path)
    if response["service"]["warning"] is not None or response["service"]["alert"] is not None:
        raise VerificationError("Amendment 02 response contains a service diagnostic")
    normalized = recover_exact_folded_output(
        response["output"], TASK_ID, str(task["source_instance_sha256"])
    )
    prefix = (
        f"KOBLITZ_MAGMA_TASK_ID={TASK_ID}\n"
        f"KOBLITZ_MAGMA_SOURCE_SHA256={task['source_instance_sha256']}\n"
    )
    if not normalized.startswith(prefix):
        raise VerificationError("normalized response identity changed")
    terminal_output = normalized[len(prefix):]
    try:
        STAGE15.require_clean_f4_output(terminal_output, int(task["seed"]))
    except STAGE15.VerificationError as error:
        raise VerificationError(str(error)) from error
    terminal = STAGE15.parse_magma_terminal(terminal_output)
    if terminal is None or terminal.get("terminal_status") != "sat":
        raise VerificationError("recovered response is not a strict SAT F4 terminal")
    if not all(
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
        and value >= 0
        for value in (terminal["cpu_seconds"], terminal["wall_seconds"])
    ):
        raise VerificationError("recovered terminal contains an invalid numeric value")

    derived = {
        "task_id": TASK_ID,
        "source_instance_sha256": task["source_instance_sha256"],
        "named_input": task["named_input"],
        "request_body": regular_record(attempt / "request-body.bin", PREDECESSOR),
        "response": regular_record(response_path, PREDECESSOR),
        "raw_output_sha256": sha256_bytes(response["output"].encode()),
        "normalized_output_sha256": sha256_bytes(normalized.encode()),
        "fold": {
            "scope": "source_identity_marker_only",
            "content_columns_before_backslash": FOLD_WIDTH,
            "first_raw_line_bytes": len(exact_folded_identity_lines(TASK_ID, task["source_instance_sha256"])[1].encode()),
            "tail_hex_characters": len(exact_folded_identity_lines(TASK_ID, task["source_instance_sha256"])[2]),
            "generic_unfolding_permitted": False,
        },
        "terminal": terminal,
        "service": response["service"],
        "classification": "sat_basis_certificate_unverified_model",
        "source_equivalent": True,
    }
    if amendment.get("derived_adjudication") != derived:
        raise VerificationError("Amendment 03 derived adjudication changed")
    return {
        "immutable_files_verified": len(inventory),
        "immutable_inventory_sha256": observed_inventory["inventory_sha256"],
        "original_receipt_outcome": receipt["outcome"],
        "original_claim_admitted": receipt["claim_admitted"],
        "original_ledger_remains_terminal": True,
        "task_retry_permitted": False,
        "derived": derived,
    }


def summarize() -> dict:
    amendment = read_json(AMENDMENT)
    if (
        amendment.get("schema") != "koblitz_magma_calculator_stage19_amendment.v1"
        or amendment.get("amendment") != 3
    ):
        raise VerificationError("unexpected Stage 19 Amendment 03 record")
    predecessor = verify_predecessor(amendment)
    if amendment.get("original_ledger_policy") != {
        "amendment02_bytes_mutable": False,
        "in_place_reclassification_permitted": False,
        "old_ledger_resume_permitted": False,
        "started_task_retry_permitted": False,
        "derived_adjudication_is_additive": True,
    }:
        raise VerificationError("Amendment 03 original-ledger policy changed")
    if amendment.get("successor_policy") != {
        "recovered_task_id": TASK_ID,
        "successor_contains_original_stage19_ordinals": list(range(2, 11)),
        "successor_local_ordinals": list(range(1, 10)),
        "predecessor_attempted_task_ids_must_be_disjoint": True,
        "future_source_digest_markers": "two independently labeled 32-hex-character halves",
        "generic_backslash_unfolding_permitted": False,
        "network_execution_in_this_amendment": False,
    }:
        raise VerificationError("Amendment 03 successor policy changed")
    return {
        "schema": SCHEMA,
        "amendment_sha256": sha256_bytes(AMENDMENT.read_bytes()),
        "immutable_amendment02": {
            "files_verified": predecessor["immutable_files_verified"],
            "inventory_sha256": predecessor["immutable_inventory_sha256"],
            "commit": PREDECESSOR_COMMIT,
        },
        "original_non_admitted_ledger": {
            "outcome": predecessor["original_receipt_outcome"],
            "claim_admitted": predecessor["original_claim_admitted"],
            "terminal": predecessor["original_ledger_remains_terminal"],
            "retry_permitted": predecessor["task_retry_permitted"],
        },
        "recovered_terminal": predecessor["derived"],
        "disposition": "verified_exact_presentation_fold_parser_false_negative",
        "public_service_f4_terminal_coverage_after_additive_adjudication": 6,
        "claim_boundary": amendment["claim_boundary"],
    }


def self_test() -> dict:
    amendment = read_json(AMENDMENT)
    predecessor = verify_predecessor(amendment)
    task = predecessor["derived"]
    source = task["source_instance_sha256"]
    folded = exact_folded_identity_lines(TASK_ID, source)
    valid = "\n".join(
        [
            *folded,
            "KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1",
            "KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse",
            "KOBLITZ_MAGMA_STATUS=SAT",
            "KOBLITZ_MAGMA_F4_DEGREES=[ 2 ]",
            "KOBLITZ_MAGMA_BASIS_SIZE=2",
            "KOBLITZ_MAGMA_CPU_SECONDS=0",
            "KOBLITZ_MAGMA_WALL_SECONDS=0",
            "",
        ]
    ) + "\n"
    recover_exact_folded_output(valid, TASK_ID, source)
    for invalid in (
        valid.replace("\\\n", "\n", 1),
        valid.replace("\\\n", "\\\n0", 1),
        valid.replace("KOBLITZ_MAGMA_SCHEMA", "KOBLITZ_MAGMA_SCHEMA\\", 1),
    ):
        try:
            recover_exact_folded_output(invalid, TASK_ID, source)
        except VerificationError:
            pass
        else:
            raise AssertionError("non-exact continuation was accepted")
    return {
        "self_test": "pass",
        "immutable_files_verified": predecessor["immutable_files_verified"],
        "exact_fold_only": True,
        "generic_unfolding_permitted": False,
        "network_requests": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        summary = self_test() if args.self_test else summarize()
    except (VerificationError, STAGE15.VerificationError) as error:
        parser.error(str(error))
    rendered = canonical_bytes(summary)
    if args.expected is not None and args.expected.read_bytes() != rendered:
        parser.error("recomputed Amendment 03 summary differs from expected")
    if args.output is not None:
        args.output.write_bytes(rendered)
    print(rendered.decode(), end="")


if __name__ == "__main__":
    main()
