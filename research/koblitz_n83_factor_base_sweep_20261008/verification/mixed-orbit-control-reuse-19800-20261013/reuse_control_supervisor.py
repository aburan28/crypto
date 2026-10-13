#!/usr/bin/env python3
"""One-shot N83 known-answer control using the retained K600 Linux worker.

This source preparation launches nothing. ``run`` requires a fresh, explicit
1,200-second grant after the audited 19,800-second local ledger. The prior
Linux build is already charged and is reused by exact digest, not rebuilt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import uuid


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD_CONTROL = HERE.parent / "mixed-orbit-known-witness-source-20261012"
sys.path.insert(0, str(OLD_CONTROL))
import control_supervisor as prior  # noqa: E402


STUDY = prior.STUDY
EVIDENCE = f"research/{STUDY}/verification/k600-mixed-orbit-first-2400-20261012"
AUDIT = f"{EVIDENCE}/audit.json"
BUILD = f"{EVIDENCE}/evidence-01/build/build.json"
BINARY = f"{EVIDENCE}/evidence-01/build/worker"
AUDIT_SHA256 = "387d393b9805f35ee1ca6328e8602133a1b9bc036088b1be1a3ba5a29bb8469c"
BUILD_SHA256 = "cd302908d34d52412e330fb987e605e288af69d5bf4b50ebbea968e7a2112804"
BINARY_SHA256 = "26dcb2b49388fcb89853eabbda01db71d148620cd465bd5458099e0606dcf606"
BINARY_SOURCE_COMMIT = "f3d27f68b68cda907e0668042ac7fa13ae6b579d"
CROSS_AUDIT = f"research/{STUDY}/verification/cross-stream-first-gate-source-20261013/audit.json"
CROSS_AUDIT_SHA256 = "dfc96f11e4a399a42559847dd39e94b753f37911a54d806e0e54e830e7735478"
CROSS_BUDGET_SHA256 = "aab1eeff9c02cd9afced71f883748c0b6fc58f7acc98631415336b57ee5d023d"
CROSS_SOURCE_COMMIT = "0764901992591c39566d720ae812c54cdc09ab75"
K600_TOTAL = 17400
PRIOR_TOTAL = 19800
VERIFY_CAP, SEARCH_CAP, CLEANUP_RESERVE = 120, 900, 180
VERIFY_WALL, SEARCH_WALL, WORKER_WALL = 90, 850, 840
TOTAL_CHARGE = VERIFY_CAP + SEARCH_CAP + CLEANUP_RESERVE
CONFLICTS, MODELS = 5000, 32
SOURCES = (
    BINARY, BUILD, AUDIT, CROSS_AUDIT,
    prior.FIXTURE, prior.SMALL_MANIFEST, prior.ENVELOPE_MANIFEST,
    f"research/{STUDY}/verification/mixed-orbit-known-witness-source-20261012/known_witness.py",
    f"research/{STUDY}/verification/mixed-orbit-known-witness-source-20261012/control_supervisor.py",
    f"research/{STUDY}/verification/mixed-orbit-control-reuse-19800-20261013/reuse_control_supervisor.py",
)


def authorization_error(record: object, output: Path) -> str | None:
    if not isinstance(record, dict):
        return "authorization is not an object"
    expected = {
        "schema": "n83.additional-local-budget/v7",
        "study": STUDY,
        "purpose": "one_known_witness_control_reuse_k600_worker_after_cross_stream",
        "prior_total_charged_seconds": str(PRIOR_TOTAL),
        "prior_k600_audit_sha256": AUDIT_SHA256,
        "prior_cross_stream_audit_sha256": CROSS_AUDIT_SHA256,
        "reuse_binary_sha256": BINARY_SHA256,
        "additional_authorized_seconds": str(TOTAL_CHARGE),
        "output_dir": str(output.resolve()),
    }
    for key, value in expected.items():
        if record.get(key) != value:
            return f"authorization {key} mismatch"
    approval = record.get("user_approval_reference")
    if not isinstance(approval, str) or not approval.strip():
        return "authorization lacks explicit approval reference"
    return None


def load_auth(path: Path, output: Path) -> tuple[dict, str]:
    raw = path.resolve(strict=True).read_bytes()
    record = json.loads(raw)
    error = authorization_error(record, output)
    if error:
        raise ValueError(error)
    return record, hashlib.sha256(raw).hexdigest()


def prior_evidence_error(checkout: Path) -> str | None:
    try:
        if (prior.sha256_file(checkout / AUDIT) != AUDIT_SHA256
                or prior.sha256_file(checkout / CROSS_AUDIT) != CROSS_AUDIT_SHA256
                or prior.sha256_file(checkout / BUILD) != BUILD_SHA256
                or prior.sha256_file(checkout / BINARY) != BINARY_SHA256
                or not prior.is_linux_elf(checkout / BINARY)):
            return "pinned prior evidence bytes differ"
        audit = json.loads((checkout / AUDIT).read_text())
        cross = json.loads((checkout / CROSS_AUDIT).read_text())
        built = json.loads((checkout / BUILD).read_text())
    except (OSError, UnicodeError, json.JSONDecodeError):
        return "pinned prior evidence missing or malformed"
    if not isinstance(audit, dict) or not isinstance(cross, dict) or not isinstance(built, dict):
        return "pinned prior evidence is not an object"
    stage = built.get("stage")
    if not isinstance(stage, dict):
        return "pinned K600 build stage missing"
    state = stage.get("container_state")
    if (audit.get("schema") != "n83.k600-mixed-orbit-first-gate-audit/v1"
            or audit.get("status") != "PASS"
            or audit.get("search_status") != "UNKNOWN_search_cap"
            or audit.get("conservative_charge_seconds") != 2400
            or audit.get("cumulative_conservative_charge_seconds") != K600_TOTAL
            or audit.get("source_commit") != BINARY_SOURCE_COMMIT
            or audit.get("binary_sha256") != BINARY_SHA256
            or audit.get("selected_best_total_runtime") is not None
            or cross.get("schema") != "n83.cross-stream-first-gate-audit/v1"
            or cross.get("status") != "PASS_receipts_with_censored_search"
            or cross.get("search_status") != "UNKNOWN_search_cap"
            or cross.get("search_models") != 0
            or cross.get("search_candidate") is not None
            or cross.get("conservative_charge_seconds") != 2400
            or cross.get("cumulative_conservative_charge_seconds") != PRIOR_TOTAL
            or cross.get("source_commit") != CROSS_SOURCE_COMMIT
            or cross.get("selected_best_total_runtime") is not None
            or not isinstance(cross.get("raw_sha256"), dict)
            or cross["raw_sha256"].get("budget.json") != CROSS_BUDGET_SHA256
            or built.get("schema") != "n83.mixed-orbit-k600-linux-build/v1"
            or built.get("status") != "PASS_build_only"
            or built.get("source_commit") != BINARY_SOURCE_COMMIT
            or built.get("binary_sha256") != BINARY_SHA256
            or built.get("container_image_id") != prior.IMAGE
            or built.get("conservative_charge_seconds") != 600
            or stage.get("cleanup_verified") is not True
            or stage.get("exit_code") != 0
            or stage.get("timed_out") is not False
            or not isinstance(state, dict)
            or state.get("ExitCode") != 0
            or state.get("OOMKilled") is not False):
        return "pinned prior audit or build disposition differs"
    return None


def run(authorization: Path, output_dir: Path, small_object: Path,
        envelope_object: Path, checkout: Path = ROOT) -> dict:
    checkout = checkout.resolve(strict=True)
    output_dir = prior.fresh_outside(output_dir, checkout)
    authorization = authorization.resolve(strict=True)
    small_object = small_object.resolve(strict=True)
    envelope_object = envelope_object.resolve(strict=True)
    _, auth_sha = load_auth(authorization, output_dir)
    commit = prior.clean_commit(checkout)
    evidence_error = prior_evidence_error(checkout)
    if evidence_error:
        raise ValueError(evidence_error)
    if (prior.sha256_file(small_object) != prior.SMALL_SHA
            or prior.sha256_file(envelope_object) != prior.ENVELOPE_SHA):
        raise ValueError("retained K1182/K2048 archive bytes differ")
    fixture = prior.collect(small_object, envelope_object)
    if fixture != json.loads((checkout / prior.FIXTURE).read_text()):
        raise ValueError("known-answer fixture differs")
    prior.image_identity()
    binary = checkout / BINARY
    inputs = {str(path): prior.sha256_file(path) for path in (
        authorization, small_object, envelope_object,
        *(checkout / relative for relative in SOURCES))}
    started = time.monotonic()
    deadline = started + TOTAL_CHARGE
    output_dir.mkdir(parents=True, exist_ok=False)
    config = {
        "schema": "n83.known-witness-reuse-supervisor-config/v1",
        "source_commit": commit,
        "binary_source_commit": BINARY_SOURCE_COMMIT,
        "reuse_binary_sha256": BINARY_SHA256,
        "prior_k600_audit_sha256": AUDIT_SHA256,
        "prior_cross_stream_audit_sha256": CROSS_AUDIT_SHA256,
        "prior_total_charged_seconds": PRIOR_TOTAL,
        "authorization_sha256": auth_sha,
        "container_image_id": prior.IMAGE,
        "input_sha256": inputs,
        "verify_charge_seconds": VERIFY_CAP,
        "search_charge_seconds": SEARCH_CAP,
        "cleanup_reserve_seconds": CLEANUP_RESERVE,
        "conservative_full_charge_seconds": TOTAL_CHARGE,
        "conflicts": CONFLICTS, "models": MODELS,
        "scope": "known_answer_control_only",
        "natural_relations": None,
        "total_index_calculus_runtime_ms": None,
    }
    prior.write_new_json(output_dir / "config.json", config)
    config_sha = prior.sha256_file(output_dir / "config.json")
    stages = []
    status = "PRODUCER_FAILURE_preflight"
    verify_passed = False
    replay_error = None
    token = uuid.uuid4().hex[:16]
    for name, wall in (("verify", VERIFY_WALL), ("control", SEARCH_WALL)):
        if time.monotonic() + wall + CLEANUP_RESERVE > deadline:
            status = "UNKNOWN_budget_before_stage"
            break
        container_name = f"n83-control-reuse-{name}-{token}"
        command = [
            "docker", "run", "--pull", "never", "--platform", "linux/arm64",
            "--name", container_name, "--network", "none", "--cpus", "1",
            "--memory", "4096m", "--memory-swap", "4096m", "--pids-limit", "64",
            "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
            "--read-only", "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m",
            "--user", f"{os.getuid()}:{os.getgid()}",
            "--mount", f"type=bind,src={binary},dst=/worker,readonly",
            "--mount", f"type=bind,src={small_object},dst=/small.jsonl.gz,readonly",
            "--mount", f"type=bind,src={envelope_object},dst=/envelope.jsonl.gz,readonly",
            "--mount", f"type=bind,src={checkout / prior.SMALL_MANIFEST},dst=/small-manifest.json,readonly",
            "--mount", f"type=bind,src={checkout / prior.ENVELOPE_MANIFEST},dst=/envelope-manifest.json,readonly",
            "--mount", f"type=bind,src={checkout / prior.FIXTURE},dst=/fixture.json,readonly",
            "--mount", f"type=bind,src={output_dir},dst=/out",
            prior.IMAGE, "/worker", "verify-control" if name == "verify" else "control-run",
            "/small-manifest.json", "/small.jsonl.gz",
            "/envelope-manifest.json", "/envelope.jsonl.gz", "/fixture.json",
        ]
        if name == "control":
            command += [str(WORKER_WALL), str(CONFLICTS), str(MODELS), "/out/control"]
        receipt = prior.stage(name=name, command=command, wall=wall,
                              output=output_dir, container_name=container_name)
        stages.append(receipt)
        state = receipt.get("container_state")
        if not receipt.get("cleanup_verified"):
            status = "PRODUCER_FAILURE_cleanup"
        elif receipt.get("timed_out"):
            status = "UNKNOWN_wall_cap"
        elif state and state.get("OOMKilled") is True:
            status = "UNKNOWN_resource_cap"
        elif receipt.get("exit_code") == 124:
            try:
                cap = json.loads((output_dir / "control/cap.json").read_text()) if name == "control" else None
            except (OSError, UnicodeError, json.JSONDecodeError):
                cap = None
            status = ("UNKNOWN_worker_wall_cap" if prior.valid_worker_cap(cap)
                      else "PRODUCER_FAILURE_worker_exit")
        elif receipt.get("exit_code") in (134, 137):
            status = "UNKNOWN_resource_cap"
        elif (receipt.get("launch_error") or receipt.get("exit_code") != 0
              or not isinstance(state, dict) or state.get("ExitCode") != 0):
            status = "PRODUCER_FAILURE_worker_exit"
        elif name == "verify":
            try:
                row = json.loads((output_dir / "verify.stdout").read_text())
            except (OSError, UnicodeError, json.JSONDecodeError):
                row = None
            if (not isinstance(row, dict) or row.get("status") != "PASS"
                    or row.get("schema") != "n83.mixed-orbit-known-witness-verify/v1"
                    or row.get("scope") != "known_answer_control_only"
                    or row.get("solver_search_executed") is not False
                    or row.get("target") != fixture["equation"]["target"]):
                status = "PRODUCER_FAILURE_verify_receipt"
            else:
                status = "PASS_verify_only"
                verify_passed = True
        else:
            try:
                worker_config = json.loads((output_dir / "control/config.json").read_text())
                worker_result = json.loads((output_dir / "control/result.json").read_text())
                replay_error = prior.control_receipt_error(
                    worker_config, worker_result, fixture, small_object, envelope_object)
            except (OSError, UnicodeError, ValueError, KeyError, TypeError,
                    EOFError, json.JSONDecodeError) as error:
                replay_error = f"control replay raised {type(error).__name__}: {error}"
            status = "PRODUCER_FAILURE_control_replay" if replay_error else worker_result["status"]
        if status != "PASS_verify_only":
            break
    try:
        frozen_changed = (prior.clean_commit(checkout) != commit
                          or any(prior.sha256_file(Path(path)) != digest
                                 for path, digest in inputs.items())
                          or prior.sha256_file(output_dir / "config.json") != config_sha)
    except (OSError, ValueError):
        frozen_changed = True
    if frozen_changed:
        status = "PRODUCER_FAILURE_frozen_input_changed"
    elapsed = time.monotonic() - started
    if elapsed > TOTAL_CHARGE:
        status = "PRODUCER_FAILURE_budget_overshoot"
    outer = {
        "schema": "n83.known-witness-reuse-supervisor-outer/v1",
        "status": status, "scope": "known_answer_control_only",
        "source_commit": commit, "binary_source_commit": BINARY_SOURCE_COMMIT,
        "prior_cross_stream_audit_sha256": CROSS_AUDIT_SHA256,
        "authorization_sha256": auth_sha,
        "config_sha256": config_sha, "stages": stages,
        "verify_passed": verify_passed, "replay_error": replay_error,
        "frozen_input_changed": frozen_changed,
        "elapsed_supervisor_seconds_L0": elapsed,
        "natural_relations": None,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    prior.write_new_json(output_dir / "outer.json", outer)
    charge = (sum(VERIFY_CAP if row["name"] == "verify" else SEARCH_CAP
                  for row in stages) + (CLEANUP_RESERVE if stages else 0))
    prior.write_new_json(output_dir / "budget.json", {
        "schema": "n83.known-witness-reuse-supervisor-budget/v1",
        "prior_total_charged_seconds": PRIOR_TOTAL,
        "prior_k600_audit_sha256": AUDIT_SHA256,
        "prior_cross_stream_audit_sha256": CROSS_AUDIT_SHA256,
        "reuse_binary_sha256": BINARY_SHA256,
        "authorization_sha256": auth_sha,
        "conservative_charge_seconds": charge,
        "cumulative_conservative_charge_seconds": PRIOR_TOTAL + charge,
        "remaining_additional_seconds": TOTAL_CHARGE - charge,
        "round_status": status,
        "selected_best_total_runtime": None,
    })
    return outer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--small-object", type=Path, required=True)
    parser.add_argument("--envelope-object", type=Path, required=True)
    args = parser.parse_args()
    row = run(args.authorization, args.output_dir, args.small_object,
              args.envelope_object)
    print(json.dumps(row, sort_keys=True))
    return 1 if row["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
