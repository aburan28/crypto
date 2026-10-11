#!/usr/bin/env python3
"""Supervise one hash-pinned N83 orbit import and one bounded SAT trial.

The 7,200-second local pilot is exhausted. This runner requires a separate,
run-specific additional-budget receipt; merely preparing this source does not
authorize or launch the gate. Every launched stage gets a no-network, one-CPU,
4 GiB, zero-swap container and leaves an outer receipt, including failures.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time
import uuid

from compact_probe_supervisor import clean_commit, sha256_file, write_new_json


STUDY = "koblitz_n83_factor_base_sweep_20261008"
ROOT = Path(__file__).resolve().parents[2]
OBJECT_SHA256 = "c5de0b26e84ad80bfb7be6a214b89c9682e356effb6079c7abf28cd1151ffb03"
OBJECT_BLAKE3 = "aea278aff2d79a592d5164d7437ddde8f8e9377fa481058dfc16afa176a31389"
MANIFEST_BLAKE3 = "d6adad9bec74350f37bf917e737006869bc0f343692f668fa4dc893afe332745"
CORPUS_BLAKE3 = "831eaa7be2387634fd610f6706103b2af9b4e274f0db59923c0280af3cbed74b"
CURVE = "icv1-f2m83-tm6151469093347-debefd74"
IMAGE = "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31"
VERIFY_WALL = 240
SEARCH_WALL = 720
WORKER_WALL = 710
TOTAL_CHARGE = 1200
BUILD_CHARGE = 600
TOTAL_ADDITIONAL_CHARGE = BUILD_CHARGE + TOTAL_CHARGE
CLEANUP_RESERVE = 120
MEMORY_MIB = 4096
TRIALS = 1
CONFLICTS = 10_000
MODELS = 32
HEX64 = re.compile(r"[0-9a-f]{64}\Z")
SOURCES = (
    "examples/koblitz_n83_structured_orbit_worker.rs",
    "examples/koblitz_n83_structured_orbit_export.rs",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/binary_semaev_chain_sat.rs",
    f"research/{STUDY}/structured-orbit-v3.json",
    f"research/{STUDY}/structured_orbit_search_supervisor.py",
    f"research/{STUDY}/verification/structured-orbit-search-gate-20261010/PROTOCOL.md",
)
MANIFEST = f"research/{STUDY}/verification/structured-d12-orbit-20261010/manifest.json"
CORPUS = f"research/{STUDY}/pilot-01/probe-corpus.json"
VALIDATION = f"research/{STUDY}/pilot-01/probe-validation.json"


def authorization_error(record: object, output: Path, build_dir: Path) -> str | None:
    if not isinstance(record, dict):
        return "additional-budget receipt must be an object"
    expected = {
        "schema": "n83.additional-local-budget/v1",
        "study": STUDY,
        "purpose": "structured_orbit_import_and_one_m5_trial",
        "prior_pilot_charged_seconds": "7200",
        "output_dir": str(output.resolve()),
        "build_dir": str(build_dir.resolve()),
    }
    for key, value in expected.items():
        if record.get(key) != value:
            return f"additional-budget {key} mismatch"
    amount = record.get("additional_authorized_seconds")
    if not isinstance(amount, str) or not amount.isdecimal():
        return "additional-budget amount must be a canonical decimal string"
    try:
        allowed = Decimal(amount)
    except (KeyError, TypeError, ValueError):
        return "additional-budget amount missing"
    if not allowed.is_finite() or allowed < TOTAL_ADDITIONAL_CHARGE:
        return "additional-budget amount below the full conservative charge"
    approval = record.get("user_approval_reference")
    if not isinstance(approval, str) or not approval.strip():
        return "additional-budget approval reference missing"
    return None


def verify_receipt_error(value: object) -> str | None:
    if not isinstance(value, dict):
        return "archive import receipt is not an object"
    expected = {
        "schema": "n83.structured-orbit-import/v1", "status": "PASS",
        "curve_slug": CURVE, "object_blake3": OBJECT_BLAKE3,
        "points": 332166, "signed_orbits": 2001,
        "x_coordinates": 166083, "solver_search_executed": False,
    }
    for key, item in expected.items():
        if type(value.get(key)) is not type(item) or value[key] != item:
            return f"archive import {key} mismatch"
    return None


def valid_worker_cap(value: object) -> bool:
    if not isinstance(value, dict):
        return False
    elapsed = value.get("elapsed_process_wall_ms")
    return (value.get("schema") == "n83.structured-orbit-cold-cap/v1"
            and value.get("status") == "UNKNOWN_budget"
            and type(value.get("budget_seconds")) is int
            and value["budget_seconds"] == WORKER_WALL
            and isinstance(elapsed, (int, float)) and not isinstance(elapsed, bool)
            and math.isfinite(elapsed) and elapsed >= 0)


def search_receipt_error(config: object, result: object) -> str | None:
    if not isinstance(config, dict) or not isinstance(result, dict):
        return "search config or result is not an object"
    expected_config = {
        "schema": "n83.structured-orbit-cold-config/v1", "study": STUDY,
        "object_blake3": OBJECT_BLAKE3, "manifest_blake3": MANIFEST_BLAKE3,
        "public_corpus_path": "/corpus.json", "validation_path": "/validation.json",
        "budget_seconds": WORKER_WALL, "max_trials": TRIALS,
        "conflict_budget": CONFLICTS, "max_models": MODELS,
        "summands": 5, "strategy": "ChainedS3FrobeniusOrbit",
        "memory_cgroup_limit_bytes": 4 * 1024**3,
        "memory_cgroup_swap_limit_bytes": 0,
    }
    for key, item in expected_config.items():
        if type(config.get(key)) is not type(item) or config[key] != item:
            return f"search config {key} mismatch"
    if not isinstance(config.get("source_blake3"), str) or not HEX64.fullmatch(config["source_blake3"]):
        return "search worker source hash missing"
    expected_result = {
        "schema": "n83.structured-orbit-cold-result/v1", "study": STUDY,
        "object_blake3": OBJECT_BLAKE3, "manifest_blake3": MANIFEST_BLAKE3,
        "factor_base_size": 332166, "orbit_count": 2001,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    for key, item in expected_result.items():
        if type(result.get(key)) is not type(item) or result[key] != item:
            return f"search result {key} mismatch"
    if result.get("public_corpus_canonical_json_blake3") != CORPUS_BLAKE3:
        return "public target corpus hash mismatch"
    peak = result.get("memory_cgroup_peak_bytes")
    if type(peak) is not int or not 0 < peak <= 4 * 1024**3:
        return "search cgroup peak is outside the declared memory limit"
    counters = (
        "trials", "relations", "independent_relations", "dependent_relations",
        "inconsistent_relations", "verification_failures", "sat_calls", "sat_unknowns",
        "sat_invalid_models", "sat_models", "sat_conflicts", "matrix_rows",
        "matrix_columns", "terminal_matrix_rank", "linear_solve_attempts",
    )
    if any(type(result.get(key)) is not int or result[key] < 0 for key in counters):
        return "search has malformed relation, rank or SAT counters"
    if (result["trials"] > TRIALS or result["relations"] > result["trials"]
            or result["independent_relations"] + result["dependent_relations"] != result["relations"]
            or result["inconsistent_relations"] or result["verification_failures"]
            or result["sat_invalid_models"]
            or result["matrix_columns"] != 2002
            or result["terminal_matrix_rank"] > result["matrix_columns"]):
        return "search relation, rank or failure counters disagree"
    if (not isinstance(result.get("attempt_records"), list)
            or not isinstance(result.get("relation_matrix"), list)
            or len(result["relation_matrix"]) != result["relations"]):
        return "search attempt or relation ledger missing"
    for relation in result["relation_matrix"]:
        if not isinstance(relation, dict):
            return "search relation row is not an object"
        if (not isinstance(relation.get("row"), list) or len(relation["row"]) != 2001
                or not all(isinstance(value, str) and value.isdecimal() for value in relation["row"])):
            return "search relation row has an invalid column vector"
        if (not isinstance(relation.get("summands"), list) or len(relation["summands"]) != 5
                or not isinstance(relation.get("summand_negated"), list)
                or len(relation["summand_negated"]) != 5):
            return "search relation has an invalid five-summand witness"
    if (not isinstance(result.get("rank_history"), list)
            or len(result["rank_history"]) != result["linear_solve_attempts"]):
        return "search rank ledger missing"
    for key in ("base_import_elapsed_ms_L0", "worker_wall_ms_L0"):
        value = result.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            return f"search {key} is not a duration"
    status = result.get("status")
    if status not in ("UNKNOWN_solver_cap", "UNKNOWN_trial_cap"):
        return "search status is not an admissible bounded outcome"
    if result.get("candidate_verified") is not False or result.get("candidate_log") is not None:
        return "one-trial gate cannot certify the 2,001-column target"
    if status == "UNKNOWN_solver_cap" and result["sat_unknowns"] == 0:
        return "solver-cap status has no unknown SAT call"
    if status == "UNKNOWN_trial_cap" and result["sat_unknowns"] != 0:
        return "trial-cap status hides an unknown SAT call"
    return None


def stage(
    *, name: str, command: list[str], wall: int, output: Path,
    container_name: str,
) -> dict:
    """Launch once, retain logs/state, and verify container removal."""
    start = time.monotonic()
    timed_out = False
    launch_error = None
    cleanup_error = None
    exit_code = None
    stdout_path, stderr_path = output / f"{name}.stdout", output / f"{name}.stderr"
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        try:
            child = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
            try:
                exit_code = child.wait(timeout=wall)
            except subprocess.TimeoutExpired:
                timed_out = True
                try:
                    killed = subprocess.run(
                        ["docker", "kill", container_name], capture_output=True,
                        text=True, timeout=15, check=False,
                    )
                    if killed.returncode != 0:
                        cleanup_error = killed.stderr.strip() or "docker kill failed"
                except (OSError, subprocess.SubprocessError) as error:
                    cleanup_error = f"docker kill failed: {error}"
                try:
                    exit_code = child.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    child.kill()
                    exit_code = child.wait()
        except (OSError, subprocess.SubprocessError) as error:
            launch_error = f"{type(error).__name__}: {error}"
    def docker(*args: str) -> subprocess.CompletedProcess[str] | None:
        try:
            return subprocess.run(["docker", *args], capture_output=True,
                                  text=True, timeout=15, check=False)
        except (OSError, subprocess.SubprocessError):
            return None
    inspected = docker("inspect", "--format", "{{json .State}}", container_name)
    state = None
    if inspected and inspected.returncode == 0:
        try:
            state = json.loads(inspected.stdout)
        except json.JSONDecodeError:
            pass
    removed = docker("rm", "-f", container_name) if state else None
    remaining = docker("ps", "-aq", "--filter", f"name=^{container_name}$")
    cleanup_ok = bool(remaining and remaining.returncode == 0 and not remaining.stdout.strip()
                      and (removed is None or removed.returncode == 0))
    if not cleanup_ok:
        cleanup_error = cleanup_error or "container removal unverified"
    return {
        "name": name, "command": command, "exit_code": exit_code,
        "timed_out": timed_out, "launch_error": launch_error,
        "cleanup_error": cleanup_error, "cleanup_verified": cleanup_ok,
        "container_state": state, "wall_elapsed_seconds_L0": time.monotonic() - start,
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_sha256": sha256_file(stderr_path),
    }


def run(
    object_path: Path, binary: Path, authorization_path: Path, output: Path,
    image: str = IMAGE, checkout: Path = ROOT,
) -> dict:
    checkout = checkout.resolve(strict=True)
    object_path, binary, authorization_path = (
        object_path.resolve(strict=True), binary.resolve(strict=True),
        authorization_path.resolve(strict=True),
    )
    output = output.resolve()
    if output.is_relative_to(checkout) or output.exists():
        raise ValueError("output must be a fresh directory outside the frozen checkout")
    auth_bytes = authorization_path.read_bytes()
    auth = json.loads(auth_bytes)
    error = authorization_error(auth, output, binary.parent)
    if error:
        raise ValueError(error)
    started = time.monotonic()
    deadline = started + TOTAL_CHARGE
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError("worker must be an executable Linux ELF binary")
    with binary.open("rb") as stream:
        magic = stream.read(4)
    if magic != b"\x7fELF":
        raise ValueError("worker must be an executable Linux ELF binary")
    if sha256_file(object_path) != OBJECT_SHA256:
        raise ValueError("downloaded S3 object SHA256 mismatch")
    manifest, corpus, validation = (checkout / MANIFEST, checkout / CORPUS, checkout / VALIDATION)
    source_paths = [checkout / relative for relative in SOURCES]
    commit = clean_commit(checkout)
    build_path = binary.parent / "build.json"
    build = json.loads(build_path.read_text())
    build_expected = {
        "schema": "n83.structured-orbit-linux-build/v1",
        "status": "PASS_build_only", "source_commit": commit,
        "binary_sha256": sha256_file(binary),
        "authorization_sha256": hashlib.sha256(auth_bytes).hexdigest(),
        "binary": str(binary), "conservative_charge_seconds": BUILD_CHARGE,
        "container_cleanup_verified": True,
        "source_clean": True, "source_same": True,
        "solver_search_executed": False,
        "exit_code": 0, "timed_out": False,
        "launch_error": None, "cleanup_error": None,
    }
    if any(type(build.get(key)) is not type(value) or build[key] != value
           for key, value in build_expected.items()):
        raise ValueError("Linux build receipt does not bind the frozen binary and authorization")
    build_state = build.get("container_state")
    build_elapsed = build.get("elapsed_seconds_L0")
    if (not isinstance(build_state, dict) or build_state.get("ExitCode") != 0
            or build_state.get("OOMKilled") is not False
            or isinstance(build_elapsed, bool) or not isinstance(build_elapsed, (int, float))
            or not math.isfinite(build_elapsed) or not 0 <= build_elapsed <= BUILD_CHARGE):
        raise ValueError("Linux build resource or wall receipt is invalid")
    inspected = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Id}} {{.Architecture}} {{.Os}}"],
        capture_output=True, text=True, check=True, timeout=15,
    )
    image_fields = inspected.stdout.strip().split()
    if (len(image_fields) != 3 or not image_fields[0].startswith("sha256:")
            or not HEX64.fullmatch(image_fields[0][7:])
            or image_fields[1:] != ["arm64", "linux"]):
        raise ValueError("local Linux arm64 image identity mismatch")
    if build.get("container_image_id") != image_fields[0]:
        raise ValueError("Linux build and search image identity differ")
    input_hashes = {str(path): sha256_file(path) for path in
                    [object_path, binary, build_path, manifest, corpus, validation, authorization_path, *source_paths]}
    if clean_commit(checkout) != commit:
        raise ValueError("source changed while freezing the gate")
    output.mkdir(parents=True, exist_ok=False)
    token = uuid.uuid4().hex[:16]
    config = {
        "schema": "n83.structured-orbit-search-config/v1", "study": STUDY,
        "source_commit": commit, "object_sha256": OBJECT_SHA256,
        "binary_sha256": input_hashes[str(binary)], "input_sha256": input_hashes,
        "build_receipt_sha256": input_hashes[str(build_path)],
        "authorization_sha256": input_hashes[str(authorization_path)],
        "container_image_id": image_fields[0], "container_architecture": image_fields[1],
        "verify_wall_seconds": VERIFY_WALL, "search_wall_seconds": SEARCH_WALL,
        "worker_wall_seconds": WORKER_WALL, "conservative_charge_seconds": TOTAL_CHARGE,
        "cleanup_reserve_seconds": CLEANUP_RESERVE,
        "memory_mib": MEMORY_MIB, "cpu_limit": 1,
        "max_trials": TRIALS, "conflict_budget": CONFLICTS, "max_models": MODELS,
        "container_token": token, "selected_best_total_runtime": None,
    }
    write_new_json(output / "config.json", config)
    config_hash = sha256_file(output / "config.json")
    stages = []
    status = "PRODUCER_FAILURE_preflight"
    import_passed = False
    search_result = None
    for name, wall in (("verify", VERIFY_WALL), ("search", SEARCH_WALL)):
        if time.monotonic() + wall + CLEANUP_RESERVE > deadline:
            status = "UNKNOWN_budget_before_stage"
            break
        container_name = f"n83-orbit-{name}-{token}"
        command = [
            "docker", "run", "--pull", "never", "--platform", "linux/arm64",
            "--name", container_name,
            "--network", "none", "--cpus", "1", "--memory", "4096m",
            "--memory-swap", "4096m", "--pids-limit", "64", "--cap-drop", "ALL",
            "--security-opt", "no-new-privileges", "--read-only",
            "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m",
            "--user", f"{os.getuid()}:{os.getgid()}",
            "--mount", f"type=bind,src={binary},dst=/worker,readonly",
            "--mount", f"type=bind,src={object_path},dst=/object.jsonl.gz,readonly",
            "--mount", f"type=bind,src={manifest},dst=/manifest.json,readonly",
            "--mount", f"type=bind,src={corpus},dst=/corpus.json,readonly",
            "--mount", f"type=bind,src={validation},dst=/validation.json,readonly",
            "--mount", f"type=bind,src={output},dst=/out",
            image_fields[0], "/worker",
        ]
        if name == "verify":
            command.extend(["verify", "/manifest.json", "/object.jsonl.gz"])
        else:
            command.extend([
                "run", "/manifest.json", "/object.jsonl.gz", "/corpus.json",
                "/validation.json", str(WORKER_WALL), str(TRIALS),
                str(CONFLICTS), str(MODELS), "/out/run",
            ])
        stage_receipt = stage(name=name, command=command, wall=wall, output=output,
                              container_name=container_name)
        stages.append(stage_receipt)
        if stage_receipt["cleanup_error"]:
            status = "PRODUCER_FAILURE_cleanup"
            break
        if stage_receipt["timed_out"]:
            status = "UNKNOWN_wall_cap"
            break
        if stage_receipt["exit_code"] == 124:
            try:
                cap = json.loads((output / "run" / "cap.json").read_text())
            except (OSError, json.JSONDecodeError):
                cap = None
            status = ("UNKNOWN_wall_cap" if valid_worker_cap(cap)
                      else "PRODUCER_FAILURE_worker_exit")
            break
        if stage_receipt["exit_code"] in (134, 137) or (
            isinstance(stage_receipt["container_state"], dict)
            and stage_receipt["container_state"].get("OOMKilled") is True
        ):
            status = "UNKNOWN_resource_cap"
            break
        if stage_receipt["launch_error"] or stage_receipt["exit_code"] != 0:
            status = "PRODUCER_FAILURE_worker_exit"
            break
        state = stage_receipt["container_state"]
        if (not isinstance(state, dict) or state.get("ExitCode") != 0
                or state.get("OOMKilled") is not False):
            status = "PRODUCER_FAILURE_container_state"
            break
        if name == "verify":
            try:
                verified = json.loads((output / "verify.stdout").read_text())
            except (OSError, json.JSONDecodeError):
                verified = None
            error = verify_receipt_error(verified)
            if error:
                status = "PRODUCER_FAILURE_import_receipt"
                break
            import_passed = True
            status = "PASS_import_only"
        else:
            try:
                worker_config = json.loads((output / "run" / "config.json").read_text())
                worker_result = json.loads((output / "run" / "result.json").read_text())
            except (OSError, json.JSONDecodeError):
                worker_config = worker_result = None
            error = search_receipt_error(worker_config, worker_result)
            if error:
                status = "PRODUCER_FAILURE_search_receipt"
            else:
                search_result = worker_result
                status = worker_result["status"]
    try:
        changed = any(sha256_file(Path(path)) != digest for path, digest in input_hashes.items())
        changed = changed or clean_commit(checkout) != commit or sha256_file(output / "config.json") != config_hash
    except (OSError, ValueError):
        changed = True
    if changed:
        status = "PRODUCER_FAILURE_frozen_input_changed"
    elapsed = time.monotonic() - started
    if elapsed > TOTAL_CHARGE:
        status = "PRODUCER_FAILURE_budget_overshoot"
    outer = {
        "schema": "n83.structured-orbit-search-outer/v1", "status": status,
        "config_sha256": config_hash, "source_commit": commit,
        "authorization_sha256": input_hashes[str(authorization_path)],
        "stages": stages, "frozen_input_changed": changed,
        "import_passed": import_passed,
        "search_stage_launched": len(stages) == 2 and stages[1]["launch_error"] is None,
        "solver_search_executed": (search_result["sat_calls"] > 0 if search_result else
                                   None if len(stages) == 2 else False),
        "elapsed_supervisor_seconds_L0": elapsed,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    write_new_json(output / "outer.json", outer)
    write_new_json(output / "budget.json", {
        "schema": "n83.structured-orbit-search-budget/v1",
        "authorization_sha256": input_hashes[str(authorization_path)],
        "prior_pilot_charged_seconds": "7200",
        "additional_authorized_seconds": auth["additional_authorized_seconds"],
        "build_conservative_charge_seconds": str(BUILD_CHARGE),
        "gate_conservative_charge_seconds": str(TOTAL_CHARGE),
        "conservative_charge_seconds": str(TOTAL_ADDITIONAL_CHARGE),
        "build_elapsed_seconds_L0": build_elapsed,
        "gate_elapsed_seconds_L0": elapsed,
        "combined_observed_seconds_L0": build_elapsed + elapsed,
        "budget_overshoot": elapsed > TOTAL_CHARGE,
        "remaining_additional_seconds": str(Decimal(auth["additional_authorized_seconds"]) - TOTAL_ADDITIONAL_CHARGE),
        "round_status": status, "selected_best_total_runtime": None,
    })
    return outer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--object", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--image", default=IMAGE)
    args = parser.parse_args()
    outer = run(args.object, args.binary, args.authorization, args.output, args.image)
    print(json.dumps(outer, sort_keys=True))
    return 1 if outer["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
