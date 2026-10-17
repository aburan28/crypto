#!/usr/bin/env python3
"""One authorized, conservatively charged N83 mixed-orbit build/import gate."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time
import uuid

from compact_probe_supervisor import clean_commit, sha256_file, write_new_json
from structured_orbit_search_supervisor import stage
from mixed_orbit_supervisor import (
    IMAGE, SMALL_MANIFEST, ENVELOPE_MANIFEST, SMALL_SHA256, ENVELOPE_SHA256,
    STUDY, image_identity, is_linux_elf, verify_receipt_error,
)


ROOT = Path(__file__).resolve().parents[2]
BUILD_WALL = 450
BUILD_CHARGE = 600
VERIFY_WALL = 240
VERIFY_CHARGE = 240
RECEIPT_RESERVE = 60
TOTAL_CHARGE = BUILD_CHARGE + VERIFY_CHARGE + RECEIPT_RESERVE
SOURCES = (
    "examples/koblitz_n83_mixed_orbit_worker.rs",
    "src/cryptanalysis/koblitz_large_prime.rs",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/binary_semaev_chain_sat.rs",
    f"research/{STUDY}/mixed_orbit_supervisor.py",
    f"research/{STUDY}/mixed_orbit_import_900_supervisor.py",
    f"research/{STUDY}/structured_orbit_search_supervisor.py",
    f"research/{STUDY}/compact_probe_supervisor.py",
    SMALL_MANIFEST, ENVELOPE_MANIFEST,
)


def authorization_error(value: object, output: Path) -> str | None:
    if not isinstance(value, dict):
        return "authorization must be an object"
    expected = {
        "schema": "n83.additional-local-budget/v1", "study": STUDY,
        "purpose": "mixed_orbit_build_import_900_gate",
        "prior_pilot_charged_seconds": "11700",
        "additional_authorized_seconds": "900",
        "output_dir": str(output.resolve()),
    }
    for key, item in expected.items():
        if value.get(key) != item:
            return f"authorization {key} mismatch"
    if not isinstance(value.get("user_approval_reference"), str) or not value["user_approval_reference"].strip():
        return "authorization user approval reference missing"
    return None


def source_hashes(checkout: Path) -> dict[str, str]:
    return {relative: sha256_file(checkout / relative) for relative in SOURCES}


def gate(authorization: Path, output: Path, small: Path, envelope: Path,
         checkout: Path = ROOT) -> dict:
    checkout = checkout.resolve(strict=True)
    output = output.resolve()
    if output.exists() or output.is_relative_to(checkout):
        raise ValueError("output must be a fresh directory outside the checkout")
    authorization = authorization.resolve(strict=True)
    small, envelope = small.resolve(strict=True), envelope.resolve(strict=True)
    raw_auth = authorization.read_bytes()
    error = authorization_error(json.loads(raw_auth), output)
    if error:
        raise ValueError(error)
    commit = clean_commit(checkout)
    sources = source_hashes(checkout)
    if sha256_file(small) != SMALL_SHA256 or sha256_file(envelope) != ENVELOPE_SHA256:
        raise ValueError("retained S3 object SHA-256 mismatch")
    image_identity()
    registry = Path(os.environ.get("CARGO_HOME", str(Path.home() / ".cargo"))) / "registry"
    if not registry.is_dir():
        raise ValueError("offline Cargo registry unavailable")
    started = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    target = output / "target"
    target.mkdir()
    auth_sha = hashlib.sha256(raw_auth).hexdigest()
    inputs = {str(path): sha256_file(path) for path in (authorization, small, envelope)}
    result = {"schema": "n83.mixed-orbit-import-900-outer/v1",
              "status": "PRODUCER_FAILURE_preflight", "source_commit": commit,
              "source_sha256": sources, "input_sha256": inputs,
              "authorization_sha256": auth_sha, "container_image_id": IMAGE,
              "stages": [], "import_passed": False,
              "total_index_calculus_runtime_ms": None,
              "selected_best_total_runtime": None}
    binary = output / "worker"
    try:
        build_name = f"n83-mixed-import-build-{uuid.uuid4().hex[:16]}"
        build_command = [
            "docker", "run", "--pull", "never", "--platform", "linux/arm64",
            "--name", build_name, "--network", "none", "--cpus", "4",
            "--memory", "8192m", "--memory-swap", "8192m", "--pids-limit", "256",
            "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
            "--user", f"{os.getuid()}:{os.getgid()}", "--tmpfs", "/cargo:rw,mode=1777",
            "--mount", f"type=bind,src={registry},dst=/cargo/registry,readonly",
            "--mount", f"type=bind,src={checkout},dst=/src,readonly",
            "--mount", f"type=bind,src={target},dst=/target",
            "--env", "CARGO_HOME=/cargo", "--env", "CARGO_TARGET_DIR=/target",
            "--workdir", "/src", IMAGE, "cargo", "build", "--release", "--offline",
            "--locked", "--no-default-features", "--example", "koblitz_n83_mixed_orbit_worker",
        ]
        built = stage(name="build", command=build_command, wall=BUILD_WALL,
                      output=output, container_name=build_name)
        result["stages"].append(built)
        state = built.get("container_state")
        artifact = target / "release/examples/koblitz_n83_mixed_orbit_worker"
        if not built.get("cleanup_verified"):
            result["status"] = "PRODUCER_FAILURE_cleanup"
        elif built.get("timed_out"):
            result["status"] = "UNKNOWN_build_wall_cap"
        elif (built.get("exit_code") != 0 or built.get("launch_error")
              or not isinstance(state, dict) or state.get("ExitCode") != 0
              or state.get("OOMKilled") is not False or not is_linux_elf(artifact)):
            result["status"] = "PRODUCER_FAILURE_build"
        else:
            shutil.copyfile(artifact, binary)
            binary.chmod(0o755)
            result["binary_sha256"] = sha256_file(binary)
            verify_name = f"n83-mixed-import-verify-{uuid.uuid4().hex[:16]}"
            verify_command = [
                "docker", "run", "--pull", "never", "--platform", "linux/arm64",
                "--name", verify_name, "--network", "none", "--cpus", "1",
                "--memory", "4096m", "--memory-swap", "4096m", "--pids-limit", "64",
                "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
                "--read-only", "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m",
                "--user", f"{os.getuid()}:{os.getgid()}",
                "--mount", f"type=bind,src={binary},dst=/worker,readonly",
                "--mount", f"type=bind,src={small},dst=/small.jsonl.gz,readonly",
                "--mount", f"type=bind,src={envelope},dst=/envelope.jsonl.gz,readonly",
                "--mount", f"type=bind,src={checkout / SMALL_MANIFEST},dst=/small-manifest.json,readonly",
                "--mount", f"type=bind,src={checkout / ENVELOPE_MANIFEST},dst=/envelope-manifest.json,readonly",
                IMAGE, "/worker", "verify", "/small-manifest.json", "/small.jsonl.gz",
                "/envelope-manifest.json", "/envelope.jsonl.gz",
            ]
            verified = stage(name="verify", command=verify_command, wall=VERIFY_WALL,
                             output=output, container_name=verify_name)
            result["stages"].append(verified)
            state = verified.get("container_state")
            if not verified.get("cleanup_verified"):
                result["status"] = "PRODUCER_FAILURE_cleanup"
            elif verified.get("timed_out"):
                result["status"] = "UNKNOWN_import_wall_cap"
            elif state and state.get("OOMKilled") is True:
                result["status"] = "UNKNOWN_import_resource_cap"
            elif (verified.get("exit_code") != 0 or verified.get("launch_error")
                  or not isinstance(state, dict) or state.get("ExitCode") != 0):
                result["status"] = "PRODUCER_FAILURE_import"
            else:
                try:
                    receipt = json.loads((output / "verify.stdout").read_text())
                except (OSError, UnicodeError, json.JSONDecodeError):
                    receipt = None
                result["import_error"] = verify_receipt_error(receipt)
                result["status"] = ("PASS_import_only" if result["import_error"] is None
                                    else "PRODUCER_FAILURE_import_receipt")
                result["import_passed"] = result["status"] == "PASS_import_only"
    except Exception as error:
        result["status"] = "PRODUCER_FAILURE_supervisor_exception"
        result["exception"] = f"{type(error).__name__}: {error}"
    try:
        frozen = (clean_commit(checkout) == commit and source_hashes(checkout) == sources
                  and all(sha256_file(Path(path)) == digest for path, digest in inputs.items()))
    except (OSError, ValueError):
        frozen = False
    result["frozen_inputs_unchanged"] = frozen
    if not frozen:
        result["status"] = "PRODUCER_FAILURE_frozen_input_changed"
        result["import_passed"] = False
    result["elapsed_supervisor_seconds_L0"] = time.monotonic() - started
    if result["elapsed_supervisor_seconds_L0"] > TOTAL_CHARGE:
        result["status"] = "PRODUCER_FAILURE_budget_overshoot"
        result["import_passed"] = False
    result["conservative_charge_seconds"] = (BUILD_CHARGE if result["stages"] else 0) + (
        VERIFY_CHARGE if len(result["stages"]) > 1 else 0) + RECEIPT_RESERVE
    write_new_json(output / "outer.json", result)
    write_new_json(output / "budget.json", {
        "schema": "n83.mixed-orbit-import-900-budget/v1",
        "authorization_sha256": auth_sha, "prior_pilot_charged_seconds": "11700",
        "additional_authorized_seconds": "900",
        "conservative_charge_seconds": result["conservative_charge_seconds"],
        "remaining_additional_seconds": str(TOTAL_CHARGE - result["conservative_charge_seconds"]),
        "status": result["status"], "selected_best_total_runtime": None,
    })
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--small-object", type=Path, required=True)
    parser.add_argument("--envelope-object", type=Path, required=True)
    args = parser.parse_args()
    result = gate(args.authorization, args.output_dir,
                  args.small_object, args.envelope_object)
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "PASS_import_only" else 1


if __name__ == "__main__":
    raise SystemExit(main())
