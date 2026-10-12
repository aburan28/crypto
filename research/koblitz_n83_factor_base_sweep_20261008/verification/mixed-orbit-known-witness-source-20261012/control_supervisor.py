#!/usr/bin/env python3
"""One-shot, budget-gated Linux supervisor for the archived N83 control.

Preparing or testing this source launches no retained SAT search. Build and
run require a separate, run-specific authorization after the 15,000-second
cumulative ledger. Every launched stage is charged its full declared cap.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
import uuid


HERE = Path(__file__).resolve().parent
STUDY_DIR = HERE.parents[1]
ROOT = STUDY_DIR.parents[1]
sys.path.insert(0, str(STUDY_DIR))
from chain_search_supervisor import parse_point  # noqa: E402
from compact_probe_supervisor import clean_commit, sha256_file, write_new_json  # noqa: E402
from mixed_orbit_supervisor import (  # noqa: E402
    ENVELOPE_MANIFEST, IMAGE, SMALL_MANIFEST, image_identity, is_linux_elf,
    replay_equation,
)
from structured_orbit_search_supervisor import stage  # noqa: E402
from known_witness import collect  # noqa: E402


STUDY = "koblitz_n83_factor_base_sweep_20261008"
FIXTURE = f"research/{STUDY}/verification/mixed-orbit-known-witness-source-20261012/fixture.json"
PRIOR_AUDIT = f"research/{STUDY}/verification/mixed-orbit-search-15000-20261012/audit.json"
PRIOR_AUDIT_SHA = "32420a62bb6ccf8f9df7f03dea9040e518606740968110dc8c0ac47cb50cafd5"
SMALL_SHA = "eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380"
ENVELOPE_SHA = "aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a"
BUILD_CAP, VERIFY_CAP, SEARCH_CAP, CLEANUP_RESERVE = 300, 120, 900, 180
BUILD_WALL, VERIFY_WALL, SEARCH_WALL, WORKER_WALL = 260, 90, 850, 840
TOTAL_CHARGE = BUILD_CAP + VERIFY_CAP + SEARCH_CAP + CLEANUP_RESERVE
RUN_CAP = VERIFY_CAP + SEARCH_CAP + CLEANUP_RESERVE
CONFLICTS, MODELS = 5000, 32
SOURCES = (
    "examples/koblitz_n83_mixed_orbit_worker.rs",
    "src/cryptanalysis/koblitz_large_prime.rs",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/binary_semaev_chain_sat.rs",
    f"research/{STUDY}/chain_search_supervisor.py",
    f"research/{STUDY}/mixed_orbit_supervisor.py",
    f"research/{STUDY}/structured_orbit_search_supervisor.py",
    f"research/{STUDY}/compact_probe_supervisor.py",
    f"research/{STUDY}/verification/mixed-orbit-known-witness-source-20261012/control_supervisor.py",
    f"research/{STUDY}/verification/mixed-orbit-known-witness-source-20261012/known_witness.py",
    FIXTURE, PRIOR_AUDIT, SMALL_MANIFEST, ENVELOPE_MANIFEST,
)


def fresh_outside(path: Path, checkout: Path) -> Path:
    path = path.resolve()
    if path.is_relative_to(checkout) or path.exists():
        raise ValueError("output must be fresh and outside the checkout")
    return path


def authorization_error(record: object, build_dir: Path, output_dir: Path) -> str | None:
    if not isinstance(record, dict):
        return "authorization is not an object"
    expected = {
        "schema": "n83.additional-local-budget/v3",
        "study": STUDY,
        "purpose": "one_mixed_orbit_known_witness_control",
        "prior_total_charged_seconds": "15000",
        "additional_authorized_seconds": str(TOTAL_CHARGE),
        "build_dir": str(build_dir.resolve()),
        "output_dir": str(output_dir.resolve()),
    }
    for key, value in expected.items():
        if record.get(key) != value:
            return f"authorization {key} mismatch"
    approval = record.get("user_approval_reference")
    if not isinstance(approval, str) or not approval.strip():
        return "authorization lacks explicit approval reference"
    return None


def load_auth(path: Path, build_dir: Path, output_dir: Path) -> tuple[dict, str]:
    raw = path.resolve(strict=True).read_bytes()
    record = json.loads(raw)
    error = authorization_error(record, build_dir, output_dir)
    if error:
        raise ValueError(error)
    return record, hashlib.sha256(raw).hexdigest()


def source_hashes(checkout: Path) -> dict[str, str]:
    return {relative: sha256_file(checkout / relative) for relative in SOURCES}


def prior_audit(checkout: Path) -> None:
    if sha256_file(checkout / PRIOR_AUDIT) != PRIOR_AUDIT_SHA:
        raise ValueError("prior 15,000-second audit hash differs")
    row = json.loads((checkout / PRIOR_AUDIT).read_text())
    if (row.get("status") != "PASS"
            or row.get("cumulative_conservative_charge_seconds") != 15000
            or row.get("small_object_sha256") != SMALL_SHA
            or row.get("envelope_object_sha256") != ENVELOPE_SHA):
        raise ValueError("prior 15,000-second audit differs")


def build(authorization: Path, build_dir: Path, output_dir: Path,
          checkout: Path = ROOT) -> dict:
    checkout = checkout.resolve(strict=True)
    build_dir = fresh_outside(build_dir, checkout)
    output_dir = fresh_outside(output_dir, checkout)
    _, auth_hash = load_auth(authorization, build_dir, output_dir)
    commit = clean_commit(checkout)
    prior_audit(checkout)
    sources = source_hashes(checkout)
    image_identity()
    registry = Path(os.environ.get("CARGO_HOME", str(Path.home() / ".cargo"))) / "registry"
    if not registry.is_dir():
        raise ValueError("offline Cargo registry unavailable")
    build_dir.mkdir(parents=True, exist_ok=False)
    target = build_dir / "target"
    target.mkdir()
    name = f"n83-control-build-{uuid.uuid4().hex[:16]}"
    command = [
        "docker", "run", "--pull", "never", "--platform", "linux/arm64",
        "--name", name, "--network", "none", "--cpus", "4",
        "--memory", "8192m", "--memory-swap", "8192m", "--pids-limit", "256",
        "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--tmpfs", "/cargo:rw,mode=1777",
        "--mount", f"type=bind,src={registry},dst=/cargo/registry,readonly",
        "--mount", f"type=bind,src={checkout},dst=/src,readonly",
        "--mount", f"type=bind,src={target},dst=/target",
        "--env", "CARGO_HOME=/cargo", "--env", "CARGO_TARGET_DIR=/target",
        "--workdir", "/src", IMAGE, "cargo", "build", "--release", "--offline",
        "--locked", "--no-default-features", "--example", "koblitz_n83_mixed_orbit_worker",
    ]
    started = time.monotonic()
    launched = stage(name="build", command=command, wall=BUILD_WALL,
                     output=build_dir, container_name=name)
    artifact = target / "release/examples/koblitz_n83_mixed_orbit_worker"
    binary = build_dir / "worker"
    state = launched.get("container_state")
    status = "PRODUCER_FAILURE_build"
    if not launched.get("cleanup_verified"):
        status = "PRODUCER_FAILURE_cleanup"
    elif launched.get("timed_out"):
        status = "UNKNOWN_build_wall_cap"
    elif (launched.get("exit_code") == 0 and not launched.get("launch_error")
          and isinstance(state, dict) and state.get("ExitCode") == 0
          and state.get("OOMKilled") is False and is_linux_elf(artifact)):
        try:
            shutil.copyfile(artifact, binary)
            binary.chmod(0o755)
            status = "PASS_build_only"
        except OSError:
            status = "PRODUCER_FAILURE_binary"
    try:
        frozen = clean_commit(checkout) == commit and source_hashes(checkout) == sources
    except (OSError, ValueError):
        frozen = False
    if not frozen:
        status = "PRODUCER_FAILURE_source_changed"
    elapsed = time.monotonic() - started
    if elapsed > BUILD_CAP:
        status = "PRODUCER_FAILURE_budget_overshoot"
    receipt = {
        "schema": "n83.known-witness-build/v1", "status": status,
        "source_commit": commit, "source_sha256": sources,
        "authorization_sha256": auth_hash, "container_image_id": IMAGE,
        "binary_sha256": sha256_file(binary) if binary.is_file() else None,
        "binary": str(binary) if binary.is_file() else None,
        "stage": launched, "elapsed_seconds_L0": elapsed,
        "conservative_charge_seconds": BUILD_CAP,
        "solver_search_executed": False,
    }
    write_new_json(build_dir / "build.json", receipt)
    return receipt


def control_receipt_error(config: object, result: object,
                          fixture: dict, small_object: Path,
                          envelope_object: Path) -> str | None:
    if not isinstance(config, dict) or not isinstance(result, dict):
        return "control receipts are not objects"
    expected_config = {
        "schema": "n83.mixed-orbit-known-witness-config/v1",
        "scope": "known_answer_control_only", "study": STUDY,
        "small_summands": 3, "residual_summands": 2,
        "conflict_budget": CONFLICTS, "max_models": MODELS,
        "budget_seconds": WORKER_WALL,
        "memory_cgroup_limit_bytes": 4 * 1024**3,
        "memory_cgroup_swap_limit_bytes": 0,
    }
    for key, value in expected_config.items():
        if type(config.get(key)) is not type(value) or config[key] != value:
            return f"control config {key} mismatch"
    if (result.get("schema") != "n83.mixed-orbit-known-witness-result/v1"
            or result.get("scope") != "known_answer_control_only"
            or result.get("target") != fixture["equation"]["target"]
            or result.get("natural_relations") is not None
            or result.get("total_index_calculus_runtime_ms") is not None
            or result.get("selected_best_total_runtime") is not None):
        return "control result scope or target mismatch"
    stats = result.get("search")
    if not isinstance(stats, dict):
        return "control search stats absent"
    for key in ("solver_calls", "conflicts", "models", "invalid_models", "patterns_attempted"):
        if type(stats.get(key)) is not int or stats[key] < 0:
            return f"control search {key} malformed"
    if stats["conflicts"] > CONFLICTS or stats["models"] > MODELS:
        return "control search exceeded cap"
    if type(stats.get("exhausted")) is not bool or type(stats.get("refuted")) is not bool:
        return "control search disposition malformed"
    status, candidate = result.get("status"), result.get("candidate")
    if status == "PASS_known_witness_solver_access":
        if (not isinstance(candidate, dict) or stats["models"] == 0
                or stats["exhausted"] or stats["refuted"]):
            return "control success lacks solver candidate"
        indices = candidate.get("small_point_indices")
        residual = candidate.get("residual_points")
        if (candidate.get("a") != "0" or candidate.get("b") != "1"
                or not isinstance(indices, list) or len(indices) != 3
                or any(type(i) is not int or not 0 <= i < 196212 for i in indices)
                or not isinstance(residual, list) or len(residual) != 2):
            return "control candidate shape differs"
        needed = set(indices)
        small = {}
        residual_set = set()
        with gzip.open(small_object, "rt") as stream:
            next(stream)
            for i, line in enumerate(stream):
                if i in needed:
                    small[i] = parse_point(json.loads(line)["point"])
        with gzip.open(envelope_object, "rt") as stream:
            next(stream)
            for i, line in enumerate(stream):
                if i >= 196212:
                    residual_set.add(parse_point(json.loads(line)["point"]))
        try:
            residual_points = [parse_point(value) for value in residual]
        except (ValueError, TypeError):
            return "control candidate residual malformed"
        if (set(small) != needed or None in small.values()
                or any(point not in residual_set for point in residual_points)
                or not replay_equation("0", "1", [small[i] for i in indices],
                                       residual_points, parse_point(fixture["equation"]["target"]))):
            return "control candidate failed independent archive replay"
    elif status == "UNKNOWN_preflight_cap":
        if (candidate is not None or stats["solver_calls"] != 0
                or stats["patterns_attempted"] != 0 or not stats["exhausted"]):
            return "preflight cap hides candidate or solver call"
    elif status == "UNKNOWN_search_cap":
        if candidate is not None or stats["solver_calls"] == 0 or not stats["exhausted"]:
            return "search cap lacks solver call or hides candidate"
    elif status == "PRODUCER_FAILURE_known_witness_refuted":
        if candidate is not None or stats.get("refuted") is not True:
            return "known-witness refutation receipt inconsistent"
    elif status == "PRODUCER_FAILURE_invalid_model":
        if candidate is not None or stats["invalid_models"] == 0:
            return "invalid-model receipt inconsistent"
    else:
        return "control status inadmissible"
    for key in ("base_import_elapsed_ms_L0", "worker_wall_ms_L0"):
        value = result.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            return f"control {key} malformed"
    peak = result.get("memory_cgroup_peak_bytes")
    if type(peak) is not int or not 0 < peak <= 4 * 1024**3:
        return "control cgroup peak malformed"
    return None


def valid_worker_cap(record: object) -> bool:
    if not isinstance(record, dict):
        return False
    elapsed = record.get("elapsed_process_wall_ms")
    return (record.get("schema") == "n83.mixed-orbit-known-witness-cap/v1"
            and record.get("status") == "UNKNOWN_budget"
            and type(record.get("budget_seconds")) is int
            and record["budget_seconds"] == WORKER_WALL
            and isinstance(elapsed, (int, float)) and not isinstance(elapsed, bool)
            and math.isfinite(elapsed) and elapsed >= 0)


def run(authorization: Path, build_dir: Path, output_dir: Path,
        small_object: Path, envelope_object: Path, checkout: Path = ROOT) -> dict:
    checkout = checkout.resolve(strict=True)
    build_dir = build_dir.resolve(strict=True)
    output_dir = fresh_outside(output_dir, checkout)
    small_object = small_object.resolve(strict=True)
    envelope_object = envelope_object.resolve(strict=True)
    _, auth_hash = load_auth(authorization, build_dir, output_dir)
    commit = clean_commit(checkout)
    prior_audit(checkout)
    built = json.loads((build_dir / "build.json").read_text())
    binary = build_dir / "worker"
    if (built.get("schema") != "n83.known-witness-build/v1"
            or built.get("status") != "PASS_build_only"
            or built.get("source_commit") != commit
            or built.get("source_sha256") != source_hashes(checkout)
            or built.get("authorization_sha256") != auth_hash
            or built.get("container_image_id") != IMAGE
            or built.get("binary_sha256") != sha256_file(binary)
            or built.get("conservative_charge_seconds") != BUILD_CAP
            or not is_linux_elf(binary)):
        raise ValueError("frozen build receipt mismatch")
    if sha256_file(small_object) != SMALL_SHA or sha256_file(envelope_object) != ENVELOPE_SHA:
        raise ValueError("retained archive bytes differ")
    fixture = collect(small_object, envelope_object)
    if fixture != json.loads((checkout / FIXTURE).read_text()):
        raise ValueError("known-answer fixture differs")
    image_identity()
    inputs = {str(path): sha256_file(path) for path in (
        authorization.resolve(strict=True), binary, small_object, envelope_object,
        build_dir / "build.json", *(checkout / relative for relative in SOURCES))}
    output_dir.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    config = {
        "schema": "n83.known-witness-supervisor-config/v1",
        "source_commit": commit, "authorization_sha256": auth_hash,
        "container_image_id": IMAGE, "input_sha256": inputs,
        "build_charge_seconds": BUILD_CAP, "verify_charge_seconds": VERIFY_CAP,
        "search_charge_seconds": SEARCH_CAP, "cleanup_reserve_seconds": CLEANUP_RESERVE,
        "conservative_full_charge_seconds": TOTAL_CHARGE,
        "conflicts": CONFLICTS, "models": MODELS,
        "scope": "known_answer_control_only",
    }
    write_new_json(output_dir / "config.json", config)
    config_sha = sha256_file(output_dir / "config.json")
    stages = []
    status = "PRODUCER_FAILURE_preflight"
    verify_passed = False
    replay_error = None
    token = uuid.uuid4().hex[:16]
    for name, wall in (("verify", VERIFY_WALL), ("control", SEARCH_WALL)):
        if time.monotonic() + wall + CLEANUP_RESERVE > started + RUN_CAP:
            status = "UNKNOWN_budget_before_stage"
            break
        container_name = f"n83-control-{name}-{token}"
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
            "--mount", f"type=bind,src={checkout / SMALL_MANIFEST},dst=/small-manifest.json,readonly",
            "--mount", f"type=bind,src={checkout / ENVELOPE_MANIFEST},dst=/envelope-manifest.json,readonly",
            "--mount", f"type=bind,src={checkout / FIXTURE},dst=/fixture.json,readonly",
            "--mount", f"type=bind,src={output_dir},dst=/out",
            IMAGE, "/worker", "verify-control" if name == "verify" else "control-run",
            "/small-manifest.json", "/small.jsonl.gz",
            "/envelope-manifest.json", "/envelope.jsonl.gz", "/fixture.json",
        ]
        if name == "control":
            command += [str(WORKER_WALL), str(CONFLICTS), str(MODELS), "/out/control"]
        receipt = stage(name=name, command=command, wall=wall,
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
            status = ("UNKNOWN_worker_wall_cap" if valid_worker_cap(cap)
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
                verify_passed = True
                status = "PASS_verify_only"
        else:
            try:
                worker_config = json.loads((output_dir / "control/config.json").read_text())
                worker_result = json.loads((output_dir / "control/result.json").read_text())
                replay_error = control_receipt_error(worker_config, worker_result, fixture,
                                                     small_object, envelope_object)
            except (OSError, UnicodeError, ValueError, KeyError, TypeError,
                    EOFError, json.JSONDecodeError) as error:
                replay_error = f"control receipt replay raised {type(error).__name__}: {error}"
            status = "PRODUCER_FAILURE_control_replay" if replay_error else worker_result["status"]
        if status not in ("PASS_verify_only", "PRODUCER_FAILURE_preflight"):
            break
    try:
        frozen_changed = (clean_commit(checkout) != commit
                          or any(sha256_file(Path(path)) != digest for path, digest in inputs.items())
                          or sha256_file(output_dir / "config.json") != config_sha)
    except (OSError, ValueError):
        frozen_changed = True
    if frozen_changed:
        status = "PRODUCER_FAILURE_frozen_input_changed"
    elapsed = time.monotonic() - started
    if elapsed > RUN_CAP:
        status = "PRODUCER_FAILURE_budget_overshoot"
    outer = {
        "schema": "n83.known-witness-supervisor-outer/v1",
        "status": status, "scope": "known_answer_control_only",
        "source_commit": commit, "authorization_sha256": auth_hash,
        "config_sha256": config_sha, "stages": stages,
        "verify_passed": verify_passed, "replay_error": replay_error,
        "frozen_input_changed": frozen_changed,
        "elapsed_supervisor_seconds_L0": elapsed,
        "natural_relations": None, "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    write_new_json(output_dir / "outer.json", outer)
    charge = BUILD_CAP + sum(VERIFY_CAP if row["name"] == "verify" else SEARCH_CAP
                             for row in stages) + (CLEANUP_RESERVE if stages else 0)
    write_new_json(output_dir / "budget.json", {
        "schema": "n83.known-witness-supervisor-budget/v1",
        "prior_total_charged_seconds": 15000,
        "authorization_sha256": auth_hash,
        "conservative_charge_seconds": charge,
        "cumulative_conservative_charge_seconds": 15000 + charge,
        "remaining_additional_seconds": TOTAL_CHARGE - charge,
        "round_status": status,
    })
    return outer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    for mode in ("build", "run"):
        sub = modes.add_parser(mode)
        sub.add_argument("--authorization", type=Path, required=True)
        sub.add_argument("--build-dir", type=Path, required=True)
        sub.add_argument("--output-dir", type=Path, required=True)
        if mode == "run":
            sub.add_argument("--small-object", type=Path, required=True)
            sub.add_argument("--envelope-object", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "build":
        row = build(args.authorization, args.build_dir, args.output_dir)
    else:
        row = run(args.authorization, args.build_dir, args.output_dir,
                  args.small_object, args.envelope_object)
    print(json.dumps(row, sort_keys=True))
    return 1 if row["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
