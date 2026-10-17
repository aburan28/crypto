#!/usr/bin/env python3
"""Opt-in, separately cold K600 Gray/hash first-trial gate; never auto-run."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
import uuid

import preflight

sys.path.insert(0, str(preflight.STUDY))
import receipt_audit as audit  # noqa: E402


HERE = Path(__file__).resolve().parent
ROOT = preflight.ROOT
PLAN = HERE / "plan.json"
PRIOR = audit.PRIOR
STAGE_CAP = preflight.CHARGE_PER_ARM
WALL = {"build": 450, "import": 300, "search": 1150}
ARM_CAP = sum(STAGE_CAP.values())
TOTAL_CAP = 2 * ARM_CAP


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def authorization_error(record: object, output: Path, plan_sha256: str) -> str | None:
    if not isinstance(record, dict):
        return "paired budget record is not an object"
    expected = {
        "schema": "n83.additional-local-budget/v7",
        "study": "koblitz_n83_factor_base_sweep_20261008",
        "purpose": "k600_gray_hash_paired_first_trial",
        "prior_conservative_charged_seconds": "23100",
        "additional_authorized_seconds": "4800",
        "plan_sha256": plan_sha256,
        "output_dir": str(output.resolve()),
    }
    for key, value in expected.items():
        if record.get(key) != value:
            return f"paired budget {key} mismatch"
    approval = record.get("user_approval_reference")
    if not isinstance(approval, str) or not approval.strip():
        return "paired budget lacks explicit approval reference"
    return None


def source_hashes(checkout: Path) -> dict[str, str]:
    paths = (
        "examples/koblitz_n83_mixed_orbit_worker.rs", "Cargo.toml", "Cargo.lock",
        str(preflight.PILOT.relative_to(ROOT) / "manifest.json"),
        str(preflight.PILOT.relative_to(ROOT) / "probe-corpus.json"),
        str(preflight.LIVE.relative_to(ROOT) / "audit.json"),
        str(preflight.LIVE.relative_to(ROOT) / "audit.py"),
        str(preflight.SCREEN.relative_to(ROOT)),
        str(preflight.CAPACITY.relative_to(ROOT)),
        str((HERE / "preflight.py").relative_to(ROOT)),
        str((HERE / "receipt_audit.py").relative_to(ROOT)),
        str((HERE / "supervisor.py").relative_to(ROOT)),
        str(PLAN.relative_to(ROOT)),
        str((preflight.STUDY / "k600_mixed_orbit_supervisor.py").relative_to(ROOT)),
        str((preflight.STUDY / "chain_search_supervisor.py").relative_to(ROOT)),
        str((preflight.STUDY / "structured_orbit_search_supervisor.py").relative_to(ROOT)),
        str((preflight.STUDY / "compact_probe_supervisor.py").relative_to(ROOT)),
    )
    sources = {relative: PRIOR.sha256_file(checkout / relative) for relative in paths}
    sources["src:git-tree"] = preflight.rust_source_tree_oid(checkout)
    return sources


def stage_command(name: str, arm: dict, checkout: Path, archives: Path,
                  build: Path, run: Path, binary: Path, registry: Path,
                  image: str, token: str) -> tuple[list[str], Path, str]:
    container = f"n83-paired-{arm['role']}-{name}-{token}"
    if name == "build":
        target = build / "target"
        command = [
            "docker", "run", "--pull", "never", "--platform", "linux/arm64",
            "--name", container, "--network", "none", "--cpus", "4",
            "--memory", "8192m", "--memory-swap", "8192m", "--pids-limit", "256",
            "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
            "--user", f"{os.getuid()}:{os.getgid()}",
            "--tmpfs", "/cargo:rw,mode=1777",
            "--mount", f"type=bind,src={registry},dst=/cargo/registry,readonly",
            "--mount", f"type=bind,src={checkout},dst=/src,readonly",
            "--mount", f"type=bind,src={target},dst=/target",
            "--env", "CARGO_HOME=/cargo", "--env", "CARGO_TARGET_DIR=/target",
            "--workdir", "/src", image, "cargo", "build", "--release", "--offline",
            "--locked", "--no-default-features", "--example", "koblitz_n83_mixed_orbit_worker",
        ]
        return command, build, container
    base = archives / f"{arm['base_object_blake3']}.jsonl.gz"
    pool = archives / f"{arm['pool_object_blake3']}.jsonl.gz"
    command = [
        "docker", "run", "--pull", "never", "--platform", "linux/arm64",
        "--name", container, "--network", "none", "--cpus", "1",
        "--memory", "4096m", "--memory-swap", "4096m", "--pids-limit", "64",
        "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
        "--read-only", "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--mount", f"type=bind,src={binary},dst=/worker,readonly",
        "--mount", f"type=bind,src={base},dst=/base.jsonl.gz,readonly",
        "--mount", f"type=bind,src={pool},dst=/pool.jsonl.gz,readonly",
        "--mount", f"type=bind,src={checkout / preflight.PILOT.relative_to(ROOT) / 'manifest.json'},dst=/pilot-manifest.json,readonly",
        "--mount", f"type=bind,src={checkout / preflight.PILOT.relative_to(ROOT) / 'probe-corpus.json'},dst=/corpus.json,readonly",
        "--mount", f"type=bind,src={run},dst=/out",
        image,
    ]
    command.extend(arm["worker_verify_argv"] if name == "import" else arm["worker_search_argv"])
    return command, run, container


def stage_issue(receipt: dict, name: str, run: Path) -> str | None:
    if not receipt.get("cleanup_verified"):
        return "PRODUCER_FAILURE_cleanup"
    if receipt.get("timed_out"):
        return "UNKNOWN_wall_cap"
    state = receipt.get("container_state")
    if isinstance(state, dict) and state.get("OOMKilled") is True:
        return "UNKNOWN_resource_cap"
    if receipt.get("exit_code") == 124 and name == "search":
        try:
            cap = json.loads((run / "search/cap.json").read_text())
        except (OSError, UnicodeError, json.JSONDecodeError):
            cap = None
        if (isinstance(cap, dict) and cap.get("schema") == "n83.cross-stream-pilot-partial-cap/v1"
                and cap.get("status") == "UNKNOWN_budget" and cap.get("budget_seconds") == 1100
                and PRIOR.finite_nonnegative(cap.get("elapsed_process_wall_ms"))):
            return "UNKNOWN_worker_wall_cap"
    if receipt.get("exit_code") in (134, 137):
        return "UNKNOWN_resource_cap"
    if receipt.get("launch_error") or receipt.get("exit_code") != 0 or not isinstance(state, dict) or state.get("ExitCode") != 0:
        return "PRODUCER_FAILURE_worker_exit"
    return None


def run(authorization: Path, output: Path, archives: Path, checkout: Path = ROOT) -> dict:
    checkout = checkout.resolve(strict=True)
    output = PRIOR.outside_fresh_checkout(output, checkout)
    authorization = authorization.resolve(strict=True)
    archives = archives.resolve(strict=True)
    saved = PLAN.read_bytes()
    plan_sha256 = hashlib.sha256(saved).hexdigest()
    record_bytes = authorization.read_bytes()
    error = authorization_error(json.loads(record_bytes), output, plan_sha256)
    if error:
        raise ValueError(error)
    plan = preflight.make_plan(archives)
    require(saved == (json.dumps(plan, indent=2, sort_keys=True) + "\n").encode(),
            "paired plan replay changed")
    commit = PRIOR.clean_commit(checkout)
    sources = source_hashes(checkout)
    image = PRIOR.image_identity()
    registry = Path(os.environ.get("CARGO_HOME", str(Path.home() / ".cargo"))) / "registry"
    require(registry.is_dir(), "offline Cargo registry unavailable")
    input_sha256 = {str(path): PRIOR.sha256_file(path) for path in
                    (authorization, PLAN, *(archives / f"{digest}.jsonl.gz" for digest in
                                            (preflight.BASE, *preflight.POOLS.values())))}
    started = time.monotonic()
    deadline = started + TOTAL_CAP
    output.mkdir(parents=True, exist_ok=False)
    PRIOR.write_new_json(output / "config.json", {
        "schema": "n83.k600-gray-hash-paired-run-config/v1", "source_commit": commit,
        "source_sha256": sources, "input_sha256": input_sha256,
        "authorization_sha256": hashlib.sha256(record_bytes).hexdigest(),
        "plan_sha256": plan_sha256, "container_image_id": image,
        "conservative_full_charge_seconds": TOTAL_CAP,
        "selected_best_total_runtime": None,
    })
    config_sha256 = PRIOR.sha256_file(output / "config.json")
    arms = []
    for arm in plan["arms"]:
        if time.monotonic() + ARM_CAP > deadline:
            break
        arm_dir = output / arm["role"]
        build, run_dir = arm_dir / "build", arm_dir / "run"
        (build / "target").mkdir(parents=True)
        run_dir.mkdir()
        binary = build / "worker"
        token = uuid.uuid4().hex[:16]
        stages = []
        status = "PRODUCER_FAILURE_preflight"
        binary_sha256 = None
        replay_error = None
        for name in ("build", "import", "search"):
            if time.monotonic() + WALL[name] + STAGE_CAP["cleanup"] > deadline:
                status = "UNKNOWN_budget_before_stage"
                break
            try:
                require(PRIOR.clean_commit(checkout) == commit and source_hashes(checkout) == sources,
                        "source changed during paired gate")
                require(all(PRIOR.sha256_file(Path(path)) == digest
                            for path, digest in input_sha256.items()),
                        "input changed during paired gate")
                if name != "build":
                    require(PRIOR.sha256_file(binary) == binary_sha256,
                            "binary changed during paired gate")
            except (OSError, ValueError) as exc:
                status = "PRODUCER_FAILURE_frozen_input_changed"
                replay_error = str(exc)
                break
            command, stage_output, container = stage_command(
                name, arm, checkout, archives, build, run_dir, binary, registry, image, token)
            receipt = PRIOR.stage(name=name, command=command, wall=WALL[name],
                                  output=stage_output, container_name=container)
            stages.append(receipt)
            issue = stage_issue(receipt, name, run_dir)
            if issue:
                status = issue
                break
            try:
                if name == "build":
                    artifact = build / "target/release/examples/koblitz_n83_mixed_orbit_worker"
                    require(PRIOR.is_linux_elf(artifact), "worker binary is not Linux ELF")
                    shutil.copyfile(artifact, binary)
                    binary.chmod(0o755)
                    binary_sha256 = PRIOR.sha256_file(binary)
                    status = "PASS_build_only"
                elif name == "import":
                    verify = json.loads((run_dir / "import.stdout").read_text())
                    audit.audit_verify(verify, arm)
                    audit.domains(archives / f"{arm['base_object_blake3']}.jsonl.gz",
                                  archives / f"{arm['pool_object_blake3']}.jsonl.gz",
                                  arm["shared_signed_orbits"])
                    status = "PASS_import_only"
                else:
                    base, residual = audit.domains(
                        archives / f"{arm['base_object_blake3']}.jsonl.gz",
                        archives / f"{arm['pool_object_blake3']}.jsonl.gz",
                        arm["shared_signed_orbits"])
                    status = audit.audit_search(
                        json.loads((run_dir / "search/config.json").read_text()),
                        json.loads((run_dir / "search/result.json").read_text()),
                        arm, base, residual)
            except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
                status = "PRODUCER_FAILURE_receipt_replay"
                replay_error = f"{type(exc).__name__}: {exc}"
                break
        outer = {
            "schema": "n83.k600-gray-hash-paired-arm-outer/v1", "role": arm["role"],
            "status": status, "stages": stages, "binary_sha256": binary_sha256,
            "replay_error": replay_error,
            "conservative_charge_seconds": ARM_CAP,
            "graph_rows": None, "matrix_rank": None,
            "total_index_calculus_runtime_ms": None,
        }
        PRIOR.write_new_json(arm_dir / "outer.json", outer)
        arms.append(outer)
        if status.startswith("PRODUCER_FAILURE"):
            break
    frozen_changed = False
    try:
        frozen_changed = (PRIOR.clean_commit(checkout) != commit
                          or source_hashes(checkout) != sources
                          or any(PRIOR.sha256_file(Path(path)) != digest
                                 for path, digest in input_sha256.items())
                          or PRIOR.sha256_file(output / "config.json") != config_sha256)
    except (OSError, ValueError):
        frozen_changed = True
    elapsed = time.monotonic() - started
    status = ("PRODUCER_FAILURE_frozen_input_changed" if frozen_changed else
              "PRODUCER_FAILURE_budget_overshoot" if elapsed > TOTAL_CAP else
              "PRODUCER_FAILURE_arm" if any(row["status"].startswith("PRODUCER_FAILURE") for row in arms) else
              "BOUNDED_UNKNOWN_incomplete_pair" if len(arms) != 2 else
              "BOUNDED_PARTIAL_paired_first_trial")
    outer = {
        "schema": "n83.k600-gray-hash-paired-outer/v1", "status": status,
        "source_commit": commit, "source_sha256": sources,
        "config_sha256": config_sha256, "input_sha256": input_sha256,
        "arms": arms, "frozen_input_changed": frozen_changed,
        "elapsed_supervisor_seconds_L0": elapsed,
        "conservative_charge_seconds": ARM_CAP * len(arms),
        "cumulative_conservative_charge_seconds": preflight.PRIOR_CHARGE + ARM_CAP * len(arms),
        "new_n83_solver_search_executed": any(
            any(stage["name"] == "search" and stage["launch_error"] is None
                and stage["container_state"] is not None
                for stage in arm["stages"]) for arm in arms),
        "graph_rows": None, "matrix_rank": None,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    PRIOR.write_new_json(output / "outer.json", outer)
    return outer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--archive-dir", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.authorization, args.output_dir, args.archive_dir)
    print(json.dumps(result, sort_keys=True))
    return 1 if result["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
