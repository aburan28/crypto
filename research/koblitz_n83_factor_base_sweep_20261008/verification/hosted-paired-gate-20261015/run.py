#!/usr/bin/env python3
"""Run the frozen K600 Gray/hash pair on one manually admitted Linux ARM runner.

This is a separate hosted protocol: it retains the mathematical inputs and
worker limits of the local plan, but records its own immutable container image
and host. It never selects a total-runtime winner from a first-trial result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import uuid


HERE = Path(__file__).resolve().parent
PAIR = HERE.parent / "k600-gray-hash-paired-protocol-20261014"
sys.path.insert(0, str(PAIR))
import preflight  # noqa: E402
import receipt_audit as audit  # noqa: E402
import supervisor  # noqa: E402


IMAGE_REF = "rust:1.98.1-bookworm@sha256:93ce27a88655056a51dbdd8f5f2d7ddc071c7b0070fb288a37b5a285fc83971e"
PLAN_SHA256 = "fbd3fb943bc0f181797d106c41ee9961b8d70768b5c65a5584002e3c1b5e488f"
TOTAL_CAP = 4800
STAGES = ("build", "import", "search")


def sha256(path: Path) -> str:
    return preflight.sha256(path)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write("\n")


def image_identity() -> tuple[str, dict]:
    result = subprocess.run(
        ["docker", "image", "inspect", IMAGE_REF,
         "--format", "{{json .}}"], capture_output=True, text=True,
        timeout=15, check=True,
    )
    row = json.loads(result.stdout)
    if row.get("Os") != "linux" or row.get("Architecture") != "arm64":
        raise ValueError("pinned image is not Linux arm64")
    image_id = row.get("Id")
    if not isinstance(image_id, str) or not image_id.startswith("sha256:"):
        raise ValueError("Docker image ID missing")
    return image_id, {
        "reference": IMAGE_REF, "image_id": image_id,
        "repo_digests": row.get("RepoDigests", []),
        "os": row["Os"], "architecture": row["Architecture"],
    }


def cgroup_probe(image_id: str) -> dict:
    command = [
        "docker", "run", "--rm", "--pull", "never", "--network", "none",
        "--cpus", "1", "--memory", "4096m", "--memory-swap", "4096m",
        image_id, "sh", "-c",
        "cat /sys/fs/cgroup/memory.max /sys/fs/cgroup/memory.swap.max",
    ]
    result = subprocess.run(command, capture_output=True, text=True,
                            timeout=30, check=True)
    if result.stdout.splitlines() != ["4294967296", "0"]:
        raise ValueError("hosted Linux cgroup does not provide 4 GiB and zero swap")
    return {"command": command, "stdout": result.stdout,
            "stderr_sha256": hashlib.sha256(result.stderr.encode()).hexdigest()}


def checked_plan(archives: Path) -> dict:
    if sha256(preflight.PLAN) != PLAN_SHA256:
        raise ValueError("frozen paired plan bytes changed")
    plan = preflight.make_plan(archives)
    expected = (json.dumps(plan, indent=2, sort_keys=True) + "\n").encode()
    if preflight.PLAN.read_bytes() != expected:
        raise ValueError("frozen paired plan replay changed")
    return plan


def checked_preparation(path: Path, source_commit: str,
                        input_sha256: dict[str, str]) -> dict:
    row = json.loads(path.read_text())
    by_name = {Path(name).name: digest for name, digest in input_sha256.items()}
    if (row.get("schema") != "n83.hosted-regenerated-archives/v1"
            or row.get("status") != "PASS_3_FROZEN_OBJECTS"
            or row.get("source_commit") != source_commit
            or row.get("workflow_dispatch_head_sha") != os.environ.get("GITHUB_SHA")
            or row.get("archive_sha256") != by_name
            or row.get("pilot_completed_base_count") != 54):
        raise ValueError("hosted preparation receipt does not bind frozen inputs")
    seconds = row.get("hosted_preparation_wall_seconds_L0")
    if not isinstance(seconds, (int, float)) or not 0 <= seconds <= 1800:
        raise ValueError("hosted preparation exceeded reserved setup envelope")
    return row


def stage_receipt(name: str, arm: dict, checkout: Path, archives: Path,
                  build: Path, run_dir: Path, binary: Path, image_id: str,
                  registry: Path, token: str) -> dict:
    command, stage_output, container = supervisor.stage_command(
        name, arm, checkout, archives, build, run_dir, binary, registry,
        image_id, token,
    )
    return supervisor.PRIOR.stage(
        name=name, command=command, wall=supervisor.WALL[name],
        output=stage_output, container_name=container,
    )


def run(archives: Path, output: Path, scratch: Path,
        preparation_receipt: Path | None = None) -> dict:
    checkout = preflight.ROOT.resolve(strict=True)
    archives = archives.resolve(strict=True)
    output = output.resolve()
    scratch = scratch.resolve()
    if output.exists() or scratch.exists() or output.is_relative_to(checkout) or scratch.is_relative_to(checkout):
        raise ValueError("output and scratch must be fresh paths outside checkout")
    if platform.system() != "Linux" or platform.machine() not in ("aarch64", "arm64"):
        raise ValueError("hosted gate requires Linux arm64")
    plan = checked_plan(archives)
    commit = supervisor.PRIOR.clean_commit(checkout)
    sources = supervisor.source_hashes(checkout)
    hosted_sources = {
        str(path.relative_to(checkout)): sha256(path)
        for path in (HERE / "run.py", HERE / "test_run.py",
                     HERE / "prepare_archives.py",
                     checkout / ".github/workflows/n83-hosted-paired-gate.yml")
    }
    image_id, image = image_identity()
    cgroup = cgroup_probe(image_id)
    registry = Path(os.environ.get("CARGO_HOME", str(Path.home() / ".cargo"))) / "registry"
    if not registry.is_dir():
        raise ValueError("offline Cargo registry missing")
    inputs = {str(archives / f"{digest}.jsonl.gz"): sha256(archives / f"{digest}.jsonl.gz")
              for digest in (preflight.BASE, *preflight.POOLS.values())}
    if preparation_receipt is None:
        raise ValueError("hosted preparation receipt required")
    preparation = checked_preparation(preparation_receipt, commit, inputs)
    started = time.monotonic()
    deadline = started + TOTAL_CAP
    output.mkdir(parents=True)
    scratch.mkdir(parents=True)
    write_json(output / "config.json", {
        "schema": "n83.hosted-gray-hash-paired-config/v1",
        "status": "HOSTED_PROTOCOL_SEPARATE_FROM_LOCAL_IMAGE",
        "source_commit": commit, "source_sha256": sources,
        "hosted_source_sha256": hosted_sources,
        "hosted_runner": os.environ.get("RUNNER_NAME"),
        "github_run_id": os.environ.get("GITHUB_RUN_ID"),
        "plan_sha256": PLAN_SHA256, "input_sha256": inputs,
        "image": image, "cgroup_probe": cgroup,
        "conservative_full_charge_seconds": TOTAL_CAP,
        "prior_local_charge_seconds": preflight.PRIOR_CHARGE,
        "preparation": preparation,
        "selected_best_total_runtime": None,
    })
    arms = []
    for arm in plan["arms"]:
        arm_out = output / arm["role"]
        run_dir = arm_out / "run"
        build = scratch / arm["role"] / "build"
        (build / "target").mkdir(parents=True)
        run_dir.mkdir(parents=True)
        binary = arm_out / "worker"
        token = uuid.uuid4().hex[:16]
        stages = []
        status = "UNKNOWN_before_stage"
        replay_error = None
        for name in STAGES:
            if time.monotonic() + supervisor.WALL[name] + preflight.CHARGE_PER_ARM["cleanup"] > deadline:
                status = "UNKNOWN_budget_before_stage"
                break
            if (supervisor.PRIOR.clean_commit(checkout) != commit
                    or supervisor.source_hashes(checkout) != sources
                    or any(sha256(checkout / path) != digest
                           for path, digest in hosted_sources.items())
                    or any(sha256(Path(path)) != digest for path, digest in inputs.items())):
                status = "PRODUCER_FAILURE_frozen_input_changed"
                break
            receipt = stage_receipt(name, arm, checkout, archives, build,
                                    run_dir, binary, image_id, registry, token)
            stages.append(receipt)
            if name == "build":
                for suffix in ("stdout", "stderr"):
                    shutil.copyfile(build / f"build.{suffix}", arm_out / f"build.{suffix}")
            write_json(arm_out / f"{name}.stage.json", receipt)
            issue = supervisor.stage_issue(receipt, name, run_dir)
            if issue:
                status = issue
                break
            try:
                if name == "build":
                    artifact = build / "target/release/examples/koblitz_n83_mixed_orbit_worker"
                    if not supervisor.PRIOR.is_linux_elf(artifact):
                        raise ValueError("built worker is not Linux ELF")
                    shutil.copyfile(artifact, binary)
                    binary.chmod(0o755)
                    status = "PASS_build_only"
                elif name == "import":
                    audit.audit_verify(json.loads((run_dir / "import.stdout").read_text()), arm)
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
            except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as error:
                status = "PRODUCER_FAILURE_receipt_replay"
                replay_error = f"{type(error).__name__}: {error}"
                break
        outer = {
            "schema": "n83.hosted-gray-hash-paired-arm/v1", "role": arm["role"],
            "status": status, "stages": stages,
            "binary_sha256": sha256(binary) if binary.is_file() else None,
            "replay_error": replay_error,
            "conservative_charge_seconds": sum(preflight.CHARGE_PER_ARM.values()),
            "graph_rows": None, "matrix_rank": None,
            "total_index_calculus_runtime_ms": None,
        }
        write_json(arm_out / "outer.json", outer)
        arms.append(outer)
        if status.startswith("PRODUCER_FAILURE"):
            break
    elapsed = time.monotonic() - started
    frozen_changed = False
    try:
        frozen_changed = (supervisor.PRIOR.clean_commit(checkout) != commit
                          or supervisor.source_hashes(checkout) != sources
                          or any(sha256(checkout / path) != digest
                                 for path, digest in hosted_sources.items())
                          or any(sha256(Path(path)) != digest for path, digest in inputs.items()))
    except (OSError, ValueError):
        frozen_changed = True
    status = ("PRODUCER_FAILURE_frozen_input_changed" if frozen_changed else
              "PRODUCER_FAILURE_budget_overshoot" if elapsed > TOTAL_CAP else
              "PRODUCER_FAILURE_arm" if any(a["status"].startswith("PRODUCER_FAILURE") for a in arms)
              else "BOUNDED_UNKNOWN_incomplete_pair" if (len(arms) != 2 or any(
                  a["status"] in ("UNKNOWN_budget_before_stage", "UNKNOWN_before_stage")
                  for a in arms))
              else "BOUNDED_PARTIAL_hosted_first_trial")
    result = {
        "schema": "n83.hosted-gray-hash-paired-outer/v1", "status": status,
        "arms": arms, "elapsed_supervisor_seconds_L0": elapsed,
        "conservative_charge_seconds": len(arms) * sum(preflight.CHARGE_PER_ARM.values()),
        "preparation_wall_seconds_L0": preparation["hosted_preparation_wall_seconds_L0"],
        "conservative_charge_including_preparation_seconds_L0": (
            len(arms) * sum(preflight.CHARGE_PER_ARM.values())
            + math.ceil(preparation["hosted_preparation_wall_seconds_L0"])),
        "frozen_input_changed": frozen_changed,
        "new_n83_solver_search_executed": any(
            any(stage["name"] == "search" and stage["launch_error"] is None
                and stage["container_state"] is not None for stage in arm["stages"])
            for arm in arms),
        "selected_best_total_runtime": None,
    }
    write_json(output / "outer.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archives", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--preparation-receipt", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.archives, args.output, args.scratch, args.preparation_receipt)
    print(json.dumps(result, sort_keys=True))
    return 1 if result["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
