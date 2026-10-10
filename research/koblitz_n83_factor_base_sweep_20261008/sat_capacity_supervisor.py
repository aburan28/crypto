#!/usr/bin/env python3
"""Bound one retained-base factored S4 model construction; do not solve it.

Docker supplies a hard cgroup memory limit with zero swap and no network.
The native Linux worker verifies those limits before importing the base.
Every launched outcome, including a cap or worker failure, leaves an outer receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import uuid


SCHEMA = "n83.factored-s4-capacity-outer/v1"
WORKER_SCHEMA = "n83.factored-s4-capacity-worker/v1"


def write_new_json(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as output:
        json.dump(value, output, indent=2, sort_keys=True)
        output.write("\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(
    panel: Path,
    columns: int,
    wall_seconds: float,
    memory_mib: int,
    binary: Path,
    output_dir: Path,
    image: str = "python:3.11-slim",
) -> dict:
    if columns not in (64, 256, 600):
        raise ValueError("capacity gate requires retained primary K=64/256/600")
    if not 0 < wall_seconds <= 86400 or not 0 < memory_mib <= 65536:
        raise ValueError("capacity wall or memory bound is outside the supported range")
    panel = panel.resolve(strict=True)
    binary = binary.resolve(strict=True)
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError("capacity worker binary must be executable")
    with binary.open("rb") as executable:
        magic = executable.read(4)
    if magic != b"\x7fELF":
        raise ValueError("capacity worker must be a Linux ELF binary")
    inspected = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Id}}"],
        capture_output=True,
        text=True,
        timeout=15,
        check=True,
    )
    image_id = inspected.stdout.strip()
    if not image_id.startswith("sha256:"):
        raise ValueError("capacity container image has no immutable image ID")
    limit_bytes = memory_mib * 1024 * 1024
    output_dir.mkdir(parents=True, exist_ok=False)
    output_dir = output_dir.resolve(strict=True)
    container_name = f"n83-sat-capacity-{uuid.uuid4().hex[:16]}"
    config = {
        "schema": "n83.factored-s4-capacity-config/v1",
        "panel_dir": str(panel),
        "columns": columns,
        "wall_seconds": wall_seconds,
        "memory_cgroup_limit_bytes": limit_bytes,
        "memory_cgroup_swap_limit_bytes": 0,
        "binary": str(binary),
        "binary_sha256": sha256_file(binary),
        "container_image": image,
        "container_image_id": image_id,
        "container_name": container_name,
        "solver_search_executed": False,
        "selected_best_total_runtime": None,
    }
    write_new_json(output_dir / "config.json", config)

    command = [
        "docker",
        "run",
        "--rm",
        "--pull",
        "never",
        "--name",
        container_name,
        "--network",
        "none",
        "--memory",
        f"{memory_mib}m",
        "--memory-swap",
        f"{memory_mib}m",
        "--pids-limit",
        "64",
        "--cap-drop",
        "ALL",
        "--security-opt",
        "no-new-privileges",
        "--read-only",
        "--tmpfs",
        "/tmp:rw,noexec,nosuid,size=64m",
        "--user",
        f"{os.getuid()}:{os.getgid()}",
        "--mount",
        f"type=bind,src={panel},dst=/panel,readonly",
        "--mount",
        f"type=bind,src={binary},dst=/worker,readonly",
        "--mount",
        f"type=bind,src={output_dir},dst=/out",
        image_id,
        "/worker",
        "primary-sat-build",
        "/panel",
        str(columns),
        str(memory_mib),
        "/out/worker.json",
    ]

    started = time.monotonic()
    timed_out = False
    exit_code: int | None = None
    launch_error: str | None = None
    cleanup_error: str | None = None
    with (output_dir / "stdout.log").open("xb") as stdout, (
        output_dir / "stderr.log"
    ).open("xb") as stderr:
        try:
            child = subprocess.Popen(
                command,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
            try:
                exit_code = child.wait(timeout=wall_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True
                try:
                    stop = subprocess.run(
                        ["docker", "kill", container_name],
                        capture_output=True,
                        text=True,
                        timeout=15,
                        check=False,
                    )
                    if stop.returncode != 0:
                        cleanup_error = stop.stderr.strip() or "docker kill failed"
                except (OSError, subprocess.SubprocessError) as error:
                    cleanup_error = f"docker kill failed: {type(error).__name__}: {error}"
                try:
                    exit_code = child.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    child.kill()
                    exit_code = child.wait()
        except (OSError, subprocess.SubprocessError) as error:
            launch_error = f"{type(error).__name__}: {error}"
    if timed_out:
        try:
            remaining = subprocess.run(
                ["docker", "ps", "-aq", "--filter", f"name=^{container_name}$"],
                capture_output=True,
                text=True,
                timeout=15,
                check=False,
            )
            if remaining.returncode != 0 or remaining.stdout.strip():
                cleanup_error = cleanup_error or "timed-out container still present"
        except (OSError, subprocess.SubprocessError) as error:
            cleanup_error = cleanup_error or f"container cleanup unverified: {error}"
    worker_path = output_dir / "worker.json"
    worker: dict | None = None
    worker_error: str | None = None
    if worker_path.exists():
        try:
            worker = json.loads(worker_path.read_text(encoding="utf-8"))
            if not isinstance(worker, dict):
                worker = None
                worker_error = "worker receipt is not a JSON object"
        except (OSError, json.JSONDecodeError) as error:
            worker_error = f"{type(error).__name__}: {error}"
    if cleanup_error is not None:
        status = "PRODUCER_FAILURE_container_cleanup"
    elif timed_out:
        status = "UNKNOWN_wall_cap"
    elif launch_error is not None:
        status = "PRODUCER_FAILURE_launch"
    elif exit_code in (134, 137):
        status = "UNKNOWN_resource_or_worker_exit"
    elif exit_code != 0:
        status = "PRODUCER_FAILURE_worker_exit"
    elif (
        worker is None
        or worker.get("schema") != WORKER_SCHEMA
        or worker.get("status") != "PASS_model_construction_only"
        or worker.get("orbit_columns") != columns
        or worker.get("memory_cgroup_limit_bytes") != limit_bytes
        or worker.get("memory_cgroup_swap_limit_bytes") != 0
        or worker.get("solver_search_executed") is not False
        or worker.get("model_lifting_executed") is not False
        or worker.get("total_index_calculus_runtime_ms") is not None
        or worker.get("selected_best_total_runtime") is not None
    ):
        status = "PRODUCER_FAILURE_worker_receipt"
    else:
        status = "PASS_model_construction_only"
    outer = {
        "schema": SCHEMA,
        "status": status,
        "config_sha256": sha256_file(output_dir / "config.json"),
        "worker_sha256": sha256_file(worker_path) if worker_path.exists() else None,
        "worker_exit_code": exit_code,
        "launch_error": launch_error,
        "cleanup_error": cleanup_error,
        "worker_error": worker_error,
        "process_wall_ms": (time.monotonic() - started) * 1000,
        "memory_cgroup_limit_bytes": limit_bytes,
        "memory_cgroup_swap_limit_bytes": 0,
        "worker_process_usage": worker.get("process_usage") if worker else None,
        "worker_cgroup_peak_bytes": worker.get("memory_cgroup_peak_bytes") if worker else None,
        "solver_search_executed": False,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    write_new_json(output_dir / "outer.json", outer)
    return outer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--columns", type=int, required=True)
    parser.add_argument("--wall-seconds", type=float, required=True)
    parser.add_argument("--memory-mib", type=int, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--image", default="python:3.11-slim")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    outer = run(
        args.panel,
        args.columns,
        args.wall_seconds,
        args.memory_mib,
        args.binary,
        args.output_dir,
        args.image,
    )
    print(json.dumps(outer, sort_keys=True))
    return 1 if outer["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
