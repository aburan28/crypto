#!/usr/bin/env python3
"""Run the frozen two-arm planted six-summand SAT access ladder under Docker limits."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FIXTURE = HERE / "fixture.json"
FIXTURE_SHA256 = "76b78ec1a68807d130fc11e54e4220b5cde619d677316e313dc247869f727a9d"
IMAGE = "rust:1.98.1-bookworm@sha256:93ce27a88655056a51dbdd8f5f2d7ddc071c7b0070fb288a37b5a285fc83971e"
STAGES = (("all_x", 600, 5000, 10), ("summand_x", 900, 10000, 6),
          ("open", 1800, 20000, 0))
ROLES = ("gray_pool", "hash_pool")


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_new(path: Path, row: dict) -> None:
    with path.open("x") as stream:
        json.dump(row, stream, sort_keys=True, indent=2)
        stream.write("\n")


def current_commit() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                   text=True).strip()


def clean_source() -> bool:
    return not subprocess.check_output(["git", "status", "--porcelain"],
                                       cwd=ROOT, text=True)


def stage_run(role: str, stage: str, wall_cap: int, conflict_cap: int,
              worker: Path, image: str, output: Path, run_id: str) -> dict:
    name = f"n83-m6-witness-{run_id}-{role}-{stage}"
    command = [
        "docker", "run", "--pull", "never", "--platform", "linux/arm64",
        "--name", name, "--network", "none", "--cpus", "1",
        "--memory", "4096m", "--memory-swap", "4096m", "--pids-limit", "64",
        "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
        "--read-only", "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--mount", f"type=bind,src={worker},dst=/worker,readonly",
        image, "/worker", role, stage, str(conflict_cap),
    ]
    started = time.monotonic()
    try:
        done = subprocess.run(command, capture_output=True, timeout=wall_cap,
                              check=False)
        exit_code, stdout, stderr = done.returncode, done.stdout, done.stderr
        timed_out = False
    except subprocess.TimeoutExpired as exc:
        exit_code, stdout, stderr = None, exc.stdout or b"", exc.stderr or b""
        timed_out = True
    elapsed = time.monotonic() - started
    (output / "solve.stdout").write_bytes(stdout)
    (output / "solve.stderr").write_bytes(stderr)
    inspected = subprocess.run(["docker", "inspect", name], capture_output=True,
                               timeout=30, check=False)
    (output / "container.inspect.json").write_bytes(inspected.stdout)
    (output / "container.inspect.stderr").write_bytes(inspected.stderr)
    cleaned = subprocess.run(["docker", "rm", "-f", name], capture_output=True,
                             timeout=30, check=False)
    (output / "container.cleanup.stdout").write_bytes(cleaned.stdout)
    (output / "container.cleanup.stderr").write_bytes(cleaned.stderr)
    if inspected.returncode != 0 or cleaned.returncode != 0:
        status = "PRODUCER_FAILURE_inspect_or_cleanup"
    elif timed_out:
        status = "UNKNOWN_outer_wall_cap"
    elif exit_code != 0:
        status = "PRODUCER_FAILURE_worker_exit"
    else:
        try:
            result = json.loads(stdout)
            require(result["schema"] == "n83.m6-planted-witness-control/v1"
                    and result["role"] == role and result["stage"] == stage
                    and result["conflict_cap"] == conflict_cap
                    and result["fixed_x_blocks"] == {"all_x": 10,
                                                      "summand_x": 6,
                                                      "open": 0}[stage]
                    and 0 <= result["conflicts"] <= conflict_cap
                    and result["fixture_blake3"] == "f07b515f10b4cfd0b0d82267a22c1d46018d278fc922de7cd7b73c77e61ed43a"
                    and result["selected_best_total_runtime"] is None,
                    "worker result identity differs")
            status = result["status"]
            allowed_sat = "SAT_verified_field_model" if stage == "open" else "SAT_verified_planted_x"
            require(status in (allowed_sat, "UNKNOWN_conflict_cap")
                    or status.startswith("PRODUCER_FAILURE"),
                    "worker status differs")
            require(not status.startswith("SAT_") or result["model_verified"] is True,
                    "SAT result lacks model verification")
        except (ValueError, KeyError, TypeError, json.JSONDecodeError):
            status = "PRODUCER_FAILURE_receipt"
    receipt = {
        "schema": "n83.m6-witness-ladder-stage/v1",
        "role": role, "stage": stage, "command": command,
        "wall_cap_seconds": wall_cap, "conflict_cap": conflict_cap,
        "elapsed_wall_seconds_L0": elapsed, "exit_code": exit_code,
        "timed_out": timed_out, "status": status,
        "stdout_sha256": sha256(output / "solve.stdout"),
        "stderr_sha256": sha256(output / "solve.stderr"),
        "inspect_exit_code": inspected.returncode,
        "cleanup_exit_code": cleaned.returncode,
    }
    write_new(output / "stage.json", receipt)
    return receipt


def execute(worker: Path, output: Path, image: str, run_id: str) -> dict:
    require(image == IMAGE and worker.is_file() and sha256(FIXTURE) == FIXTURE_SHA256
            and run_id.isdecimal() and clean_source(),
            "frozen worker, fixture, image, run ID or source differs")
    fixture = json.loads(FIXTURE.read_text())
    require(fixture["schema"] == "n83.m6-planted-witness-ladder/v1"
            and [arm["role"] for arm in fixture["arms"]] == list(ROLES),
            "frozen fixture arms differ")
    commit = current_commit()
    output.mkdir(parents=True, exist_ok=False)
    write_new(output / "config.json", {
        "schema": "n83.m6-witness-ladder-config/v1",
        "source_commit": commit, "worker_sha256": sha256(worker),
        "fixture_sha256": FIXTURE_SHA256, "image": image, "run_id": run_id,
        "conservative_full_charge_seconds": 10800,
        "selected_best_total_runtime": None,
    })
    started = time.monotonic()
    receipts = []
    for role in ROLES:
        for stage, wall_cap, conflict_cap, _ in STAGES:
            stage_dir = output / role / stage
            stage_dir.mkdir(parents=True)
            receipt = stage_run(role, stage, wall_cap, conflict_cap, worker,
                                image, stage_dir, run_id)
            receipts.append({"role": role, "stage": stage, "status": receipt["status"],
                             "stage_receipt": str(stage_dir.relative_to(output) / "stage.json")})
            if receipt["status"].startswith("PRODUCER_FAILURE"):
                break
        if receipts[-1]["status"].startswith("PRODUCER_FAILURE"):
            break
    changed = current_commit() != commit or not clean_source()
    status = ("PRODUCER_FAILURE_source_changed" if changed else
              "PRODUCER_FAILURE_stage" if any(r["status"].startswith("PRODUCER_FAILURE")
                                              for r in receipts) else
              "BOUNDED_UNKNOWN_incomplete_ladder" if len(receipts) != 6 else
              "BOUNDED_PARTIAL_six_stage_access_ladder")
    outer = {
        "schema": "n83.m6-witness-ladder-outer/v1",
        "status": status, "source_commit": commit, "run_id": run_id,
        "stages": receipts, "elapsed_supervisor_wall_seconds_L0": time.monotonic() - started,
        "conservative_full_charge_seconds": 10800,
        "natural_relations": None, "matrix_rank": None,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    write_new(output / "outer.json", outer)
    return outer


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--image", default=IMAGE)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    result = execute(args.worker.resolve(strict=True), args.output,
                     args.image, args.run_id)
    print(json.dumps({"status": result["status"],
                      "stages": [row["status"] for row in result["stages"]]}, sort_keys=True))
    if result["status"].startswith("PRODUCER_FAILURE"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
