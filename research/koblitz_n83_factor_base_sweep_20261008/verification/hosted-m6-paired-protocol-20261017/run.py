#!/usr/bin/env python3
"""Guarded networkless 4+2 Gray/hash search using three exact pilot archives."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import plan

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
MANIFEST = ROOT / "research/koblitz_n83_factor_base_sweep_20261008/pilot-01/manifest.json"
CORPUS = ROOT / "research/koblitz_n83_factor_base_sweep_20261008/pilot-01/probe-corpus.json"
IMAGE = "rust:1.98.1-bookworm@sha256:93ce27a88655056a51dbdd8f5f2d7ddc071c7b0070fb288a37b5a285fc83971e"


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_new(path: Path, value: dict) -> None:
    with path.open("x") as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write("\n")


def container_stage(stage: str, arm: dict, archives: Path, worker: Path,
                    output: Path, image: str, run_id: str, wall: int) -> dict:
    role = arm["role"]
    name = f"n83-m6-{run_id}-{role}-{stage}"
    base = archives / f"{arm['base_object_blake3']}.jsonl.gz"
    pool = archives / f"{arm['pool_object_blake3']}.jsonl.gz"
    require(digest(base) == arm["base_compressed_sha256"]
            and digest(pool) == arm["pool_compressed_sha256"],
            "frozen archive bytes changed before stage")
    argv = [
        "docker", "run", "--pull", "never", "--platform", "linux/arm64",
        "--name", name, "--network", "none", "--cpus", "1",
        "--memory", "4096m", "--memory-swap", "4096m", "--pids-limit", "64",
        "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
        "--read-only", "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--mount", f"type=bind,src={worker},dst=/worker,readonly",
        "--mount", f"type=bind,src={base},dst=/base.jsonl.gz,readonly",
        "--mount", f"type=bind,src={pool},dst=/pool.jsonl.gz,readonly",
        "--mount", f"type=bind,src={MANIFEST},dst=/pilot-manifest.json,readonly",
        "--mount", f"type=bind,src={CORPUS},dst=/corpus.json,readonly",
        "--mount", f"type=bind,src={output},dst=/out",
        image,
        *(arm["worker_verify_argv"] if stage == "import" else arm["worker_search_argv"]),
    ]
    started = time.monotonic()
    timed_out = False
    try:
        result = subprocess.run(argv, capture_output=True, timeout=wall, check=False)
        code = result.returncode
        stdout, stderr = result.stdout, result.stderr
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        code = None
        stdout, stderr = exc.stdout or b"", exc.stderr or b""
    elapsed = time.monotonic() - started
    stdout_path = output / f"{stage}.stdout"
    stderr_path = output / f"{stage}.stderr"
    stdout_path.write_bytes(stdout)
    stderr_path.write_bytes(stderr)
    inspected = subprocess.run(["docker", "inspect", name], capture_output=True,
                               timeout=30, check=False)
    (output / f"{stage}.inspect.json").write_bytes(inspected.stdout)
    (output / f"{stage}.inspect.stderr").write_bytes(inspected.stderr)
    cleaned = subprocess.run(["docker", "rm", "-f", name], capture_output=True,
                             timeout=30, check=False)
    (output / f"{stage}.cleanup.stdout").write_bytes(cleaned.stdout)
    (output / f"{stage}.cleanup.stderr").write_bytes(cleaned.stderr)
    receipt = {
        "schema": "n83.hosted-m6-paired-stage/v1", "stage": stage,
        "role": role, "command": argv, "wall_cap_seconds": wall,
        "elapsed_wall_seconds_L0": elapsed, "exit_code": code,
        "timed_out": timed_out,
        "stdout_sha256": digest(stdout_path), "stderr_sha256": digest(stderr_path),
        "inspect_exit_code": inspected.returncode,
        "cleanup_exit_code": cleaned.returncode,
        "cleanup_stdout": cleaned.stdout.decode(errors="replace").strip(),
    }
    write_new(output / f"{stage}.stage.json", receipt)
    return receipt


def execute(archives: Path, worker: Path, output: Path, image: str,
            run_id: str) -> dict:
    expected = plan.build()
    require(plan.SAVED.read_text() == json.dumps(expected, indent=2, sort_keys=True) + "\n",
            "frozen 4+2 plan differs")
    require(image == IMAGE and worker.is_file() and MANIFEST.is_file()
            and CORPUS.is_file() and run_id.isdecimal(),
            "worker, public inputs or image differ")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                     text=True).strip()
    require(not subprocess.check_output(["git", "status", "--porcelain"],
                                        cwd=ROOT, text=True), "source checkout is dirty")
    for arm in expected["arms"]:
        for kind in ("base", "pool"):
            archive = archives / f"{arm[f'{kind}_object_blake3']}.jsonl.gz"
            require(digest(archive) == arm[f"{kind}_compressed_sha256"],
                    "frozen pilot archive differs")
    output.mkdir(parents=True, exist_ok=False)
    write_new(output / "config.json", {
        "schema": "n83.hosted-m6-paired-config/v1",
        "source_commit": commit, "plan_sha256": digest(plan.SAVED),
        "worker_sha256": digest(worker), "manifest_sha256": digest(MANIFEST),
        "corpus_sha256": digest(CORPUS), "image": image,
        "run_id": run_id, "conservative_full_charge_seconds": 7500,
        "selected_best_total_runtime": None,
    })
    started = time.monotonic()
    arms = []
    for arm in expected["arms"]:
        arm_dir = output / arm["role"]
        arm_dir.mkdir()
        imported = container_stage("import", arm, archives, worker, arm_dir,
                                   image, run_id, 300)
        if (imported["timed_out"] or imported["exit_code"] != 0
                or imported["cleanup_exit_code"] != 0):
            arms.append({"role": arm["role"], "status": "PRODUCER_FAILURE_import",
                         "import": imported, "search": None})
            break
        try:
            verify = json.loads((arm_dir / "import.stdout").read_text())
            require(verify["schema"] == "n83.cross-stream-pilot-import/v1"
                    and verify["status"] == "PASS"
                    and verify["small_object_blake3"] == arm["base_object_blake3"]
                    and verify["pool_object_blake3"] == arm["pool_object_blake3"]
                    and verify["shared_signed_orbits"] == arm["shared_signed_orbits"]
                    and verify["residual_signed_orbits"] == arm["residual_signed_orbits"]
                    and verify["solver_search_executed"] is False,
                    "native exact import did not match frozen domains")
        except (ValueError, KeyError, json.JSONDecodeError):
            arms.append({"role": arm["role"], "status": "PRODUCER_FAILURE_import_receipt",
                         "import": imported, "search": None})
            break
        searched = container_stage("search", arm, archives, worker, arm_dir,
                                   image, run_id, 2400)
        if searched["timed_out"]:
            status = "UNKNOWN_outer_wall_cap"
        elif searched["exit_code"] == 124:
            status = "UNKNOWN_worker_wall_cap"
        elif searched["exit_code"] != 0 or searched["cleanup_exit_code"] != 0:
            status = "PRODUCER_FAILURE_search"
        else:
            try:
                status = json.loads((arm_dir / "search/result.json").read_text())["status"]
            except (OSError, KeyError, json.JSONDecodeError):
                status = "PRODUCER_FAILURE_result_receipt"
        arms.append({"role": arm["role"], "status": status,
                     "import": imported, "search": searched})
        write_new(arm_dir / "outer.json", arms[-1])
        if status.startswith("PRODUCER_FAILURE"):
            break
    changed = (subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                       text=True).strip() != commit
               or bool(subprocess.check_output(["git", "status", "--porcelain"],
                                               cwd=ROOT, text=True)))
    status = ("PRODUCER_FAILURE_frozen_source_changed" if changed else
              "PRODUCER_FAILURE_arm" if any(a["status"].startswith("PRODUCER_FAILURE")
                                            for a in arms) else
              "BOUNDED_UNKNOWN_incomplete_pair" if len(arms) != 2 else
              "BOUNDED_PARTIAL_one_query_pair")
    result = {
        "schema": "n83.hosted-m6-paired-outer/v1",
        "status": status,
        "source_commit": commit, "run_id": run_id,
        "arms": arms, "elapsed_supervisor_wall_seconds_L0": time.monotonic() - started,
        "conservative_full_charge_seconds": 7500,
        "natural_partial_yield": None, "matrix_rank": None,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    write_new(output / "outer.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archives", type=Path, required=True)
    parser.add_argument("--worker", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--image", default=IMAGE)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    result = execute(args.archives.resolve(strict=True), args.worker.resolve(strict=True),
                     args.output, args.image, args.run_id)
    print(json.dumps({"status": result["status"], "arms": [a["status"] for a in result["arms"]]},
                     sort_keys=True))
    if result["status"].startswith("PRODUCER_FAILURE"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
