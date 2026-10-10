#!/usr/bin/env python3
"""Bound construction of one declared N83 v2 factor-base object.

This is a construction gate. Replay and S3 publication are separate commands.
The Linux exporter runs with a hard Docker memory limit, no swap or network,
and an outer wall deadline. Every launched outcome receives an outer receipt.
"""

from __future__ import annotations

import argparse
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
SCHEMA = "n83.v2-object-construction-outer/v1"
SIZES = (1182, 2048, 4096, 8192, 16627)
POLICIES = ("public_x_sequential", "public_x_hash", "public_x_gray_prefix")
SEEDS = (2026100801, 2026100802, 2026100803)
HEX64 = re.compile(r"[0-9a-f]{64}\Z")


def manifest_error(manifest: object, config: dict, object_dir: Path) -> str | None:
    if not isinstance(manifest, dict):
        return "manifest is not an object"
    expected = {
        "schema": "n83.factor-base-panel/v2-size-frontier",
        "study": STUDY,
        "status": "completed_factor_base_object",
        "source_commit": config["source_commit"],
        "budget_seconds": config["worker_wall_seconds"],
        "expected_base_count": 1,
        "completed_base_count": 1,
        "selected_best_total_runtime": None,
    }
    for key, value in expected.items():
        if key not in manifest or type(manifest[key]) is not type(value) or manifest[key] != value:
            return f"manifest {key} disagrees with frozen configuration"
    for key in ("source_blake3", "size_design_blake3"):
        digest = manifest.get(key)
        if not isinstance(digest, str) or HEX64.fullmatch(digest) is None:
            return f"manifest {key} is not a digest"
    rows = manifest.get("bases")
    if not isinstance(rows, list) or len(rows) != 1 or not isinstance(rows[0], dict):
        return "manifest needs exactly one base row"
    row = rows[0]
    for key, value in {
        "a": config["curve_a"], "policy": config["policy"],
        "columns": config["columns"], "seed": config["seed"],
        "points": 166 * config["columns"],
        "runtime_rank_eligible": False,
        "total_index_calculus_runtime_ms": None,
        "relation_yield": None, "matrix_rank": None, "verified_dlp": None,
    }.items():
        if key not in row or type(row[key]) is not type(value) or row[key] != value:
            return f"base row {key} disagrees with frozen configuration"
    for key in ("compressed_blake3", "plain_blake3", "point_set_blake3"):
        digest = row.get(key)
        if not isinstance(digest, str) or HEX64.fullmatch(digest) is None:
            return f"base row {key} is not a digest"
    object_name = f"objects/{row['compressed_blake3']}.jsonl.gz"
    if row.get("object") != object_name:
        return "object name is not bound to the compressed digest"
    uri = ("s3://crypto-autoresearcher/factor-bases/icv1/etc/"
           f"{STUDY}/v2-size-frontier/a{config['curve_a']}/{object_name}")
    if row.get("s3_uri") != uri:
        return "object destination is outside the declared v2 prefix"
    object_path = object_dir / object_name
    size = row.get("bytes")
    if type(size) is not int or size <= 0 or not object_path.is_file() or object_path.stat().st_size != size:
        return "compressed object is missing or has a wrong size"
    return None


def docker_command(config: dict, checkout: Path, repository_root: Path,
                   binary: Path, output_dir: Path) -> list[str]:
    memory_mib = config["memory_mib"]
    return [
        "docker", "run", "--rm", "--pull", "never", "--name", config["container_name"],
        "--network", "none", "--memory", f"{memory_mib}m",
        "--memory-swap", f"{memory_mib}m", "--cpus", "1", "--pids-limit", "64",
        "--cap-drop", "ALL", "--security-opt", "no-new-privileges", "--read-only",
        "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--env", "GIT_CONFIG_COUNT=1",
        "--env", "GIT_CONFIG_KEY_0=safe.directory",
        "--env", f"GIT_CONFIG_VALUE_0={checkout}",
        "--mount", f"type=bind,src={repository_root},dst={repository_root},readonly",
        "--mount", f"type=bind,src={binary},dst=/worker,readonly",
        "--mount", f"type=bind,src={output_dir},dst=/out",
        "--workdir", str(checkout), config["container_image_id"],
        "/worker", "v2-construct-one", "/out/object", str(config["curve_a"]),
        config["policy"], str(config["columns"]), str(config["seed"]),
        str(config["worker_wall_seconds"]),
    ]


def run(checkout: Path, repository_root: Path, binary: Path, output_dir: Path,
        curve_a: int, policy: str, columns: int, seed: int,
        wall_seconds: float, memory_mib: int,
        image: str = "python:3.11-slim") -> dict:
    if curve_a not in (0, 1) or policy not in POLICIES or columns not in SIZES or seed not in SEEDS:
        raise ValueError("object selection is outside the declared v2 grid")
    if not math.isfinite(wall_seconds) or not 1 <= wall_seconds <= 7200:
        raise ValueError("wall cap must be finite and within 1..7200 seconds")
    if not 128 <= memory_mib <= 65536:
        raise ValueError("memory cap must be within 128..65536 MiB")
    checkout = checkout.resolve(strict=True)
    repository_root = repository_root.resolve(strict=True)
    binary = binary.resolve(strict=True)
    if not checkout.is_relative_to(repository_root) or checkout == repository_root:
        raise ValueError("checkout must be inside the mounted repository root")
    if output_dir.resolve().is_relative_to(repository_root):
        raise ValueError("output must be outside the frozen repository")
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError("worker must be an executable Linux ELF binary")
    with binary.open("rb") as executable:
        if executable.read(4) != b"\x7fELF":
            raise ValueError("worker must be an executable Linux ELF binary")
    git_common = subprocess.run(
        ["git", "-C", str(checkout), "rev-parse", "--path-format=absolute", "--git-common-dir"],
        capture_output=True, text=True, check=True, timeout=15,
    ).stdout.strip()
    if not Path(git_common).resolve(strict=True).is_relative_to(repository_root):
        raise ValueError("git common directory is outside the read-only repository mount")
    commit = clean_commit(checkout)
    inspected = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Id}}"],
        capture_output=True, text=True, check=True, timeout=15,
    )
    image_id = inspected.stdout.strip()
    if not image_id.startswith("sha256:") or HEX64.fullmatch(image_id[7:]) is None:
        raise ValueError("container image has no immutable local ID")
    output_dir.mkdir(parents=True, exist_ok=False)
    output_dir = output_dir.resolve(strict=True)
    config = {
        "schema": "n83.v2-object-construction-config/v1", "study": STUDY,
        "source_commit": commit, "checkout": str(checkout),
        "repository_root": str(repository_root),
        "curve_a": curve_a, "policy": policy, "columns": columns, "seed": seed,
        "wall_seconds": wall_seconds, "worker_wall_seconds": math.ceil(wall_seconds),
        "memory_mib": memory_mib, "memory_cgroup_limit_bytes": memory_mib * 1024 * 1024,
        "memory_cgroup_swap_limit_bytes": 0, "cpu_limit": 1.0,
        "binary_sha256": sha256_file(binary), "container_image_id": image_id,
        "container_name": f"n83-v2-object-{uuid.uuid4().hex[:16]}",
        "s3_upload_executed": False, "selected_best_total_runtime": None,
    }
    write_new_json(output_dir / "config.json", config)
    command = docker_command(config, checkout, repository_root, binary, output_dir)
    started = time.monotonic()
    timed_out = False
    exit_code: int | None = None
    launch_error: str | None = None
    cleanup_error: str | None = None
    with (output_dir / "stdout.log").open("xb") as stdout, (output_dir / "stderr.log").open("xb") as stderr:
        try:
            child = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
            try:
                exit_code = child.wait(timeout=wall_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True
                try:
                    killed = subprocess.run(
                        ["docker", "kill", config["container_name"]],
                        capture_output=True, text=True, check=False, timeout=15,
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
    if timed_out:
        try:
            remaining = subprocess.run(
                ["docker", "ps", "-aq", "--filter", f"name=^{config['container_name']}$"],
                capture_output=True, text=True, check=False, timeout=15,
            )
            if remaining.returncode != 0 or remaining.stdout.strip():
                cleanup_error = cleanup_error or "timed-out container still present"
        except (OSError, subprocess.SubprocessError) as error:
            cleanup_error = cleanup_error or f"container cleanup unverified: {error}"
    object_dir = output_dir / "object"
    manifest_path = object_dir / "manifest.json"
    cap_path = object_dir / "cap.json"
    manifest: object = None
    cap: object = None
    parse_error: str | None = None
    for path, name in ((manifest_path, "manifest"), (cap_path, "cap")):
        if path.exists():
            try:
                loaded = json.loads(path.read_text(encoding="utf-8"))
                if name == "manifest":
                    manifest = loaded
                else:
                    cap = loaded
            except (OSError, json.JSONDecodeError) as error:
                parse_error = f"{name}: {error}"
    source_changed = False
    try:
        source_changed = clean_commit(checkout) != commit or sha256_file(binary) != config["binary_sha256"]
    except (OSError, ValueError, subprocess.SubprocessError):
        source_changed = True
    valid_cap = (isinstance(cap, dict) and
                 cap.get("schema") == "n83.factor-base-object-cap/v2-size-frontier" and
                 cap.get("status") == "UNKNOWN_budget" and
                 cap.get("curve_a") == curve_a and cap.get("policy") == policy and
                 cap.get("orbit_columns") == columns and cap.get("seed") == seed and
                 cap.get("budget_seconds") == config["worker_wall_seconds"])
    manifest_problem = manifest_error(manifest, config, object_dir)
    if cleanup_error:
        status = "PRODUCER_FAILURE_container_cleanup"
    elif source_changed:
        status = "PRODUCER_FAILURE_frozen_input_changed"
    elif timed_out or exit_code == 124 and valid_cap:
        status = "UNKNOWN_wall_cap"
    elif launch_error:
        status = "PRODUCER_FAILURE_launch"
    elif exit_code in (134, 137):
        status = "UNKNOWN_resource_or_worker_exit"
    elif exit_code != 0:
        status = "PRODUCER_FAILURE_worker_exit"
    elif parse_error or manifest_problem or cap_path.exists():
        status = "PRODUCER_FAILURE_manifest"
    else:
        status = "PASS_construction_only"
    outer = {
        "schema": SCHEMA, "status": status, "source_commit": commit,
        "config_sha256": sha256_file(output_dir / "config.json"),
        "manifest_sha256": sha256_file(manifest_path) if manifest_path.exists() else None,
        "stdout_sha256": sha256_file(output_dir / "stdout.log"),
        "stderr_sha256": sha256_file(output_dir / "stderr.log"),
        "worker_exit_code": exit_code, "launch_error": launch_error,
        "cleanup_error": cleanup_error, "parse_error": parse_error,
        "manifest_error": manifest_problem, "source_changed": source_changed,
        "process_wall_ms": (time.monotonic() - started) * 1000,
        "memory_cgroup_limit_bytes": config["memory_cgroup_limit_bytes"],
        "memory_cgroup_swap_limit_bytes": 0,
        "replay_executed": False, "s3_upload_executed": False,
        "total_index_calculus_runtime_ms": None, "selected_best_total_runtime": None,
    }
    write_new_json(output_dir / "outer.json", outer)
    return outer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--curve-a", type=int, required=True)
    parser.add_argument("--policy", choices=POLICIES, required=True)
    parser.add_argument("--columns", type=int, choices=SIZES, required=True)
    parser.add_argument("--seed", type=int, choices=SEEDS, required=True)
    parser.add_argument("--wall-seconds", type=float, required=True)
    parser.add_argument("--memory-mib", type=int, required=True)
    parser.add_argument("--image", default="python:3.11-slim")
    args = parser.parse_args()
    outer = run(args.checkout, args.repository_root, args.binary, args.output_dir,
                args.curve_a, args.policy, args.columns, args.seed,
                args.wall_seconds, args.memory_mib, args.image)
    print(json.dumps(outer, sort_keys=True))
    return 1 if outer["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
