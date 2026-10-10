#!/usr/bin/env python3
"""Bound one retained N83 compact point probe and retain every launched outcome.

The Linux worker checks its cgroup and compiled source against a small
supervisor-made snapshot. Docker supplies the hard memory and process-wall
limits. A PASS here certifies only a bounded diagnostic probe, never a
complete index-calculus runtime or a factor-base winner.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import uuid


STUDY = "koblitz_n83_factor_base_sweep_20261008"
SCHEMA = "n83.primary-compact-probe-outer/v1"
WORKER_SCHEMA = "n83.primary-compact-s3-probe/v1"
ATTESTATION_SCHEMA = "n83.primary-probe-source-attestation/v1"
POLICIES = ("public_x_sequential", "public_x_hash", "public_x_gray_prefix")
SEEDS = (2026100801, 2026100802, 2026100803)
SOURCES = (
    "examples/koblitz_n83_factor_base_export.rs",
    f"research/{STUDY}/primary_adapter.rs",
    f"research/{STUDY}/compact_cold.rs",
)
ORDER = "2417851639230796216685689"
HEX40 = re.compile(r"[0-9a-f]{40}\Z")
HEX64 = re.compile(r"[0-9a-f]{64}\Z")


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


def clean_commit(checkout: Path) -> str:
    status = subprocess.run(
        ["git", "-C", str(checkout), "status", "--porcelain", "--untracked-files=normal"],
        capture_output=True, text=True, timeout=15, check=True,
    )
    if status.stdout:
        raise ValueError("compact probe requires a clean frozen checkout")
    head = subprocess.run(
        ["git", "-C", str(checkout), "rev-parse", "HEAD"],
        capture_output=True, text=True, timeout=15, check=True,
    ).stdout.strip()
    if HEX40.fullmatch(head) is None:
        raise ValueError("compact probe has no full source commit")
    return head


def candidate_states(columns: int, pair_mode: str) -> int:
    pairs = columns * columns if pair_mode == "ordered" else columns * (columns + 1) // 2
    return 83 * pairs


def receipt_error(worker: object, config: dict) -> str | None:
    if not isinstance(worker, dict):
        return "worker receipt is not an object"
    expected = {
        "schema": WORKER_SCHEMA, "study": STUDY, "curve_a": 0, "fixture": 0,
        "orbit_columns": config["columns"], "policy": config["policy"],
        "seed": config["seed"], "max_candidate_states": config["max_candidate_states"],
        "memory_cgroup_limit_bytes": config["memory_cgroup_limit_bytes"],
        "memory_cgroup_swap_limit_bytes": 0, "source_commit": config["source_commit"],
        "source_attestation_mode": "supervised_snapshot", "rank_stage_executed": False,
        "column_log_verification": False, "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    for key, value in expected.items():
        if worker.get(key) != value:
            return f"worker {key} disagrees with frozen configuration"
    if worker.get("status") not in ("HIT", "MISS"):
        return "worker status is neither HIT nor MISS"
    for key in ("source_exporter_blake3", "source_adapter_blake3", "source_compact_blake3",
                "panel_manifest_blake3", "point_set_blake3", "public_corpus_canonical_json_blake3"):
        if not isinstance(worker.get(key), str) or HEX64.fullmatch(worker[key]) is None:
            return f"worker {key} is not a digest"
    for key in ("base_import_ms", "target_validation_ms", "point_probe_ms",
                "relation_row_ms", "process_wall_ms"):
        number = worker.get(key)
        if isinstance(number, bool) or not isinstance(number, (int, float)) or not math.isfinite(number) or number < 0:
            return f"worker {key} is not a nonnegative duration"
    peak = worker.get("memory_cgroup_peak_bytes")
    if type(peak) is not int or peak < 0:
        return "worker cgroup peak is invalid"
    probe = worker.get("probe")
    if not isinstance(probe, dict):
        return "worker probe is absent"
    for key, value in {
        "schema": "compact-four-sum-point-probe/v1", "status": worker["status"],
        "n": 83, "curve_a": 0, "orbit_columns": config["columns"],
        "unordered_pairs": config["pair_mode"] == "unordered",
        "candidate_states": config["max_candidate_states"],
        "rank_stage_executed": False, "total_index_calculus_runtime_ms": None,
    }.items():
        if probe.get(key) != value:
            return f"nested probe {key} disagrees with configuration"
    hit = worker["status"] == "HIT"
    if worker.get("relation_row_stage_executed") is not hit:
        return "relation-row stage disagrees with hit status"
    if hit:
        indices = probe.get("point_indices")
        if (not isinstance(indices, list) or len(indices) != 4 or
                any(type(index) is not int or not 0 <= index < 166 * config["columns"] for index in indices)):
            return "hit has invalid point indices"
        row = worker.get("relation_row")
        if (probe.get("group_verified") is not True or not isinstance(row, dict) or
                row.get("schema") != "n83.primary-wide-relation-row/v1" or
                row.get("modulus_decimal") != ORDER or
                row.get("orbit_columns") != config["columns"] or
                row.get("group_verified") is not True or
                not isinstance(row.get("dense_decimal_blake3"), str) or
                HEX64.fullmatch(row["dense_decimal_blake3"]) is None):
            return "hit lacks a checked full-width relation row"
        nonzero = row.get("nonzero")
        if not isinstance(nonzero, list) or not 1 <= len(nonzero) <= 4:
            return "hit row has invalid sparsity"
        prior = -1
        for entry in nonzero:
            if not isinstance(entry, dict):
                return "hit row entry is not an object"
            column, coefficient = entry.get("column"), entry.get("coefficient_decimal")
            if (type(column) is not int or not prior < column < config["columns"] or
                    not isinstance(coefficient, str) or not coefficient.isdecimal() or
                    not 0 < int(coefficient) < int(ORDER)):
                return "hit row has invalid column or coefficient"
            prior = column
    elif worker.get("relation_row") is not None or probe.get("point_indices") is not None:
        return "miss unexpectedly carries a relation"
    return None


def run(
    panel: Path, columns: int, policy: str, seed: int, pair_mode: str,
    wall_seconds: float, memory_mib: int, binary: Path, output_dir: Path,
    image: str = "python:3.11-slim", checkout: Path | None = None,
) -> dict:
    if columns not in (64, 256, 600) or policy not in POLICIES or seed not in SEEDS:
        raise ValueError("compact probe selection is outside the retained primary panel")
    if pair_mode not in ("ordered", "unordered"):
        raise ValueError("pair mode must be ordered or unordered")
    if not math.isfinite(wall_seconds) or not 0 < wall_seconds <= 86400:
        raise ValueError("wall cap must be finite and within 1..86400 seconds")
    if not 0 < memory_mib <= 65536:
        raise ValueError("memory cap must be within 1..65536 MiB")
    panel, binary = panel.resolve(strict=True), binary.resolve(strict=True)
    checkout = (checkout or Path(__file__).resolve().parents[2]).resolve(strict=True)
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError("compact probe worker must be an executable Linux ELF binary")
    with binary.open("rb") as executable:
        magic = executable.read(4)
    if magic != b"\x7fELF":
        raise ValueError("compact probe worker must be an executable Linux ELF binary")
    if output_dir.resolve().is_relative_to(checkout):
        raise ValueError("probe output directory must be outside the frozen checkout")
    commit = clean_commit(checkout)
    image_result = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Id}}"],
        capture_output=True, text=True, timeout=15, check=True,
    )
    image_id = image_result.stdout.strip()
    if not image_id.startswith("sha256:") or HEX64.fullmatch(image_id[7:]) is None:
        raise ValueError("container image has no immutable local ID")
    limit_bytes = memory_mib * 1024 * 1024
    states = candidate_states(columns, pair_mode)
    started = time.monotonic()
    output_dir.mkdir(parents=True, exist_ok=False)
    output_dir = output_dir.resolve(strict=True)
    source_dir = output_dir / "source"
    for relative in SOURCES:
        destination = source_dir / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(checkout / relative, destination)
    attestation = {"schema": ATTESTATION_SCHEMA, "source_commit": commit, "status_clean": True}
    write_new_json(source_dir / "attestation.json", attestation)
    source_hashes = {relative: sha256_file(source_dir / relative) for relative in SOURCES}
    source_hashes["attestation.json"] = sha256_file(source_dir / "attestation.json")
    if clean_commit(checkout) != commit:
        raise ValueError("checkout changed while freezing compact probe sources")
    config = {
        "schema": "n83.primary-compact-probe-config/v1", "study": STUDY,
        "panel_dir": str(panel), "columns": columns, "policy": policy, "seed": seed,
        "pair_mode": pair_mode, "max_candidate_states": states,
        "wall_seconds": wall_seconds, "memory_cgroup_limit_bytes": limit_bytes,
        "memory_cgroup_swap_limit_bytes": 0, "cpu_limit": 1.0,
        "binary": str(binary), "binary_sha256": sha256_file(binary),
        "source_commit": commit, "source_sha256": source_hashes,
        "container_image": image, "container_image_id": image_id,
        "container_name": f"n83-compact-probe-{uuid.uuid4().hex[:16]}",
        "rank_stage_executed": False, "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    write_new_json(output_dir / "config.json", config)
    config_hash = sha256_file(output_dir / "config.json")
    command = [
        "docker", "run", "--rm", "--pull", "never", "--name", config["container_name"],
        "--network", "none", "--memory", f"{memory_mib}m",
        "--memory-swap", f"{memory_mib}m", "--cpus", "1", "--pids-limit", "64",
        "--cap-drop", "ALL", "--security-opt", "no-new-privileges", "--read-only",
        "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m", "--user", f"{os.getuid()}:{os.getgid()}",
        "--mount", f"type=bind,src={panel},dst=/panel,readonly",
        "--mount", f"type=bind,src={binary},dst=/worker,readonly",
        "--mount", f"type=bind,src={source_dir},dst=/source,readonly",
        "--mount", f"type=bind,src={output_dir},dst=/out",
        "--env", "ICV1_FROZEN_SOURCE_DIR=/source", image_id,
        "/worker", "primary-s3-probe", "/panel", str(columns), policy, str(seed),
        pair_mode, str(states), str(memory_mib), "/out/worker.json",
    ]
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
                    stopped = subprocess.run(
                        ["docker", "kill", config["container_name"]],
                        capture_output=True, text=True, timeout=15, check=False,
                    )
                    if stopped.returncode != 0:
                        cleanup_error = stopped.stderr.strip() or "docker kill failed"
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
                capture_output=True, text=True, timeout=15, check=False,
            )
            if remaining.returncode != 0 or remaining.stdout.strip():
                cleanup_error = cleanup_error or "timed-out container still present"
        except (OSError, subprocess.SubprocessError) as error:
            cleanup_error = cleanup_error or f"container cleanup unverified: {error}"
    worker_path = output_dir / "worker.json"
    worker: object = None
    parse_error: str | None = None
    if worker_path.exists():
        try:
            worker = json.loads(worker_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            parse_error = f"{type(error).__name__}: {error}"
    try:
        source_changed = any(
            sha256_file(source_dir / relative) != digest
            for relative, digest in source_hashes.items()
        )
    except OSError:
        source_changed = True
    try:
        config_changed = sha256_file(output_dir / "config.json") != config_hash
    except OSError:
        config_changed = True
    try:
        binary_changed = sha256_file(binary) != config["binary_sha256"]
    except OSError:
        binary_changed = True
    if cleanup_error:
        status = "PRODUCER_FAILURE_container_cleanup"
    elif source_changed or config_changed or binary_changed:
        status = "PRODUCER_FAILURE_frozen_input_changed"
    elif timed_out:
        status = "UNKNOWN_wall_cap"
    elif launch_error:
        status = "PRODUCER_FAILURE_launch"
    elif exit_code in (134, 137):
        status = "UNKNOWN_resource_or_worker_exit"
    elif exit_code != 0:
        status = "PRODUCER_FAILURE_worker_exit"
    else:
        parse_error = parse_error or receipt_error(worker, config)
        status = "PASS_probe_hit" if isinstance(worker, dict) and worker.get("status") == "HIT" else "PASS_probe_miss"
        if parse_error:
            status = "PRODUCER_FAILURE_worker_receipt"
    outer = {
        "schema": SCHEMA, "status": status,
        "config_sha256": config_hash,
        "worker_sha256": sha256_file(worker_path) if worker_path.exists() else None,
        "stdout_sha256": sha256_file(output_dir / "stdout.log"),
        "stderr_sha256": sha256_file(output_dir / "stderr.log"),
        "worker_exit_code": exit_code, "launch_error": launch_error,
        "cleanup_error": cleanup_error, "worker_error": parse_error,
        "source_changed": source_changed,
        "config_changed": config_changed,
        "binary_changed": binary_changed,
        "process_wall_ms": (time.monotonic() - started) * 1000,
        "memory_cgroup_limit_bytes": limit_bytes, "memory_cgroup_swap_limit_bytes": 0,
        "worker_cgroup_peak_bytes": worker.get("memory_cgroup_peak_bytes") if isinstance(worker, dict) else None,
        "rank_stage_executed": False, "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    write_new_json(output_dir / "outer.json", outer)
    return outer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--columns", type=int, required=True)
    parser.add_argument("--policy", choices=POLICIES, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--pair-mode", choices=("ordered", "unordered"), required=True)
    parser.add_argument("--wall-seconds", type=float, required=True)
    parser.add_argument("--memory-mib", type=int, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--image", default="python:3.11-slim")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    outer = run(args.panel, args.columns, args.policy, args.seed, args.pair_mode,
                args.wall_seconds, args.memory_mib, args.binary, args.output_dir, args.image)
    print(json.dumps(outer, sort_keys=True))
    return 1 if outer["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
