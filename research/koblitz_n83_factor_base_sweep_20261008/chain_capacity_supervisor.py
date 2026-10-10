#!/usr/bin/env python3
"""Bound one exact retained-base chained-S3 model construction.

The worker checks a hard cgroup memory limit and its compiled source snapshot.
Docker supplies a wall cap, zero swap, and no network. Every launched outcome
is retained; a construction receipt never implies SAT search or a relation.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import uuid

from compact_probe_supervisor import clean_commit, sha256_file, write_new_json


STUDY = "koblitz_n83_factor_base_sweep_20261008"
SCHEMA = "n83.chain-s3-capacity-outer/v1"
WORKER_SCHEMA = "n83.chain-s3-capacity-worker/v1"
ATTESTATION_SCHEMA = "n83.chain-s3-source-attestation/v1"
POLICIES = ("public_x_sequential", "public_x_hash", "public_x_gray_prefix")
SEEDS = (2026100801, 2026100802, 2026100803)
V2_SIZES = (1182, 2048, 4096, 8192, 16627)
SOURCES = (
    "examples/koblitz_n83_factor_base_export.rs",
    f"research/{STUDY}/primary_adapter.rs",
    "src/cryptanalysis/binary_semaev_chain_sat.rs",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    f"research/{STUDY}/size-frontier-v2.json",
)
HEX64 = re.compile(r"[0-9a-f]{64}\Z")


def receipt_error(worker: object, config: dict) -> str | None:
    if not isinstance(worker, dict):
        return "worker receipt is not an object"
    expected = {
        "schema": WORKER_SCHEMA, "study": STUDY, "curve_a": 0, "fixture": 0,
        "orbit_columns": config["columns"], "policy": config["policy"],
        "seed": config["seed"], "summands": config["summands"],
        "identity_mask": config["identity_mask"],
        "max_variables": config["max_variables"],
        "max_domain_clauses": config["max_domain_clauses"],
        "memory_cgroup_limit_bytes": config["memory_cgroup_limit_bytes"],
        "memory_cgroup_swap_limit_bytes": 0,
        "source_commit": config["source_commit"],
        "source_attestation_mode": "supervised_snapshot",
        "solver_search_executed": False, "relation_stage_executed": False,
        "rank_stage_executed": False, "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    for key, value in expected.items():
        if worker.get(key) != value:
            return f"worker {key} disagrees with frozen configuration"
    status = worker.get("status")
    if status not in ("PASS_model_construction_only", "UNKNOWN_variable_cap", "UNKNOWN_domain_clause_cap"):
        return "worker status is outside the capacity contract"
    for key in (
        "source_exporter_blake3", "source_adapter_blake3", "source_chain_blake3",
        "source_index_calculus_blake3", "panel_manifest_blake3",
        "point_set_blake3", "public_corpus_canonical_json_blake3",
    ):
        digest = worker.get(key)
        if not isinstance(digest, str) or HEX64.fullmatch(digest) is None:
            return f"worker {key} is not a digest"
    for key in ("base_import_ms", "target_validation_ms", "model_construction_ms", "process_wall_ms"):
        value = worker.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            return f"worker {key} is not a nonnegative duration"
    for key in ("memory_cgroup_peak_bytes", "required_max_variables", "required_domain_clauses"):
        value = worker.get(key)
        if type(value) is not int or value < 0:
            return f"worker {key} is not a nonnegative count"
    if worker["required_max_variables"] != {5: 84743, 6: 105908}[config["summands"]]:
        return "worker source variable count disagrees with the pinned arity"
    if worker.get("legal_x_coordinates") != 83 * config["columns"] or worker["required_domain_clauses"] == 0:
        return "worker coordinate domain disagrees with retained base"
    if not isinstance(worker.get("object"), str) or not worker["object"]:
        return "worker has no selected base object"
    counts = ("sat_variables", "sat_clauses", "sat_xor_rows", "sat_and_gates", "sat_s3_nodes")
    if status == "UNKNOWN_variable_cap":
        if worker["required_max_variables"] <= config["max_variables"]:
            return "variable cap is not exceeded"
    elif status == "UNKNOWN_domain_clause_cap":
        if (worker["required_max_variables"] > config["max_variables"] or
                worker["required_domain_clauses"] <= config["max_domain_clauses"]):
            return "domain cap is not the first exceeded cap"
    else:
        if (worker["required_max_variables"] > config["max_variables"] or
                worker["required_domain_clauses"] > config["max_domain_clauses"]):
            return "successful model exceeded a declared preflight cap"
        for key in counts:
            if type(worker.get(key)) is not int or worker[key] < 0:
                return f"successful model has invalid {key}"
        if (worker["sat_variables"] > config["max_variables"] or
                worker["sat_clauses"] < 3 * worker["sat_and_gates"] or
                not 0 <= worker["sat_s3_nodes"] <= config["summands"] - 1):
            return "successful model counts are inconsistent"
    if status != "PASS_model_construction_only" and any(worker.get(key) is not None for key in counts):
        return "capped model unexpectedly has installed SAT counts"
    return None


def run(
    panel: Path, columns: int, policy: str, seed: int, summands: int,
    identity_mask: int, max_variables: int, max_domain_clauses: int,
    wall_seconds: float, memory_mib: int, binary: Path, output_dir: Path,
    image: str = "python:3.11-slim", checkout: Path | None = None,
) -> dict:
    if columns not in (64, 256, 600, *V2_SIZES) or policy not in POLICIES or seed not in SEEDS:
        raise ValueError("chain selection is outside retained primary panel")
    if summands not in (5, 6) or not 0 <= identity_mask < 1 << (summands - 2):
        raise ValueError("chain arity or identity pattern is unsupported")
    if not 0 < max_variables <= 2**31 - 1 or not 0 < max_domain_clauses <= 2**31 - 1:
        raise ValueError("chain variable or domain cap is invalid")
    if not math.isfinite(wall_seconds) or not 0 < wall_seconds <= 86400:
        raise ValueError("wall cap must be finite and within 1..86400 seconds")
    if not 0 < memory_mib <= 65536:
        raise ValueError("memory cap must be within 1..65536 MiB")
    panel, binary = panel.resolve(strict=True), binary.resolve(strict=True)
    checkout = (checkout or Path(__file__).resolve().parents[2]).resolve(strict=True)
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError("worker must be an executable Linux ELF binary")
    with binary.open("rb") as executable:
        if executable.read(4) != b"\x7fELF":
            raise ValueError("worker must be a Linux ELF binary")
    if output_dir.resolve().is_relative_to(checkout):
        raise ValueError("capacity output must be outside the frozen checkout")
    commit = clean_commit(checkout)
    inspected = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Id}}"],
        capture_output=True, text=True, timeout=15, check=True,
    )
    image_id = inspected.stdout.strip()
    if not image_id.startswith("sha256:") or HEX64.fullmatch(image_id[7:]) is None:
        raise ValueError("container image has no immutable local ID")
    started = time.monotonic()
    output_dir.mkdir(parents=True, exist_ok=False)
    output_dir = output_dir.resolve(strict=True)
    source_dir = output_dir / "source"
    for relative in SOURCES:
        destination = source_dir / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(checkout / relative, destination)
    write_new_json(source_dir / "attestation.json", {
        "schema": ATTESTATION_SCHEMA, "source_commit": commit, "status_clean": True,
    })
    source_hashes = {relative: sha256_file(source_dir / relative) for relative in SOURCES}
    source_hashes["attestation.json"] = sha256_file(source_dir / "attestation.json")
    if clean_commit(checkout) != commit:
        raise ValueError("checkout changed while freezing chained-S3 sources")
    config = {
        "schema": "n83.chain-s3-capacity-config/v1", "study": STUDY,
        "panel_dir": str(panel), "columns": columns, "policy": policy, "seed": seed,
        "summands": summands, "identity_mask": identity_mask,
        "max_variables": max_variables, "max_domain_clauses": max_domain_clauses,
        "wall_seconds": wall_seconds, "memory_cgroup_limit_bytes": memory_mib * 1024 * 1024,
        "memory_cgroup_swap_limit_bytes": 0, "cpu_limit": 1.0,
        "binary": str(binary), "binary_sha256": sha256_file(binary),
        "source_commit": commit, "source_sha256": source_hashes,
        "container_image": image, "container_image_id": image_id,
        "container_name": f"n83-chain-capacity-{uuid.uuid4().hex[:16]}",
        "solver_search_executed": False, "selected_best_total_runtime": None,
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
        "/worker", "primary-chain-build", "/panel", str(columns), policy, str(seed),
        str(summands), str(identity_mask), str(max_variables), str(max_domain_clauses),
        str(memory_mib), "/out/worker.json",
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
    worker_error: str | None = None
    if worker_path.exists():
        try:
            worker = json.loads(worker_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            worker_error = f"{type(error).__name__}: {error}"
    try:
        source_changed = any(sha256_file(source_dir / path) != digest for path, digest in source_hashes.items())
        config_changed = sha256_file(output_dir / "config.json") != config_hash
        binary_changed = sha256_file(binary) != config["binary_sha256"]
    except OSError:
        source_changed = config_changed = binary_changed = True
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
        worker_error = worker_error or receipt_error(worker, config)
        status = worker["status"] if isinstance(worker, dict) and worker_error is None else "PRODUCER_FAILURE_worker_receipt"
    outer = {
        "schema": SCHEMA, "status": status, "config_sha256": config_hash,
        "worker_sha256": sha256_file(worker_path) if worker_path.exists() else None,
        "stdout_sha256": sha256_file(output_dir / "stdout.log"),
        "stderr_sha256": sha256_file(output_dir / "stderr.log"),
        "worker_exit_code": exit_code, "launch_error": launch_error,
        "cleanup_error": cleanup_error, "worker_error": worker_error,
        "source_changed": source_changed, "config_changed": config_changed,
        "binary_changed": binary_changed,
        "process_wall_ms": (time.monotonic() - started) * 1000,
        "memory_cgroup_limit_bytes": config["memory_cgroup_limit_bytes"],
        "memory_cgroup_swap_limit_bytes": 0,
        "worker_cgroup_peak_bytes": worker.get("memory_cgroup_peak_bytes") if isinstance(worker, dict) else None,
        "solver_search_executed": False, "total_index_calculus_runtime_ms": None,
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
    parser.add_argument("--summands", type=int, choices=(5, 6), required=True)
    parser.add_argument("--identity-mask", type=int, required=True)
    parser.add_argument("--max-variables", type=int, required=True)
    parser.add_argument("--max-domain-clauses", type=int, required=True)
    parser.add_argument("--wall-seconds", type=float, required=True)
    parser.add_argument("--memory-mib", type=int, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--image", default="python:3.11-slim")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    outer = run(
        args.panel, args.columns, args.policy, args.seed, args.summands,
        args.identity_mask, args.max_variables, args.max_domain_clauses,
        args.wall_seconds, args.memory_mib, args.binary, args.output_dir, args.image,
    )
    print(json.dumps(outer, sort_keys=True))
    return 1 if outer["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
