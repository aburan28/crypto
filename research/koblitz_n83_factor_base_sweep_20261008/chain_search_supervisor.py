#!/usr/bin/env python3
"""Run one retained primary N83 chained-S3 search under hard resource caps.

This produces a stage receipt, not a fully charged cold-runtime comparison.
The worker checks its cgroup and source snapshot; the outer process limits
wall time, disables network, and retains every launched outcome.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
import uuid

from chain_capacity_supervisor import (
    ATTESTATION_SCHEMA, HEX64, POLICIES, SEEDS, SOURCES,
    clean_commit, sha256_file, write_new_json,
)


STUDY = "koblitz_n83_factor_base_sweep_20261008"
SCHEMA = "n83.chain-s3-search-outer/v1"
WORKER_SCHEMA = "n83.primary-cold-result/v1"
WORKER_CONFIG_SCHEMA = "n83.primary-cold-config/v1"
SIZES = (64, 256, 600, 1182, 2048, 4096, 8192, 16627)
COUNTERS = (
    "factor_base_size", "orbit_count", "trials", "relations",
    "independent_relations", "dependent_relations", "inconsistent_relations",
    "verification_failures", "sat_unknowns", "sat_invalid_models",
    "sat_group_rejected_models", "sat_calls", "sat_refutations",
    "sat_models", "sat_conflicts", "linear_solve_attempts",
    "direct_relations_skipped",
)


def receipt_error(worker: object, worker_config: object, config: dict) -> str | None:
    if not isinstance(worker, dict) or not isinstance(worker_config, dict):
        return "worker summary or config is not a JSON object"
    limits = {
        "max_variables": config["max_variables"],
        "max_domain_clauses": config["max_domain_clauses"],
        "max_models": config["max_models"],
        "conflict_budget": config["conflict_budget"],
    }
    expected = {
        "schema": WORKER_CONFIG_SCHEMA, "study": STUDY,
        "curve_a": 0, "fixture": 0, "orbit_columns": config["columns"],
        "summands": config["summands"], "strategy": "chain-s3",
        "policy": config["policy"], "seed": config["seed"],
        "max_trials": config["max_trials"],
        "budget_seconds": config["worker_wall_seconds"],
        "allow_direct_relation": False,
        "source_commit": config["source_commit"],
        "source_attestation_mode": "supervised_snapshot",
        "fixture_dir": "/fixtures",
        "chain_s3_limits": limits,
        "memory_cgroup_limit_bytes": config["memory_cgroup_limit_bytes"],
        "memory_cgroup_swap_limit_bytes": 0,
    }
    for key, value in expected.items():
        if key not in worker_config or type(worker_config[key]) is not type(value) or worker_config[key] != value:
            return f"worker config {key} disagrees with frozen inputs"
    adapter_digest = worker_config.get("source_adapter_blake3")
    if not isinstance(adapter_digest, str) or HEX64.fullmatch(adapter_digest) is None:
        return "worker config adapter source is not a digest"
    expected["schema"] = WORKER_SCHEMA
    expected.pop("allow_direct_relation")
    expected.pop("memory_cgroup_swap_limit_bytes")
    for key, value in expected.items():
        if key not in worker or type(worker[key]) is not type(value) or worker[key] != value:
            return f"worker summary {key} disagrees with frozen inputs"
    if worker.get("memory_cgroup_swap_limit_bytes") != 0:
        return "worker did not report zero swap"
    status = worker.get("status")
    if status not in ("PREFLIGHT_ONLY", "UNKNOWN_solver_cap", "UNKNOWN_trial_cap",
                      "INADMISSIBLE_cofactor_class", "PASS_verified_target_only"):
        return "worker status is outside the search receipt contract"
    for key in ("point_set_blake3", "panel_manifest_blake3",
                "public_corpus_canonical_json_blake3", "source_chain_blake3",
                "source_index_calculus_blake3"):
        digest = worker.get(key)
        if not isinstance(digest, str) or HEX64.fullmatch(digest) is None:
            return f"worker {key} is not a digest"
    if not isinstance(worker.get("object"), str) or not worker["object"].startswith("objects/"):
        return "worker selected object is invalid"
    peak = worker.get("memory_cgroup_peak_bytes")
    if type(peak) is not int or not 0 <= peak <= config["memory_cgroup_limit_bytes"]:
        return "worker cgroup peak is invalid"
    for key in ("base_import_ms", "target_validation_ms", "solver_ms",
                "post_solver_validation_ms", "pre_summary_process_wall_ms"):
        value = worker.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            return f"worker {key} is not a nonnegative duration"
    if (worker.get("column_log_verification") is not False or
            worker.get("total_index_calculus_runtime_ms") is not None or
            worker.get("selected_best_total_runtime") is not None):
        return "worker claims unmeasured completion fields"
    report = worker.get("report")
    if not isinstance(report, dict):
        return "worker report is absent"
    for key in COUNTERS:
        if type(report.get(key)) is not int or report[key] < 0:
            return f"worker report {key} is not a nonnegative counter"
    if (report["trials"] > config["max_trials"] or
            report["orbit_count"] != config["columns"] or
            report["factor_base_size"] != 166 * config["columns"] or
            report["independent_relations"] > report["orbit_count"] or
            report["independent_relations"] + report["dependent_relations"] != report["relations"] or
            report["sat_refutations"] > report["sat_calls"] or
            report["sat_unknowns"] > report["sat_calls"] or
            report["sat_models"] > report["sat_calls"] or
            report["sat_group_rejected_models"] > report["sat_models"] or
            report["inconsistent_relations"] != 0 or
            report["verification_failures"] != 0 or
            report["sat_invalid_models"] != 0):
        return "worker rank, relation, or verification counts are inconsistent"
    for key in ("relation_collection_ns", "linear_algebra_ns"):
        value = report.get(key)
        if not isinstance(value, str) or not value.isdecimal():
            return f"worker {key} is not an integer nanosecond count"
    if worker.get("solver_stage_executed") is not (report["sat_calls"] > 0):
        return "worker SAT-stage flag disagrees with solver calls"
    if status == "PREFLIGHT_ONLY":
        if (config["max_trials"] != 0 or report["trials"] != 0 or
                report["sat_calls"] != 0 or report["relations"] != 0 or
                report["linear_solve_attempts"] != 0):
            return "preflight unexpectedly searched"
    elif config["max_trials"] == 0:
        return "zero-trial run has a search disposition"
    if status == "UNKNOWN_solver_cap" and report["sat_unknowns"] == 0:
        return "solver-cap disposition has no inconclusive SAT call"
    if status == "UNKNOWN_trial_cap" and report["sat_unknowns"] != 0:
        return "trial-cap disposition hides an inconclusive SAT call"
    verified = status == "PASS_verified_target_only"
    log = worker.get("verified_log")
    if verified:
        if (not isinstance(log, str) or not log.isdecimal() or
                report["independent_relations"] == 0 or report["linear_solve_attempts"] == 0 or
                not worker["solver_stage_executed"]):
            return "verified target lacks rank, linear solve, or scalar receipt"
    elif log is not None:
        return "nonverified run contains a target scalar"
    return None


def run(panel: Path, fixtures: Path, columns: int, policy: str, seed: int, summands: int,
        max_trials: int, max_variables: int, max_domain_clauses: int,
        max_models: int, conflict_budget: int, wall_seconds: float,
        memory_mib: int, binary: Path, output_dir: Path,
        image: str = "python:3.11-slim", checkout: Path | None = None) -> dict:
    if columns not in SIZES or policy not in POLICIES or seed not in SEEDS:
        raise ValueError("search selection is outside the retained primary panel")
    if summands not in (5, 6) or not 0 <= max_trials <= 1000:
        raise ValueError("search arity or trial cap is unsupported")
    if (not 0 < max_variables <= 2**31 - 1 or
            not 0 < max_domain_clauses <= 2**31 - 1 or
            not 0 < max_models <= 1_000_000 or
            not 0 < conflict_budget <= 10_000_000_000):
        raise ValueError("chain SAT cap is outside the worker contract")
    if not math.isfinite(wall_seconds) or not 1 <= wall_seconds <= 7200:
        raise ValueError("wall cap must be finite and within 1..7200 seconds")
    if not 128 <= memory_mib <= 65536:
        raise ValueError("memory cap must be within 128..65536 MiB")
    panel, fixtures, binary = (
        panel.resolve(strict=True), fixtures.resolve(strict=True), binary.resolve(strict=True)
    )
    panel_files = ("manifest.json", "replay.json", "upload-receipt.json")
    fixture_files = ("probe-corpus.json", "probe-validation.json")
    if not all((panel / name).is_file() for name in panel_files):
        raise ValueError("panel needs manifest, replay, and S3 round-trip receipts")
    if not all((fixtures / name).is_file() for name in fixture_files):
        raise ValueError("fixture directory needs the public corpus and validation sidecar")
    checkout = (checkout or Path(__file__).resolve().parents[2]).resolve(strict=True)
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError("search worker must be an executable Linux ELF binary")
    with binary.open("rb") as executable:
        if executable.read(4) != b"\x7fELF":
            raise ValueError("search worker must be a Linux ELF binary")
    if output_dir.resolve().is_relative_to(checkout):
        raise ValueError("search output must be outside the frozen checkout")
    commit = clean_commit(checkout)
    inspected = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Id}}"],
        capture_output=True, text=True, check=True, timeout=15,
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
        raise ValueError("checkout changed while freezing search sources")
    config = {
        "schema": "n83.chain-s3-search-config/v1", "study": STUDY,
        "panel_dir": str(panel), "fixture_dir": str(fixtures),
        "panel_sha256": {name: sha256_file(panel / name) for name in panel_files},
        "fixture_sha256": {name: sha256_file(fixtures / name) for name in fixture_files},
        "columns": columns, "policy": policy,
        "seed": seed, "summands": summands, "max_trials": max_trials,
        "max_variables": max_variables, "max_domain_clauses": max_domain_clauses,
        "max_models": max_models, "conflict_budget": conflict_budget,
        "wall_seconds": wall_seconds, "worker_wall_seconds": math.ceil(wall_seconds),
        "memory_cgroup_limit_bytes": memory_mib * 1024 * 1024,
        "memory_cgroup_swap_limit_bytes": 0, "cpu_limit": 1.0,
        "binary": str(binary), "binary_sha256": sha256_file(binary),
        "source_commit": commit, "source_sha256": source_hashes,
        "container_image": image, "container_image_id": image_id,
        "container_name": f"n83-chain-search-{uuid.uuid4().hex[:16]}",
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
        "--mount", f"type=bind,src={fixtures},dst=/fixtures,readonly",
        "--mount", f"type=bind,src={binary},dst=/worker,readonly",
        "--mount", f"type=bind,src={source_dir},dst=/source,readonly",
        "--mount", f"type=bind,src={output_dir},dst=/out",
        "--env", "ICV1_FROZEN_SOURCE_DIR=/source", image_id,
        "/worker", "primary-chain-cold", "/panel", "/fixtures", str(columns), policy, str(seed),
        str(summands), str(max_trials), str(max_variables), str(max_domain_clauses),
        str(max_models), str(conflict_budget), str(config["worker_wall_seconds"]),
        str(memory_mib), "/out/run",
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
                        capture_output=True, text=True, check=False, timeout=15,
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
                capture_output=True, text=True, check=False, timeout=15,
            )
            if remaining.returncode != 0 or remaining.stdout.strip():
                cleanup_error = cleanup_error or "timed-out container still present"
        except (OSError, subprocess.SubprocessError) as error:
            cleanup_error = cleanup_error or f"container cleanup unverified: {error}"
    run_dir = output_dir / "run"
    worker_config_path = run_dir / "config.json"
    worker_path = run_dir / "summary.json"
    cap_path = run_dir / "cap.json"
    parse_error: str | None = None
    parsed: dict[str, object] = {}
    for key, path in (("worker_config", worker_config_path), ("worker", worker_path), ("cap", cap_path)):
        if path.exists():
            try:
                parsed[key] = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as error:
                parse_error = f"{key}: {error}"
    worker = parsed.get("worker")
    worker_config = parsed.get("worker_config")
    cap = parsed.get("cap")
    try:
        source_changed = any(sha256_file(source_dir / path) != digest for path, digest in source_hashes.items())
        config_changed = sha256_file(output_dir / "config.json") != config_hash
        binary_changed = sha256_file(binary) != config["binary_sha256"]
        input_changed = (
            any(sha256_file(panel / name) != digest for name, digest in config["panel_sha256"].items()) or
            any(sha256_file(fixtures / name) != digest for name, digest in config["fixture_sha256"].items())
        )
    except OSError:
        source_changed = config_changed = binary_changed = input_changed = True
    valid_cap = (isinstance(cap, dict) and cap.get("schema") == "n83.primary-cold-cap/v1" and
                 cap.get("status") == "UNKNOWN_budget" and
                 cap.get("orbit_columns") == columns and cap.get("summands") == summands and
                 cap.get("strategy") == "chain-s3" and cap.get("max_trials") == max_trials and
                 cap.get("budget_seconds") == config["worker_wall_seconds"])
    worker_error = receipt_error(worker, worker_config, config)
    if cleanup_error:
        status = "PRODUCER_FAILURE_container_cleanup"
    elif source_changed or config_changed or binary_changed or input_changed:
        status = "PRODUCER_FAILURE_frozen_input_changed"
    elif timed_out or exit_code == 124 and valid_cap:
        status = "UNKNOWN_wall_cap"
    elif launch_error:
        status = "PRODUCER_FAILURE_launch"
    elif exit_code in (134, 137):
        status = "UNKNOWN_resource_or_worker_exit"
    elif exit_code != 0:
        status = "PRODUCER_FAILURE_worker_exit"
    elif parse_error or cap_path.exists() or worker_error:
        status = "PRODUCER_FAILURE_worker_receipt"
    else:
        status = worker["status"]
    outer = {
        "schema": SCHEMA, "status": status, "config_sha256": config_hash,
        "worker_config_sha256": sha256_file(worker_config_path) if worker_config_path.exists() else None,
        "worker_sha256": sha256_file(worker_path) if worker_path.exists() else None,
        "stdout_sha256": sha256_file(output_dir / "stdout.log"),
        "stderr_sha256": sha256_file(output_dir / "stderr.log"),
        "worker_exit_code": exit_code, "launch_error": launch_error,
        "cleanup_error": cleanup_error, "parse_error": parse_error,
        "worker_error": worker_error, "source_changed": source_changed,
        "config_changed": config_changed, "binary_changed": binary_changed,
        "input_changed": input_changed,
        "process_wall_ms": (time.monotonic() - started) * 1000,
        "memory_cgroup_limit_bytes": config["memory_cgroup_limit_bytes"],
        "memory_cgroup_swap_limit_bytes": 0,
        "worker_cgroup_peak_bytes": worker.get("memory_cgroup_peak_bytes") if isinstance(worker, dict) else None,
        "report": worker.get("report") if isinstance(worker, dict) and worker_error is None else None,
        "total_index_calculus_runtime_ms": None, "selected_best_total_runtime": None,
    }
    write_new_json(output_dir / "outer.json", outer)
    return outer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--columns", type=int, choices=SIZES, required=True)
    parser.add_argument("--policy", choices=POLICIES, required=True)
    parser.add_argument("--seed", type=int, choices=SEEDS, required=True)
    parser.add_argument("--summands", type=int, choices=(5, 6), required=True)
    parser.add_argument("--max-trials", type=int, required=True)
    parser.add_argument("--max-variables", type=int, required=True)
    parser.add_argument("--max-domain-clauses", type=int, required=True)
    parser.add_argument("--max-models", type=int, required=True)
    parser.add_argument("--conflict-budget", type=int, required=True)
    parser.add_argument("--wall-seconds", type=float, required=True)
    parser.add_argument("--memory-mib", type=int, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--image", default="python:3.11-slim")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    outer = run(args.panel, args.fixtures, args.columns, args.policy, args.seed, args.summands,
                args.max_trials, args.max_variables, args.max_domain_clauses,
                args.max_models, args.conflict_budget, args.wall_seconds,
                args.memory_mib, args.binary, args.output_dir, args.image)
    print(json.dumps(outer, sort_keys=True))
    return 1 if outer["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
