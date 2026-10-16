#!/usr/bin/env python3
"""Independently classify one retained v7 known-answer control receipt set.

This reads existing files only. It never launches the worker or a solver.
An unvalidated producer result is preserved as UNVALIDATED, never promoted.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STUDY = "koblitz_n83_factor_base_sweep_20261008"
sys.path.insert(0, str(ROOT / "research" / STUDY))
from chain_search_supervisor import parse_point, point_add  # noqa: E402


K600_AUDIT = "387d393b9805f35ee1ca6328e8602133a1b9bc036088b1be1a3ba5a29bb8469c"
CROSS_AUDIT = "dfc96f11e4a399a42559847dd39e94b753f37911a54d806e0e54e830e7735478"
BINARY = "26dcb2b49388fcb89853eabbda01db71d148620cd465bd5458099e0606dcf606"
BUILD = "cd302908d34d52412e330fb987e605e288af69d5bf4b50ebbea968e7a2112804"
SMALL = "eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380"
ENVELOPE = "aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a"
FIXTURE_BLAKE3 = "7e257cf3f278051094e27bd6df30ac8e78a0f55739b11f8b6837f23e5637703d"
FIXTURE_SHA256 = "fc203db28eed0164e47b0c755c3d9a09f80df1acd5e4c7e09dfe172fa5f0a538"
IMAGE = "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31"
FIXTURE = (ROOT / "research" / STUDY / "verification" /
           "mixed-orbit-known-witness-source-20261012/fixture.json")
KNOWN_STATUSES = {"PASS_known_witness_solver_access", "UNKNOWN_preflight_cap",
                  "UNKNOWN_search_cap", "PRODUCER_FAILURE_known_witness_refuted",
                  "PRODUCER_FAILURE_invalid_model"}


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def record(path: Path) -> dict:
    row = json.loads(path.read_text())
    require(isinstance(row, dict), f"{path.name} is not an object")
    return row


def exact(row: dict, expected: dict, label: str) -> None:
    for key, value in expected.items():
        require(type(row.get(key)) is type(value) and row[key] == value,
                f"{label}.{key} differs")


def option(command: list[str], flag: str) -> str:
    require(flag in command and command.index(flag) + 1 < len(command),
            f"missing Docker {flag}")
    return command[command.index(flag) + 1]


def candidate_sum_matches(points: list[tuple[int, int]], target: tuple[int, int]) -> bool:
    result = None
    for point in points:
        result = point_add(result, point)
    return result == target


def candidate_from_archives(candidate: dict, fixture: dict, small: Path,
                            envelope: Path) -> bool:
    indices = candidate.get("small_point_indices")
    residual = candidate.get("residual_points")
    if (candidate.get("a") != "0" or candidate.get("b") != "1"
            or not isinstance(indices, list) or len(indices) != 3
            or any(type(index) is not int or not 0 <= index < 196212 for index in indices)
            or not isinstance(residual, list) or len(residual) != 2):
        return False
    try:
        needed = set(indices)
        small_points: dict[int, tuple[int, int]] = {}
        with gzip.open(small, "rt", encoding="utf-8") as stream:
            next(stream)
            for index, line in enumerate(stream):
                if index in needed:
                    point = parse_point(json.loads(line)["point"])
                    require(point is not None, "small point is infinity")
                    small_points[index] = point
        domain = set()
        with gzip.open(envelope, "rt", encoding="utf-8") as stream:
            next(stream)
            for index, line in enumerate(stream):
                if index >= 196212:
                    point = parse_point(json.loads(line)["point"])
                    require(point is not None, "residual point is infinity")
                    domain.add(point)
        points = [parse_point(value) for value in residual]
        target = parse_point(fixture["equation"]["target"])
        return (set(small_points) == needed and target is not None
                and all(point is not None and point in domain for point in points)
                and candidate_sum_matches([small_points[index] for index in indices]
                                          + points, target))
    except (OSError, EOFError, UnicodeError, StopIteration, KeyError, TypeError,
            ValueError, json.JSONDecodeError):
        return False


def classify_worker(config: dict, result: dict, fixture: dict,
                    small: Path, envelope: Path) -> str:
    exact(config, {
        "schema": "n83.mixed-orbit-known-witness-config/v1",
        "scope": "known_answer_control_only", "study": STUDY,
        "fixture_blake3": FIXTURE_BLAKE3,
        "small_object_blake3": fixture["archive"]["small"]["s3_uri"].split("/")[-1][:-9],
        "envelope_object_blake3": fixture["archive"]["envelope"]["s3_uri"].split("/")[-1][:-9],
        "small_summands": 3, "residual_summands": 2,
        "budget_seconds": 840, "conflict_budget": 5000, "max_models": 32,
        "memory_cgroup_limit_bytes": 4294967296,
        "memory_cgroup_swap_limit_bytes": 0,
    }, "worker config")
    exact(result, {"schema": "n83.mixed-orbit-known-witness-result/v1",
                   "scope": "known_answer_control_only",
                   "target": fixture["equation"]["target"],
                   "natural_relations": None,
                   "total_index_calculus_runtime_ms": None,
                   "selected_best_total_runtime": None}, "worker result")
    status, stats = result.get("status"), result.get("search")
    require(status in KNOWN_STATUSES and isinstance(stats, dict), "worker status or search missing")
    for key, cap in (("conflicts", 5000), ("models", 32),
                     ("invalid_models", 256), ("group_rejected_models", 256),
                     ("solver_calls", 8), ("patterns_attempted", 8)):
        value = stats.get(key)
        require(type(value) is int and 0 <= value <= cap,
                f"worker {key} outside cap")
    exact(stats, {"domain_clauses": 465738}, "worker search")
    require(type(stats.get("exhausted")) is bool and type(stats.get("refuted")) is bool,
            "worker disposition malformed")
    for key, cap in (("peak_variables", 101636), ("peak_xors", 4980)):
        value = stats.get(key)
        require(type(value) is int and 0 <= value <= cap, f"worker {key} malformed")
    require(type(stats.get("peak_clauses")) is int and stats["peak_clauses"] >= 0,
            "worker peak clauses malformed")
    for key in ("base_import_elapsed_ms_L0", "worker_wall_ms_L0"):
        value = result.get(key)
        require(type(value) in (int, float) and math.isfinite(value)
                and 0 <= value <= 840000, f"worker {key} malformed")
    require(result["base_import_elapsed_ms_L0"] <= result["worker_wall_ms_L0"],
            "import exceeds worker wall")
    peak = result.get("memory_cgroup_peak_bytes")
    require(type(peak) is int and 0 < peak <= 4294967296,
            "worker memory peak malformed")
    candidate = result.get("candidate")
    if status == "PASS_known_witness_solver_access":
        require(isinstance(candidate, dict) and stats["models"] > 0
                and not stats["exhausted"] and not stats["refuted"]
                and candidate_from_archives(candidate, fixture, small, envelope),
                "candidate lacks independent archive/group replay")
        return "VALIDATED_CONTROL_ACCESS"
    require(candidate is None, "non-success status hides candidate")
    if status == "UNKNOWN_preflight_cap":
        require(stats["solver_calls"] == stats["patterns_attempted"] == 0
                and stats["exhausted"], "preflight cap inconsistent")
    elif status == "UNKNOWN_search_cap":
        require(stats["solver_calls"] > 0 and stats["exhausted"],
                "search cap inconsistent")
    elif status == "PRODUCER_FAILURE_known_witness_refuted":
        require(stats["refuted"], "known witness refutation lacks flag")
    else:
        require(stats["invalid_models"] > 0, "invalid model lacks counter")
    return "VALIDATED_CENSORED" if status.startswith("UNKNOWN_") else "VALIDATED_PRODUCER_FAILURE"


def check_stage(stage: dict, name: str, run: Path) -> None:
    exact(stage, {"name": name, "exit_code": 0, "timed_out": False,
                  "launch_error": None, "cleanup_error": None,
                  "cleanup_verified": True}, "stage")
    state = stage.get("container_state")
    require(isinstance(state, dict), "container state missing")
    exact(state, {"ExitCode": 0, "OOMKilled": False,
                  "Running": False, "Status": "exited"}, "container")
    wall = stage.get("wall_elapsed_seconds_L0")
    require(type(wall) in (int, float) and math.isfinite(wall)
            and 0 <= wall <= (90 if name == "verify" else 850),
            "stage wall cap differs")
    command = stage.get("command")
    require(isinstance(command, list) and command[:2] == ["docker", "run"],
            "stage Docker command missing")
    for flag, expected in (("--pull", "never"), ("--platform", "linux/arm64"),
                           ("--network", "none"), ("--cpus", "1"),
                           ("--memory", "4096m"), ("--memory-swap", "4096m"),
                           ("--pids-limit", "64"),
                           ("--security-opt", "no-new-privileges")):
        require(option(command, flag) == expected, f"stage {flag} differs")
    require("--read-only" in command
            and option(command, "--cap-drop") == "ALL"
            and IMAGE in command
            and all(any(f"dst=/{item},readonly" in argument for argument in command)
                    for item in ("worker", "small.jsonl.gz", "envelope.jsonl.gz",
                                 "small-manifest.json", "envelope-manifest.json",
                                 "fixture.json")), "stage readonly mounts differ")
    tail = ["/worker", "verify-control" if name == "verify" else "control-run",
            "/small-manifest.json", "/small.jsonl.gz",
            "/envelope-manifest.json", "/envelope.jsonl.gz", "/fixture.json"]
    if name == "control":
        tail += ["840", "5000", "32", "/out/control"]
    require(command[-len(tail):] == tail, "stage worker arguments differ")
    for suffix in ("stdout", "stderr"):
        require(sha(run / f"{name}.{suffix}") == stage.get(f"{suffix}_sha256"),
                "stage log hash differs")


def audit(run: Path, authorization: Path, small: Path, envelope: Path) -> dict:
    fixture = record(FIXTURE)
    require(sha(FIXTURE) == FIXTURE_SHA256, "known witness fixture differs")
    require(sha(small) == SMALL and sha(envelope) == ENVELOPE,
            "archive compressed hashes differ")
    config, outer, budget = (record(run / name) for name in
                             ("config.json", "outer.json", "budget.json"))
    auth = record(authorization)
    exact(auth, {"schema": "n83.additional-local-budget/v7",
                 "study": STUDY,
                 "purpose": "one_known_witness_control_reuse_k600_worker_after_cross_stream",
                 "prior_total_charged_seconds": "19800",
                 "prior_k600_audit_sha256": K600_AUDIT,
                 "prior_cross_stream_audit_sha256": CROSS_AUDIT,
                 "reuse_binary_sha256": BINARY,
                 "additional_authorized_seconds": "1200",
                 "output_dir": str(run.resolve())}, "authorization")
    require(isinstance(auth.get("user_approval_reference"), str)
            and auth["user_approval_reference"].strip(), "approval reference missing")
    auth_hash = sha(authorization)
    exact(config, {"schema": "n83.known-witness-reuse-supervisor-config/v1",
                   "binary_source_commit": "f3d27f68b68cda907e0668042ac7fa13ae6b579d",
                   "reuse_binary_sha256": BINARY,
                   "prior_k600_audit_sha256": K600_AUDIT,
                   "prior_cross_stream_audit_sha256": CROSS_AUDIT,
                   "prior_total_charged_seconds": 19800,
                   "authorization_sha256": auth_hash,
                   "container_image_id": IMAGE,
                   "verify_charge_seconds": 120,
                   "search_charge_seconds": 900,
                   "cleanup_reserve_seconds": 180,
                   "conservative_full_charge_seconds": 1200,
                   "conflicts": 5000, "models": 32,
                   "scope": "known_answer_control_only",
                   "natural_relations": None,
                   "total_index_calculus_runtime_ms": None}, "supervisor config")
    commit = config.get("source_commit")
    require(isinstance(commit, str) and len(commit) == 40
            and all(char in "0123456789abcdef" for char in commit),
            "frozen source commit malformed")
    inputs = config.get("input_sha256")
    require(isinstance(inputs, dict) and len(inputs) >= 10,
            "frozen input inventory missing")
    for path, digest in inputs.items():
        source = Path(path)
        require(isinstance(digest, str) and sha(source) == digest,
                "frozen input bytes differ")
        if source.is_relative_to(ROOT):
            relative = source.relative_to(ROOT)
            blob = subprocess.run(["git", "show", f"{commit}:{relative}"],
                                  cwd=ROOT, capture_output=True, check=True).stdout
            require(hashlib.sha256(blob).hexdigest() == digest,
                    f"frozen Git blob differs: {relative}")
    require(inputs.get(str(authorization.resolve())) == auth_hash
            and inputs.get(str(small.resolve())) == SMALL
            and inputs.get(str(envelope.resolve())) == ENVELOPE,
            "input path bindings differ")
    fixed_paths = {
        ROOT / "research" / STUDY /
        "verification/k600-mixed-orbit-first-2400-20261012/audit.json": K600_AUDIT,
        ROOT / "research" / STUDY /
        "verification/cross-stream-first-gate-source-20261013/audit.json": CROSS_AUDIT,
        ROOT / "research" / STUDY /
        "verification/k600-mixed-orbit-first-2400-20261012/evidence-01/build/build.json": BUILD,
        ROOT / "research" / STUDY /
        "verification/k600-mixed-orbit-first-2400-20261012/evidence-01/build/worker": BINARY,
        FIXTURE: FIXTURE_SHA256,
    }
    require(all(inputs.get(str(path)) == digest for path, digest in fixed_paths.items()),
            "pinned prior evidence path or hash differs")
    exact(outer, {"schema": "n83.known-witness-reuse-supervisor-outer/v1",
                  "scope": "known_answer_control_only",
                  "source_commit": commit,
                  "binary_source_commit": config["binary_source_commit"],
                  "prior_cross_stream_audit_sha256": CROSS_AUDIT,
                  "authorization_sha256": auth_hash,
                  "config_sha256": sha(run / "config.json"),
                  "verify_passed": True, "replay_error": None,
                  "frozen_input_changed": False,
                  "natural_relations": None,
                  "total_index_calculus_runtime_ms": None,
                  "selected_best_total_runtime": None}, "outer")
    status = outer.get("status")
    require(status in KNOWN_STATUSES, "outer status is not a completed worker result")
    elapsed = outer.get("elapsed_supervisor_seconds_L0")
    require(type(elapsed) in (int, float) and math.isfinite(elapsed)
            and 0 <= elapsed <= 1200, "outer wall cap differs")
    stages = outer.get("stages")
    require(isinstance(stages, list) and len(stages) == 2,
            "verify/control stage sequence missing")
    check_stage(stages[0], "verify", run)
    check_stage(stages[1], "control", run)
    mounted = {
        "worker": ROOT / "research" / STUDY /
                  "verification/k600-mixed-orbit-first-2400-20261012/evidence-01/build/worker",
        "small.jsonl.gz": small.resolve(), "envelope.jsonl.gz": envelope.resolve(),
        "small-manifest.json": ROOT / "research" / STUDY /
                               "verification/matched-k1182-k2048-900s-20261011/evidence-01/panel-k1182/manifest.json",
        "envelope-manifest.json": ROOT / "research" / STUDY /
                                  "verification/remaining-budget-k2048-20261010/manifest.json",
        "fixture.json": FIXTURE,
    }
    for stage in stages:
        command = stage["command"]
        for target, source in mounted.items():
            require(f"type=bind,src={source},dst=/{target},readonly" in command,
                    f"{stage['name']} {target} input mount differs")
        require(f"type=bind,src={run.resolve()},dst=/out" in command,
                f"{stage['name']} output mount differs")
    verify = record(run / "verify.stdout")
    exact(verify, {"schema": "n83.mixed-orbit-known-witness-verify/v1",
                   "status": "PASS", "scope": "known_answer_control_only",
                   "target": fixture["equation"]["target"],
                   "solver_search_executed": False}, "verify")
    worker_config, result = (record(run / name) for name in
                             ("control/config.json", "control/result.json"))
    classification = classify_worker(worker_config, result, fixture, small, envelope)
    require(result["status"] == status, "outer and worker status differ")
    exact(budget, {"schema": "n83.known-witness-reuse-supervisor-budget/v1",
                   "prior_total_charged_seconds": 19800,
                   "prior_k600_audit_sha256": K600_AUDIT,
                   "prior_cross_stream_audit_sha256": CROSS_AUDIT,
                   "reuse_binary_sha256": BINARY,
                   "authorization_sha256": auth_hash,
                   "conservative_charge_seconds": 1200,
                   "cumulative_conservative_charge_seconds": 21000,
                   "remaining_additional_seconds": 0,
                   "round_status": status,
                   "selected_best_total_runtime": None}, "budget")
    return {"schema": "n83.known-witness-v7-independent-audit/v1",
            "status": classification, "worker_status": status,
            "source_commit": commit,
            "conservative_charge_seconds": 1200,
            "cumulative_conservative_charge_seconds": 21000,
            "natural_relations": None,
            "total_index_calculus_runtime_ms": None,
            "selected_best_total_runtime": None}


def raw_hashes(run: Path) -> dict[str, str]:
    return {str(path.relative_to(run)): sha(path) for path in sorted(run.rglob("*"))
            if path.is_file()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--authorization", required=True, type=Path)
    parser.add_argument("--small-object", required=True, type=Path)
    parser.add_argument("--envelope-object", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        row = audit(args.run_dir, args.authorization,
                    args.small_object, args.envelope_object)
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError,
            json.JSONDecodeError) as error:
        row = {"schema": "n83.known-witness-v7-independent-audit/v1",
               "status": "UNVALIDATED", "reason": f"{type(error).__name__}: {error}",
               "selected_best_total_runtime": None}
    row["raw_sha256"] = raw_hashes(args.run_dir) if args.run_dir.is_dir() else {}
    require(not args.output.exists(), "audit output already exists")
    args.output.write_text(json.dumps(row, sort_keys=True, indent=2) + "\n")
    print(json.dumps({key: value for key, value in row.items() if key != "raw_sha256"},
                     sort_keys=True))
    return 0 if row["status"] == "VALIDATED_CONTROL_ACCESS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
