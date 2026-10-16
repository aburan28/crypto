#!/usr/bin/env python3
"""Independent raw-receipt audit of the first retained N83 cross-stream gate."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess


BASE = "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e"
POOL = "e5600aa76a88246a4c3f73fce83ba73479e79c34e88d6205fffdb91ce09e0a99"
COMMIT = "0764901992591c39566d720ae812c54cdc09ab75"
BASE_SHA256 = "3e85c9e419fcddd3c9faed3c84470bf96fcad9178e48e1d73629369335841ecf"
POOL_SHA256 = "74841a1e290e84904f898a9c952e91a795e34498f48e4529ca1c621a503b62cc"
AUDIT_PRIOR = "387d393b9805f35ee1ca6328e8602133a1b9bc036088b1be1a3ba5a29bb8469c"


def must(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def document(path: Path) -> dict:
    value = json.loads(path.read_text())
    must(isinstance(value, dict), f"{path.name} is not a JSON object")
    return value


def exact(row: dict, expected: dict, label: str) -> None:
    for key, value in expected.items():
        must(type(row.get(key)) is type(value) and row[key] == value,
             f"{label}.{key} differs")


def archive_set(path: Path, policy: str) -> set[tuple[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        header = json.loads(next(stream))
        exact(header, {"schema": "n83.factor-base/v1",
                       "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
                       "policy": policy, "seed": 2026100801,
                       "closure": "signed_frobenius", "orbit_columns": 600,
                       "point_count": 99600}, "archive header")
        points: set[tuple[str, str]] = set()
        count = 0
        for count, line in enumerate(stream, 1):
            row = json.loads(line)
            point = row.get("point")
            must(isinstance(point, list) and len(point) == 2
                 and all(isinstance(part, str) and part.isascii() and part.isdecimal()
                         for part in point), "archive point malformed")
            exact(row, {"column": (count - 1) // 166,
                        "phase": ((count - 1) % 166) // 2,
                        "negative": (count - 1) % 2 == 1}, "archive row")
            points.add((point[0], point[1]))
        must(count == 99600 and len(points) == count, "archive has missing or repeated points")
        return points


def audit(run: Path, base: Path, pool: Path) -> dict:
    must(sha256(base) == BASE_SHA256 and sha256(pool) == POOL_SHA256,
         "downloaded compressed archive hash differs")
    must(base.stat().st_size == 3025798 and pool.stat().st_size == 3020387,
         "compressed archive length differs")
    base_points = archive_set(base, "public_x_hash")
    pool_points = archive_set(pool, "public_x_sequential")
    must(base_points.isdisjoint(pool_points), "point domains overlap")
    must({x for x, _ in base_points}.isdisjoint({x for x, _ in pool_points}),
         "x domains overlap")

    outer, budget = document(run / "outer.json"), document(run / "budget.json")
    config = document(run / "config.json")
    verify = document(run / "run/verify.stdout")
    worker_config = document(run / "run/search/config.json")
    result = document(run / "run/search/result.json")
    exact(outer, {"schema": "n83.cross-stream-first-gate-outer/v1",
                  "status": "UNKNOWN_search_cap", "source_commit": COMMIT,
                  "frozen_input_changed": False, "replay_error": None,
                  "total_index_calculus_runtime_ms": None,
                  "selected_best_total_runtime": None,
                  "overlap": {"shared_signed_orbits": 0,
                              "residual_signed_orbits": 600}}, "outer")
    exact(config, {"schema": "n83.cross-stream-first-gate-config/v1",
                   "source_commit": COMMIT,
                   "base_object_blake3": BASE, "pool_object_blake3": POOL,
                   "prior_total_charged_seconds": 17400,
                   "conservative_full_charge_seconds": 2400,
                   "selected_best_total_runtime": None}, "supervisor config")
    must(outer["config_sha256"] == sha256(run / "config.json"),
         "supervisor config digest differs")
    must(outer["authorization_sha256"] == config["authorization_sha256"]
         == budget["authorization_sha256"], "authorization binding differs")
    must(outer["binary_sha256"] == sha256(run / "build/worker"),
         "Linux binary digest differs")
    source_hashes = outer.get("source_sha256")
    must(isinstance(source_hashes, dict) and source_hashes == config.get("source_sha256")
         and len(source_hashes) >= 10, "frozen source hash map differs")
    checkout = Path(__file__).resolve().parents[3]
    for relative, digest in source_hashes.items():
        must(isinstance(relative, str) and isinstance(digest, str),
             "frozen source hash entry malformed")
        blob = subprocess.run(["git", "show", f"{COMMIT}:{relative}"],
                              cwd=checkout, capture_output=True, check=True).stdout
        must(hashlib.sha256(blob).hexdigest() == digest,
             f"frozen source differs at {relative}")
    inputs = outer.get("input_sha256")
    must(isinstance(inputs, dict) and inputs == config.get("input_sha256")
         and len(inputs) == 5, "frozen input hash map differs")
    for path, digest in inputs.items():
        must(sha256(Path(path)) == digest, "frozen input bytes differ")
    must(inputs.get(str(base.resolve())) == BASE_SHA256
         and inputs.get(str(pool.resolve())) == POOL_SHA256,
         "archive input binding differs")
    authorizations = [Path(path) for path, digest in inputs.items()
                      if digest == outer["authorization_sha256"]]
    must(len(authorizations) == 1, "authorization input missing")
    authorization = document(authorizations[0])
    exact(authorization, {
        "schema": "n83.additional-local-budget/v6",
        "purpose": "one_k600_hash_vs_sequential_cross_stream_partial",
        "prior_total_charged_seconds": "17400",
        "prior_k600_audit_sha256": AUDIT_PRIOR,
        "base_object_blake3": BASE, "pool_object_blake3": POOL,
        "additional_authorized_seconds": "2400",
        "output_dir": str(run.resolve()),
    }, "authorization")
    must(isinstance(authorization.get("user_approval_reference"), str)
         and authorization["user_approval_reference"].strip(),
         "authorization approval reference missing")
    exact(budget, {"schema": "n83.cross-stream-first-gate-budget/v1",
                   "prior_total_charged_seconds": 17400,
                   "prior_k600_audit_sha256": AUDIT_PRIOR,
                   "conservative_charge_seconds": 2400,
                   "cumulative_conservative_charge_seconds": 19800,
                   "remaining_additional_seconds": 0,
                   "round_status": "UNKNOWN_search_cap",
                   "selected_best_total_runtime": None}, "budget")
    must(type(outer.get("elapsed_supervisor_seconds_L0")) in (int, float)
         and math.isfinite(outer["elapsed_supervisor_seconds_L0"])
         and 0 < outer["elapsed_supervisor_seconds_L0"] <= 2400,
         "supervisor wall time differs")
    exact(verify, {"schema": "n83.cross-stream-pilot-import/v1", "status": "PASS",
                   "small_object_blake3": BASE, "pool_object_blake3": POOL,
                   "small_policy": "public_x_hash", "pool_policy": "public_x_sequential",
                   "small_seed": 2026100801, "pool_seed": 2026100801,
                   "small_columns": 600, "pool_columns": 600,
                   "small_points": 99600, "pool_points": 99600,
                   "shared_signed_orbits": 0, "small_only_signed_orbits": 600,
                   "residual_signed_orbits": 600, "residual_points": 99600,
                   "solver_search_executed": False}, "verify")
    exact(worker_config, {"schema": "n83.cross-stream-pilot-partial-config/v1",
                          "base_object_blake3": BASE, "pool_object_blake3": POOL,
                          "budget_seconds": 1100, "conflict_budget": 5000,
                          "max_models": 32, "base_summands": 3,
                          "residual_summands": 2, "order_first_pair": True,
                          "max_variables": 150000,
                          "max_domain_clauses": 12000000,
                          "memory_cgroup_limit_bytes": 4294967296,
                          "memory_cgroup_swap_limit_bytes": 0}, "worker config")
    exact(result, {"schema": "n83.cross-stream-pilot-partial-result/v1",
                   "status": "UNKNOWN_search_cap", "base_object_blake3": BASE,
                   "pool_object_blake3": POOL,
                   "base_policy": "public_x_hash", "pool_policy": "public_x_sequential",
                   "base_signed_orbits": 600, "pool_signed_orbits": 600,
                   "shared_signed_orbits": 0, "residual_signed_orbits": 600,
                   "trial": 0, "candidate": None, "graph_rows": None,
                   "matrix_rank": None, "total_index_calculus_runtime_ms": None,
                   "selected_best_total_runtime": None}, "search result")
    exact(result["search"], {"conflicts": 5000, "models": 0,
                             "invalid_models": 0, "group_rejected_models": 0,
                             "solver_calls": 1, "patterns_attempted": 1,
                             "patterns_refuted": 0,
                             "exhausted": True, "refuted": False,
                             "domain_clauses": 269830, "peak_variables": 97080,
                             "peak_xors": 4980, "peak_clauses": 518164},
          "search counters")
    must(0 < result["memory_cgroup_peak_bytes"] <= 4294967296,
         "memory peak invalid")
    must(0 < result["base_import_elapsed_ms_L0"]
         < result["worker_wall_ms_L0"] <= 1100000,
         "worker wall or import timer invalid")
    stages = outer.get("stages")
    must(isinstance(stages, list) and [row.get("name") for row in stages]
         == ["build", "verify", "search"], "stage order differs")
    for stage, cap in zip(stages, (450, 300, 1150)):
        state = stage.get("container_state")
        must(isinstance(state, dict) and state.get("ExitCode") == 0
             and state.get("OOMKilled") is False
             and stage.get("exit_code") == 0
             and stage.get("timed_out") is False
             and stage.get("launch_error") is None
             and stage.get("cleanup_error") is None
             and stage.get("cleanup_verified") is True
             and 0 < stage.get("wall_elapsed_seconds_L0", -1) <= cap,
             f"{stage.get('name')} stage failed or exceeded wall cap")
        command = stage.get("command")
        must(isinstance(command, list) and command[:2] == ["docker", "run"]
             and "--network" in command and command[command.index("--network") + 1] == "none",
             "stage network isolation differs")
        folder = run / ("build" if stage["name"] == "build" else "run")
        must(stage["stdout_sha256"] == sha256(folder / f"{stage['name']}.stdout")
             and stage["stderr_sha256"] == sha256(folder / f"{stage['name']}.stderr"),
             "stage log digest differs")
    files = ("outer.json", "budget.json", "config.json", "build/worker",
             "build/build.stdout", "build/build.stderr", "run/verify.stdout",
             "run/verify.stderr", "run/search.stdout", "run/search.stderr",
             "run/search/config.json", "run/search/result.json")
    return {
        "schema": "n83.cross-stream-first-gate-audit/v1",
        "status": "PASS_receipts_with_censored_search",
        "source_commit": COMMIT,
        "base_object_blake3": BASE, "pool_object_blake3": POOL,
        "base_object_sha256": BASE_SHA256, "pool_object_sha256": POOL_SHA256,
        "shared_signed_orbits": 0, "residual_signed_orbits": 600,
        "search_status": "UNKNOWN_search_cap", "search_conflicts": 5000,
        "search_models": 0, "search_candidate": None,
        "conservative_charge_seconds": 2400,
        "cumulative_conservative_charge_seconds": 19800,
        "elapsed_supervisor_seconds_L0": outer["elapsed_supervisor_seconds_L0"],
        "selected_best_total_runtime": None,
        "raw_sha256": {name: sha256(run / name) for name in files},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--pool", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    receipt = audit(args.run, args.base, args.pool)
    if args.output.exists():
        must(document(args.output) == receipt, "existing audit differs")
    else:
        args.output.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({key: value for key, value in receipt.items() if key != "raw_sha256"},
                     sort_keys=True))


if __name__ == "__main__":
    main()
