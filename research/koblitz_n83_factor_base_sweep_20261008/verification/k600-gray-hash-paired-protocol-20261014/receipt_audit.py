#!/usr/bin/env python3
"""Independent fail-closed replay of one K600 paired-arm worker receipt."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
import sys

import blake3
import preflight


STUDY = preflight.STUDY
sys.path.insert(0, str(STUDY))
PRIOR = preflight.load_module("n83_paired_prior", STUDY / "k600_mixed_orbit_supervisor.py")
CAPACITY = preflight.load_module("n83_paired_capacity_audit", preflight.CAPACITY)
POINTS_PER_ORBIT = 166
PUBLIC_TARGET = (1834913121178197242856238, 8797964190732846213981377)


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def archive_points(path: Path) -> list[tuple[int, int]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        header = json.loads(next(stream))
        require(header.get("schema") == "n83.factor-base/v1"
                and header.get("point_count") == 600 * POINTS_PER_ORBIT,
                "archive header changed")
        points = []
        for index, line in enumerate(stream):
            row = json.loads(line)
            point = PRIOR.parse_point(row.get("point"))
            require(point is not None and row.get("column") == index // POINTS_PER_ORBIT
                    and row.get("phase") == (index % POINTS_PER_ORBIT) // 2
                    and row.get("negative") is (index % 2 == 1),
                    "archive row or orbit labeling changed")
            points.append(point)
    require(len(points) == 600 * POINTS_PER_ORBIT and len(set(points)) == len(points),
            "archive point count or uniqueness changed")
    return points


def domains(base: Path, pool: Path, shared_expected: int) -> tuple[list, set]:
    base_points, pool_points = archive_points(base), archive_points(pool)
    base_set = set(base_points)
    base_x = {point[0] for point in base_points}
    shared = 0
    residual = set()
    for offset in range(0, len(pool_points), POINTS_PER_ORBIT):
        orbit = pool_points[offset:offset + POINTS_PER_ORBIT]
        overlap = sum(point in base_set for point in orbit)
        require(overlap in (0, POINTS_PER_ORBIT), "partial signed-orbit overlap")
        if overlap:
            shared += 1
        else:
            require(not any(point[0] in base_x for point in orbit),
                    "residual x domain overlaps base")
            residual.update(orbit)
    require(shared == shared_expected and len(residual) == (600 - shared) * POINTS_PER_ORBIT,
            "signed-orbit overlap changed")
    return base_points, residual


def audit_verify(record: object, arm: dict) -> None:
    require(isinstance(record, dict), "verify record missing")
    expected = {
        "schema": "n83.cross-stream-pilot-import/v1",
        "status": "PASS", "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
        "small_object_blake3": arm["base_object_blake3"],
        "small_policy": "public_x_sequential", "small_seed": 2026100802,
        "small_columns": 600, "small_points": 99600,
        "pool_object_blake3": arm["pool_object_blake3"],
        "pool_policy": "public_x_gray_prefix" if arm["role"] == "gray_pool" else "public_x_hash",
        "pool_seed": 2026100801, "pool_columns": 600, "pool_points": 99600,
        "shared_signed_orbits": arm["shared_signed_orbits"],
        "small_only_signed_orbits": 600 - arm["shared_signed_orbits"],
        "residual_signed_orbits": arm["residual_signed_orbits"],
        "residual_points": arm["residual_signed_orbits"] * POINTS_PER_ORBIT,
        "solver_search_executed": False,
    }
    for key, value in expected.items():
        require(type(record.get(key)) is type(value) and record[key] == value,
                f"verify {key} changed")


def audit_search(config: object, result: object, arm: dict,
                 base_points: list, residual: set) -> str:
    require(isinstance(config, dict) and isinstance(result, dict), "search receipts missing")
    limits = {
        "schema": "n83.cross-stream-pilot-partial-config/v1",
        "study": "koblitz_n83_factor_base_sweep_20261008",
        "base_object_blake3": arm["base_object_blake3"],
        "pool_object_blake3": arm["pool_object_blake3"],
        "pilot_manifest_blake3": blake3.blake3(
            preflight.PINNED_SHA256["manifest"][0].read_bytes()).hexdigest(),
        "source_blake3": blake3.blake3(
            preflight.PINNED_SHA256["worker_source"][0].read_bytes()).hexdigest(),
        "public_corpus_path": "/corpus.json",
        "public_corpus_canonical_json_blake3": PRIOR.CORPUS_BLAKE3,
        "budget_seconds": 1100, "conflict_budget": 5000,
        "max_models": 32, "base_summands": 3, "residual_summands": 2,
        "order_first_pair": True, "max_variables": 150000,
        "max_domain_clauses": 12000000,
        "memory_cgroup_limit_bytes": 4294967296,
        "memory_cgroup_swap_limit_bytes": 0,
    }
    for key, value in limits.items():
        require(type(config.get(key)) is type(value) and config[key] == value,
                f"search config {key} changed")
    expected_result = {
        "schema": "n83.cross-stream-pilot-partial-result/v1",
        "study": "koblitz_n83_factor_base_sweep_20261008",
        "base_object_blake3": arm["base_object_blake3"],
        "pool_object_blake3": arm["pool_object_blake3"],
        "base_policy": "public_x_sequential", "base_seed": 2026100802,
        "pool_policy": "public_x_gray_prefix" if arm["role"] == "gray_pool" else "public_x_hash",
        "pool_seed": 2026100801,
        "public_corpus_canonical_json_blake3": PRIOR.CORPUS_BLAKE3,
        "base_signed_orbits": 600, "pool_signed_orbits": 600,
        "shared_signed_orbits": arm["shared_signed_orbits"],
        "residual_signed_orbits": arm["residual_signed_orbits"],
        "trial": 0, "a": "1833078656608906180176350",
        "b": "163089659398503559150141",
        "graph_rows": None, "matrix_rank": None,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    for key, value in expected_result.items():
        require(type(result.get(key)) is type(value) and result[key] == value,
                f"search result {key} changed")
    require(PRIOR.finite_nonnegative(result.get("base_import_elapsed_ms_L0"))
            and PRIOR.finite_nonnegative(result.get("worker_wall_ms_L0"))
            and result["base_import_elapsed_ms_L0"] <= result["worker_wall_ms_L0"] <= 1100000,
            "worker runtime counters changed")
    peak = result.get("memory_cgroup_peak_bytes")
    require(type(peak) is int and 0 < peak <= 4294967296, "worker memory counter changed")
    stats = result.get("search")
    require(isinstance(stats, dict), "search counters missing")
    shape = CAPACITY.requirements(600, 600 + arm["residual_signed_orbits"], 5, 2)
    require(stats.get("domain_clauses") == shape["domain_clauses"]
            and type(stats.get("peak_variables")) is int
            and 0 <= stats["peak_variables"] <= shape["total_variables"]
            and type(stats.get("peak_xors")) is int
            and 0 <= stats["peak_xors"] <= shape["xor_rows"]
            and type(stats.get("peak_clauses")) is int
            and stats["peak_clauses"] >= 0,
            "search capacity counters changed")
    for key, cap in (("conflicts", 5000), ("models", 256),
                     ("invalid_models", 256), ("group_rejected_models", 256),
                     ("solver_calls", 8), ("patterns_attempted", 8)):
        value = stats.get(key)
        require(type(value) is int and 0 <= value <= cap, f"search {key} cap changed")
    require(type(stats.get("exhausted")) is bool and type(stats.get("refuted")) is bool,
            "search disposition counters missing")
    status = result.get("status")
    require(status in {"PASS_one_bounded_partial", "PASS_one_trial_refuted",
                       "UNKNOWN_search_cap", "UNKNOWN_preflight_cap",
                       "PRODUCER_FAILURE_invalid_model"}, "search status unclassified")
    candidate = result.get("candidate")
    require((candidate is not None) == (status == "PASS_one_bounded_partial"),
            "candidate/status mismatch")
    require((stats["invalid_models"] > 0) == (status == "PRODUCER_FAILURE_invalid_model"),
            "invalid-model disposition mismatch")
    if status == "UNKNOWN_search_cap":
        require(stats["exhausted"], "search cap lacks exhausted counter")
    elif status == "UNKNOWN_preflight_cap":
        require(stats["exhausted"] and stats["patterns_attempted"] == 0
                and stats["solver_calls"] == 0, "preflight cap has solver calls")
    elif status == "PASS_one_trial_refuted":
        require(stats["refuted"] and not stats["exhausted"],
                "refutation counter changed")
    elif status == "PRODUCER_FAILURE_invalid_model":
        require(stats["invalid_models"] > 0, "invalid-model counter missing")
    if candidate is not None:
        require(isinstance(candidate, dict) and candidate.get("source_id") == 1
                and (candidate.get("a"), candidate.get("b")) ==
                (expected_result["a"], expected_result["b"]),
                "candidate source changed")
        indices = candidate.get("base_point_indices")
        points = [PRIOR.parse_point(point) for point in candidate.get("residual_points", [])]
        require(isinstance(indices, list) and len(indices) == 3
                and all(type(index) is int and 0 <= index < len(base_points) for index in indices)
                and len(points) == 2 and all(point in residual for point in points),
                "candidate outside exact domains")
        require(PRIOR.replay_equation(candidate["a"], candidate["b"],
                                      [base_points[index] for index in indices],
                                      points, PUBLIC_TARGET),
                "candidate group equation fails independent replay")
    return status
