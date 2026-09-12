#!/usr/bin/env python3
"""Compose and verify the hosted Stage 44 n=53 parallel-query result."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage44_n53_parallel as stage44


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
ARCHIVE = STAGE / "koblitz-stage44-n53-parallel-evidence-34707644672.tar.gz"
ARCHIVE_SHA = STAGE / "koblitz-stage44-n53-parallel-evidence-34707644672.tar.gz.sha256"
WORKFLOW = STAGE / "stage-44-workflow-34707644672.json"
ARTIFACTS = STAGE / "stage-44-artifacts-34707644672.json"
DEFAULT_OUTPUT = STAGE / "stage-44-n53-parallel-result-20260912"
ARCHIVE_ROOT = "koblitz-stage44-n53-parallel-34707644672"
RUN_ID = 34707644672
RUN_COMMIT = "14ceb064a93b1b3c4c6cf3f4318f656311c98055"
ARTIFACT_ID = 10301749682
ARTIFACT_DIGEST = "sha256:33a19620b5e4bb632b946f850c77f7814868c9dec6a49e1840bd0e30e79100fa"
SCHEMA = "koblitz_stage44_hosted_result.v1"
SEAL_SCHEMA = "koblitz_stage44_hosted_result_seal.v1"


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def timestamp(value: str) -> datetime:
    require(isinstance(value, str) and value.endswith("Z"), "workflow timestamp is invalid")
    return datetime.fromisoformat(value[:-1] + "+00:00")


def archive_identity() -> dict[str, Any]:
    fields = ARCHIVE_SHA.read_text().split()
    relative = str(ARCHIVE.relative_to(REPO))
    require(len(fields) == 2 and fields[1] == relative, "Stage-44 archive sidecar changed")
    digest = custody.sha256_file(ARCHIVE, "Stage-44 archive")
    require(fields[0] == digest, "Stage-44 archive hash changed")
    return {"path": relative, "bytes": ARCHIVE.stat().st_size, "sha256": digest}


def extract(destination: Path) -> Path:
    with tarfile.open(ARCHIVE, "r:gz") as source:
        members = source.getmembers()
        require(bool(members), "Stage-44 archive is empty")
        roots: set[str] = set()
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, f"unsafe archive path: {member.name}")
            require(member.isfile() or member.isdir(), f"unsupported archive member: {member.name}")
            if path.parts:
                roots.add(path.parts[0])
        require(roots == {ARCHIVE_ROOT}, "Stage-44 archive root changed")
        source.extractall(destination, filter="data")
    return destination / ARCHIVE_ROOT


def workflow_accounting() -> dict[str, Any]:
    value = load(WORKFLOW, "Stage-44 workflow metadata")
    require(value.get("databaseId") == RUN_ID and value.get("headSha") == RUN_COMMIT, "Stage-44 workflow identity changed")
    require(value.get("status") == "completed" and value.get("conclusion") == "success", "Stage-44 workflow failed")
    jobs = value["jobs"]
    require({row["name"] for row in jobs} == {"validate-control", "n53-parallel-query"} and all(row["conclusion"] == "success" for row in jobs), "Stage-44 workflow jobs changed")
    production = next(row for row in jobs if row["name"] == "n53-parallel-query")
    return {
        "run_id": RUN_ID,
        "head_sha": RUN_COMMIT,
        "url": value["url"],
        "workflow_wall_seconds": (timestamp(value["updatedAt"]) - timestamp(value["createdAt"])).total_seconds(),
        "production_job_id": production["databaseId"],
        "production_job_url": production["url"],
        "production_job_wall_seconds": (timestamp(production["completedAt"]) - timestamp(production["startedAt"])).total_seconds(),
    }


def artifact_accounting() -> dict[str, Any]:
    value = load(ARTIFACTS, "Stage-44 artifact metadata")
    rows = value["artifacts"]
    require(value.get("total_count") == 1 and len(rows) == 1, "Stage-44 artifact count changed")
    row = rows[0]
    require(row.get("id") == ARTIFACT_ID and row.get("digest") == ARTIFACT_DIGEST and row.get("expired") is False, "Stage-44 artifact identity changed")
    require(row.get("workflow_run", {}).get("id") == RUN_ID and row["workflow_run"].get("head_sha") == RUN_COMMIT, "Stage-44 artifact source changed")
    return {key: row[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}


def compact_summary(summary: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "query_mode", "query_parallel_threads", "query_parallel_waves",
        "query_parallel_chunks", "factor_base_points", "orbit_columns",
        "matrix_columns", "admitted_relations", "target_trials",
        "full_rank_at_relation", "setup_ms", "collection_ms",
        "charged_total_ms", "support_queries", "query_batch_inversions",
        "query_x_filter_rejections", "query_exact_table_misses",
        "published_q", "published_fixture_scalar", "recovered_fixture_scalar",
        "linear_solution_verified", "all_relations_group_verified",
    )
    return {key: summary.get(key) for key in keys}


def result() -> dict[str, Any]:
    archive = archive_identity()
    with tempfile.TemporaryDirectory() as directory:
        root = extract(Path(directory))
        verification = stage44.verify(root / "stage44-build", root / "stage44-run")
        embedded = load(root / "stage44-verification.json", "embedded Stage-44 verification")
        require(embedded == verification, "Stage-44 embedded verification changed")
        run = load(root / "stage44-run/result.json", "Stage-44 run result")
        build_seal = custody.sha256_file(root / "stage44-build/result-seal.json", "Stage-44 build seal")
        run_seal = custody.sha256_file(root / "stage44-run/result-seal.json", "Stage-44 run seal")
    comparison = run["comparison"]
    require(comparison["parallel_direct_wall_speedup"] > 1 and comparison["parallel_whole_process_crossover"] is False, "Stage-44 result boundary changed")
    return {
        "schema": SCHEMA,
        "status": "hosted_n53_parallel_query_verified",
        "source_commit": RUN_COMMIT,
        "archive": archive,
        "workflow": workflow_accounting(),
        "artifact": artifact_accounting(),
        "build_seal_sha256": build_seal,
        "run_seal_sha256": run_seal,
        "verification": verification,
        "comparison": comparison,
        "sequential": compact_summary(run["sequential_summary"]),
        "parallel": compact_summary(run["parallel_summary"]),
        "processes": {
            "sequential_direct": run["sequential_process"]["metrics"],
            "parallel_direct": run["parallel_process"]["metrics"],
            "rho": run["rho_process"]["metrics"],
        },
        "resources": {
            "build_outer": run["build_binding"]["outer_resources"],
            "science_outer": run["outer_resources"],
            "charged_total_core_seconds_available": run["charged_total_core_seconds_available"],
            "charged_sequential_wall_seconds_available": run["charged_sequential_wall_seconds_available"],
            "maximum_sampled_process_tree_rss_bytes": run["maximum_sampled_process_tree_rss_bytes"],
        },
        "same_public_target": True,
        "exact_relation_hash_equality": True,
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "fixture_scalar_used_by_collector": False,
        "sat_conflicts": None,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "remaining_blockers": [
            "parallel n=53 direct remains 5.974650 times same-target rho wall",
            "fresh build plus parallel direct remains 25.281709 times rho wall",
            "parallel direct spends 1.482930 times sequential process CPU",
            "unknown-scalar n=53, licensed Magma, and unaffiliated review remain absent",
            "parallel query changes a finite wall constant, not the asymptotic exponent",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-44 result output must be new")
    output.mkdir()
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "hosted_result_frozen", "workflow_run_id": RUN_ID, "workflow_commit": RUN_COMMIT, "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-44 result")
    require(committed == result(), "committed Stage-44 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-44 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-44 result seal changed")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-44 result inventory changed")
    return committed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    compose = sub.add_parser("compose")
    compose.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    check = sub.add_parser("verify")
    check.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    try:
        value = freeze(args.output.resolve()) if args.command == "compose" else verify(args.output.resolve())
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, stage44.Stage44Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage44-result: {error}")


if __name__ == "__main__":
    main()
