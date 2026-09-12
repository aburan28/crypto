#!/usr/bin/env python3
"""Compose and verify the hosted n=53 Stages 46 through 48."""

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
import run_koblitz_stage47_n53_width as stage47
import run_koblitz_stage48_n53_width as stage48


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
DEFAULT_OUTPUT = STAGE / "stage-49-n53-width-results-20260912"
SCHEMA = "koblitz_stage49_width_results.v1"
SEAL_SCHEMA = "koblitz_stage49_width_results_seal.v1"

CONFIG = {
    46: {
        "run_id": 34709425290,
        "commit": "6dcad390a3022d71731d4b04d6e3e16a826d00a9",
        "artifact_id": 10302812323,
        "artifact_digest": "sha256:dafb19ce53c47ca73e4c571afbbf7052de5c2fdb814ec64515790e213039c915",
        "archive": "koblitz-stage46-n53-scratch-evidence-34709425290.tar.gz",
        "archive_root": "koblitz-stage46-n53-scratch-34709425290",
        "prefix": "stage44",
        "job": "n53-parallel-query",
        "runner": stage44,
    },
    47: {
        "run_id": 34709970211,
        "commit": "03ee86d87bae3f3ded460dab552220fea1bd0ab2",
        "artifact_id": 10303530259,
        "artifact_digest": "sha256:889b11bc7237865400c06df81b2ad63ada70df4b9585e00eb904d3977f88d621",
        "archive": "koblitz-stage47-n53-width-evidence-34709970211.tar.gz",
        "archive_root": "koblitz-stage47-n53-width-34709970211",
        "prefix": "stage47",
        "job": "n53-parallel-width",
        "runner": stage47,
    },
    48: {
        "run_id": 34710311549,
        "commit": "e84de09979c7656d4045952f75ce832bac257d2e",
        "artifact_id": 10302004328,
        "artifact_digest": "sha256:c66f7122f31c00d0a8a0b1bb8bef82034e560bf345bf691e53a49a548c63c0e5",
        "archive": "koblitz-stage48-n53-width-evidence-34710311549.tar.gz",
        "archive_root": "koblitz-stage48-n53-width-34710311549",
        "prefix": "stage48",
        "job": "n53-parallel-width",
        "runner": stage48,
    },
}


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


def archive_identity(stage: int, config: dict[str, Any]) -> dict[str, Any]:
    archive = STAGE / config["archive"]
    sidecar = archive.with_name(archive.name + ".sha256")
    fields = sidecar.read_text().split()
    relative = str(archive.relative_to(REPO))
    require(len(fields) == 2 and fields[1] == relative, f"Stage-{stage} archive sidecar changed")
    digest = custody.sha256_file(archive, f"Stage-{stage} archive")
    require(fields[0] == digest, f"Stage-{stage} archive hash changed")
    return {"path": relative, "bytes": archive.stat().st_size, "sha256": digest}


def extract(stage: int, config: dict[str, Any], destination: Path) -> Path:
    archive = STAGE / config["archive"]
    with tarfile.open(archive, "r:gz") as source:
        members = source.getmembers()
        require(bool(members), f"Stage-{stage} archive is empty")
        roots: set[str] = set()
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, f"unsafe archive path: {member.name}")
            require(member.isfile() or member.isdir(), f"unsupported archive member: {member.name}")
            if path.parts:
                roots.add(path.parts[0])
        require(roots == {config["archive_root"]}, f"Stage-{stage} archive root changed")
        source.extractall(destination, filter="data")
    return destination / config["archive_root"]


def workflow_accounting(stage: int, config: dict[str, Any]) -> dict[str, Any]:
    run_id = config["run_id"]
    value = load(STAGE / f"stage-{stage}-workflow-{run_id}.json", f"Stage-{stage} workflow")
    require(value.get("databaseId") == run_id and value.get("headSha") == config["commit"], f"Stage-{stage} workflow identity changed")
    require(value.get("status") == "completed" and value.get("conclusion") == "success", f"Stage-{stage} workflow failed")
    jobs = value["jobs"]
    require({row["name"] for row in jobs} == {"validate-control", config["job"]}, f"Stage-{stage} workflow jobs changed")
    require(all(row["status"] == "completed" and row["conclusion"] == "success" for row in jobs), f"Stage-{stage} job failed")
    production = next(row for row in jobs if row["name"] == config["job"])
    return {
        "run_id": run_id,
        "head_sha": config["commit"],
        "url": value["url"],
        "workflow_wall_seconds": (timestamp(value["updatedAt"]) - timestamp(value["createdAt"])).total_seconds(),
        "production_job_id": production["databaseId"],
        "production_job_url": production["url"],
        "production_job_wall_seconds": (timestamp(production["completedAt"]) - timestamp(production["startedAt"])).total_seconds(),
    }


def artifact_accounting(stage: int, config: dict[str, Any]) -> dict[str, Any]:
    run_id = config["run_id"]
    value = load(STAGE / f"stage-{stage}-artifacts-{run_id}.json", f"Stage-{stage} artifacts")
    require(value.get("total_count") == 1 and len(value.get("artifacts", [])) == 1, f"Stage-{stage} artifact count changed")
    row = value["artifacts"][0]
    require(row.get("id") == config["artifact_id"] and row.get("digest") == config["artifact_digest"], f"Stage-{stage} artifact identity changed")
    require(row.get("expired") is False and row.get("workflow_run", {}).get("id") == run_id, f"Stage-{stage} artifact source changed")
    require(row["workflow_run"].get("head_sha") == config["commit"], f"Stage-{stage} artifact commit changed")
    return {key: row[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}


def compact_summary(summary: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "query_mode", "query_parallel_threads", "query_parallel_waves", "query_parallel_chunks",
        "factor_base_points", "orbit_columns", "matrix_columns", "admitted_relations", "target_trials",
        "full_rank_at_relation", "setup_ms", "collection_ms", "charged_total_ms", "support_queries",
        "query_batch_inversions", "query_x_filter_rejections", "query_exact_table_misses", "published_q",
        "published_fixture_scalar", "recovered_fixture_scalar", "linear_solution_verified",
        "all_relations_group_verified",
    )
    return {key: summary.get(key) for key in keys}


def one_result(stage: int, config: dict[str, Any]) -> dict[str, Any]:
    with tempfile.TemporaryDirectory() as directory:
        root = extract(stage, config, Path(directory))
        prefix = config["prefix"]
        verification = config["runner"].verify(root / f"{prefix}-build", root / f"{prefix}-run")
        embedded = load(root / f"{prefix}-verification.json", f"embedded Stage-{stage} verification")
        require(verification == embedded, f"Stage-{stage} embedded verification changed")
        run = load(root / f"{prefix}-run/result.json", f"Stage-{stage} run")
        build_seal = custody.sha256_file(root / f"{prefix}-build/result-seal.json", f"Stage-{stage} build seal")
        run_seal = custody.sha256_file(root / f"{prefix}-run/result-seal.json", f"Stage-{stage} run seal")
    if stage == 46:
        baseline_summary = run["sequential_summary"]
        candidate_summary = run["parallel_summary"]
        baseline_process = run["sequential_process"]
        candidate_process = run["parallel_process"]
    else:
        baseline_summary = run["baseline_summary"]
        candidate_summary = run["candidate_summary"]
        baseline_process = run["baseline_process"]
        candidate_process = run["candidate_process"]
    return {
        "source_commit": config["commit"],
        "archive": archive_identity(stage, config),
        "workflow": workflow_accounting(stage, config),
        "artifact": artifact_accounting(stage, config),
        "build_seal_sha256": build_seal,
        "run_seal_sha256": run_seal,
        "verification": verification,
        "comparison": run["comparison"],
        "baseline": compact_summary(baseline_summary),
        "candidate": compact_summary(candidate_summary),
        "processes": {
            "baseline": baseline_process["metrics"],
            "candidate": candidate_process["metrics"],
            "rho": run["rho_process"]["metrics"],
        },
        "resources": {
            "build_outer": run["build_binding"]["outer_resources"],
            "science_outer": run["outer_resources"],
            "charged_total_core_seconds_available": run["charged_total_core_seconds_available"],
            "charged_sequential_wall_seconds_available": run["charged_sequential_wall_seconds_available"],
            "maximum_sampled_process_tree_rss_bytes": run["maximum_sampled_process_tree_rss_bytes"],
        },
    }


def result() -> dict[str, Any]:
    stages = {str(stage): one_result(stage, config) for stage, config in CONFIG.items()}
    stage47_result = stages["47"]
    stage48_result = stages["48"]
    require(stage47_result["comparison"]["candidate_direct_wall_speedup"] > 1, "Stage-47 candidate did not improve wall")
    require(stage47_result["comparison"]["candidate_direct_core_seconds_ratio"] < 1, "Stage-47 candidate did not improve CPU")
    require(stage48_result["comparison"]["candidate_direct_wall_speedup"] > 1, "Stage-48 candidate did not improve wall")
    require(stage48_result["comparison"]["candidate_direct_core_seconds_ratio"] < 1, "Stage-48 candidate did not improve CPU")
    latest = stage48_result["verification"]
    return {
        "schema": SCHEMA,
        "status": "hosted_n53_scratch_and_width_results_verified",
        "stages": stages,
        "preferred_query_mode": "pair_pair_parallel_4096",
        "latest_same_target_direct_over_rho_wall_ratio": latest["candidate_over_rho_wall_ratio"],
        "latest_fresh_build_plus_direct_over_rho_wall_ratio": latest["fresh_build_plus_candidate_over_rho_wall_ratio"],
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
            f"preferred n=53 direct remains {latest['candidate_over_rho_wall_ratio']:.6f} times same-target rho wall",
            f"fresh build plus preferred direct remains {latest['fresh_build_plus_candidate_over_rho_wall_ratio']:.6f} times rho wall",
            "unknown-scalar n=53, licensed same-instance Magma, and unaffiliated review remain absent",
            "scratch reuse and batch width change finite constants, not the asymptotic exponent",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-49 output must be new")
    output.mkdir()
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "hosted_results_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-49 result")
    require(committed == result(), "committed Stage-49 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-49 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-49 result seal changed")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-49 result inventory changed")
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
    except (OSError, ValueError, KeyError, VerificationError, stage44.Stage44Error, stage47.Stage47Error, stage48.Stage48Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage49-results: {error}")


if __name__ == "__main__":
    main()
