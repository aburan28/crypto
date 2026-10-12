#!/usr/bin/env python3
"""Compose and verify the Stage 41 tuning and Stage 42 hosted n=53 result."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage42_n53_same_target as stage42


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
TUNING_ARCHIVE = STAGE / "koblitz-stage41-n53-tuning-evidence-20260912.tar.gz"
TUNING_SHA = STAGE / "koblitz-stage41-n53-tuning-evidence-20260912.tar.gz.sha256"
HOSTED_ARCHIVE = STAGE / "koblitz-stage42-n53-same-target-evidence-34705118094.tar.gz"
HOSTED_SHA = STAGE / "koblitz-stage42-n53-same-target-evidence-34705118094.tar.gz.sha256"
WORKFLOW = STAGE / "stage-42-workflow-34705118094.json"
ARTIFACTS = STAGE / "stage-42-artifacts-34705118094.json"
DEFAULT_OUTPUT = STAGE / "stage-42-n53-same-target-result-20260912"

TUNING_ROOT = "koblitz-stage41-n53-tuning-20260912"
HOSTED_ROOT = "koblitz-stage42-n53-same-target-34705118094"
RUN_ID = 34705118094
RUN_COMMIT = "0654c7dd16b8dfd7a61b8770583ac39316b5ff6b"
ARTIFACT_ID = 10300903216
ARTIFACT_DIGEST = "sha256:6b753a15ce81f1ffebb0298d3044e88b8eecb40c55c91091098012d2ffa42467"
SCHEMA = "koblitz_stage42_hosted_result.v1"
SEAL_SCHEMA = "koblitz_stage42_hosted_result_seal.v1"


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


def archive_identity(path: Path, sidecar: Path, context: str) -> dict[str, Any]:
    fields = sidecar.read_text().split()
    relative = str(path.relative_to(REPO))
    require(len(fields) == 2 and fields[1] == relative, f"{context} sidecar is invalid")
    digest = custody.sha256_file(path, context)
    require(fields[0] == digest, f"{context} hash changed")
    return {"path": relative, "bytes": path.stat().st_size, "sha256": digest}


def extract(archive: Path, expected_root: str, destination: Path) -> Path:
    with tarfile.open(archive, "r:gz") as source:
        members = source.getmembers()
        require(bool(members), f"{archive.name} is empty")
        roots: set[str] = set()
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, f"unsafe archive path: {member.name}")
            require(member.isfile() or member.isdir(), f"unsupported archive member: {member.name}")
            if path.parts:
                roots.add(path.parts[0])
        require(roots == {expected_root}, f"{archive.name} root changed")
        source.extractall(destination, filter="data")
    return destination / expected_root


def json_lines(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def candidate(root: Path, directory: str, label: str) -> dict[str, Any]:
    path = root / directory
    metrics = load(path / "direct.metrics.json", f"{label} metrics")
    rows = json_lines(path / "direct.stdout.jsonl")
    base = next(row for row in rows if row.get("kind") == "point_defined_factor_base")
    summary = next(row for row in rows if row.get("kind") == "relation_rank_summary")
    require(metrics.get("returncode") == 0 and metrics.get("timed_out") is False and metrics.get("orphan_group_terminated") is False, f"{label} process failed")
    require(summary.get("status") == "RANK_PLUS_32" and summary.get("linear_solution_verified") is True, f"{label} rank endpoint failed")
    require(summary.get("all_relations_group_verified") is True and base.get("selection_uses_scalar_labels") is False, f"{label} evidence boundary failed")
    return {
        "label": label,
        "eta": summary["eta"],
        "pair_mode": summary["pair_index_mode"],
        "query_mode": summary["query_mode"],
        "decomposition_arity": summary["decomposition_arity"],
        "factor_base_points": summary["factor_base_points"],
        "orbit_columns": summary["orbit_columns"],
        "relations": summary["admitted_relations"],
        "trials": summary["target_trials"],
        "setup_seconds": summary["setup_ms"] / 1000.0,
        "collection_seconds": summary["collection_ms"] / 1000.0,
        "charged_seconds": summary["charged_total_ms"] / 1000.0,
        "whole_process_wall_seconds": metrics["metrics"]["wall_seconds"],
        "total_core_seconds": metrics["metrics"]["total_core_seconds"],
        "peak_rss_bytes": metrics["metrics"]["peak_rss_bytes"],
        "support_table_allocated_bytes": base["support_table_allocated_bytes"],
        "support_index_entries": base["support_index_entries"],
        "support_queries": summary["support_queries"],
        "query_canonicalization_maps": summary["query_canonicalization_maps"],
        "query_x_filter_rejections": summary["query_x_filter_rejections"],
        "published_fixture_scalar": summary["published_fixture_scalar"],
        "recovered_fixture_scalar": summary["recovered_fixture_scalar"],
        "base_hash": summary["base_hash"],
    }


def pointwise(root: Path) -> dict[str, Any]:
    run = root / "autolab-pointwise"
    state = load(run / "state.json", "Stage-41 autolab state")
    claim = load(run / "artifacts/claim_draft.json", "Stage-41 autolab claim")
    check = load(run / "artifacts/claim_check.json", "Stage-41 autolab claim check")
    manifest = load(run / "artifacts/review_manifest.json", "Stage-41 autolab manifest")
    mismatches = [name for name, digest in manifest["files"].items() if not (run / name).is_file() or custody.sha256_file(run / name, name) != digest]
    require(not mismatches, f"Stage-41 autolab manifest changed: {mismatches}")
    require(state.get("status") == "PENDING_INDEPENDENT_VALIDATION" and check.get("status") == "PASS", "Stage-41 autolab status changed")
    rows = json_lines(run / "logs/direct.stdout.jsonl")
    base = next(row for row in rows if row.get("kind") == "point_defined_factor_base")
    summary = next(row for row in rows if row.get("kind") == "relation_rank_summary")
    receipt = load(run / "receipts/direct.resource.json", "Stage-41 direct receipt")
    return {
        "label": "pointwise-autolab",
        "status": state["status"],
        "claim_check": check["status"],
        "eta": summary["eta"],
        "pair_mode": summary["pair_index_mode"],
        "query_mode": summary["query_mode"],
        "factor_base_points": summary["factor_base_points"],
        "orbit_columns": summary["orbit_columns"],
        "relations": summary["admitted_relations"],
        "trials": summary["target_trials"],
        "setup_seconds": summary["setup_ms"] / 1000.0,
        "collection_seconds": summary["collection_ms"] / 1000.0,
        "charged_seconds": summary["charged_total_ms"] / 1000.0,
        "whole_process_wall_seconds": receipt["whole_process_wall_ms"] / 1000.0,
        "total_core_seconds": (receipt["children_cpu_user_ms"] + receipt["children_cpu_system_ms"]) / 1000.0,
        "peak_rss_bytes": None,
        "resource_cap_enforced": False,
        "support_table_allocated_bytes": base["support_table_allocated_bytes"],
        "query_canonicalization_maps": summary["query_canonicalization_maps"],
        "rho_wall_seconds": claim["whole_process_wall_ms"]["rho"] / 1000.0,
        "direct_over_rho_wall_ratio": claim["ic_cost"] / claim["rho_cost"],
    }


def tuning_result(root: Path) -> dict[str, Any]:
    rows = [
        candidate(root, "signed-quotient-fiber64", "signed-quotient fiber64"),
        candidate(root, "signed-quotient-batch-inverse", "signed-quotient batch inverse"),
        candidate(root, "signed-expanded-fiber64", "signed-expanded fiber64"),
        candidate(root, "eta16-expanded-pairpair256", "eta16 pair-pair256"),
        candidate(root, "eta96-expanded-pairpair256", "eta96 pair-pair256"),
        candidate(root, "eta128-expanded-pairpair256", "eta128 pair-pair256"),
        candidate(root, "eta192-expanded-pairpair256", "eta192 pair-pair256"),
        candidate(root, "eta1024-expanded-pairpair256", "eta1024 pair-pair256"),
        candidate(root, "eta128-pairpair1024-rejected", "eta128 pair-pair1024"),
    ]
    eta_rows = [row for row in rows if row["label"].startswith("eta") and row["query_mode"] == "pair_pair_256"]
    winner = min(eta_rows, key=lambda row: row["whole_process_wall_seconds"])
    require(winner["label"] == "eta128 pair-pair256", "Stage-41 eta winner changed")
    width_1024 = next(row for row in rows if row["label"] == "eta128 pair-pair1024")
    require(width_1024["whole_process_wall_seconds"] > winner["whole_process_wall_seconds"], "Stage-41 rejected width no longer loses")
    local = root / "local-same-target"
    direct_metrics = load(local / "direct.metrics.json", "local same-target direct metrics")
    rho_metrics = load(local / "rho.metrics.json", "local same-target rho metrics")
    direct = next(row for row in json_lines(local / "direct.stdout.jsonl") if row.get("kind") == "relation_rank_summary")
    rho = json_lines(local / "rho.stdout.jsonl")[0]
    require(direct["published_q"] == rho["published_q"] and direct["recovered_fixture_scalar"] == rho["recovered_fixture_scalar"], "local same-target evidence changed")
    same_target = {
        "published_q": direct["published_q"],
        "published_fixture_scalar": direct["published_fixture_scalar"],
        "direct_wall_seconds": direct_metrics["metrics"]["wall_seconds"],
        "rho_wall_seconds": rho_metrics["metrics"]["wall_seconds"],
        "direct_over_rho_wall_ratio": direct_metrics["metrics"]["wall_seconds"] / rho_metrics["metrics"]["wall_seconds"],
        "direct_peak_rss_bytes": direct_metrics["metrics"]["peak_rss_bytes"],
        "rho_peak_rss_bytes": rho_metrics["metrics"]["peak_rss_bytes"],
    }
    return {"pointwise": pointwise(root), "rows": rows, "winner": winner, "local_same_target": same_target}


def workflow_accounting() -> dict[str, Any]:
    value = load(WORKFLOW, "Stage-42 workflow metadata")
    require(value.get("databaseId") == RUN_ID and value.get("headSha") == RUN_COMMIT, "Stage-42 workflow identity changed")
    require(value.get("status") == "completed" and value.get("conclusion") == "success", "Stage-42 workflow is incomplete")
    jobs = value["jobs"]
    require({row["name"] for row in jobs} == {"validate-control", "n53-same-target"} and all(row["conclusion"] == "success" for row in jobs), "Stage-42 workflow jobs changed")
    production = next(row for row in jobs if row["name"] == "n53-same-target")
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
    value = load(ARTIFACTS, "Stage-42 artifact metadata")
    rows = value["artifacts"]
    require(value.get("total_count") == 1 and len(rows) == 1, "Stage-42 artifact count changed")
    row = rows[0]
    require(row.get("id") == ARTIFACT_ID and row.get("digest") == ARTIFACT_DIGEST and row.get("expired") is False, "Stage-42 artifact changed")
    require(row.get("workflow_run", {}).get("id") == RUN_ID and row["workflow_run"].get("head_sha") == RUN_COMMIT, "Stage-42 artifact source changed")
    return {key: row[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}


def result() -> dict[str, Any]:
    tuning_archive = archive_identity(TUNING_ARCHIVE, TUNING_SHA, "Stage-41 tuning archive")
    hosted_archive = archive_identity(HOSTED_ARCHIVE, HOSTED_SHA, "Stage-42 hosted archive")
    with tempfile.TemporaryDirectory() as directory:
        destination = Path(directory)
        tuning = tuning_result(extract(TUNING_ARCHIVE, TUNING_ROOT, destination / "tuning"))
        hosted = extract(HOSTED_ARCHIVE, HOSTED_ROOT, destination / "hosted")
        verification = stage42.verify(hosted / "stage42-build", hosted / "stage42-run")
        embedded = load(hosted / "stage42-verification.json", "embedded Stage-42 verification")
        require(embedded == verification, "Stage-42 embedded verification changed")
        run = load(hosted / "stage42-run/result.json", "Stage-42 run result")
        build_seal = custody.sha256_file(hosted / "stage42-build/result-seal.json", "Stage-42 build seal")
        run_seal = custody.sha256_file(hosted / "stage42-run/result-seal.json", "Stage-42 run seal")
    observation = run["observation"]
    base_keys = (
        "n", "a", "eta", "orbit_columns", "factor_base_points", "base_hash",
        "point_selection", "field_x_values_scanned", "base_construction_ms",
        "support_index_ms", "pair_index_mode", "support_index_entries",
        "support_table_allocated_bytes", "selection_uses_scalar_labels",
    )
    direct_keys = (
        "status", "decomposition_arity", "orbit_columns", "factor_base_points",
        "matrix_columns", "admitted_relations", "target_trials", "setup_ms",
        "collection_ms", "charged_total_ms", "query_mode", "support_queries",
        "query_batch_inversions", "query_canonicalization_maps",
        "query_x_filter_rejections", "full_rank_at_relation", "published_q",
        "published_fixture_scalar", "fixture_scalar_source",
        "recovered_fixture_scalar", "linear_solution_verified",
        "all_relations_group_verified",
    )
    rho_keys = (
        "n", "a", "quotient_mode", "arithmetic_backend", "automorphism_size",
        "published_q", "published_fixture_scalar", "fixture_scalar_source",
        "recovered_fixture_scalar", "walk_steps", "restarts", "setup_ms",
        "walk_ms", "validation_ms", "total_ms", "verified",
        "reference_group_validation",
    )
    compact_observation = {
        "base": {key: observation["base"][key] for key in base_keys},
        "direct": {key: observation["direct"][key] for key in direct_keys},
        "rho": {key: observation["rho"][key] for key in rho_keys},
    }
    return {
        "schema": SCHEMA,
        "status": "hosted_n53_same_target_verified",
        "source_commit": RUN_COMMIT,
        "tuning_archive": tuning_archive,
        "hosted_archive": hosted_archive,
        "workflow": workflow_accounting(),
        "artifact": artifact_accounting(),
        "build_seal_sha256": build_seal,
        "run_seal_sha256": run_seal,
        "verification": verification,
        "tuning": tuning,
        "hosted": {
            "observation": compact_observation,
            "comparison": run["comparison"],
            "direct_process": run["direct_process"],
            "rho_process": run["rho_process"],
            "outer_resources": run["outer_resources"],
            "build_outer_resources": run["build_binding"]["outer_resources"],
            "charged_total_core_seconds_available": run["charged_total_core_seconds_available"],
            "charged_sequential_wall_seconds_available": run["charged_sequential_wall_seconds_available"],
            "maximum_sampled_process_tree_rss_bytes": run["maximum_sampled_process_tree_rss_bytes"],
        },
        "same_public_target": True,
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "fixture_scalar_used_by_collector": False,
        "target_scalar_retained_for_validation": True,
        "pointwise_resource_cap_enforced": False,
        "sat_conflicts": None,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "remaining_blockers": [
            "hosted same-target direct wall remains 12.13 times signed-Frobenius rho",
            "fresh build plus direct remains 37.70 times rho wall",
            "the original autolab pointwise arm did not enforce its memory cap or sample peak RSS",
            "licensed same-instance Magma and unaffiliated reproduction remain absent",
            "the finite n=53 optimization does not change the asymptotic exponent",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-42 result output must be new")
    output.mkdir()
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "hosted_result_frozen", "workflow_run_id": RUN_ID, "workflow_commit": RUN_COMMIT, "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-42 result")
    require(committed == result(), "committed Stage-42 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-42 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(claimed == custody.canonical_sha256(payload) and seal.get("schema") == SEAL_SCHEMA, "Stage-42 result seal changed")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal["inventory"] and custody.canonical_sha256(inventory) == seal["inventory_sha256"], "Stage-42 result inventory changed")
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
    except (OSError, ValueError, KeyError, VerificationError, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage42-result: {error}")


if __name__ == "__main__":
    main()
