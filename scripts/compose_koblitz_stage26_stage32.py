#!/usr/bin/env python3
"""Compose the immutable Stage-26 matrix with the two-row Stage-32 correction."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import copy
from datetime import datetime
import json
import math
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage32_wdsat_correction as stage32_tool
import score_koblitz_stage26_affinity as stage26_score_tool
import verify_koblitz_stage26_result as stage26_verifier


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
BASE_SCORE_ROOT = STAGE / "stage-26-affinity-matrix-result-20260911"
BASE_ARCHIVE = STAGE / "koblitz-stage26-terminal-evidence-successor-01-20260911.tar.gz"
STAGE32_ARCHIVE = STAGE / "koblitz-stage32-wdsat-capacity-evidence-34649735091.tar.gz"
STAGE32_ARCHIVE_SHA = STAGE / "koblitz-stage32-wdsat-capacity-evidence-34649735091.tar.gz.sha256"
STAGE32_WORKFLOW = STAGE / "stage-32-wdsat-workflow-34649735091.json"
STAGE32_ARTIFACTS = STAGE / "stage-32-wdsat-artifacts-34649735091.json"
DEFAULT_OUTPUT = STAGE / "stage-26-affinity-matrix-stage32-corrected-result-20260911"

SCHEMA = "koblitz_stage26_stage32_composed_score.v1"
SEAL_SCHEMA = "koblitz_stage26_stage32_composed_score_seal.v1"
STAGE32_RUN = 34649735091
STAGE32_COMMIT = "8ad01b0c8eb7bd90cb83976a080cea5cd518f0fb"
STAGE32_ARTIFACT_ID = 10283915853
STAGE32_ARTIFACT_DIGEST = "sha256:9eacd0472395b475129e652658ba4917843f60b9b2d17386903f244e61f9843c"
STAGE32_ARCHIVE_ROOT = "koblitz-stage32-run-34649735091-download-20260911"
CORRECTED_CELL = "n59-l9-m3-standard-a1-f0"
CORRECTED_IDS = stage32_tool.FAILED_IDS
REQUIRED_BUFFERS = dict(zip(CORRECTED_IDS, stage32_tool.REQUIRED_BUFFERS, strict=True))


class CompositionError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CompositionError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def timestamp(value: str) -> datetime:
    require(isinstance(value, str) and value.endswith("Z"), "workflow timestamp is invalid")
    return datetime.fromisoformat(value[:-1] + "+00:00")


def archive_identity(archive: Path, sidecar: Path) -> dict[str, Any]:
    fields = sidecar.read_text().split()
    require(
        len(fields) == 2 and fields[1] == str(archive.relative_to(REPO)),
        "Stage-32 archive sidecar is invalid",
    )
    digest = custody.sha256_file(archive, "Stage-32 evidence archive")
    require(fields[0] == digest, "Stage-32 archive hash changed")
    return {
        "path": str(archive.relative_to(REPO)),
        "bytes": archive.stat().st_size,
        "sha256": digest,
    }


def extract_safe(archive: Path, destination: Path) -> Path:
    with tarfile.open(archive, "r:gz") as source:
        members = source.getmembers()
        require(bool(members), f"empty archive: {archive}")
        roots: set[str] = set()
        for member in members:
            name = PurePosixPath(member.name)
            require(not name.is_absolute() and ".." not in name.parts, f"unsafe archive path: {member.name}")
            require(member.isfile() or member.isdir(), f"unsupported archive entry: {member.name}")
            if name.parts:
                roots.add(name.parts[0])
        require(len(roots) == 1, f"archive must have one top-level directory: {archive}")
        source.extractall(destination, filter="data")
    return destination / next(iter(roots))


def workflow_accounting() -> dict[str, Any]:
    workflow = load(STAGE32_WORKFLOW, "Stage-32 workflow metadata")
    require(workflow.get("databaseId") == STAGE32_RUN, "Stage-32 workflow ID changed")
    require(workflow.get("headSha") == STAGE32_COMMIT, "Stage-32 workflow commit changed")
    require(workflow.get("status") == "completed" and workflow.get("conclusion") == "success", "Stage-32 workflow is incomplete")
    jobs = workflow.get("jobs")
    require(isinstance(jobs, list) and {row.get("name") for row in jobs} == {"validate-control", "correction"}, "Stage-32 job inventory changed")
    for row in jobs:
        require(row.get("status") == "completed" and row.get("conclusion") == "success", f"Stage-32 job failed: {row.get('name')}")
    created, updated = timestamp(workflow["createdAt"]), timestamp(workflow["updatedAt"])
    correction = next(row for row in jobs if row["name"] == "correction")
    return {
        "run_id": STAGE32_RUN,
        "head_sha": STAGE32_COMMIT,
        "url": workflow.get("url"),
        "workflow_wall_seconds": (updated - created).total_seconds(),
        "correction_job_id": correction.get("databaseId"),
        "correction_job_url": correction.get("url"),
        "correction_job_wall_seconds": (timestamp(correction["completedAt"]) - timestamp(correction["startedAt"])).total_seconds(),
    }


def artifact_accounting() -> dict[str, Any]:
    value = load(STAGE32_ARTIFACTS, "Stage-32 artifact metadata")
    rows = value.get("artifacts")
    require(value.get("total_count") == 1 and isinstance(rows, list) and len(rows) == 1, "Stage-32 artifact count changed")
    row = rows[0]
    require(row.get("id") == STAGE32_ARTIFACT_ID, "Stage-32 artifact ID changed")
    require(row.get("name") == f"koblitz-stage32-wdsat-capacity-{STAGE32_RUN}", "Stage-32 artifact name changed")
    require(row.get("digest") == STAGE32_ARTIFACT_DIGEST and row.get("expired") is False, "Stage-32 artifact digest or state changed")
    source = row.get("workflow_run", {})
    require(source.get("id") == STAGE32_RUN and source.get("head_sha") == STAGE32_COMMIT, "Stage-32 artifact source changed")
    require(isinstance(row.get("size_in_bytes"), int) and row["size_in_bytes"] > 0, "Stage-32 artifact size is invalid")
    return {key: row[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}


def unique_task(root: Path, blind_id: str) -> tuple[Path, dict[str, Any]]:
    paths = list(root.rglob(f"*{blind_id}/task-result.json"))
    require(len(paths) == 1, f"expected one task for {blind_id}, found {len(paths)}")
    return paths[0], load(paths[0], f"task {blind_id}")


def process_metrics(row: dict[str, Any], context: str) -> dict[str, Any]:
    process = row.get("process")
    require(isinstance(process, dict), f"{context} lacks a process receipt")
    metrics = process.get("metrics")
    require(isinstance(metrics, dict), f"{context} lacks process metrics")
    for name in ("total_core_seconds", "wall_seconds"):
        require(isinstance(metrics.get(name), (int, float)) and metrics[name] >= 0, f"{context} has invalid {name}")
    require(isinstance(metrics.get("peak_rss_bytes"), int) and metrics["peak_rss_bytes"] > 0, f"{context} has invalid peak RSS")
    return metrics


def correction_rows(base_archive_root: Path, stage32_root: Path, base: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    replacements: list[dict[str, Any]] = []
    corrected = copy.deepcopy(base["rows"])
    for blind_id in CORRECTED_IDS:
        original_path, original = unique_task(base_archive_root, blind_id)
        corrected_path, successor = unique_task(stage32_root / "stage32-correction", blind_id)
        base_rows = [row for row in corrected if row["blind_instance_id"] == blind_id and row["backend"] == "wdsat"]
        require(len(base_rows) == 1, f"base score WDSat row missing: {blind_id}")
        score_row = base_rows[0]
        require(score_row["cell_id"] == CORRECTED_CELL and score_row["solver_status"] == "solver_error", f"base score row is not the expected error: {blind_id}")
        original_backends = [row for row in original.get("backends", []) if row.get("solver") == "wdsat"]
        require(len(original_backends) == 1, f"original WDSat task row missing: {blind_id}")
        original_row = original_backends[0]
        require(original_row.get("status") == "solver_error" and original_row.get("returncode") == -6 and original_row.get("timed_out") is False, f"original failure changed: {blind_id}")
        require(successor.get("blind_instance_id") == blind_id and successor.get("cell_id") == CORRECTED_CELL, f"Stage-32 task identity changed: {blind_id}")
        require(successor.get("source_instance_id") == original.get("source_instance_id") and successor.get("target") == original.get("target"), f"Stage-32 task source changed: {blind_id}")
        require(successor.get("source_before") == original.get("source_artifacts_before") == successor.get("source_after"), f"Stage-32 frozen source differs: {blind_id}")
        require(successor.get("requirements", {}).get("max_buffer_size") == REQUIRED_BUFFERS[blind_id], f"Stage-32 buffer requirement changed: {blind_id}")
        new_row = successor.get("result", {})
        require(new_row.get("solver") == "wdsat" and new_row.get("status") == "timeout_inconclusive", f"Stage-32 terminal status changed: {blind_id}")
        require(new_row.get("timed_out") is True and new_row.get("returncode") == -15 and new_row.get("conflicts") is None, f"Stage-32 timeout receipt changed: {blind_id}")
        require(score_row["classification"] == "inconclusive", f"original classification changed: {blind_id}")
        score_row["solver_status"] = "timeout_inconclusive"
        replacements.append({
            "blind_instance_id": blind_id,
            "cell_id": CORRECTED_CELL,
            "source_instance_id": original["source_instance_id"],
            "target": original["target"],
            "original_task_path_in_archive": str(original_path.relative_to(base_archive_root)),
            "corrected_task_path_in_archive": str(corrected_path.relative_to(stage32_root)),
            "required_max_buffer_size": REQUIRED_BUFFERS[blind_id],
            "original_status": "solver_error",
            "corrected_status": "timeout_inconclusive",
            "classification_before_and_after": "inconclusive",
            "original_process_metrics": process_metrics(original_row, f"original WDSat {blind_id}"),
            "corrected_process_metrics": process_metrics(new_row, f"corrected WDSat {blind_id}"),
        })
    return corrected, replacements


def aggregate_rows(rows: list[dict[str, Any]]) -> tuple[dict[str, int], list[dict[str, Any]]]:
    grouped: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    for row in rows:
        key = (row["cell_id"], row["backend"], row["target_class"])
        grouped[key][row["classification"]] += 1
        grouped[key][f"status:{row['solver_status']}"] += 1
    per_cell = [
        {"cell_id": key[0], "backend": key[1], "target_class": key[2], "counts": dict(sorted(value.items()))}
        for key, value in sorted(grouped.items())
    ]
    return dict(sorted(Counter(row["classification"] for row in rows).items())), per_cell


def corrected_backend_resources(base: dict[str, Any], replacements: list[dict[str, Any]]) -> dict[str, Any]:
    resources = copy.deepcopy(base["backend_resources"])
    wdsat = resources["wdsat"]
    old = [row["original_process_metrics"] for row in replacements]
    new = [row["corrected_process_metrics"] for row in replacements]
    wdsat["total_core_seconds"] = round(wdsat["total_core_seconds"] - math.fsum(row["total_core_seconds"] for row in old) + math.fsum(row["total_core_seconds"] for row in new), 12)
    wdsat["summed_process_wall_seconds"] = round(wdsat["summed_process_wall_seconds"] - math.fsum(row["wall_seconds"] for row in old) + math.fsum(row["wall_seconds"] for row in new), 12)
    wdsat["peak_rss_bytes"] = max(wdsat["peak_rss_bytes"], *(row["peak_rss_bytes"] for row in new))
    return resources


def compose_result() -> dict[str, Any]:
    base_verification = stage26_verifier.verify()
    base = load(BASE_SCORE_ROOT / "score.json", "Stage-26 score")
    base_seal = load(BASE_SCORE_ROOT / "score-seal.json", "Stage-26 score seal")
    require(base_verification["score_sha256"] == custody.sha256_file(BASE_SCORE_ROOT / "score.json", "Stage-26 score"), "Stage-26 score binding changed")
    archive = archive_identity(STAGE32_ARCHIVE, STAGE32_ARCHIVE_SHA)
    workflow = workflow_accounting()
    artifact = artifact_accounting()
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory)
        base_root = extract_safe(BASE_ARCHIVE, temp / "stage26")
        stage32_root = extract_safe(STAGE32_ARCHIVE, temp / "stage32")
        require(stage32_root.name == STAGE32_ARCHIVE_ROOT, "Stage-32 archive root changed")
        correction_verification = stage32_tool.verify(stage32_root / "stage32-correction")
        correction = load(stage32_root / "stage32-correction/result.json", "Stage-32 correction result")
        require(correction.get("packet_verification") == base.get("packet_verification"), "Stage-32 packet differs from Stage-26")
        require(correction.get("inputs") == CORRECTED_IDS and correction.get("corrected_statuses") == {"timeout_inconclusive": 2}, "Stage-32 correction set changed")
        require(correction.get("affinity", {}).get("child_effective") == [0] and correction.get("affinity", {}).get("parent_effective") == [0], "Stage-32 one-CPU affinity changed")
        build = correction.get("wdsat_build", {})
        require(build.get("config_sha256") == "4a73c3b5a14ded98f749d282b597594ae3ac05355a39cadcb9345a90bf735352", "Stage-32 config changed")
        require(build.get("limits", {}).get("max_buffer_size") == 32804, "Stage-32 build capacity changed")
        require(build.get("implementation_state", {}).get("commit") == STAGE32_COMMIT and build.get("implementation_state", {}).get("dirty") is False, "Stage-32 implementation state changed")
        corrected_rows, replacements = correction_rows(base_root, stage32_root, base)
        acquisition = stage26_score_tool.parse_gnu_time(stage32_root / "stage32-acquisition/wdsat-clone.time")
        correction_seal_sha256 = custody.sha256_file(
            stage32_root / "stage32-correction/result-seal.json", "Stage-32 result seal"
        )
    classification_counts, per_cell = aggregate_rows(corrected_rows)
    require(classification_counts == base["classification_counts"] == {"inconclusive": 219, "true_negative": 120, "true_positive": 141}, "classification totals changed")
    original_error_rows = sum(row["solver_status"] == "solver_error" for row in base["rows"])
    corrected_error_rows = sum(row["solver_status"] == "solver_error" for row in corrected_rows)
    require(original_error_rows == 2 and corrected_error_rows == 0, "solver-error replacement count changed")
    outer = correction["outer_resources"]
    build_resources = correction["wdsat_build"]["resources"]
    increment_core = math.fsum((outer["total_core_seconds"], build_resources["total_core_seconds"], acquisition["total_core_seconds"]))
    result = copy.deepcopy(base)
    result.update({
        "schema": SCHEMA,
        "status": "complete_verified_stage32_corrected_four_cell_panel",
        "base_stage26_score": {
            "path": str((BASE_SCORE_ROOT / "score.json").relative_to(REPO)),
            "bytes": (BASE_SCORE_ROOT / "score.json").stat().st_size,
            "sha256": custody.sha256_file(BASE_SCORE_ROOT / "score.json", "Stage-26 score"),
            "score_inventory_sha256": base_seal["inventory_sha256"],
            "verification": base_verification,
            "immutable_original_solver_error_rows": 2,
        },
        "rows": corrected_rows,
        "per_cell_backend_class": per_cell,
        "classification_counts": classification_counts,
        "backend_resources": corrected_backend_resources(base, replacements),
        "stage32_wdsat_correction": {
            "archive": archive,
            "artifact": artifact,
            "workflow": workflow,
            "verification": correction_verification,
            "result_seal_sha256": correction_seal_sha256,
            "replacements": replacements,
            "corrected_rows": 2,
            "solver_error_rows_before": 2,
            "solver_error_rows_after": 0,
            "classification_counts_changed": False,
            "acquisition": acquisition,
            "build_resources": build_resources,
            "outer_resources": outer,
            "charged_increment_total_core_seconds": increment_core,
            "unmetered": [
                "Git checkout after the metered WDSat clone",
                "Stage-26 backend artifact download and extraction",
                "workflow checkout and artifact upload CPU and memory",
            ],
        },
        "campaign_outer_resources": {
            "stage26_summed_single_core_elapsed_seconds": base["cell_outer_resources"]["summed_single_core_elapsed_seconds"],
            "stage32_single_core_elapsed_seconds": outer["single_core_elapsed_seconds"],
            "summed_single_core_elapsed_seconds": base["cell_outer_resources"]["summed_single_core_elapsed_seconds"] + outer["single_core_elapsed_seconds"],
            "stage26_summed_total_core_seconds": base["cell_outer_resources"]["summed_total_core_seconds"],
            "stage32_total_core_seconds": outer["total_core_seconds"],
            "summed_outer_total_core_seconds": base["cell_outer_resources"]["summed_total_core_seconds"] + outer["total_core_seconds"],
            "maximum_sampled_process_tree_rss_bytes": max(base["cell_outer_resources"]["maximum_sampled_process_tree_rss_bytes"], outer["sampled_peak_process_tree_rss_bytes"]),
            "stage26_workflow_wall_seconds": base["workflow"]["workflow_wall_seconds"],
            "stage32_workflow_wall_seconds": workflow["workflow_wall_seconds"],
            "summed_workflow_wall_seconds": base["workflow"]["workflow_wall_seconds"] + workflow["workflow_wall_seconds"],
        },
        "maximum_recorded_individual_or_tree_rss_bytes": max(
            base["tool_accounting"]["build_maximum_individual_process_rss_bytes"],
            base["tool_accounting"]["acquisition_maximum_process_rss_bytes"],
            base["cell_outer_resources"]["maximum_sampled_process_tree_rss_bytes"],
            base["matched_direct_mitm"]["resources"]["maximum_sampled_process_tree_rss_bytes"],
            outer["sampled_peak_process_tree_rss_bytes"],
            build_resources["largest_child_peak_rss_bytes"],
            acquisition["peak_rss_bytes"],
        ),
        "charged_total_core_seconds_available": base["charged_total_core_seconds_available"] + increment_core,
        "charged_scope": base["charged_scope"] + "; plus the metered Stage-32 WDSat source clone, fresh corrected build, source verification, and two one-CPU correction attempts",
        "corrected_panel_solver_error_rows": 0,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    })
    audit = copy.deepcopy(result["completion_gate_audit"])
    audit["2_same_instance_backend_matrix"]["proved"].append("the two n=59 WDSat capacity failures were rerun to clean timeout terminals with an exact-capacity fresh build")
    audit["3_resource_fields"]["proved"].append("Stage-32 correction acquisition, build, one-CPU outer envelope, process-tree memory, process wall, and workflow wall")
    result["completion_gate_audit"] = audit
    result["remaining_blockers"] = [
        "licensed Magma F4 has not executed the same 160-input packet",
        "the four-cell PDP panel is not a complete end-to-end index-calculus run",
        "factor-base discovery, relation collection, linear algebra, and rho are bound in separate controls rather than one same-instance full-cost experiment",
        "some workflow checkout, artifact transfer, package, and compiler aggregate CPU or memory costs remain unmetered",
        "unaffiliated reproduction and novelty review remain absent",
    ]
    return result


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "composition output must be new")
    result = compose_result()
    output.mkdir(parents=False)
    custody.write_json_new(output / "score.json", result)
    inventory = custody.all_regular_inventory(output, {"score-seal.json"})
    payload = {
        "schema": SEAL_SCHEMA,
        "status": "composition_frozen",
        "stage26_run_id": stage26_score_tool.EXPECTED_RUN,
        "stage32_run_id": STAGE32_RUN,
        "inventory": inventory,
        "inventory_sha256": custody.canonical_sha256(inventory),
    }
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "score-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    score = load(output / "score.json", "composed score")
    expected = compose_result()
    require(score == expected, "committed composed score differs from source evidence")
    seal = load(output / "score-seal.json", "composed score seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "composed score seal is invalid")
    inventory = custody.all_regular_inventory(output, {"score-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), "composed score inventory changed")
    return {
        "schema": "koblitz_stage26_stage32_composed_verification.v1",
        "status": "base_matrix_and_two_row_capacity_correction_verified",
        "stage26_run_id": stage26_score_tool.EXPECTED_RUN,
        "stage32_run_id": STAGE32_RUN,
        "instances": score["instances"],
        "backend_rows": score["backend_rows"],
        "corrected_rows": score["stage32_wdsat_correction"]["corrected_rows"],
        "solver_error_rows_after": score["corrected_panel_solver_error_rows"],
        "classification_counts": score["classification_counts"],
        "wdsat_total_core_seconds_corrected_panel": score["backend_resources"]["wdsat"]["total_core_seconds"],
        "stage32_charged_increment_total_core_seconds": score["stage32_wdsat_correction"]["charged_increment_total_core_seconds"],
        "charged_total_core_seconds_available": score["charged_total_core_seconds_available"],
        "summed_workflow_wall_seconds": score["campaign_outer_resources"]["summed_workflow_wall_seconds"],
        "maximum_recorded_individual_or_tree_rss_bytes": score["maximum_recorded_individual_or_tree_rss_bytes"],
        "licensed_magma_complete": score["licensed_magma_f4_same_instance_panel_complete"],
        "independent_external_reproduction_satisfied": score["independent_external_reproduction_satisfied"],
        "full_cost_gate_passed": score["full_cost_gate_passed"],
        "koblitz_index_calculus_sota": score["koblitz_index_calculus_sota"],
        "score_sha256": custody.sha256_file(output / "score.json", "composed score"),
        "score_inventory_sha256": seal["inventory_sha256"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    compose_parser = sub.add_parser("compose")
    compose_parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    verify_parser = sub.add_parser("verify")
    verify_parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    try:
        value = freeze(args.output.resolve()) if args.command == "compose" else verify(args.output.resolve())
        print(json.dumps(value, indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        CompositionError,
        custody.PhaseBError,
        stage32_tool.Stage32Error,
        stage26_score_tool.Stage26ScoreError,
        stage26_verifier.VerificationError,
    ) as error:
        raise SystemExit(f"stage26-stage32-composition: {error}")


if __name__ == "__main__":
    main()
