#!/usr/bin/env python3
"""Verify, score, and compose four sealed Stage-26 affinity cells."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime
import json
import math
from pathlib import Path
import statistics
from typing import Any

import koblitz_stage26_affinity_inputs as packet_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_stage26_affinity_cell as cell_tool
import run_koblitz_stage25_single_core as stage25_tool
import score_koblitz_blind_pdp_phase_b as phase_b_score
import verify_koblitz_stage27_result as stage27_verifier


REPO = Path(__file__).resolve().parents[1]
DEFAULT_PROTOCOL = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-balanced-pdp-phase-b-protocol.json"
SCHEMA = "koblitz_stage26_affinity_score.v1"
SEAL_SCHEMA = "koblitz_stage26_affinity_score_seal.v1"
EXPECTED_RUN = 34632018379
EXPECTED_COMMIT = "03968a5da2a511abb723652529de34dc60e10942"
DIRECT_MITM_SCORE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-27-direct-mitm-result-20260911/score.json"
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"


class Stage26ScoreError(RuntimeError):
    pass


def read_json(path: Path, context: str) -> dict[str, Any]:
    value, _ = phase_b.read_json(path, context)
    if not isinstance(value, dict):
        raise Stage26ScoreError(f"{context} must be an object")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage26ScoreError(message)


def parse_timestamp(value: str) -> datetime:
    require(isinstance(value, str) and value.endswith("Z"), "workflow timestamp is invalid")
    return datetime.fromisoformat(value[:-1] + "+00:00")


def elapsed_time(text: str) -> float:
    parts = text.strip().split(":")
    require(2 <= len(parts) <= 3, "GNU time elapsed value is invalid")
    values = [float(value) for value in parts]
    return values[0] * 3600 + values[1] * 60 + values[2] if len(values) == 3 else values[0] * 60 + values[1]


def parse_gnu_time(path: Path) -> dict[str, Any]:
    fields: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith("User time (seconds):"):
            fields["user"] = line.rsplit(":", 1)[1].strip()
        elif line.startswith("System time (seconds):"):
            fields["system"] = line.rsplit(":", 1)[1].strip()
        elif line.startswith("Elapsed (wall clock) time"):
            fields["elapsed"] = line.split("): ", 1)[1].strip()
        elif line.startswith("Maximum resident set size (kbytes):"):
            fields["rss"] = line.rsplit(":", 1)[1].strip()
        elif line.startswith("Exit status:"):
            fields["exit"] = line.rsplit(":", 1)[1].strip()
    require(set(fields) == {"user", "system", "elapsed", "rss", "exit"}, f"incomplete GNU time receipt: {path}")
    user, system = float(fields["user"]), float(fields["system"])
    require(fields["exit"] == "0" and min(user, system) >= 0, f"failed acquisition receipt: {path}")
    return {
        "path": path.name,
        "sha256": phase_b.sha256_file(path, "acquisition time receipt"),
        "user_seconds": user,
        "system_seconds": system,
        "total_core_seconds": user + system,
        "wall_seconds": elapsed_time(fields["elapsed"]),
        "peak_rss_bytes": int(fields["rss"]) * 1024,
        "exit_status": 0,
    }


def tool_costs(tools_root: Path) -> dict[str, Any]:
    acquisition_paths = sorted((tools_root / "acquisition").glob("*.time"))
    require([path.name for path in acquisition_paths] == [
        "cadiback.time", "cadical.time", "cms.time", "rust-crates.time", "wdsat.time"
    ], "tool acquisition receipt inventory changed")
    acquisition = [parse_gnu_time(path) for path in acquisition_paths]
    build_paths = {
        "wdsat": tools_root / "build-evidence/wdsat/receipt.json",
        "rust": tools_root / "build-evidence/rust/receipt.json",
        "cryptominisat": tools_root / "build-evidence/cryptominisat/receipt.json",
    }
    builds: dict[str, Any] = {}
    for name, path in build_paths.items():
        receipt = read_json(path, f"{name} build receipt")
        require(receipt.get("status") == "completed", f"{name} build did not complete")
        resources = receipt.get("resources")
        require(isinstance(resources, dict), f"{name} build lacks resources")
        builds[name] = {
            "receipt_sha256": phase_b.sha256_file(path, f"{name} build receipt"),
            "resources": resources,
        }
    return {
        "acquisition": acquisition,
        "acquisition_total_core_seconds": math.fsum(row["total_core_seconds"] for row in acquisition),
        "acquisition_summed_wall_seconds": math.fsum(row["wall_seconds"] for row in acquisition),
        "acquisition_maximum_process_rss_bytes": max(row["peak_rss_bytes"] for row in acquisition),
        "builds": builds,
        "build_total_core_seconds": math.fsum(row["resources"]["total_core_seconds"] for row in builds.values()),
        "build_summed_process_wall_seconds": math.fsum(row["resources"]["summed_process_wall_seconds"] for row in builds.values()),
        "build_maximum_individual_process_rss_bytes": max(
            row["resources"].get("largest_child_peak_rss_bytes", row["resources"].get("peak_rss_bytes", 0))
            for row in builds.values()
        ),
        "excluded": [
            "preinstalled hosted-runner operating system and toolchain acquisition",
            "apt package installation CPU and memory beyond hosted job timestamps",
            "simultaneous aggregate memory of parallel compiler children",
            "checkout commands adjacent to metered source clones",
            "the workflow-level packet precheck before the runner's separately charged packet verification",
        ],
    }


def workflow_accounting(path: Path, expected_cells: set[str]) -> dict[str, Any]:
    workflow = read_json(path, "workflow metadata")
    require(workflow.get("databaseId") == EXPECTED_RUN, "workflow run ID changed")
    require(workflow.get("headSha") == EXPECTED_COMMIT, "workflow source commit changed")
    require(workflow.get("status") == "completed" and workflow.get("conclusion") == "success", "workflow run is not successfully complete")
    jobs = workflow.get("jobs")
    require(isinstance(jobs, list), "workflow metadata lacks jobs")
    cell_jobs: dict[str, dict[str, Any]] = {}
    for job in jobs:
        name = job.get("name", "")
        if name.startswith("one-cpu-cell (") and name.endswith(")"):
            cell_jobs[name[len("one-cpu-cell ("):-1]] = job
    require(set(cell_jobs) == expected_cells, "workflow cell job inventory changed")
    for cell, job in cell_jobs.items():
        require(job.get("status") == "completed" and job.get("conclusion") == "success", f"workflow cell {cell} is incomplete")
    created = parse_timestamp(workflow["createdAt"])
    updated = parse_timestamp(workflow["updatedAt"])
    starts = [parse_timestamp(job["startedAt"]) for job in cell_jobs.values()]
    ends = [parse_timestamp(job["completedAt"]) for job in cell_jobs.values()]
    return {
        "run_id": EXPECTED_RUN,
        "head_sha": EXPECTED_COMMIT,
        "url": workflow.get("url"),
        "workflow_wall_seconds": (updated - created).total_seconds(),
        "parallel_cell_job_span_seconds": (max(ends) - min(starts)).total_seconds(),
        "cell_jobs": {
            cell: {
                "job_id": job.get("databaseId"),
                "started_at": job["startedAt"],
                "completed_at": job["completedAt"],
                "job_wall_seconds": (parse_timestamp(job["completedAt"]) - parse_timestamp(job["startedAt"])).total_seconds(),
                "url": job.get("url"),
            }
            for cell, job in sorted(cell_jobs.items())
        },
    }


def matched_direct_mitm() -> dict[str, Any]:
    verification = stage27_verifier.verify()
    score = read_json(DIRECT_MITM_SCORE, "matched direct-MITM score")
    require(score.get("status") == "complete_verified_truth_scored_direct_mitm_panel", "matched direct-MITM score is incomplete")
    require(score.get("instances") == 160 and score.get("outcomes") == 160, "matched direct-MITM counts changed")
    require(score.get("classification_counts") == {"true_negative": 80, "true_positive": 80}, "matched direct-MITM classifications changed")
    require(score.get("full_cost_gate_passed") is False and score.get("koblitz_index_calculus_sota") is False, "matched direct-MITM score widened its claim")
    return {
        "verification": verification,
        "score": {
            "path": str(DIRECT_MITM_SCORE.relative_to(REPO)),
            "bytes": DIRECT_MITM_SCORE.stat().st_size,
            "sha256": phase_b.sha256_file(DIRECT_MITM_SCORE, "matched direct-MITM score"),
        },
        "workflow": score["workflow"],
        "instances": score["instances"],
        "outcomes": score["outcomes"],
        "classification_counts": score["classification_counts"],
        "resources": score["resources"],
        "operations": score["operations"],
        "conflicts": score["conflicts"],
        "conflict_semantics": score["conflict_semantics"],
    }


def unknown_scalar_control() -> dict[str, Any]:
    names = {
        "result": "stage-25-single-core-result.json",
        "verification": "stage-25-single-core-verification.json",
        "artifact": "stage-25-single-core-artifact.json",
        "math_replay": "stage-25-single-core-math-replay.json",
        "affinity": "stage-25-single-core-affinity.json",
        "result_seal": "stage-25-single-core-result-seal.json",
    }
    values = {key: read_json(STAGE / name, f"Stage-25 {key}") for key, name in names.items()}
    result = stage25_tool.validate_result(values["result"])
    affinity = stage25_tool.validate_affinity(values["affinity"])
    require(result["affinity"] == affinity, "Stage-25 result and affinity differ")
    artifact = values["artifact"]
    require(artifact.get("status") == "github_artifact_downloaded_and_portable_stage23_verified", "Stage-25 artifact is incomplete")
    for name, identity in artifact.get("committed_files", {}).items():
        path = STAGE / name
        require(path.is_file() and path.stat().st_size == identity["bytes"] and phase_b.sha256_file(path, "Stage-25 committed file") == identity["sha256"], "Stage-25 committed file changed")
    verification = values["verification"]
    require(verification.get("status") == "single_cpu_affinity_receipt_verified" and verification.get("completed_rows") == 5, "Stage-25 verification is incomplete")
    require(verification.get("single_core_elapsed_seconds") == result["measurements"]["single_core_elapsed_seconds"], "Stage-25 single-core time changed")
    replay = values["math_replay"]
    require(replay.get("status") == "PASS" and replay.get("check_count") == 1251, "Stage-25 retained-math replay is incomplete")
    require(replay.get("retained_mathematical_witness_replay_completed") is True and replay.get("independent_mathematical_payload_replay_completed") is False, "Stage-25 replay boundary changed")
    require(replay.get("attempt_totals") == {"attempts": 437, "conflicts": 28422672, "models": 252, "relation_found": 252, "solver_calls": 437, "unknown": 185}, "Stage-25 relation totals changed")
    return {
        "file_identities": {key: {"path": name, "bytes": (STAGE / name).stat().st_size, "sha256": phase_b.sha256_file(STAGE / name, f"Stage-25 {key}")} for key, name in names.items()},
        "status": result["status"],
        "profile": result["profile"],
        "completed_unknown_scalar_targets": 5,
        "factor_base_logs_known_by_construction": False,
        "target_scalars_known_by_construction": False,
        "measurements": result["measurements"],
        "attempt_totals": replay["attempt_totals"],
        "factor_base": replay["factor_base"],
        "discovery": replay["discovery"],
        "rho_totals": replay["rho_totals"],
        "rho_charge_totals": replay["rho_charge_totals"],
        "ratios": replay["ratios"],
        "retained_math_checks": replay["check_count"],
        "retained_mathematical_witness_replay_completed": True,
        "independent_mathematical_payload_replay_completed": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def natural_relation_yield_control() -> dict[str, Any]:
    summary_path = STAGE / "stage-21-relation-yield-result-summary-20260910.json"
    seal_path = STAGE / "stage-21-relation-yield-result-seal-20260910.json"
    summary = read_json(summary_path, "Stage-21 relation-yield summary")
    seal = read_json(seal_path, "Stage-21 relation-yield seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(claimed == phase_b.canonical_sha256(payload), "Stage-21 relation-yield seal is invalid")
    repository_summary = seal.get("repository_artifacts", {}).get("summary", {})
    require(repository_summary.get("path") == str(summary_path.relative_to(REPO)), "Stage-21 summary path changed")
    require(repository_summary.get("bytes") == summary_path.stat().st_size and repository_summary.get("sha256") == phase_b.sha256_file(summary_path, "Stage-21 summary"), "Stage-21 summary identity changed")
    factor_base = summary.get("factor_base", {})
    require(
        factor_base.get("selection_used_target") is False
        and factor_base.get("selection_used_relation_yield") is False
        and factor_base.get("factor_base_discrete_log_labels_constructed") is False
        and factor_base.get("target_subgroup_enumerated_for_selection") is False,
        "Stage-21 algebraic factor-base boundary changed",
    )
    natural = summary.get("measurement", {}).get("natural", {})
    require(natural.get("targets") == 256 and natural.get("hits") == 163 and natural.get("misses") == 93, "Stage-21 natural-yield totals changed")
    admission = summary.get("verification_and_admission", {})
    require(admission.get("measurement_admission_status") == "pending_independent_payload_replay" and admission.get("scientific_measurement_admitted") is False, "Stage-21 admission boundary changed")
    return {
        "summary": {"path": str(summary_path.relative_to(REPO)), "bytes": summary_path.stat().st_size, "sha256": phase_b.sha256_file(summary_path, "Stage-21 summary")},
        "seal": {"path": str(seal_path.relative_to(REPO)), "bytes": seal_path.stat().st_size, "sha256": phase_b.sha256_file(seal_path, "Stage-21 seal")},
        "status": summary["status"],
        "factor_base": factor_base,
        "natural": natural,
        "planted_sat": summary["measurement"]["planted_sat"],
        "proven_unsat": summary["measurement"]["proven_unsat"],
        "internal_timing_seconds": summary["internal_timing_seconds"],
        "operation_counts": summary["operation_counts"],
        "resource_totals": summary["resource_totals"],
        "verification_and_admission": admission,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def completion_gate_audit() -> dict[str, Any]:
    return {
        "1_full_cost_accounting": {
            "status": "partial",
            "proved": [
                "Stage-21 algebraic factor-base and exact pair-oracle process resources",
                "Stage-25 factor-base discovery, relation collection, linear solve, IC, rho, and one-CPU outer envelope",
                "Stage-26 and Stage-27 matched PDP cell, build, acquisition, conflict, operation, wall, and memory records",
            ],
            "missing": [
                "one same-instance end-to-end index-calculus experiment through relation collection and linear algebra at n=31, n=41, and the larger PDP regime",
                "hosted package-installation CPU and memory",
                "simultaneous aggregate memory for parallel compiler processes",
            ],
        },
        "2_same_instance_backend_matrix": {
            "status": "partial",
            "proved": ["native XOR SAT, WDSat, CryptoMiniSat, and direct MITM on the same 160 inputs", "standard and GGMP factor-base cells at n=31"],
            "missing": ["licensed Magma F4 execution on all 160 inputs"],
        },
        "3_resource_fields": {
            "status": "partial",
            "proved": ["one-CPU elapsed, total core-seconds, solver conflicts where exposed, process RSS, sampled process-tree RSS, process wall, and workflow wall"],
            "missing": ["Magma process resources", "CPU and memory for hosted package installation and adjacent checkout commands"],
        },
        "4_scaling": {
            "status": "finite_pdp_complete_not_end_to_end_ic",
            "proved": ["balanced PDP cells at n=31, n=41, and n=59"],
            "missing": ["end-to-end index-calculus scaling at those degrees"],
        },
        "5_unknown_scalar": {
            "status": "finite_degree23_complete",
            "proved": ["five public unknown-scalar targets completed without constructed target scalars or factor-base log labels"],
            "missing": ["larger unknown-scalar end-to-end regimes"],
        },
        "6_pollard_rho": {
            "status": "finite_degree23_complete_full_scope_partial",
            "proved": ["same-target signed-Frobenius/negation rho controls for all five unknown-scalar targets"],
            "missing": ["one fully unified cost comparison at the n=31, n=41, and larger end-to-end regimes"],
        },
        "7_external_review": {
            "status": "missing",
            "proved": ["fresh hosted execution of project-authored workflows"],
            "missing": ["unaffiliated reproduction", "source-pinned novelty verdict"],
        },
        "overall": "strong internal engineering and public toy-research improvement; not a Koblitz index-calculus SOTA",
    }


def validate_task_source(task_root: Path, packet: Path, packet_record: dict[str, Any]) -> None:
    expected = {Path(row["path"]).name: row for row in packet_record["files"]}
    instance = task_root / "instance"
    require(set(path.name for path in instance.iterdir()) == set(expected), "task source inventory changed")
    for name, identity in expected.items():
        data = phase_b.regular_file_bytes(instance / name, "cell task source")
        require(len(data) == identity["bytes"] and phase_b.sha256_bytes(data) == identity["sha256"], "cell task source differs from packet")


def backend_metrics(row: dict[str, Any]) -> list[dict[str, Any]]:
    values = [row["process"]["metrics"]]
    validation = row.get("point_witness_validation")
    if isinstance(validation, dict):
        values.append(validation["process"]["metrics"])
    return values


def score(args: argparse.Namespace) -> dict[str, Any]:
    if args.output.exists() or args.output.is_symlink():
        raise Stage26ScoreError("score output must be new")
    packet_check = packet_tool.verify(args.packet)
    packet_manifest = read_json(args.packet / "packet-manifest.json", "packet manifest")
    blind_bundle = read_json(args.packet / "source/blind-bundle.json", "blind bundle")
    protocol, _ = phase_b.read_json(args.protocol, "Phase-B protocol")
    phase_b.validate_protocol(protocol)
    raw_seal, oracle = phase_b_score.load_authenticated_oracle(args.phase_a_seal, args.oracle_ledger, protocol)
    truth = phase_b_score.oracle_index(oracle, blind_bundle, protocol)
    cell_roots: dict[str, Path] = {}
    for root in args.cell_root:
        result = read_json(root / "result.json", "cell result")
        cell_id = result.get("cell_id")
        require(cell_id in packet_manifest["cells"] and cell_id not in cell_roots, "duplicate or unexpected cell root")
        cell_roots[cell_id] = root
    require(set(cell_roots) == set(packet_manifest["cells"]), "four-cell result set is incomplete")
    packet_by_id = {row["blind_instance_id"]: row for row in packet_manifest["instances"]}
    blind_by_id = {row["blind_instance_id"]: row for row in blind_bundle["instances"]}
    rows: list[dict[str, Any]] = []
    grouped: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    backend_resources: dict[str, dict[str, Any]] = {
        name: {"total_core_seconds": 0.0, "summed_process_wall_seconds": 0.0, "peak_rss_bytes": 0, "conflicts": []}
        for name in phase_b.BACKENDS
    }
    cells: dict[str, Any] = {}
    for cell_id, root in sorted(cell_roots.items()):
        verification = cell_tool.verify(root)
        result = read_json(root / "result.json", "cell result")
        require(result.get("cell_id") == cell_id and result.get("instances") == 40 and result.get("backend_outcomes") == 120, "cell result counts changed")
        require(result.get("packet_verification") == packet_check, "cell packet verification changed")
        cell_rows = []
        for ordinal, blind_id in enumerate(packet_manifest["cells"][cell_id]):
            task_root = root / "tasks" / f"{ordinal:06d}-{blind_id}"
            task = read_json(task_root / "task-result.json", "cell task result")
            packet_record, blind = packet_by_id[blind_id], blind_by_id[blind_id]
            require(task.get("schema") == "koblitz_stage26_frozen_export_task.v1", "cell task schema changed")
            require(task.get("ordinal") == ordinal and task.get("blind_instance_id") == blind_id, "cell task identity changed")
            require(task.get("cell_id") == cell_id and task.get("target") == blind["target"], "cell task target changed")
            require(task.get("source_instance_id") == packet_record["source_instance_id"], "cell source ID changed")
            require(task.get("status") == "terminal_outputs_recorded" and task.get("source_artifacts_unchanged") is True, "cell task is incomplete")
            require(task.get("source_artifacts_before") == task.get("source_artifacts_after"), "solver changed cell source")
            validate_task_source(task_root, args.packet, packet_record)
            task_manifest = read_json(task_root / "instance/manifest.json", "cell task manifest")
            source_variables = task_manifest.get("source_variables")
            require(isinstance(source_variables, int) and source_variables > 0, "cell task source-variable count is invalid")
            admission_task = dict(task)
            admission_task["wdsat_requirements"] = {"max_anf_id": source_variables + 1}
            source = task.get("source_verification", {})
            require(source.get("status") == "verified" and source.get("source_instance_id") == packet_record["source_instance_id"], "cell source verification changed")
            phase_b_score.validate_process_summary(source.get("process"), f"{blind_id} source verification")
            backends = task.get("backends")
            require(isinstance(backends, list) and [row.get("solver") for row in backends] == list(phase_b.BACKENDS), "cell backend order changed")
            expected = truth[blind_id]["target_class"]
            for backend_row in backends:
                backend = backend_row["solver"]
                require(backend_row.get("status") in phase_b_score.KNOWN_SOLVER_STATUSES, "unknown solver status")
                phase_b_score.validate_backend_admission(backend_row, blind_id, admission_task, task_root)
                outcome = phase_b_score.classification(expected, backend_row["status"])
                grouped[(cell_id, backend, expected)][outcome] += 1
                grouped[(cell_id, backend, expected)][f"status:{backend_row['status']}"] += 1
                for metrics in backend_metrics(backend_row):
                    backend_resources[backend]["total_core_seconds"] += metrics["total_core_seconds"]
                    backend_resources[backend]["summed_process_wall_seconds"] += metrics["wall_seconds"]
                    backend_resources[backend]["peak_rss_bytes"] = max(backend_resources[backend]["peak_rss_bytes"], metrics["peak_rss_bytes"])
                conflict = backend_row.get("conflicts")
                if isinstance(conflict, int) and conflict >= 0:
                    backend_resources[backend]["conflicts"].append(conflict)
                row = {
                    "blind_instance_id": blind_id,
                    "source_system_id": blind["source_system_id"],
                    "cell_id": cell_id,
                    "target_class": expected,
                    "backend": backend,
                    "solver_status": backend_row["status"],
                    "classification": outcome,
                    "conflicts": conflict,
                }
                rows.append(row)
                cell_rows.append(row)
        cells[cell_id] = {
            "verification": verification,
            "result_seal_sha256": phase_b.sha256_file(root / "result-seal.json", "cell result seal"),
            "outer_resources": result["outer_resources"],
            "process_resources": result["process_resources"],
            "classification_counts": dict(sorted(Counter(row["classification"] for row in cell_rows).items())),
        }
    require(len(rows) == 480, "scored row count is not 480")
    for backend, resource_row in backend_resources.items():
        values = sorted(resource_row.pop("conflicts"))
        resource_row.update({
            "total_core_seconds": round(resource_row["total_core_seconds"], 12),
            "summed_process_wall_seconds": round(resource_row["summed_process_wall_seconds"], 12),
            "conflicts_reported": len(values),
            "conflicts_sum": sum(values),
            "conflicts_min": min(values) if values else None,
            "conflicts_median": statistics.median(values) if values else None,
            "conflicts_max": max(values) if values else None,
            "conflict_values_sha256": phase_b.canonical_sha256(values),
        })
    tool_accounting = tool_costs(args.tools_root)
    workflow = workflow_accounting(args.workflow_metadata, set(cell_roots))
    direct_mitm = matched_direct_mitm()
    unknown_scalar = unknown_scalar_control()
    relation_yield = natural_relation_yield_control()
    outer_core = math.fsum(row["outer_resources"]["total_core_seconds"] for row in cells.values())
    outer_elapsed = math.fsum(row["outer_resources"]["single_core_elapsed_seconds"] for row in cells.values())
    result = {
        "schema": SCHEMA,
        "status": "complete_verified_truth_scored_four_cell_panel",
        "packet_verification": packet_check,
        "phase_a_seal_sha256": phase_b.sha256_file(args.phase_a_seal, "Phase-A seal"),
        "phase_a_oracle_sha256": raw_seal["oracle_ledger_sha256"],
        "workflow": workflow,
        "cells": cells,
        "instances": 160,
        "backend_rows": 480,
        "per_cell_backend_class": [
            {"cell_id": key[0], "backend": key[1], "target_class": key[2], "counts": dict(sorted(value.items()))}
            for key, value in sorted(grouped.items())
        ],
        "classification_counts": dict(sorted(Counter(row["classification"] for row in rows).items())),
        "backend_resources": backend_resources,
        "cell_outer_resources": {
            "summed_single_core_elapsed_seconds": outer_elapsed,
            "summed_total_core_seconds": outer_core,
            "maximum_sampled_process_tree_rss_bytes": max(row["outer_resources"]["sampled_peak_process_tree_rss_bytes"] for row in cells.values()),
            "parallel_cell_job_span_seconds": workflow["parallel_cell_job_span_seconds"],
        },
        "tool_accounting": tool_accounting,
        "matched_direct_mitm": direct_mitm,
        "unknown_scalar_index_calculus_and_rho": unknown_scalar,
        "natural_relation_yield": relation_yield,
        "charged_total_core_seconds_available": outer_core + direct_mitm["resources"]["summed_outer_total_core_seconds"] + tool_accounting["acquisition_total_core_seconds"] + tool_accounting["build_total_core_seconds"],
        "charged_scope": "tool source acquisition, exact tool builds, eight one-CPU packet/cell envelopes, 480 SAT-backend outcomes, and 160 matched direct-MITM outcomes",
        "rows": rows,
        "factor_base_algebraic_without_target_subgroup_enumeration": True,
        "factor_base_logs_known_by_construction": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_f4_same_instance_panel_complete": False,
        "full_end_to_end_index_calculus_cost_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "remaining_blockers": [
            "licensed Magma F4 has not executed the same 160-input packet",
            "the four-cell PDP panel is not a complete end-to-end index-calculus run",
            "factor-base discovery, relation collection, linear algebra, and rho are bound in separate controls rather than one same-instance full-cost experiment",
            "hosted package installation lacks CPU and memory receipts",
            "unaffiliated reproduction and novelty review remain absent",
        ],
        "completion_gate_audit": completion_gate_audit(),
    }
    args.output.mkdir(parents=False)
    phase_b.write_json_new(args.output / "score.json", result)
    inventory = phase_b.all_regular_inventory(args.output, {"score-seal.json"})
    payload = {
        "schema": SEAL_SCHEMA,
        "status": "score_frozen",
        "workflow_run_id": EXPECTED_RUN,
        "workflow_commit": EXPECTED_COMMIT,
        "inventory": inventory,
        "inventory_sha256": phase_b.canonical_sha256(inventory),
    }
    seal = {**payload, "seal_payload_sha256": phase_b.canonical_sha256(payload)}
    phase_b.write_json_new(args.output / "score-seal.json", seal)
    return seal


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--packet", type=Path, required=True)
    value.add_argument("--cell-root", type=Path, action="append", required=True)
    value.add_argument("--tools-root", type=Path, required=True)
    value.add_argument("--workflow-metadata", type=Path, required=True)
    value.add_argument("--phase-a-seal", type=Path, required=True)
    value.add_argument("--oracle-ledger", type=Path, required=True)
    value.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    value.add_argument("--output", type=Path, required=True)
    return value


def main() -> None:
    args = parser().parse_args()
    for name in ("packet", "tools_root", "workflow_metadata", "phase_a_seal", "oracle_ledger", "protocol"):
        setattr(args, name, getattr(args, name).resolve())
    args.cell_root = [root.resolve() for root in args.cell_root]
    args.output = args.output.resolve()
    try:
        print(json.dumps(score(args), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, Stage26ScoreError, phase_b.PhaseBError, cell_tool.Stage26Error, packet_tool.Stage26InputError) as error:
        raise SystemExit(f"stage26-score: {error}")


if __name__ == "__main__":
    main()
