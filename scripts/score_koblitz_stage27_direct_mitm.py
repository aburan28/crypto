#!/usr/bin/env python3
"""Verify and truth-score the four-cell Stage-27 direct-MITM matrix."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
import json
import math
from pathlib import Path
from typing import Any

import koblitz_stage26_affinity_inputs as packet_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_stage27_direct_mitm as stage27
import score_koblitz_blind_pdp_phase_b as phase_b_score


SCHEMA = "koblitz_stage27_direct_mitm_score.v1"
SEAL_SCHEMA = "koblitz_stage27_direct_mitm_score_seal.v1"
RUN_ID = 34633920325
RUN_COMMIT = "02225ed61e7331c11a1e4665b953a8384b0b08de"


class Stage27ScoreError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage27ScoreError(message)


def read_json(path: Path, context: str) -> dict[str, Any]:
    value, _ = phase_b.read_json(path, context)
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def timestamp(value: str) -> datetime:
    require(isinstance(value, str) and value.endswith("Z"), "workflow timestamp is invalid")
    return datetime.fromisoformat(value[:-1] + "+00:00")


def workflow_receipt(path: Path, cells: set[str]) -> dict[str, Any]:
    value = read_json(path, "Stage-27 workflow metadata")
    require(value.get("databaseId") == RUN_ID, "Stage-27 workflow ID changed")
    require(value.get("headSha") == RUN_COMMIT, "Stage-27 workflow commit changed")
    require(value.get("status") == "completed" and value.get("conclusion") == "success", "Stage-27 workflow is incomplete")
    jobs = value.get("jobs")
    require(isinstance(jobs, list), "Stage-27 workflow lacks jobs")
    cell_jobs = {}
    for job in jobs:
        name = job.get("name", "")
        if name.startswith("direct-mitm-cell (") and name.endswith(")"):
            cell_jobs[name[len("direct-mitm-cell ("):-1]] = job
    require(set(cell_jobs) == cells, "Stage-27 workflow cell inventory changed")
    for cell, job in cell_jobs.items():
        require(job.get("status") == "completed" and job.get("conclusion") == "success", f"Stage-27 workflow cell {cell} failed")
    return {
        "run_id": RUN_ID,
        "commit": RUN_COMMIT,
        "url": value.get("url"),
        "workflow_wall_seconds": (timestamp(value["updatedAt"]) - timestamp(value["createdAt"])).total_seconds(),
        "parallel_cell_job_span_seconds": (
            max(timestamp(job["completedAt"]) for job in cell_jobs.values())
            - min(timestamp(job["startedAt"]) for job in cell_jobs.values())
        ).total_seconds(),
        "cell_jobs": {
            cell: {
                "job_id": job.get("databaseId"),
                "started_at": job["startedAt"],
                "completed_at": job["completedAt"],
                "wall_seconds": (timestamp(job["completedAt"]) - timestamp(job["startedAt"])).total_seconds(),
                "url": job.get("url"),
            }
            for cell, job in sorted(cell_jobs.items())
        },
    }


def artifact_receipt(path: Path, cells: set[str]) -> dict[str, Any]:
    value = read_json(path, "Stage-27 artifact metadata")
    artifacts = value.get("artifacts")
    require(value.get("total_count") == 4 and isinstance(artifacts, list) and len(artifacts) == 4, "Stage-27 artifact count changed")
    expected_names = {f"koblitz-stage27-{cell}-{RUN_ID}" for cell in cells}
    require({row.get("name") for row in artifacts} == expected_names, "Stage-27 artifact names changed")
    records = {}
    for artifact in artifacts:
        name = artifact["name"]
        require(artifact.get("expired") is False, f"Stage-27 artifact expired: {name}")
        require(isinstance(artifact.get("id"), int) and artifact["id"] > 0, f"Stage-27 artifact ID invalid: {name}")
        phase_b.require_hex64(artifact.get("digest", "").removeprefix("sha256:"), f"Stage-27 artifact digest {name}")
        require(isinstance(artifact.get("size_in_bytes"), int) and artifact["size_in_bytes"] > 0, f"Stage-27 artifact size invalid: {name}")
        cell = name[len("koblitz-stage27-") : -len(f"-{RUN_ID}")]
        records[cell] = {
            "id": artifact["id"],
            "name": name,
            "size_in_bytes": artifact["size_in_bytes"],
            "digest": artifact["digest"],
            "expires_at": artifact["expires_at"],
        }
    require(set(records) == cells, "Stage-27 artifact-to-cell mapping changed")
    return records


def validate_source(task_root: Path, packet_record: dict[str, Any]) -> None:
    expected = {Path(row["path"]).name: row for row in packet_record["files"]}
    instance = task_root / "instance"
    require({path.name for path in instance.iterdir()} == set(expected), "Stage-27 task source inventory changed")
    for name, identity in expected.items():
        data = phase_b.regular_file_bytes(instance / name, "Stage-27 task source")
        require(len(data) == identity["bytes"] and phase_b.sha256_bytes(data) == identity["sha256"], "Stage-27 task source differs from packet")


def score(args: argparse.Namespace) -> dict[str, Any]:
    require(not args.output.exists() and not args.output.is_symlink(), "score output must be new")
    packet_check = packet_tool.verify(args.packet)
    packet = read_json(args.packet / "packet-manifest.json", "Stage-27 packet manifest")
    blind = read_json(args.packet / "source/blind-bundle.json", "Stage-27 blind bundle")
    protocol, _ = phase_b.read_json(stage27.PROTOCOL, "Stage-27 Phase-B protocol")
    phase_b.validate_protocol(protocol)
    raw_seal, oracle = phase_b_score.load_authenticated_oracle(args.phase_a_seal, args.oracle_ledger, protocol)
    truth = phase_b_score.oracle_index(oracle, blind, protocol)
    packet_by_id = {row["blind_instance_id"]: row for row in packet["instances"]}
    cell_roots = {}
    for root in args.cell_root:
        result = read_json(root / "result.json", "Stage-27 cell result")
        cell = result.get("cell_id")
        require(cell in packet["cells"] and cell not in cell_roots, "duplicate or unexpected Stage-27 cell root")
        cell_roots[cell] = root
    require(set(cell_roots) == set(packet["cells"]), "Stage-27 four-cell result is incomplete")
    workflow = workflow_receipt(args.workflow_metadata, set(cell_roots))
    artifacts = artifact_receipt(args.artifact_metadata, set(cell_roots))
    classifications = Counter()
    cells = {}
    rows = []
    for cell_id, root in sorted(cell_roots.items()):
        verification = stage27.verify(root)
        result = read_json(root / "result.json", "Stage-27 cell result")
        require(result.get("cell_id") == cell_id and result.get("instances") == 40 and result.get("outcomes") == 40, "Stage-27 cell counts changed")
        require(result.get("status_counts") == {"sat": 20, "unsat": 20}, "Stage-27 cell did not return the balanced terminal panel")
        cell_classification = Counter()
        additions = 0
        pairs = 0
        for ordinal, blind_id in enumerate(packet["cells"][cell_id]):
            task_root = root / "tasks" / f"{ordinal:06d}-{blind_id}"
            task = read_json(task_root / "task-result.json", "Stage-27 task result")
            record = packet_by_id[blind_id]
            require(task.get("schema") == "koblitz_stage27_direct_mitm_task.v1", "Stage-27 task schema changed")
            require(task.get("ordinal") == ordinal and task.get("blind_instance_id") == blind_id, "Stage-27 task identity changed")
            require(task.get("cell_id") == cell_id and task.get("target") == record["target"], "Stage-27 task target changed")
            require(task.get("source_instance_id") == record["source_instance_id"], "Stage-27 task source ID changed")
            require(task.get("status") == "terminal_direct_mitm_recorded" and task.get("source_before") == task.get("source_after"), "Stage-27 task source custody failed")
            validate_source(task_root, record)
            row = task.get("result")
            require(isinstance(row, dict), "Stage-27 task lacks its result")
            stage27.validate_terminal(row, blind_id)
            phase_b_score.validate_process_summary(row.get("process"), f"{blind_id} direct MITM")
            report = row["backend_report"]
            require(report.get("source_instance_id") == record["source_instance_id"], "Stage-27 backend source ID changed")
            expected = truth[blind_id]["target_class"]
            classification = phase_b_score.classification(expected, row["status"])
            require(classification in {"true_positive", "true_negative"}, "Stage-27 direct MITM produced a false or inconclusive outcome")
            classifications[classification] += 1
            cell_classification[classification] += 1
            additions += report["group_additions"]
            pairs += report["pair_entries"]
            rows.append({
                "blind_instance_id": blind_id,
                "cell_id": cell_id,
                "target_class": expected,
                "solver_status": row["status"],
                "classification": classification,
                "factor_points": report["factor_points"],
                "pair_entries": report["pair_entries"],
                "group_additions": report["group_additions"],
            })
        require(cell_classification == Counter({"true_positive": 20, "true_negative": 20}), "Stage-27 cell truth scoring changed")
        cells[cell_id] = {
            "artifact": artifacts[cell_id],
            "verification": verification,
            "result_seal_sha256": phase_b.sha256_file(root / "result-seal.json", "Stage-27 cell seal"),
            "classifications": dict(sorted(cell_classification.items())),
            "group_additions": additions,
            "pair_entries": pairs,
            "factor_points": result["factor_points"],
            "process_resources": result["process_resources"],
            "outer_resources": result["outer_resources"],
        }
    require(len(rows) == 160 and classifications == Counter({"true_positive": 80, "true_negative": 80}), "Stage-27 global truth scoring changed")
    outer_core = math.fsum(row["outer_resources"]["total_core_seconds"] for row in cells.values())
    outer_elapsed = math.fsum(row["outer_resources"]["single_core_elapsed_seconds"] for row in cells.values())
    process_core = math.fsum(row["process_resources"]["total_core_seconds"] for row in cells.values())
    result = {
        "schema": SCHEMA,
        "status": "complete_verified_truth_scored_direct_mitm_panel",
        "packet_verification": packet_check,
        "phase_a_seal_sha256": phase_b.sha256_file(args.phase_a_seal, "Stage-27 Phase-A seal"),
        "phase_a_oracle_sha256": raw_seal["oracle_ledger_sha256"],
        "workflow": workflow,
        "cells": cells,
        "instances": 160,
        "outcomes": 160,
        "classification_counts": dict(sorted(classifications.items())),
        "resources": {
            "summed_single_core_elapsed_seconds": outer_elapsed,
            "summed_outer_total_core_seconds": outer_core,
            "summed_direct_mitm_process_core_seconds": process_core,
            "maximum_sampled_process_tree_rss_bytes": max(row["outer_resources"]["sampled_peak_process_tree_rss_bytes"] for row in cells.values()),
            "maximum_individual_direct_mitm_process_rss_bytes": max(row["process_resources"]["maximum_individual_process_rss_bytes"] for row in cells.values()),
            "parallel_cell_job_span_seconds": workflow["parallel_cell_job_span_seconds"],
            "workflow_wall_seconds": workflow["workflow_wall_seconds"],
        },
        "operations": {
            "group_additions": sum(row["group_additions"] for row in rows),
            "pair_entries": sum(row["pair_entries"] for row in rows),
        },
        "conflicts": None,
        "conflict_semantics": "direct MITM exposes group additions and pair-table entries rather than SAT conflicts",
        "rows": rows,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "claim_boundary": "matched public synthetic PDP direct-decomposition control only",
    }
    args.output.mkdir(parents=False)
    phase_b.write_json_new(args.output / "score.json", result)
    inventory = phase_b.all_regular_inventory(args.output, {"score-seal.json"})
    payload = {
        "schema": SEAL_SCHEMA,
        "status": "direct_mitm_score_frozen",
        "workflow_run_id": RUN_ID,
        "workflow_commit": RUN_COMMIT,
        "inventory": inventory,
        "inventory_sha256": phase_b.canonical_sha256(inventory),
    }
    seal = {**payload, "seal_payload_sha256": phase_b.canonical_sha256(payload)}
    phase_b.write_json_new(args.output / "score-seal.json", seal)
    return seal


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet", type=Path, required=True)
    parser.add_argument("--cell-root", type=Path, action="append", required=True)
    parser.add_argument("--phase-a-seal", type=Path, required=True)
    parser.add_argument("--oracle-ledger", type=Path, required=True)
    parser.add_argument("--workflow-metadata", type=Path, required=True)
    parser.add_argument("--artifact-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for name in ("packet", "phase_a_seal", "oracle_ledger", "workflow_metadata", "artifact_metadata", "output"):
        setattr(args, name, getattr(args, name).resolve())
    args.cell_root = [root.resolve() for root in args.cell_root]
    try:
        print(json.dumps(score(args), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, Stage27ScoreError, stage27.Stage27Error, phase_b.PhaseBError, packet_tool.Stage26InputError) as error:
        raise SystemExit(f"stage27-score: {error}")


if __name__ == "__main__":
    main()
