#!/usr/bin/env python3
"""Run one exact Stage-26 cell through direct MITM on one Linux CPU."""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
from typing import Any

import koblitz_stage26_affinity_inputs as packet_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_pdp_matrix as matrix
import run_koblitz_stage26_affinity_cell as stage26


REPO = Path(__file__).resolve().parents[1]
PROTOCOL = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-balanced-pdp-phase-b-protocol.json"
SCHEMA = "koblitz_stage27_direct_mitm_cell_result.v1"
SEAL_SCHEMA = "koblitz_stage27_direct_mitm_cell_seal.v1"
PACKET_SHA256 = "c937afb9b172d114768b0b96a4b5ccf66e91fbf278b77c58ed49f0fd74af37f7"
BACKEND_SHA256 = "667fc00681e114acaaf662a8c00c214354e1a08e236a6974fb08e52504711ea0"
STAGE26_TOOL_RUN = 34632018379


class Stage27Error(RuntimeError):
    pass


def read_json(path: Path, context: str) -> dict[str, Any]:
    value, _ = phase_b.read_json(path, context)
    if not isinstance(value, dict):
        raise Stage27Error(f"{context} must be an object")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage27Error(message)


def validate_terminal(row: dict[str, Any], blind_id: str) -> None:
    require(row.get("solver") == "direct-mitm", f"{blind_id} used the wrong backend")
    require(row.get("status") in {"sat", "unsat"}, f"{blind_id} direct MITM is not terminal")
    process = row.get("process")
    require(isinstance(process, dict), f"{blind_id} lacks its MITM process")
    require(process.get("returncode") == 0 and process.get("timed_out") is False, f"{blind_id} MITM did not complete cleanly")
    report = row.get("backend_report")
    require(isinstance(report, dict), f"{blind_id} lacks its backend report")
    require(
        report.get("schema") == "koblitz_pdp_isolated_backend.v1"
        and report.get("backend") == "direct-mitm"
        and report.get("source_instance_verified") is True
        and report.get("regenerated_source_exact") is True
        and report.get("exhaustive") is True,
        f"{blind_id} MITM report violates the source or exhaustive-search contract",
    )
    require(isinstance(report.get("factor_points"), int) and report["factor_points"] > 0, f"{blind_id} factor-point count is invalid")
    require(isinstance(report.get("pair_entries"), int) and report["pair_entries"] > 0, f"{blind_id} pair-table size is invalid")
    require(isinstance(report.get("group_additions"), int) and report["group_additions"] > 0, f"{blind_id} group-addition count is invalid")
    if row["status"] == "sat":
        require(row.get("source_witness_valid") is True, f"{blind_id} SAT lacks its exact point witness")
        witness = report.get("witness_indices")
        require(isinstance(witness, list) and len(witness) == 3 and all(type(value) is int and value >= 0 for value in witness), f"{blind_id} SAT witness indices are invalid")
    else:
        require(row.get("source_witness_valid") is None and report.get("witness_indices") is None, f"{blind_id} UNSAT carries a witness")


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(not args.output.exists() and not args.output.is_symlink(), "output must be new")
    affinity = stage26.singleton_affinity(args.cpu)
    self_before = resource.getrusage(resource.RUSAGE_SELF)
    children_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    sampler = stage26.TreeSampler()
    sampler.__enter__()
    packet_check = packet_tool.verify(args.packet)
    require(packet_check["packet_inventory_sha256"] == PACKET_SHA256, "packet identity changed")
    packet = read_json(args.packet / "packet-manifest.json", "packet manifest")
    require(args.cell in packet["cells"] and len(packet["cells"][args.cell]) == 40, "cell is missing or incomplete")
    backend = args.backend.resolve(strict=True)
    require(phase_b.sha256_file(backend, "direct MITM backend") == BACKEND_SHA256, "backend differs from the frozen Stage-26 tool")
    protocol, _ = phase_b.read_json(PROTOCOL, "Phase-B protocol")
    phase_b.validate_protocol(protocol)
    environment = phase_b.safe_child_environment()
    markers = phase_b.forbidden_material(protocol)
    identity = phase_b.executable_identity(backend, "direct MITM backend")
    records = {row["blind_instance_id"]: row for row in packet["instances"]}
    args.output.mkdir(parents=True)
    (args.output / "tasks").mkdir()
    phase_b.write_json_new(args.output / "affinity.json", affinity)
    plan = {
        "schema": "koblitz_stage27_direct_mitm_plan.v1",
        "cell_id": args.cell,
        "instance_count": 40,
        "packet_verification": packet_check,
        "backend_identity": identity,
        "backend_origin": {
            "stage26_workflow_run": STAGE26_TOOL_RUN,
            "artifact_name": f"koblitz-stage26-tools-{STAGE26_TOOL_RUN}",
            "artifact_digest": "sha256:d1628f8464942aaeb3af95f1d118140a770cf923ae6ad502d8541201fd3c291f",
            "build_cost_charged_in_stage26": True,
        },
        "affinity": affinity,
        "truth_scoring_status": "withheld_until_all_four_direct_mitm_cells_are_frozen",
    }
    phase_b.write_json_new(args.output / "execution-plan.json", plan)
    results = []
    timeout = protocol["execution"]["per_process_watchdog_seconds"]
    for ordinal, blind_id in enumerate(packet["cells"][args.cell]):
        record = records[blind_id]
        task_root = args.output / "tasks" / f"{ordinal:06d}-{blind_id}"
        task_root.mkdir()
        instance_root = task_root / "instance"
        stage26.copy_instance(args.packet, record, instance_root)
        manifest_path = instance_root / "manifest.json"
        manifest = read_json(manifest_path, "direct MITM manifest")
        require(manifest.get("blind_instance_id") == blind_id, "task manifest blind ID changed")
        before = phase_b.frozen_source_snapshot(instance_root, manifest)
        process = phase_b.run_metered(
            role="direct-mitm",
            command=[str(backend), "direct-mitm", str(manifest_path)],
            cwd=instance_root,
            task_root=task_root,
            input_paths=[manifest_path, *phase_b.manifest_source_paths(instance_root, manifest)],
            timeout=timeout,
            meter=args.meter.resolve(strict=True),
            environment=environment,
            markers=markers,
            expected_executable_sha256=BACKEND_SHA256,
        )
        phase_b.require_source_unchanged(instance_root, manifest, before, "direct MITM")
        row = matrix.isolated_backend_status(phase_b.process_record_for_matrix(process), "direct-mitm", manifest)
        row["process"] = phase_b.compact_process(process)
        validate_terminal(row, blind_id)
        task = {
            "schema": "koblitz_stage27_direct_mitm_task.v1",
            "ordinal": ordinal,
            "blind_instance_id": blind_id,
            "cell_id": args.cell,
            "source_instance_id": record["source_instance_id"],
            "target": record["target"],
            "source_before": before,
            "source_after": phase_b.frozen_source_snapshot(instance_root, manifest),
            "result": row,
            "status": "terminal_direct_mitm_recorded",
        }
        require(task["source_before"] == task["source_after"], "direct MITM changed source bytes")
        phase_b.write_json_new(task_root / "task-result.json", task)
        results.append(task)
    sampler.__exit__(None, None, None)
    elapsed = time.monotonic() - started
    self_after = resource.getrusage(resource.RUSAGE_SELF)
    children_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    total_core = (
        self_after.ru_utime - self_before.ru_utime + self_after.ru_stime - self_before.ru_stime
        + children_after.ru_utime - children_before.ru_utime
        + children_after.ru_stime - children_before.ru_stime
    )
    metrics = [task["result"]["process"]["metrics"] for task in results]
    reports = [task["result"]["backend_report"] for task in results]
    statuses = Counter(task["result"]["status"] for task in results)
    result = {
        "schema": SCHEMA,
        "status": "complete_truth_free_single_cpu_direct_mitm_cell",
        "cell_id": args.cell,
        "instances": 40,
        "outcomes": 40,
        "status_counts": dict(sorted(statuses.items())),
        "group_additions": {
            "sum": sum(row["group_additions"] for row in reports),
            "min": min(row["group_additions"] for row in reports),
            "median": statistics_median(row["group_additions"] for row in reports),
            "max": max(row["group_additions"] for row in reports),
        },
        "pair_entries": {
            "sum": sum(row["pair_entries"] for row in reports),
            "min": min(row["pair_entries"] for row in reports),
            "median": statistics_median(row["pair_entries"] for row in reports),
            "max": max(row["pair_entries"] for row in reports),
        },
        "factor_points": sorted(set(row["factor_points"] for row in reports)),
        "process_resources": {
            "process_receipts": len(metrics),
            "total_core_seconds": round(math.fsum(row["total_core_seconds"] for row in metrics), 12),
            "summed_process_wall_seconds": round(math.fsum(row["wall_seconds"] for row in metrics), 12),
            "maximum_individual_process_rss_bytes": max(row["peak_rss_bytes"] for row in metrics),
        },
        "outer_resources": {
            "single_core_elapsed_seconds": elapsed,
            "total_core_seconds": total_core,
            "one_cpu_utilization_percent": 100.0 * total_core / elapsed,
            "sampled_peak_process_tree_rss_bytes": sampler.peak,
            "process_tree_rss_samples": sampler.samples,
            "sample_interval_seconds": 0.02,
            "scope": "packet integrity verification, setup, exact source verification, factor-point materialisation, exhaustive pair table, and target search",
        },
        "conflicts": None,
        "conflict_semantics": "direct MITM exposes group additions and pair entries rather than SAT conflicts",
        "truth_labels_present": False,
        "known_witnesses_present": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    phase_b.write_json_new(args.output / "result.json", result)
    inventory = phase_b.all_regular_inventory(args.output, {"result-seal.json"})
    payload = {
        "schema": SEAL_SCHEMA,
        "status": "direct_mitm_cell_frozen",
        "cell_id": args.cell,
        "packet_inventory_sha256": PACKET_SHA256,
        "inventory": inventory,
        "inventory_sha256": phase_b.canonical_sha256(inventory),
    }
    seal = {**payload, "seal_payload_sha256": phase_b.canonical_sha256(payload)}
    phase_b.write_json_new(args.output / "result-seal.json", seal)
    return seal


def statistics_median(values: Any) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    return float(ordered[middle]) if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2


def verify(output: Path) -> dict[str, Any]:
    seal = read_json(output / "result-seal.json", "direct MITM seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == phase_b.canonical_sha256(payload), "direct MITM seal is invalid")
    inventory = phase_b.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal.get("inventory") and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"), "direct MITM inventory changed")
    result = read_json(output / "result.json", "direct MITM result")
    require(result.get("schema") == SCHEMA and result.get("status") == "complete_truth_free_single_cpu_direct_mitm_cell", "direct MITM result is incomplete")
    require(result.get("instances") == 40 and result.get("outcomes") == 40, "direct MITM result counts changed")
    require(sum(result.get("status_counts", {}).values()) == 40, "direct MITM terminal counts changed")
    require(result.get("truth_labels_present") is False and result.get("known_witnesses_present") is False, "direct MITM execution was not truth-free")
    require(result.get("koblitz_index_calculus_sota") is False, "direct MITM result overclaims SOTA")
    return {"schema": "koblitz_stage27_direct_mitm_verification.v1", "status": "verified", "cell_id": result["cell_id"], "inventory_sha256": seal["inventory_sha256"]}


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest="command", required=True)
    run_p = sub.add_parser("run")
    run_p.add_argument("--packet", type=Path, required=True)
    run_p.add_argument("--cell", required=True)
    run_p.add_argument("--output", type=Path, required=True)
    run_p.add_argument("--backend", type=Path, required=True)
    run_p.add_argument("--meter", type=Path, default=phase_b.DEFAULT_METER)
    run_p.add_argument("--cpu", type=int)
    verify_p = sub.add_parser("verify")
    verify_p.add_argument("--output", type=Path, required=True)
    return root


def main() -> None:
    args = parser().parse_args()
    try:
        if args.command == "run":
            args.packet, args.output = args.packet.resolve(), args.output.resolve()
            value = run(args)
        else:
            value = verify(args.output.resolve())
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, Stage27Error, phase_b.PhaseBError, packet_tool.Stage26InputError) as error:
        raise SystemExit(f"stage27-direct-mitm: {error}")


if __name__ == "__main__":
    main()
