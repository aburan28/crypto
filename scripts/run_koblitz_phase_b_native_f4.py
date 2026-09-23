#!/usr/bin/env python3
"""Run and score the repository's native F4 on the frozen Phase-B PDP packet.

The ``run`` subcommand receives only the authenticated blind bundle.  It
exports each selected public instance, verifies the historical solver-source
bytes, then runs ``koblitz_pdp_backend native-f4`` in a fresh, single-threaded
metered process.  The separate ``score`` subcommand may open the already
published Phase-B score only after the F4 run has an immutable inventory seal.

This adds a native F4 comparison arm.  It does not stand in for licensed Magma
F4, and it does not turn a PDP-stage measurement into an end-to-end index-
calculus result.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import resource
import time
from typing import Any

import build_koblitz_phase_b_native_f4 as tool_builder
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_pdp_matrix as matrix


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
DEFAULT_PROTOCOL = GATES / "stage-20-balanced-pdp-phase-b-protocol.json"
DEFAULT_BUNDLE = (
    GATES
    / "stage-20-phase-b-terminal-evidence-successor-04-20260910"
    / "inputs/blind-bundle.json"
)
DEFAULT_TRUTH = (
    GATES
    / "stage-20-phase-b-terminal-evidence-successor-04-20260910"
    / "score/score.json"
)
RUN_PLAN_SCHEMA = "koblitz_phase_b_native_f4_plan.v1"
TASK_SCHEMA = "koblitz_phase_b_native_f4_task.v1"
RUN_SUMMARY_SCHEMA = "koblitz_phase_b_native_f4_summary.v1"
RUN_SEAL_SCHEMA = "koblitz_phase_b_native_f4_run_seal.v1"
SCORE_SCHEMA = "koblitz_phase_b_native_f4_score.v1"
SCORE_SEAL_SCHEMA = "koblitz_phase_b_native_f4_score_seal.v1"
BACKEND = "native-f4"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise phase_b.PhaseBError(message)


def read_inputs(protocol_path: Path, bundle_path: Path) -> tuple[dict[str, Any], bytes, dict[str, Any], bytes]:
    protocol, protocol_bytes = phase_b.read_json(protocol_path, "Phase-B protocol")
    phase_b.validate_protocol(protocol)
    bundle, bundle_bytes = phase_b.read_json(bundle_path, "Phase-B blind bundle")
    phase_b.validate_blind_bundle(bundle, protocol)
    require(
        phase_b.sha256_bytes(bundle_bytes)
        == protocol["phase_a_binding"]["blind_bundle_sha256"],
        "blind bundle differs from the frozen Phase-A commitment",
    )
    return protocol, protocol_bytes, bundle, bundle_bytes


def select_instances(bundle: dict[str, Any], args: argparse.Namespace) -> list[dict[str, Any]]:
    instances = bundle["instances"]
    if args.blind_instance_id:
        wanted = set(args.blind_instance_id)
        selected = [row for row in instances if row["blind_instance_id"] in wanted]
        require(len(selected) == len(wanted), "one or more requested blind instance ids are absent")
    elif args.cell:
        selected = [row for row in instances if row["cell_id"] == args.cell]
        require(selected, f"cell {args.cell!r} is absent from the blind bundle")
    elif args.all:
        selected = list(instances)
    else:
        raise phase_b.PhaseBError("select --blind-instance-id, --cell, or --all")
    if args.max_instances is not None:
        require(args.max_instances > 0, "--max-instances must be positive")
        selected = selected[: args.max_instances]
    require(selected, "instance selection is empty")
    return selected


def process_resources(records: list[dict[str, Any]]) -> dict[str, Any]:
    metrics = [record["metrics"] for record in records]
    return {
        "process_receipts": len(records),
        "summed_process_wall_seconds": round(
            math.fsum(row["wall_seconds"] for row in metrics), 12
        ),
        "total_core_seconds": round(
            math.fsum(row["total_core_seconds"] for row in metrics), 12
        ),
        "single_core_seconds": round(
            math.fsum(row["single_core_seconds"] for row in metrics), 12
        ),
        "maximum_individual_process_rss_bytes": max(
            (row["peak_rss_bytes"] for row in metrics), default=0
        ),
        "single_core_semantics": (
            "sum of user plus system CPU from sequential processes whose child environment "
            "requests one thread; wall time is reported separately"
        ),
        "peak_rss_semantics": (
            "maximum fresh-process high-water RSS; processes run sequentially, so this is "
            "also the child-process panel peak"
        ),
    }


def compact_tool_identity(identity: dict[str, Any]) -> dict[str, Any]:
    return {key: identity[key] for key in ("path", "bytes", "sha256")}


def validate_f4_result(row: dict[str, Any], manifest: dict[str, Any]) -> None:
    require(row["solver"] == BACKEND, "native F4 parser returned the wrong backend")
    require(
        row["source_instance_id"] == manifest["source_instance"]["id_blake3"],
        "native F4 changed the authenticated source identity",
    )
    require(
        row["status"]
        in {
            "sat",
            "unsat",
            "unknown_inconclusive",
            "timeout_inconclusive",
            "not_run_resource_cap",
        },
        f"native F4 returned a rejected terminal {row['status']!r}",
    )
    report = row["backend_report"]
    if row["status"] == "sat":
        require(
            report.get("source_model_valid") is True
            and report.get("source_witness_valid") is True,
            "native F4 SAT lacks exact polynomial and curve-group validation",
        )
    if row["status"] == "unsat":
        require(
            report.get("exhaustive") is True,
            "native F4 UNSAT lacks a complete basis and root-extraction terminal",
        )
    require(report.get("conflicts") is None, "native F4 invented a SAT conflict count")


def run_one(
    *,
    ordinal: int,
    instance: dict[str, Any],
    output: Path,
    protocol: dict[str, Any],
    tools: dict[str, Path],
    identities: dict[str, dict[str, Any]],
    environment: dict[str, str],
    f4_budget_seconds: int,
    watchdog_seconds: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    prepared = phase_b.export_one_task(
        index=ordinal,
        instance=instance,
        output=output,
        protocol=protocol,
        tools=tools,
        identities=identities,
        environment=environment,
    )
    task_root = prepared["task_root"]
    instance_root = prepared["instance_root"]
    manifest_path = prepared["manifest_path"]
    manifest = prepared["manifest"]
    source_before = prepared["source_before"]
    markers = phase_b.forbidden_material(protocol)
    record = phase_b.run_metered(
        role=BACKEND,
        command=[
            str(tools["backend"].resolve()),
            BACKEND,
            str(manifest_path.resolve()),
            str(f4_budget_seconds),
        ],
        cwd=instance_root,
        task_root=task_root,
        input_paths=[
            manifest_path,
            *phase_b.manifest_source_paths(instance_root, manifest),
        ],
        timeout=watchdog_seconds,
        meter=tools["meter"],
        environment=environment,
        markers=markers,
        expected_executable_sha256=identities["backend"]["sha256"],
    )
    phase_b.require_source_unchanged(instance_root, manifest, source_before, BACKEND)
    row = matrix.isolated_backend_status(
        phase_b.process_record_for_matrix(record), BACKEND, manifest
    )
    validate_f4_result(row, manifest)
    task = {
        "schema": TASK_SCHEMA,
        "ordinal": ordinal,
        "blind_instance_id": instance["blind_instance_id"],
        "source_system_id": instance["source_system_id"],
        "cell_id": instance["cell_id"],
        "n": instance["n"],
        "ell": instance["ell"],
        "m": instance["m"],
        "basis": instance["basis"],
        "curve_a": instance["curve_a"],
        "factor_index": instance["factor_index"],
        "target": instance["target"],
        "source_instance_id": manifest["source_instance"]["id_blake3"],
        "source_artifacts_before": source_before,
        "source_artifacts_after": phase_b.frozen_source_snapshot(instance_root, manifest),
        "export": prepared["export_result"]["export"],
        "source_verification": prepared["export_result"]["source_verification"],
        "result": row,
        "truth_labels_present": False,
        "known_witnesses_present": False,
        "status": "terminal_native_f4_recorded",
    }
    require(
        task["source_artifacts_before"] == task["source_artifacts_after"],
        "native F4 changed a frozen source artifact",
    )
    phase_b.write_json_new(task_root / "task-result.json", task)
    return task, [
        prepared["export_result"]["export"],
        prepared["export_result"]["source_verification"]["process"],
        phase_b.compact_process(record),
    ]


def run_panel(args: argparse.Namespace) -> dict[str, Any]:
    protocol_path = args.protocol.resolve(strict=True)
    bundle_path = args.bundle.resolve(strict=True)
    protocol, protocol_bytes, bundle, bundle_bytes = read_inputs(protocol_path, bundle_path)
    selected = select_instances(bundle, args)
    require(not args.output.exists() and not args.output.is_symlink(), "run output must be new")
    require(args.f4_budget_seconds > 0, "--f4-budget-seconds must be positive")
    watchdog = args.watchdog_seconds or protocol["execution"]["per_process_watchdog_seconds"]
    require(
        watchdog > args.f4_budget_seconds,
        "the outer watchdog must exceed the internal F4 budget",
    )
    state = phase_b.git_state()
    require(not state["dirty"], "native F4 evidence requires a clean source checkout")
    phase_b.require_phase_a_ancestor(protocol["phase_a_binding"]["source_revision"])

    tools = {
        "exporter": args.exporter.resolve(strict=True),
        "backend": args.backend.resolve(strict=True),
        "meter": args.meter.resolve(strict=True),
    }
    identities = {
        name: phase_b.executable_identity(path, name, executable=name != "meter")
        for name, path in tools.items()
    }
    source_objects = tool_builder.source_objects(REPO, state["commit"])
    build = tool_builder.validate_receipt(
        args.rust_build_receipt.resolve(strict=True),
        {name: identities[name] for name in ("exporter", "backend")},
        source_objects,
        state,
    )
    environment = phase_b.safe_child_environment()
    require(environment.get("RAYON_NUM_THREADS") == "1", "F4 must be bound to one Rayon thread")
    plan = {
        "schema": RUN_PLAN_SCHEMA,
        "created_at": phase_b.now(),
        "claim_boundary": (
            "Native F4 same-instance PDP-stage evidence only. It is distinct from licensed "
            "Magma F4 and excludes relation collection, linear algebra and scalar recovery."
        ),
        "protocol_path": str(protocol_path),
        "protocol_sha256": phase_b.sha256_bytes(protocol_bytes),
        "blind_bundle_path": str(bundle_path),
        "blind_bundle_sha256": phase_b.sha256_bytes(bundle_bytes),
        "source_revision": state,
        "selected_instance_count": len(selected),
        "full_instance_count": len(bundle["instances"]),
        "selected_blind_instance_ids": [row["blind_instance_id"] for row in selected],
        "selection_order": "authenticated blind-bundle presentation order",
        "tool_identities": {
            name: compact_tool_identity(identity) for name, identity in identities.items()
        },
        "rust_build": build,
        "child_environment": environment,
        "backend": BACKEND,
        "solver": "f4-f2",
        "f4_budget_seconds": args.f4_budget_seconds,
        "watchdog_seconds": watchdog,
        "one_process_at_a_time": True,
        "single_thread_requested": True,
        "conflicts": None,
        "truth_labels_present": False,
        "known_witnesses_present": False,
    }
    if args.plan:
        return plan

    self_before = resource.getrusage(resource.RUSAGE_SELF)
    children_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    outer_started = time.monotonic()
    args.output.mkdir(parents=True)
    (args.output / "tasks").mkdir()
    phase_b.write_json_new(args.output / "execution-plan.json", plan)
    tool_builder.copy_capsule(
        args.rust_build_receipt.resolve(strict=True), args.output / "tool-builds" / "rust"
    )

    tasks: list[dict[str, Any]] = []
    processes: list[dict[str, Any]] = []
    for ordinal, instance in enumerate(selected):
        task, task_processes = run_one(
            ordinal=ordinal,
            instance=instance,
            output=args.output,
            protocol=protocol,
            tools=tools,
            identities=identities,
            environment=environment,
            f4_budget_seconds=args.f4_budget_seconds,
            watchdog_seconds=watchdog,
        )
        tasks.append(task)
        processes.extend(task_processes)

    outer_wall = time.monotonic() - outer_started
    self_after = resource.getrusage(resource.RUSAGE_SELF)
    children_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    outer_core = (
        self_after.ru_utime
        - self_before.ru_utime
        + self_after.ru_stime
        - self_before.ru_stime
        + children_after.ru_utime
        - children_before.ru_utime
        + children_after.ru_stime
        - children_before.ru_stime
    )
    statuses = Counter(task["result"]["status"] for task in tasks)
    f4_reports = [task["result"]["backend_report"] for task in tasks]
    f4_costs = [report.get("cost") for report in f4_reports if isinstance(report.get("cost"), dict)]
    f4_processes = [processes[index] for index in range(2, len(processes), 3)]
    summary = {
        "schema": RUN_SUMMARY_SCHEMA,
        "status": "complete_blind_native_f4_selection",
        "claim_boundary": plan["claim_boundary"],
        "selected_instances": len(tasks),
        "status_counts": dict(sorted(statuses.items())),
        "cells": dict(sorted(Counter(task["cell_id"] for task in tasks).items())),
        "backend": BACKEND,
        "solver": "f4-f2",
        "solver_operation_unit": "word XORs (elimination only)",
        "solver_word_xors": sum(cost["ops"] for cost in f4_costs),
        "solver_peak_internal_matrix_bytes": max(
            (cost["peak_bytes"] for cost in f4_costs), default=0
        ),
        "solver_conflicts": None,
        "solver_conflict_semantics": "F4 exposes algebraic and matrix counters, not SAT conflicts",
        "f4_process_resources": process_resources(f4_processes),
        "all_task_process_resources": process_resources(processes),
        "outer_resources": {
            "wall_seconds": outer_wall,
            "total_core_seconds": outer_core,
            "scope": (
                "Python custody driver plus all waited export, source-verification and native-F4 "
                "children; separately charged clean build is in tool-builds/rust"
            ),
        },
        "separately_charged_rust_build": build,
        "truth_labels_present": False,
        "known_witnesses_present": False,
        "licensed_magma_f4_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    phase_b.write_json_new(args.output / "run-summary.json", summary)
    inventory = phase_b.all_regular_inventory(args.output, {"run-seal.json"})
    seal = {
        "schema": RUN_SEAL_SCHEMA,
        "status": "native_f4_run_frozen_before_truth_scoring",
        "protocol_sha256": phase_b.sha256_bytes(protocol_bytes),
        "blind_bundle_sha256": phase_b.sha256_bytes(bundle_bytes),
        "selected_blind_instance_ids": plan["selected_blind_instance_ids"],
        "inventory": inventory,
        "inventory_sha256": phase_b.canonical_sha256(inventory),
    }
    seal["seal_payload_sha256"] = phase_b.canonical_sha256(seal)
    phase_b.write_json_new(args.output / "run-seal.json", seal)
    return summary


def validate_run_seal(run_root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    seal, _ = phase_b.read_json(run_root / "run-seal.json", "native F4 run seal")
    require(seal.get("schema") == RUN_SEAL_SCHEMA, "native F4 run seal schema changed")
    payload = dict(seal)
    payload_hash = payload.pop("seal_payload_sha256", None)
    require(
        phase_b.canonical_sha256(payload) == payload_hash,
        "native F4 run seal self-hash is invalid",
    )
    inventory = phase_b.all_regular_inventory(run_root, {"run-seal.json"})
    require(inventory == seal.get("inventory"), "native F4 run inventory changed")
    require(
        phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"),
        "native F4 run inventory hash changed",
    )
    plan, _ = phase_b.read_json(run_root / "execution-plan.json", "native F4 execution plan")
    require(plan.get("schema") == RUN_PLAN_SCHEMA, "native F4 execution plan schema changed")
    return seal, plan


def truth_map(score: dict[str, Any]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in score.get("rows", []):
        blind_id = row.get("blind_instance_id")
        truth = {
            "target_class": row.get("target_class"),
            "cell_id": row.get("cell_id"),
            "source_system_id": row.get("source_system_id"),
        }
        if blind_id in result:
            require(result[blind_id] == truth, f"Phase-B truth rows disagree for {blind_id}")
        else:
            result[blind_id] = truth
    return result


def score_panel(args: argparse.Namespace) -> dict[str, Any]:
    run_root = args.run.resolve(strict=True)
    seal, plan = validate_run_seal(run_root)
    truth, truth_bytes = phase_b.read_json(args.truth.resolve(strict=True), "Phase-B score")
    require(
        truth.get("schema") == "koblitz_pdp_phase_b_score.v1"
        and truth.get("full_panel_scored") is True,
        "truth source is not the completed frozen Phase-B score",
    )
    mapping = truth_map(truth)
    rows = []
    classifications = Counter()
    for ordinal, blind_id in enumerate(plan["selected_blind_instance_ids"]):
        task_path = phase_b.task_directory(run_root, ordinal, blind_id) / "task-result.json"
        task, _ = phase_b.read_json(task_path, "native F4 task result")
        require(task.get("schema") == TASK_SCHEMA, f"task schema changed for {blind_id}")
        require(task.get("blind_instance_id") == blind_id, f"task identity changed for {blind_id}")
        expected = mapping.get(blind_id)
        require(expected is not None, f"truth source lacks {blind_id}")
        require(
            task["cell_id"] == expected["cell_id"]
            and task["source_system_id"] == expected["source_system_id"],
            f"truth join changed source identity for {blind_id}",
        )
        status = task["result"]["status"]
        target_class = expected["target_class"]
        if status in {"unknown_inconclusive", "timeout_inconclusive", "not_run_resource_cap"}:
            classification = "inconclusive"
        elif status == "sat" and target_class == "decomposable":
            classification = "true_positive"
        elif status == "unsat" and target_class == "nondecomposable":
            classification = "true_negative"
        elif status == "sat":
            classification = "false_positive"
        else:
            classification = "false_negative"
        classifications[classification] += 1
        rows.append(
            {
                "blind_instance_id": blind_id,
                "source_system_id": task["source_system_id"],
                "cell_id": task["cell_id"],
                "target_class": target_class,
                "solver_status": status,
                "classification": classification,
                "conflicts": None,
            }
        )
    require(
        classifications["false_positive"] == 0 and classifications["false_negative"] == 0,
        "native F4 produced a false classification",
    )
    run_summary, _ = phase_b.read_json(run_root / "run-summary.json", "native F4 summary")
    score = {
        "schema": SCORE_SCHEMA,
        "status": "native_f4_selection_scored",
        "created_at": phase_b.now(),
        "run_inventory_sha256": seal["inventory_sha256"],
        "truth_source_path": str(args.truth.resolve()),
        "truth_source_sha256": phase_b.sha256_bytes(truth_bytes),
        "selected_instances": len(rows),
        "classification_counts": dict(sorted(classifications.items())),
        "rows": rows,
        "resources": run_summary,
        "claim_boundary": (
            "Internal native-F4 same-instance PDP-stage evidence. Licensed Magma F4, full "
            "index-calculus cost, rho crossover and independent reproduction remain open."
        ),
        "licensed_magma_f4_complete": False,
        "full_cost_gate_passed": False,
        "independent_external_reproduction_satisfied": False,
        "koblitz_index_calculus_sota": False,
    }
    require(not args.output.exists() and not args.output.is_symlink(), "score output must be new")
    args.output.mkdir(parents=True)
    phase_b.write_json_new(args.output / "score.json", score)
    inventory = phase_b.all_regular_inventory(args.output, {"score-seal.json"})
    score_seal = {
        "schema": SCORE_SEAL_SCHEMA,
        "status": "native_f4_score_frozen",
        "run_inventory_sha256": seal["inventory_sha256"],
        "inventory": inventory,
        "inventory_sha256": phase_b.canonical_sha256(inventory),
    }
    score_seal["seal_payload_sha256"] = phase_b.canonical_sha256(score_seal)
    phase_b.write_json_new(args.output / "score-seal.json", score_seal)
    return score


def parser() -> argparse.ArgumentParser:
    out = argparse.ArgumentParser(description=__doc__)
    sub = out.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="run native F4 without truth labels")
    run.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    run.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--exporter", type=Path, required=True)
    run.add_argument("--backend", type=Path, required=True)
    run.add_argument("--rust-build-receipt", type=Path, required=True)
    run.add_argument("--meter", type=Path, default=phase_b.DEFAULT_METER)
    selection = run.add_mutually_exclusive_group(required=True)
    selection.add_argument("--blind-instance-id", action="append")
    selection.add_argument("--cell")
    selection.add_argument("--all", action="store_true")
    run.add_argument("--max-instances", type=int)
    run.add_argument("--f4-budget-seconds", type=int, default=115)
    run.add_argument("--watchdog-seconds", type=int)
    run.add_argument("--plan", action="store_true")
    score = sub.add_parser("score", help="score an already sealed native-F4 run")
    score.add_argument("--run", type=Path, required=True)
    score.add_argument("--truth", type=Path, default=DEFAULT_TRUTH)
    score.add_argument("--output", type=Path, required=True)
    return out


def main() -> None:
    args = parser().parse_args()
    try:
        result = run_panel(args) if args.command == "run" else score_panel(args)
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, phase_b.PhaseBError) as error:
        parser().exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
