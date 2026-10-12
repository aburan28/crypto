#!/usr/bin/env python3
"""Compose the sealed Stage-159 single-target native-F4 result."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_phase_b_native_f4 as native_f4
import run_koblitz_stage26_affinity_cell as stage26
import run_koblitz_stage27_direct_mitm as stage27


REPO = Path(__file__).resolve().parents[1]
GATES = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
STAGE = GATES / "stage-159-native-f4-single-target-20260922"
BLIND_ID = "b-23f77043e135d4b3b9d8ea41443fe2ba16b23932b1f37b7f78499ce6fca75f30"
CELL = "n59-l9-m3-standard-a1-f0"
SOURCE_ID = "fb0ea6bb30be3dea841c609a2668939128b4b799d4c0e051b0b52dc6da05f0fa"
RESULT_SCHEMA = "koblitz_stage159_native_f4_single_target.v1"
SEAL_SCHEMA = "koblitz_stage159_native_f4_single_target_seal.v1"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise phase_b.PhaseBError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value, _ = phase_b.read_json(path, context)
    require(isinstance(value, dict), f"{context} must be an object")
    return value


def only_task(root: Path) -> dict[str, Any]:
    paths = sorted((root / "tasks").glob("*/task-result.json"))
    require(len(paths) == 1, f"{root.name} must contain one task")
    task = load(paths[0], f"{root.name} task")
    require(
        task.get("blind_instance_id") == BLIND_ID
        and task.get("cell_id") == CELL
        and task.get("source_instance_id") == SOURCE_ID,
        f"{root.name} task identity changed",
    )
    return task


def verify_run(name: str) -> tuple[dict[str, Any], dict[str, Any]]:
    root = STAGE / name
    seal, plan = native_f4.validate_run_seal(root)
    require(plan["selected_blind_instance_ids"] == [BLIND_ID], f"{name} selection changed")
    return only_task(root), seal


def verify_score(name: str, expected: str) -> dict[str, Any]:
    root = STAGE / name
    seal = load(root / "score-seal.json", f"{name} score seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == native_f4.SCORE_SEAL_SCHEMA
        and claimed == phase_b.canonical_sha256(payload),
        f"{name} score seal changed",
    )
    inventory = phase_b.all_regular_inventory(root, {"score-seal.json"})
    require(
        inventory == seal.get("inventory")
        and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"),
        f"{name} score inventory changed",
    )
    score = load(root / "score.json", f"{name} score")
    require(
        score.get("classification_counts") == {expected: 1}
        and score.get("selected_instances") == 1,
        f"{name} classification changed",
    )
    return score


def verify_failed_attempt(name: str) -> dict[str, Any]:
    root = STAGE / name
    seal = load(root / "attempt-seal.json", f"{name} seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == "koblitz_phase_b_native_f4_failed_attempt_seal.v1"
        and claimed == phase_b.canonical_sha256(payload),
        f"{name} seal changed",
    )
    inventory = phase_b.all_regular_inventory(root, {"attempt-seal.json"})
    require(inventory == seal.get("inventory"), f"{name} inventory changed")
    summary = load(root / "failure-summary.json", f"{name} failure summary")
    require(
        summary.get("status") == "watchdog_timeout_inconclusive"
        and summary.get("classification") == "inconclusive"
        and summary.get("blind_instance_id") == BLIND_ID,
        f"{name} failure boundary changed",
    )
    return summary


def compact_metrics(row: dict[str, Any]) -> dict[str, Any]:
    return {
        key: row.get(key)
        for key in (
            "wall_seconds",
            "single_core_seconds",
            "total_core_seconds",
            "peak_rss_bytes",
            "meter",
        )
    }


def ratio(candidate: float, reference: float) -> float:
    return candidate / reference


def metric_records(root: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(root.rglob("*.metrics.json")):
        value = load(path, str(path))
        metrics = value.get("metrics")
        if isinstance(metrics, dict) and "total_core_seconds" in metrics:
            rows.append(metrics)
    return rows


def add_resources(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "resource_components": len(rows),
        "summed_wall_seconds": round(math.fsum(row["wall_seconds"] for row in rows), 12),
        "total_core_seconds": round(math.fsum(row["total_core_seconds"] for row in rows), 12),
        "peak_rss_bytes": max((row["peak_rss_bytes"] for row in rows), default=0),
    }


def run_build_receipt(run_name: str) -> dict[str, Any]:
    return load(STAGE / run_name / "tool-builds/rust/receipt.json", f"{run_name} build")


def compose(stage26_root: Path, stage27_root: Path) -> dict[str, Any]:
    failed_1 = verify_failed_attempt("attempt-01-watchdog-timeout")
    failed_2 = verify_failed_attempt("attempt-02-budgeted-terminal")
    task3, seal3 = verify_run("attempt-03-fixed-x1")
    task4, seal4 = verify_run("attempt-04-fixed-x1-long")
    task5, seal5 = verify_run("attempt-05-rational-x1")
    task6, seal6 = verify_run("attempt-06-trace")
    verify_score("attempt-03-score", "inconclusive")
    verify_score("attempt-04-score", "true_positive")
    verify_score("attempt-05-score", "true_positive")
    verify_score("attempt-06-score", "true_positive")

    verification26 = stage26.verify(stage26_root)
    verification27 = stage27.verify(stage27_root)
    require(verification26["status"] == "verified", "Stage-26 control failed verification")
    require(verification27["status"] == "verified", "Stage-27 control failed verification")
    control26_paths = sorted(stage26_root.glob(f"tasks/*{BLIND_ID}/task-result.json"))
    control27_paths = sorted(stage27_root.glob(f"tasks/*{BLIND_ID}/task-result.json"))
    require(len(control26_paths) == len(control27_paths) == 1, "same-target control task missing")
    control26 = load(control26_paths[0], "Stage-26 same-target task")
    control27 = load(control27_paths[0], "Stage-27 same-target task")
    require(
        control26.get("blind_instance_id") == control27.get("blind_instance_id") == BLIND_ID
        and control26.get("source_instance_id") == control27.get("source_instance_id") == SOURCE_ID,
        "control source identity changed",
    )
    controls = STAGE / "controls"
    require(not controls.exists(), "Stage-159 controls output already exists")
    controls.mkdir()
    phase_b.write_json_new(controls / "stage26-task-result.json", control26)
    phase_b.write_json_new(controls / "stage27-task-result.json", control27)

    stage26_score = load(GATES / "stage-26-affinity-matrix-result-20260911/score.json", "Stage-26 score")
    stage27_score = load(GATES / "stage-27-direct-mitm-result-20260911/score.json", "Stage-27 score")
    artifact26 = stage26_score["artifacts"]["koblitz-stage26-n59-l9-m3-standard-a1-f0-34632018379"]
    artifact27 = stage27_score["cells"][CELL]["artifact"]
    provenance = {
        "schema": "koblitz_stage159_same_target_control_provenance.v1",
        "blind_instance_id": BLIND_ID,
        "source_instance_id": SOURCE_ID,
        "stage26": {
            "artifact": artifact26,
            "verification": verification26,
            "artifact_inventory_sha256": phase_b.canonical_sha256(
                phase_b.all_regular_inventory(stage26_root)
            ),
            "task_sha256": phase_b.sha256_file(control26_paths[0], "Stage-26 task"),
        },
        "stage27": {
            "artifact": artifact27,
            "verification": verification27,
            "artifact_inventory_sha256": phase_b.canonical_sha256(
                phase_b.all_regular_inventory(stage27_root)
            ),
            "task_sha256": phase_b.sha256_file(control27_paths[0], "Stage-27 task"),
        },
    }
    phase_b.write_json_new(controls / "provenance.json", provenance)

    local_root = STAGE / "local-same-host-controls"
    local_seal = load(local_root / "result-seal.json", "local control seal")
    local_payload = dict(local_seal)
    local_claimed = local_payload.pop("seal_payload_sha256", None)
    require(
        local_claimed == phase_b.canonical_sha256(local_payload)
        and phase_b.all_regular_inventory(local_root, {"result-seal.json"})
        == local_seal["inventory"],
        "local same-host controls changed",
    )
    local = load(local_root / "summary.json", "local same-host controls")
    require(
        local["native_xor"]["status"] == "unknown_inconclusive"
        and local["direct_mitm"]["status"] == "sat"
        and local["direct_mitm"]["source_witness_valid"] is True,
        "local control terminals changed",
    )

    native26 = next(row for row in control26["backends"] if row["solver"] == "native-xor")
    wdsat26 = next(row for row in control26["backends"] if row["solver"] == "wdsat")
    cms26 = next(row for row in control26["backends"] if row["solver"] == "cryptominisat")
    mitm27 = control27["result"]
    f4_3 = task3["result"]
    f4_4 = task4["result"]
    f4_5 = task5["result"]
    f4_6 = task6["result"]
    for row, status in (
        (native26, "unknown_inconclusive"),
        (wdsat26, "timeout_inconclusive"),
        (cms26, "timeout_inconclusive"),
        (mitm27, "sat"),
        (f4_5, "sat"),
    ):
        require(row["status"] == status, f"same-target terminal changed for {row['solver']}")

    table = [
        {
            "backend": "native-f4-selected",
            "formulation": "fixed-X1 direct S4; non-rational X1 masks skipped",
            "host": "local arm64",
            "status": f4_5["status"],
            "classification": "true_positive",
            "wall_seconds": f4_5["metrics"]["wall_seconds"],
            "total_core_seconds": f4_5["metrics"]["total_core_seconds"],
            "peak_rss_bytes": f4_5["metrics"]["peak_rss_bytes"],
            "conflicts": None,
            "operations": {
                "word_xors_elimination_only": f4_5["backend_report"]["cost"]["ops"],
                "f4_calls": f4_5["backend_report"]["cost"]["extra"]["f4_calls"],
            },
        },
        {
            "backend": "native-xor",
            "formulation": "frozen S4 correspondence system",
            "host": "local arm64",
            "status": local["native_xor"]["status"],
            "classification": "inconclusive",
            "wall_seconds": local["native_xor"]["metrics"]["wall_seconds"],
            "total_core_seconds": local["native_xor"]["metrics"]["total_core_seconds"],
            "peak_rss_bytes": local["native_xor"]["metrics"]["peak_rss_bytes"],
            "conflicts": local["native_xor"]["conflicts"],
            "operations": None,
        },
        {
            "backend": "direct-mitm",
            "formulation": "exact factor-point pair table",
            "host": "local arm64",
            "status": local["direct_mitm"]["status"],
            "classification": "true_positive",
            "wall_seconds": local["direct_mitm"]["metrics"]["wall_seconds"],
            "total_core_seconds": local["direct_mitm"]["metrics"]["total_core_seconds"],
            "peak_rss_bytes": local["direct_mitm"]["metrics"]["peak_rss_bytes"],
            "conflicts": None,
            "operations": {
                "group_additions": local["direct_mitm"]["backend_report"]["group_additions"],
                "pair_entries": local["direct_mitm"]["backend_report"]["pair_entries"],
            },
        },
        {
            "backend": "wdsat",
            "formulation": "frozen S4 ANF",
            "host": "GitHub Linux x86_64",
            "status": wdsat26["status"],
            "classification": "inconclusive",
            "wall_seconds": wdsat26["metrics"]["wall_seconds"],
            "total_core_seconds": wdsat26["metrics"]["total_core_seconds"],
            "peak_rss_bytes": wdsat26["metrics"]["peak_rss_bytes"],
            "conflicts": None,
            "operations": None,
        },
        {
            "backend": "cryptominisat",
            "formulation": "frozen S4 CNF plus XOR",
            "host": "GitHub Linux x86_64",
            "status": cms26["status"],
            "classification": "inconclusive",
            "wall_seconds": cms26["metrics"]["wall_seconds"],
            "total_core_seconds": cms26["metrics"]["total_core_seconds"],
            "peak_rss_bytes": cms26["metrics"]["peak_rss_bytes"],
            "conflicts": None,
            "operations": None,
        },
        {
            "backend": "magma-f4",
            "formulation": "frozen direct Boolean F4 source",
            "host": None,
            "status": "not_run_licensed_tool_missing",
            "classification": "not_run",
            "wall_seconds": None,
            "total_core_seconds": None,
            "peak_rss_bytes": None,
            "conflicts": None,
            "operations": None,
        },
        {
            "backend": "ggmp",
            "formulation": "GGMP invariant-kernel base exists only in the frozen n31 cell",
            "host": None,
            "status": "not_same_instance",
            "classification": "not_applicable",
            "wall_seconds": None,
            "total_core_seconds": None,
            "peak_rss_bytes": None,
            "conflicts": None,
            "operations": None,
        },
    ]

    selected_metrics = f4_5["metrics"]
    baseline_metrics = f4_4["metrics"]
    trace_metrics = f4_6["metrics"]
    mitm_metrics = local["direct_mitm"]["metrics"]
    build5 = run_build_receipt("attempt-05-rational-x1")
    build_resources_selected = build5["resources"]
    comparison = {
        "selected_over_full_fixed_x1": {
            "wall_ratio": ratio(selected_metrics["wall_seconds"], baseline_metrics["wall_seconds"]),
            "core_ratio": ratio(selected_metrics["total_core_seconds"], baseline_metrics["total_core_seconds"]),
            "rss_ratio": ratio(selected_metrics["peak_rss_bytes"], baseline_metrics["peak_rss_bytes"]),
            "word_xor_ratio": ratio(
                f4_5["backend_report"]["cost"]["ops"],
                f4_4["backend_report"]["cost"]["ops"],
            ),
            "same_verified_witness_x": [
                point["x"] for point in f4_5["backend_report"]["witness_points"]
            ]
            == [point["x"] for point in f4_4["backend_report"]["witness_points"]],
        },
        "rejected_trace_over_selected": {
            "wall_ratio": ratio(trace_metrics["wall_seconds"], selected_metrics["wall_seconds"]),
            "core_ratio": ratio(trace_metrics["total_core_seconds"], selected_metrics["total_core_seconds"]),
            "rss_ratio": ratio(trace_metrics["peak_rss_bytes"], selected_metrics["peak_rss_bytes"]),
            "word_xor_ratio": ratio(
                f4_6["backend_report"]["cost"]["ops"],
                f4_5["backend_report"]["cost"]["ops"],
            ),
        },
        "selected_over_same_host_direct_mitm": {
            "wall_ratio": ratio(selected_metrics["wall_seconds"], mitm_metrics["wall_seconds"]),
            "core_ratio": ratio(selected_metrics["total_core_seconds"], mitm_metrics["total_core_seconds"]),
            "rss_ratio": ratio(selected_metrics["peak_rss_bytes"], mitm_metrics["peak_rss_bytes"]),
        },
        "selected_cold_build_plus_run": {
            "wall_seconds": build_resources_selected["summed_process_wall_seconds"]
            + load(STAGE / "attempt-05-rational-x1/run-summary.json", "selected summary")[
                "outer_resources"
            ]["wall_seconds"],
            "total_core_seconds": build_resources_selected["total_core_seconds"]
            + load(STAGE / "attempt-05-rational-x1/run-summary.json", "selected summary")[
                "outer_resources"
            ]["total_core_seconds"],
            "peak_rss_bytes": max(
                build_resources_selected["largest_child_peak_rss_bytes"],
                selected_metrics["peak_rss_bytes"],
            ),
        },
    }

    # Non-overlapping measured campaign lower bound.  Build receipts are
    # deduplicated; completed runs use their inclusive outer receipt, while the
    # two early runner failures have only their preserved child receipts.
    build_receipts = {}
    for name in (
        "attempt-01-watchdog-timeout",
        "attempt-02-budgeted-terminal",
        "attempt-03-fixed-x1",
        "attempt-05-rational-x1",
        "attempt-06-trace",
    ):
        receipt = run_build_receipt(name)
        build_receipts[receipt["receipt_payload_sha256"]] = receipt
    build_rows = [receipt["resources"] for receipt in build_receipts.values()]
    science_rows = []
    for name in ("attempt-01-watchdog-timeout", "attempt-02-budgeted-terminal"):
        task_root = next((STAGE / name / "tasks").iterdir())
        science_rows.extend(metric_records(task_root))
    for name in (
        "attempt-03-fixed-x1",
        "attempt-04-fixed-x1-long",
        "attempt-05-rational-x1",
        "attempt-06-trace",
    ):
        outer = load(STAGE / name / "run-summary.json", f"{name} summary")["outer_resources"]
        science_rows.append(
            {
                "wall_seconds": outer["wall_seconds"],
                "total_core_seconds": outer["total_core_seconds"],
                "peak_rss_bytes": load(
                    STAGE / name / "run-summary.json", f"{name} summary"
                )["all_task_process_resources"]["maximum_individual_process_rss_bytes"],
            }
        )
    science_rows.append(
        {
            "wall_seconds": local["outer_resources"]["wall_seconds"],
            "total_core_seconds": local["outer_resources"]["total_core_seconds"],
            "peak_rss_bytes": max(
                local["native_xor"]["metrics"]["peak_rss_bytes"],
                local["direct_mitm"]["metrics"]["peak_rss_bytes"],
            ),
        }
    )
    development_rows = metric_records(STAGE / "development")
    measured = [
        {
            "wall_seconds": row["summed_process_wall_seconds"],
            "total_core_seconds": row["total_core_seconds"],
            "peak_rss_bytes": row["largest_child_peak_rss_bytes"],
        }
        for row in build_rows
    ] + science_rows + development_rows
    campaign = {
        "measured_lower_bound": add_resources(measured),
        "distinct_successful_clean_builds": len(build_receipts),
        "failed_build_and_dependency_fetch_records": len(development_rows),
        "outer_driver_missing_for_attempts": [1, 2],
        "unmetered_development_commands_present": True,
        "unmetered_scope": [
            "incremental cargo check and focused cargo tests",
            "Python unit tests and evidence composition",
            "GitHub artifact downloads and local file copies",
        ],
        "complete_campaign_cost": None,
        "reason_complete_is_null": "some development commands predated an outer meter; measured totals are a lower bound",
    }

    result = {
        "schema": RESULT_SCHEMA,
        "status": "complete_single_target_native_f4_true_positive",
        "date": "2026-09-22",
        "claim_boundary": (
            "Strong internal engineering and one finite public toy-PDP native-F4 result. "
            "It is not licensed Magma F4, a full 160-input native-F4 panel, an end-to-end "
            "index-calculus speedup, or a Koblitz index-calculus SOTA."
        ),
        "instance": {
            "curve": "K_1 over GF(2^59)",
            "n": 59,
            "ell": 9,
            "m": 3,
            "basis": "standard polynomial subspace",
            "blind_instance_id": BLIND_ID,
            "source_instance_id": SOURCE_ID,
            "target": task5["target"],
            "target_class_opened_after_run_seal": "decomposable",
        },
        "factor_base": {
            "definition": "span(1,z,...,z^8), with exact rational-lift membership checks",
            "factor_base_logs_known_by_construction": False,
            "target_subgroup_enumerated": False,
            "discrete_log_labels_used": False,
        },
        "selected_native_f4": {
            "source_revision": load(
                STAGE / "attempt-05-rational-x1/execution-plan.json", "selected plan"
            )["source_revision"],
            "run_inventory_sha256": seal5["inventory_sha256"],
            "status": f4_5["status"],
            "classification": "true_positive",
            "source_model_valid": f4_5["source_model_valid"],
            "source_witness_valid": f4_5["source_witness_valid"],
            "metrics": f4_5["metrics"],
            "report": f4_5["backend_report"],
        },
        "same_target_table": table,
        "comparisons": comparison,
        "attempt_history": {
            "attempt_01": failed_1,
            "attempt_02": failed_2,
            "attempt_03_run_inventory_sha256": seal3["inventory_sha256"],
            "attempt_04_run_inventory_sha256": seal4["inventory_sha256"],
            "attempt_05_run_inventory_sha256": seal5["inventory_sha256"],
            "attempt_06_run_inventory_sha256": seal6["inventory_sha256"],
        },
        "control_provenance": provenance,
        "campaign_accounting": campaign,
        "gates": {
            "all_seven_gates_passed": False,
            "1_all_stage_resource_charging": "partial_measured_lower_bound_complete_total_null",
            "2_same_instance_solver_matrix": "partial_native_f4_one_n59_target_added_licensed_magma_and_full_native_panel_missing",
            "3_single_core_core_memory_conflicts_wall": "partial_native_f4_complete_conflicts_null_magma_missing",
            "4_n31_n41_larger_pdp": "inherited_satisfied_finite_coverage",
            "5_unknown_scalar_without_known_base_logs": "inherited_satisfied_finite_controls",
            "6_full_cost_vs_automorphism_rho": "unchanged_not_tested_by_pdp_stage_and_current_full_cost_gate_false",
            "7_independent_external_reproduction": "missing",
        },
        "licensed_magma_f4_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    phase_b.write_json_new(STAGE / "result.json", result)
    return result


def write_markdown(result: dict[str, Any]) -> None:
    rows = result["same_target_table"]
    lines = [
        "# Stage 159: native F4 on one frozen n=59 Phase-B target",
        "",
        result["claim_boundary"],
        "",
        "The selected fixed-X1 F4 arm is class-blind during execution. It verifies the historical source export, defines the base as the public polynomial subspace, skips only x-values with no rational curve lift, and accepts SAT only after an exact three-point group check. Truth scoring occurs after the run seal.",
        "",
        "| Backend | Host | Terminal | Wall s | Core s | Peak RSS B | Conflicts / native work |",
        "|:--|:--|:--|--:|--:|--:|:--|",
    ]
    for row in rows:
        work = (
            "null"
            if row["operations"] is None and row["conflicts"] is None
            else str(row["operations"] if row["operations"] is not None else row["conflicts"])
        )
        lines.append(
            "| {backend} | {host} | {status} | {wall} | {core} | {rss} | {work} |".format(
                backend=row["backend"],
                host=row["host"] or "not run",
                status=row["status"],
                wall="null" if row["wall_seconds"] is None else f"{row['wall_seconds']:.6f}",
                core="null" if row["total_core_seconds"] is None else f"{row['total_core_seconds']:.6f}",
                rss="null" if row["peak_rss_bytes"] is None else row["peak_rss_bytes"],
                work=work,
            )
        )
    c = result["comparisons"]
    lines += [
        "",
        f"Rational-X1 membership reduced F4 wall to {c['selected_over_full_fixed_x1']['wall_ratio']:.3f} of the complete fixed-X1 baseline, core to {c['selected_over_full_fixed_x1']['core_ratio']:.3f}, and elimination word XORs to {c['selected_over_full_fixed_x1']['word_xor_ratio']:.3f}; the exact witness is unchanged.",
        "",
        f"On the same local host and target, selected native F4 is {c['selected_over_same_host_direct_mitm']['wall_ratio']:.2f} times direct MITM wall, {c['selected_over_same_host_direct_mitm']['core_ratio']:.2f} times its CPU, and {c['selected_over_same_host_direct_mitm']['rss_ratio']:.2f} times its RSS. The trace constraint is rejected for speed: {c['rejected_trace_over_selected']['wall_ratio']:.3f} wall and {c['rejected_trace_over_selected']['word_xor_ratio']:.3f} word XORs, despite lower RSS.",
        "",
        "The WDSat and CryptoMiniSat rows are exact-target frozen Linux receipts, but their host differs from the local F4 host; their seconds are descriptive and are not used as paired speed ratios. Licensed Magma remains unexecuted. GGMP has no same n=59 cell. Conflict count is null for F4 and MITM because neither exposes SAT conflicts.",
        "",
        "The measured campaign total is a lower bound because some early developer checks lacked an outer meter. Full end-to-end IC cost and the ratio to automorphism-optimized rho remain unchanged and false. This does not establish a SOTA.",
        "",
    ]
    phase_b.write_new(STAGE / "RESULTS.md", ("\n".join(lines)).encode())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage26-cell", type=Path, required=True)
    parser.add_argument("--stage27-cell", type=Path, required=True)
    args = parser.parse_args()
    try:
        for path in (STAGE / "result.json", STAGE / "RESULTS.md", STAGE / "result-seal.json"):
            require(not path.exists(), f"refusing to overwrite {path}")
        result = compose(args.stage26_cell.resolve(strict=True), args.stage27_cell.resolve(strict=True))
        write_markdown(result)
        inventory = phase_b.all_regular_inventory(STAGE, {"result-seal.json"})
        seal = {
            "schema": SEAL_SCHEMA,
            "status": "stage159_native_f4_result_frozen",
            "result_sha256": phase_b.sha256_file(STAGE / "result.json", "Stage-159 result"),
            "inventory": inventory,
            "inventory_sha256": phase_b.canonical_sha256(inventory),
        }
        seal["seal_payload_sha256"] = phase_b.canonical_sha256(seal)
        phase_b.write_json_new(STAGE / "result-seal.json", seal)
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, StopIteration, phase_b.PhaseBError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
