#!/usr/bin/env python3
"""Compose the current seven-gate audit through Stage 186."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


STAGE = Path(__file__).resolve().parent
GATES = STAGE.parent


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stage_dir(number: int) -> Path:
    matches = list(GATES.glob(f"stage-{number}-*"))
    if len(matches) != 1:
        raise SystemExit(f"expected one Stage {number}, found {len(matches)}")
    return matches[0]


records = {}
for number in range(174, 187):
    root = stage_dir(number)
    result_path = root / "result.json"
    verification_path = root / "verification.json"
    result = json.loads(result_path.read_text())
    verification = json.loads(verification_path.read_text())
    if verification.get("failed"):
        raise SystemExit(f"Stage {number} verification has failures")
    records[number] = {
        "directory": root.name,
        "result": result,
        "result_sha256": sha256(result_path),
        "verification_sha256": sha256(verification_path),
        "verification_checks": verification["checks"],
    }

base = records[174]["result"]["campaign_accounting"]["measured_lower_bound_through_stage174"]
incremental = []
for number, field in (
    (175, "campaign_charge"),
    (176, "campaign_charge"),
    (177, "campaign_charge"),
    (178, "campaign_charge"),
    (179, "campaign_charge"),
):
    incremental.append((number, records[number]["result"][field]))

# Stage 180's result carries one inherited Stage 179 build. Count only its
# fresh four-column build and four screen processes.
stage180 = records[180]["result"]
stage180_processes = [
    stage180["builds"]["four"]["process"],
    *(run["process"] for run in stage180["runs"]),
]
incremental.append(
    (
        180,
        {
            "components": len(stage180_processes),
            "wall_seconds_sum": sum(x["metrics"]["wall_seconds"] for x in stage180_processes),
            "total_core_seconds_sum": sum(
                x["metrics"]["total_core_seconds"] for x in stage180_processes
            ),
            "peak_rss_bytes_max": max(x["metrics"]["peak_rss_bytes"] for x in stage180_processes),
        },
    )
)
for number, field in (
    (181, "campaign_charge"),
    (182, "new_campaign_charge"),
    (183, "campaign_charge"),
    (184, "new_campaign_charge"),
    (185, "campaign_charge"),
    (186, "new_campaign_charge"),
):
    incremental.append((number, records[number]["result"][field]))

increment = {
    "resource_components": sum(item["components"] for _, item in incremental),
    "summed_wall_seconds": sum(item["wall_seconds_sum"] for _, item in incremental),
    "total_core_seconds": sum(item["total_core_seconds_sum"] for _, item in incremental),
    "peak_rss_bytes": max(item["peak_rss_bytes_max"] for _, item in incremental),
}
cumulative = {
    "resource_components": base["resource_components"] + increment["resource_components"],
    "summed_wall_seconds": base["summed_wall_seconds"] + increment["summed_wall_seconds"],
    "total_core_seconds": base["total_core_seconds"] + increment["total_core_seconds"],
    "peak_rss_bytes": max(base["peak_rss_bytes"], increment["peak_rss_bytes"]),
    "complete_campaign_cost": None,
    "reason_complete_is_null": "interactive incremental compiles and tests outside process meters remain; the value is a measured lower bound with inherited setup de-duplicated",
}

comments_path = STAGE / "external-comments.json"
comments = json.loads(comments_path.read_text())
authors = sorted({comment["user"]["login"] for comment in comments})
external = [comment for comment in comments if comment["user"]["login"] != "aburan28"]

decisions = {
    number: records[number]["result"].get("decision", {}).get("status")
    or records[number]["result"].get("status")
    for number in range(175, 187)
}
expected_decisions = {
    175: "REJECTED_REGRESSION",
    176: "REJECTED_NO_JOINT_SCREEN_WIN",
    177: "REJECTED_PAIRED_TIMING_REGRESSION",
    178: "REJECTED_CPU_THRESHOLD",
    179: "REJECTED_WIDER_TABLE_REGRESSION",
    180: "REJECTED_SCREEN_WORK_REGRESSION",
    181: "RETAIN_OPT_IN_CONFIRMATION_WALL_FAIL",
    182: "REJECTED_SINGLE_CORE_TIMING_REGRESSION",
    183: "ACCEPTED_FOR_FULL_M4RI_STACK",
    184: "ACCEPTED_FOR_REPOSITORY_DEFAULT",
    185: "SELECTED_DEFAULT_REPLAY_PASS",
    186: "REJECTED_SCREEN_CPU_REGRESSION",
}
if decisions != expected_decisions:
    raise SystemExit(f"decision chain mismatch: {decisions}")

gates = {
    "1_all_stage_resource_charging": "partial_measured_lower_bound_extended_through_stage186_complete_campaign_total_null",
    "2_same_instance_solver_matrix": "partial_native_f4_wdsat_cryptominisat_direct_mitm_present_ggmp_same_cell_proven_inapplicable_licensed_magma_missing",
    "3_single_core_core_memory_conflicts_wall": "partial_native_f4_parallel_and_single_core_complete_conflicts_null_magma_missing",
    "4_n31_n41_larger_pdp": "inherited_satisfied_finite_coverage",
    "5_unknown_scalar_without_known_base_logs": "inherited_satisfied_finite_controls",
    "6_full_cost_vs_automorphism_rho": "unchanged_current_full_cost_gate_false",
    "7_independent_external_reproduction": "requested_but_missing_zero_unaffiliated_comments",
    "all_seven_gates_passed": False,
}
result = {
    "schema": "koblitz_stage187_current_gate_audit.v1",
    "predecessor": "stage-186-dense-pair-elimination-policy-20261001",
    "stage_records": {
        str(number): {key: value for key, value in record.items() if key != "result"}
        for number, record in records.items()
    },
    "decision_chain": {str(k): v for k, v in decisions.items()},
    "selected_implementation": {
        "f4_engine": "current repository Boolean F4",
        "pair_selector": "dense exact LCM groups and submask cover lookup through 20 variables",
        "pair_selector_default": True,
        "quadratic_control": "F4_F2_DENSE_PAIR_SELECT=0",
        "elimination": "current five-column BlockTables",
        "full_m4ri_control": "F4_F2_FULL_M4RI=1",
        "factor_base": "algebraic polynomial subspace with exact rational lift checks",
        "target_subgroup_enumerated": False,
        "discrete_log_labels_used": False,
    },
    "stage175_through_186_increment": increment,
    "measured_campaign_lower_bound_through_stage186": cumulative,
    "external_review_snapshot": {
        "issue": "https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks/issues/1",
        "comments": len(comments),
        "authors": authors,
        "unaffiliated_comments": len(external),
        "snapshot_sha256": sha256(comments_path),
        "snapshot_bytes": comments_path.stat().st_size,
    },
    "gates": gates,
    "full_cost_gate_passed": False,
    "independent_reproduction_and_novelty_review_passed": False,
    "sota_claim": False,
    "claim_boundary": "Dense pair selection is a reproducible one-target current-F4 engineering improvement. Full attack cost, licensed Magma, matched automorphism-rho crossover, and unaffiliated reproduction/novelty remain open; this is not Koblitz index-calculus SOTA.",
}
(STAGE / "audit.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 187: current seven-gate audit",
    "",
    "Stages 175–186 reconcile PR #660 with current `main`, run the frozen single target through current repository F4, preserve four rejected tuning families, add an opt-in full-matrix M4RI control, and select dense exact pair updates as the repository default after two independent same-binary panels plus an exact-commit replay.",
    "",
    f"The unique Stage 175–186 increment is {increment['summed_wall_seconds']:.6f} wall seconds, {increment['total_core_seconds']:.6f} core-seconds, {increment['peak_rss_bytes']} bytes maximum RSS, and {increment['resource_components']} process components. The measured campaign lower bound through Stage 186 is {cumulative['summed_wall_seconds']:.6f} wall seconds, {cumulative['total_core_seconds']:.6f} core-seconds, {cumulative['peak_rss_bytes']} bytes maximum RSS, and {cumulative['resource_components']} components. Complete cost remains null because some interactive compile/test work was not outer-metered.",
    "",
    "The selected Phase B implementation is current five-column BlockTables plus dense exact pair selection. `F4_F2_DENSE_PAIR_SELECT=0` retains the quadratic control; `F4_F2_FULL_M4RI=1` retains full M4RI as a research control. The algebraic factor base still enumerates neither the target subgroup nor discrete-log labels.",
    "",
    f"The public reproduction issue snapshot contains {len(comments)} comments from {', '.join(authors)} and {len(external)} unaffiliated comments. Independent reproduction and novelty review therefore remain missing.",
    "",
    "Gates 1–3 remain partial, gates 4–5 retain finite-control coverage, gate 6 remains false, and gate 7 remains open. All seven gates and the SOTA claim remain false.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"increment": increment, "cumulative": cumulative, "external": result["external_review_snapshot"], "gates": gates}, indent=2))
