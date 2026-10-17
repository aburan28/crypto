#!/usr/bin/env python3
"""Agent-facing autolab control plane for the IC boundary ledger.

Reads docs/ic/boundary_targets.json (schema_version 2), fail-closed validates
measurement reports, and launches public-synthetic Koblitz producer runs
against the next admitted beat targets.

This runner is local to the crypto repository. It does not impersonate any
external Autoresearcher dispatcher or Coordinator.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import shlex
import json
import math
import os
import platform
import resource
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
RESEARCH_DIR = HERE.parent
REPO = HERE.parents[2]
PROTOCOL_PATH = HERE / "protocol.json"
RUNS_DIR = HERE / "runs"
CURRENT_PATH = RUNS_DIR / "current.json"
LOCK_PATH = RUNS_DIR / "autolab.lock"
TASK_ID = "TASK-IC-BOUNDARY-AUTOLAB-20260910"

# Abstract measurement_schema keys -> accepted concrete aliases in reports.
FIELD_ALIASES: dict[str, tuple[str, ...]] = {
    "n_or_bits": ("n_or_bits", "n", "bits"),
    "factor_base_size_F": ("factor_base_size_F", "factor_base_size", "F", "|F|"),
    "orbit_count_K": ("orbit_count_K", "orbit_columns", "K", "orbit_columns_K"),
    "dimension_l_or_dim": ("dimension_l_or_dim", "dimension", "l", "ell", "dim"),
    "construction_method": ("construction_method",),
    "materialized": ("materialized",),
    "construction_wall_ms": ("construction_wall_ms",),
    "retained_bytes": ("retained_bytes",),
    "m_summands": ("m_summands", "m", "summands"),
    "unknowns": ("unknowns",),
    "system_degree": ("system_degree",),
    "eq_var_ratio": ("eq_var_ratio",),
    "ffd_or_degree_of_regularity": (
        "ffd_or_degree_of_regularity",
        "ffd",
        "degree_of_regularity",
        "DoR",
        "dor",
    ),
    "oracle_class": ("oracle_class",),
    "median_ms_per_target": ("median_ms_per_target",),
    "largest_solvable": ("largest_solvable",),
    "base_id_or_hash": ("base_id_or_hash", "base_id", "base_hash"),
    "eta_or_coverage_policy": ("eta_or_coverage_policy", "eta", "coverage_policy"),
    "pr_decomposition_or_hit_rate_with_ci": (
        "pr_decomposition_or_hit_rate_with_ci",
        "pr_decomposition",
        "hit_rate",
        "hit_rate_with_ci",
    ),
    "trials_per_relation": ("trials_per_relation",),
    "target_mix": ("target_mix",),
    "orbit_columns_K": ("orbit_columns_K", "orbit_columns", "K"),
    "relations_collected": ("relations_collected",),
    "relations_needed": ("relations_needed",),
    "surplus": ("surplus",),
    "matrix_dims": ("matrix_dims",),
    "sparse_or_dense": ("sparse_or_dense",),
    "rank_accumulation": ("rank_accumulation",),
    "la_wall_ms_or_la_charged_ms": (
        "la_wall_ms_or_la_charged_ms",
        "la_wall_ms",
        "la_charged_ms",
    ),
    "recovered_d_verified": ("recovered_d_verified", "recovered_d", "d_verified"),
    "stage_timers": ("stage_timers",),
    "claim_boundary": ("claim_boundary",),
    "timing_class": ("timing_class",),
    "target_count": ("target_count",),
    "ic_target_hash": ("ic_target_hash",),
    "rho_target_hash": ("rho_target_hash",),
    "ic_online_wall_ms": ("ic_online_wall_ms",),
    "rho_online_wall_ms": ("rho_online_wall_ms",),
    "online_speedup": ("online_speedup",),
    "online_interval": ("online_interval",),
    "same_resource_envelope": ("same_resource_envelope",),
    "ic_scalar_verified": ("ic_scalar_verified",),
    "rho_scalar_verified": ("rho_scalar_verified",),
    "independent_validation": ("independent_validation",),
    "ic_replay_certificate_sha256": ("ic_replay_certificate_sha256",),
    "rho_replay_certificate_sha256": ("rho_replay_certificate_sha256",),
    "ic_resource_envelope": ("ic_resource_envelope",),
    "rho_resource_envelope": ("rho_resource_envelope",),
    "rho_policy": ("rho_policy",),
    "ic_cost": ("ic_cost",),
    "rho_cost": ("rho_cost",),
    "automorphism_discount": ("automorphism_discount",),
    "all_stages_charged_same_series": ("all_stages_charged_same_series",),
    "verdict": ("verdict",),
    "independent_replay_pointer": ("independent_replay_pointer",),
    "fixture_hash": ("fixture_hash",),
    "executable_or_source_hash": (
        "executable_or_source_hash",
        "executable_hash",
        "source_hash",
    ),
    "host_id": ("host_id",),
    "resource_caps": ("resource_caps",),
    "seeds": ("seeds", "seed"),
    "claim_boundary_non_claims": (
        "claim_boundary_non_claims",
        "non_claims",
        "claim_boundary",
    ),
    "operation_accounting": ("operation_accounting",),
}

ONLINE_REQUIRED_STAGES: dict[str, tuple[str, ...]] = {
    "ic_included_stages": (
        "target_query",
        "target_PDP",
        "target_relation_check",
        "target_descent",
        "target_recovery_check",
    ),
    "rho_included_stages": ("walk", "collision", "recovery_check"),
}


class AutolabError(RuntimeError):
    """Operator-facing failure that must not be treated as a ledger beat."""


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def independent_replay_pointer(run: Path) -> str:
    """Return the reserved path for this run's independent replay receipt."""
    replay_path = run / "validation" / "independent_replay.json"
    try:
        return replay_path.relative_to(REPO).as_posix()
    except ValueError:
        return str(replay_path)


def sha256(path: Path | str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AutolabError(message)


def load_protocol() -> dict[str, Any]:
    protocol = read_json(PROTOCOL_PATH)
    require(protocol.get("task_id") == TASK_ID, "protocol task_id mismatch")
    return protocol


def ledger_path(protocol: dict[str, Any]) -> Path:
    return REPO / protocol["ledger"]["path"]


def load_ledger(protocol: dict[str, Any]) -> dict[str, Any]:
    path = ledger_path(protocol)
    require(path.is_file(), f"boundary ledger missing: {path}")
    ledger = read_json(path)
    required = protocol["ledger"]["required_schema_version"]
    require(
        ledger.get("schema_version") == required,
        f"ledger schema_version must be {required}, got {ledger.get('schema_version')}",
    )
    require(
        isinstance(ledger.get("measurement_schema"), dict),
        "ledger missing measurement_schema (fail closed)",
    )
    require(
        ledger["measurement_schema"].get("fail_closed") is True,
        "measurement_schema.fail_closed must be true",
    )
    return ledger


def field_present(report: dict[str, Any], abstract_key: str) -> bool:
    aliases = FIELD_ALIASES.get(abstract_key, (abstract_key,))
    for alias in aliases:
        if alias not in report:
            continue
        value = report[alias]
        if value is None:
            continue
        if isinstance(value, str) and not value.strip():
            continue
        return True
    return False


def missing_fields(report: dict[str, Any], keys: list[str]) -> list[str]:
    return [key for key in keys if not field_present(report, key)]


def operation_accounting_errors(accounting: Any) -> list[str]:
    """Check the operation-count block that sits next to the wall ratio.

    Native IC probes and rho steps have different costs.  Keep both counts
    and their units, while requiring a calibration receipt before reporting
    an operation speedup.
    """
    if not isinstance(accounting, dict):
        return ["operation_accounting must be an object"]
    errors = []
    ic_ops = accounting.get("ic_online_operations")
    rho_ops = accounting.get("rho_online_operations")
    for name, value in (("ic_online_operations", ic_ops), ("rho_online_operations", rho_ops)):
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            errors.append(f"operation_accounting.{name} must be nonnegative and finite")
    units = accounting.get("operation_units")
    if not isinstance(units, dict) or not all(
        isinstance(units.get(arm), str) and units[arm].strip() for arm in ("ic", "rho")
    ):
        errors.append("operation_accounting.operation_units must name the ic and rho units")
    assumption = accounting.get("unit_assumption")
    if not isinstance(assumption, str) or not assumption.strip():
        errors.append("operation_accounting.unit_assumption must state how the units compare")
    status = accounting.get("comparison_status")
    if status not in ("native_counters_only", "calibrated_common_unit"):
        errors.append("operation_accounting.comparison_status must name the calibration state")
    raw_quotient = accounting.get("rho_per_ic_native_counter")
    if not errors:
        if ic_ops == 0:
            if raw_quotient is not None:
                errors.append("operation_accounting.rho_per_ic_native_counter must be null for zero IC probes")
        elif raw_quotient is not None:
            expected = float(rho_ops) / float(ic_ops)
            if type(raw_quotient) not in (int, float) or not math.isfinite(raw_quotient) or raw_quotient < 0 or not math.isclose(
                float(raw_quotient), expected, rel_tol=1e-8, abs_tol=1e-12
            ):
                errors.append("operation_accounting.rho_per_ic_native_counter does not match the raw counts")
    ratio = accounting.get("ops_speedup_online")
    if ratio is not None:
        calibrated = (
            status == "calibrated_common_unit"
            and isinstance(units, dict)
            and units.get("ic") == units.get("rho")
            and isinstance(accounting.get("calibration_receipt"), str)
            and bool(accounting["calibration_receipt"].strip())
        )
        if not calibrated:
            errors.append("operation_accounting.ops_speedup_online requires a calibrated common unit and receipt")
        elif positive_cost(ic_ops) is None or positive_cost(rho_ops) is None:
            errors.append("operation_accounting calibrated counts must be positive")
        elif positive_cost(ratio) is None:
            errors.append("operation_accounting.ops_speedup_online must be positive and finite")
        elif not errors:
            expected = float(rho_ops) / float(ic_ops)
            if not math.isclose(float(ratio), expected, rel_tol=1e-8, abs_tol=1e-12):
                errors.append("operation_accounting.ops_speedup_online does not match calibrated counts")
    elif status == "calibrated_common_unit":
        errors.append("operation_accounting.ops_speedup_online is required for calibrated counts")
    return errors


def hardware_accounting_errors(accounting: Any) -> list[str]:
    """Check the optional common retired-instruction counter for a paired target."""
    if not isinstance(accounting, dict):
        return ["hardware_accounting must be an object"]
    errors = []
    status = accounting.get("status")
    if status not in ("measured_common_counter", "counter_unavailable"):
        errors.append("hardware_accounting.status is invalid")
    if accounting.get("unit") != "calling-thread user instructions retired":
        errors.append("hardware_accounting.unit must name the measured unit")
    if not isinstance(accounting.get("method"), str) or not accounting["method"].strip():
        errors.append("hardware_accounting.method must describe the counter boundary")
    ic_count = accounting.get("ic_online_instructions")
    rho_count = accounting.get("rho_online_instructions")
    for arm, count in (("ic", ic_count), ("rho", rho_count)):
        if count is not None and (type(count) is not int or count <= 0):
            errors.append(f"hardware_accounting.{arm}_online_instructions must be positive when present")
    ratio = accounting.get("rho_per_ic_instructions")
    if status == "measured_common_counter":
        for arm in ("ic", "rho"):
            enabled = accounting.get(f"{arm}_time_enabled_ns")
            running = accounting.get(f"{arm}_time_running_ns")
            if (type(enabled) is not int or enabled <= 0 or
                    type(running) is not int or running != enabled):
                errors.append(f"hardware_accounting.{arm} counter must run throughout its enabled interval")
            if accounting.get(f"{arm}_counter_error"):
                errors.append(f"hardware_accounting.{arm} counter error must be null")
        if type(ic_count) is not int or ic_count <= 0 or type(rho_count) is not int or rho_count <= 0:
            errors.append("hardware_accounting measured status requires both online instruction counts")
        elif positive_cost(ratio) is None or not math.isclose(
            float(ratio), rho_count / ic_count, rel_tol=1e-9, abs_tol=1e-12
        ):
            errors.append("hardware_accounting.rho_per_ic_instructions must match the counts")
    elif status == "counter_unavailable" and ratio is not None:
        errors.append("hardware_accounting.rho_per_ic_instructions must be null when unavailable")
    return errors


def validate_claim(
    report: dict[str, Any],
    *,
    stage: str,
    ledger: dict[str, Any],
) -> dict[str, Any]:
    schema = ledger["measurement_schema"]
    require(stage in schema, f"unknown stage for measurement schema: {stage}")
    stage_schema = schema[stage]
    required = list(stage_schema.get("required", []))
    global_required = list(schema.get("global_provenance_required", []))
    missing_stage = missing_fields(report, required)
    missing_global = missing_fields(report, global_required)
    pairing_errors: list[str] = []
    if stage == "vs_rho":
        paired = report.get("paired_target")
        if report.get("target_count") != 1:
            pairing_errors.append("target_count must be exactly one")
        if not isinstance(paired, dict):
            pairing_errors.append("paired_target must be an object")
        else:
            ic_q = paired.get("ic_public_q")
            rho_q = paired.get("rho_public_q")
            if not isinstance(ic_q, list) or not isinstance(rho_q, list):
                pairing_errors.append("both arms must record public target coordinates")
            elif ic_q != rho_q:
                pairing_errors.append("IC and rho public targets differ")
            if paired.get("same_public_point") is not True:
                pairing_errors.append("same_public_point is not verified")
        if report.get("ic_verified") is not True:
            pairing_errors.append("IC target recovery is not verified")
        if report.get("rho_verified") is not True:
            pairing_errors.append("rho target recovery is not verified")
        if positive_cost(report.get("ic_online_ms")) is None:
            pairing_errors.append("positive ic_online_ms is required")
        if positive_cost(report.get("rho_online_ms")) is None:
            pairing_errors.append("positive rho_online_ms is required")
        if positive_cost(report.get("online_speedup")) is None:
            pairing_errors.append("positive same-target online_speedup is required")
        elif positive_cost(report.get("ic_online_ms")) is not None and positive_cost(report.get("rho_online_ms")) is not None:
            expected = float(report["rho_online_ms"]) / float(report["ic_online_ms"])
            if not math.isclose(float(report["online_speedup"]), expected, rel_tol=1e-8, abs_tol=1e-9):
                pairing_errors.append("online_speedup does not equal rho_online_ms / ic_online_ms")
        for arm in ("ic", "rho"):
            phases = report.get(f"{arm}_online_phase_ms")
            total = report.get(f"{arm}_online_ms")
            if not isinstance(phases, dict) or not phases:
                pairing_errors.append(f"{arm}_online_phase_ms must record exclusive phases")
            elif any(type(value) not in (int, float) or not math.isfinite(value) or value < 0 for value in phases.values()):
                pairing_errors.append(f"{arm}_online_phase_ms must contain finite nonnegative costs")
            elif positive_cost(total) is not None:
                phase_sum = math.fsum(float(value) for value in phases.values())
                if not math.isclose(phase_sum, float(total), rel_tol=1e-8, abs_tol=0.02):
                    pairing_errors.append(f"{arm} online phase costs do not sum to online time")
            else:
                pairing_errors.append(f"{arm}_online_ms is missing")
        if report.get("all_stages_charged_same_series") is not True:
            pairing_errors.append("same-target charged intervals are not confirmed")
        if field_present(report, "operation_accounting"):
            pairing_errors.extend(operation_accounting_errors(report["operation_accounting"]))
    validation_errors: list[str] = []

    if stage == "vs_rho":
        paired = report.get("paired_target")
        if paired is not None:
            if not isinstance(paired, dict):
                validation_errors.append("paired_target must be an object")
            elif (
                paired.get("same_public_point") is not True
                or paired.get("ic_public_q") != paired.get("rho_public_q")
            ):
                validation_errors.append("paired_target must contain the same public point in both arms")
        if "all_stages_charged_same_series" in report and report["all_stages_charged_same_series"] is not True:
            validation_errors.append("all_stages_charged_same_series must be true when present")
        if field_present(report, "operation_accounting"):
            validation_errors.extend(operation_accounting_errors(report["operation_accounting"]))
        for key in ("candidate_id", "workload_id", "run_id"):
            value = report.get(key)
            if not isinstance(value, str) or not value.strip():
                validation_errors.append(f"{key} must be a nonempty string")
        candidate_id = report.get("candidate_id")
        workload_id = report.get("workload_id")
        run_id = report.get("run_id")
        candidate_manifest_hash = report.get("candidate_manifest_sha256")
        workload_manifest_hash = report.get("workload_manifest_sha256")
        if not isinstance(candidate_manifest_hash, str) or len(candidate_manifest_hash) != 64:
            validation_errors.append("candidate_manifest_sha256 must be a full SHA-256 hex digest")
        elif any(ch not in "0123456789abcdef" for ch in candidate_manifest_hash):
            validation_errors.append("candidate_manifest_sha256 must be lowercase hex")
        if not isinstance(workload_manifest_hash, str) or len(workload_manifest_hash) != 64:
            validation_errors.append("workload_manifest_sha256 must be a full SHA-256 hex digest")
        elif any(ch not in "0123456789abcdef" for ch in workload_manifest_hash):
            validation_errors.append("workload_manifest_sha256 must be lowercase hex")
        if isinstance(candidate_id, str) and isinstance(candidate_manifest_hash, str):
            import re
            match = re.search(r"h([0-9a-f]{12,64})$", candidate_id)
            if not candidate_id.startswith("IC1") or not match:
                validation_errors.append("candidate_id must use the IC1 identity format with a digest suffix")
            elif not candidate_manifest_hash.startswith(match.group(1)):
                validation_errors.append("candidate_id digest must match candidate_manifest_sha256")
        if isinstance(workload_id, str) and isinstance(workload_manifest_hash, str):
            if len(workload_id) != 12 or any(ch not in "0123456789abcdef" for ch in workload_id):
                validation_errors.append("workload_id must be 12 lowercase hex digits")
            elif not workload_manifest_hash.startswith(workload_id):
                validation_errors.append("workload_id must match workload_manifest_sha256")
        if all(isinstance(value, str) and value for value in (candidate_id, workload_id, run_id)):
            import re
            if not re.fullmatch(re.escape(candidate_id) + r"W" + re.escape(workload_id) + r"R[1-9][0-9]*", run_id):
                validation_errors.append("run_id must be <candidate_id>W<workload_id>R<run-number>")
        if type(report.get("target_count")) is not int or report["target_count"] != 1:
            validation_errors.append("target_count must equal 1")
        target_hashes = (report.get("ic_target_hash"), report.get("rho_target_hash"))
        if not all(isinstance(value, str) and value.strip() for value in target_hashes):
            validation_errors.append("IC and rho target hashes must be nonempty strings")
        elif target_hashes[0] != target_hashes[1]:
            validation_errors.append("IC and rho target hashes must match")
        timing_classes = stage_schema.get("timing_class_enum", [])
        if "primary_ic_online_phase_keys" in stage_schema:
            if report.get("timing_class") != "single_target_online":
                validation_errors.append("timing_class must be single_target_online for a primary vs_rho claim")
        elif report.get("timing_class") not in timing_classes:
            validation_errors.append(f"timing_class must be one of {timing_classes}")
        record_class = report.get("record_class")
        if record_class not in stage_schema.get("record_class_enum", []):
            validation_errors.append("record_class must identify exploratory or controlled wall evidence")
        if "hardware_accounting" in report:
            validation_errors.extend(hardware_accounting_errors(report["hardware_accounting"]))
        controlled_speedup = report.get("controlled_online_speedup")
        if controlled_speedup is None:
            if record_class == "verified_answer_controlled_wall":
                validation_errors.append("controlled wall record requires controlled_online_speedup")
        else:
            if record_class != "verified_answer_controlled_wall":
                validation_errors.append("controlled_online_speedup requires a controlled wall record")
            measured_speedup = positive_cost(report.get("online_speedup"))
            if positive_cost(controlled_speedup) is None or measured_speedup is None or not math.isclose(
                float(controlled_speedup), measured_speedup, rel_tol=1e-9, abs_tol=1e-12
            ):
                validation_errors.append("controlled_online_speedup must equal the paired online ratio")
            receipt = report.get("host_isolation_receipt")
            if not isinstance(receipt, dict) or receipt.get("status") != "PASS" or not isinstance(
                receipt.get("path"), str
            ) or not receipt["path"].strip():
                validation_errors.append("controlled_online_speedup requires a host-isolation receipt")
        if report.get("same_resource_envelope") is not True:
            validation_errors.append("same_resource_envelope must be true")
        for key in ("ic_scalar_verified", "rho_scalar_verified"):
            if report.get(key) is not True:
                validation_errors.append(f"{key} must be true")
        if report.get("independent_validation") is not True:
            validation_errors.append("independent_validation must be true")
        for key in ("ic_replay_certificate_sha256", "rho_replay_certificate_sha256"):
            digest = report.get(key)
            if not isinstance(digest, str) or len(digest) != 64 or any(
                ch not in "0123456789abcdef" for ch in digest
            ):
                validation_errors.append(f"{key} must be a full lowercase SHA-256 digest")
        ic_resources = report.get("ic_resource_envelope")
        rho_resources = report.get("rho_resource_envelope")
        if not isinstance(ic_resources, dict) or not ic_resources:
            validation_errors.append("ic_resource_envelope must be a nonempty object")
        if not isinstance(rho_resources, dict) or not rho_resources:
            validation_errors.append("rho_resource_envelope must be a nonempty object")
        if isinstance(ic_resources, dict) and isinstance(rho_resources, dict):
            if ic_resources != rho_resources:
                validation_errors.append("IC and rho resource envelopes must match exactly")

        ic_ms = positive_cost(report.get("ic_online_wall_ms"))
        rho_ms = positive_cost(report.get("rho_online_wall_ms"))
        speedup = positive_cost(report.get("online_speedup"))
        if ic_ms is None or rho_ms is None or speedup is None:
            validation_errors.append("online times and speedup must be positive finite numbers")
        elif not math.isclose(speedup, rho_ms / ic_ms, rel_tol=1e-9, abs_tol=1e-12):
            validation_errors.append("online_speedup must equal rho_online_wall_ms / ic_online_wall_ms")
        for arm, wall_ms in (("ic", ic_ms), ("rho", rho_ms)):
            for alias in (f"{arm}_cost", f"{arm}_online_ms"):
                alias_ms = positive_cost(report.get(alias))
                if alias_ms is None or wall_ms is None or not math.isclose(
                    alias_ms, wall_ms, rel_tol=1e-9, abs_tol=1e-12
                ):
                    validation_errors.append(f"{alias} must equal {arm}_online_wall_ms")

        discount = report.get("automorphism_discount")
        policy = report.get("rho_policy")
        if not isinstance(discount, dict) or type(discount.get("A")) is not int or discount["A"] < 1:
            validation_errors.append("automorphism_discount.A must be a positive integer")
        elif isinstance(policy, dict) and discount["A"] != policy.get("automorphism_size"):
            validation_errors.append("automorphism_discount.A must match rho_policy.automorphism_size")

        phase_costs = report.get("ic_online_phase_ms")
        phase_fields = stage_schema.get("ic_online_phase_fields")
        if phase_fields is None:
            primary_keys = stage_schema.get("primary_ic_online_phase_keys", [])
            phase_fields = [f"T_{key}_ms" for key in primary_keys]
        if not phase_fields:
            validation_errors.append("ledger must specify IC online phase fields")
        if not isinstance(phase_costs, dict):
            validation_errors.append("ic_online_phase_ms must be an object")
        else:
            unexpected = set(phase_costs) - set(phase_fields)
            if unexpected:
                validation_errors.append(
                    f"ic_online_phase_ms has unexpected fields: {sorted(unexpected)}"
                )
            phase_values = []
            for key in phase_fields:
                value = phase_costs.get(key)
                if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                    validation_errors.append(f"ic_online_phase_ms.{key} must be a finite nonnegative number")
                else:
                    phase_values.append(float(value))
            if len(phase_values) == len(phase_fields) and ic_ms is not None:
                if not math.isclose(math.fsum(phase_values), ic_ms, rel_tol=1e-6, abs_tol=1e-3):
                    validation_errors.append("IC exclusive phase costs must sum to ic_online_wall_ms")

        rho_phases = report.get("rho_online_phase_ms")
        rho_phase_fields = stage_schema.get("rho_online_phase_fields", [])
        if not isinstance(rho_phases, dict):
            validation_errors.append("rho_online_phase_ms must be an object")
        elif rho_phase_fields:
            if set(rho_phases) != set(rho_phase_fields):
                validation_errors.append("rho_online_phase_ms must contain exactly the declared phases")
            elif all(type(rho_phases[key]) in (int, float) and math.isfinite(rho_phases[key])
                     and rho_phases[key] >= 0 for key in rho_phase_fields) and rho_ms is not None:
                if not math.isclose(math.fsum(float(rho_phases[key]) for key in rho_phase_fields),
                                    rho_ms, rel_tol=1e-6, abs_tol=1e-3):
                    validation_errors.append("rho exclusive phase costs must sum to rho_online_wall_ms")
            else:
                validation_errors.append("rho_online_phase_ms must contain finite nonnegative costs")

        interval = report.get("online_interval")
        interval_fields = (
            stage_schema.get("online_interval_event_fields", [])
            + stage_schema.get("online_interval_stage_fields", [])
        )
        if not isinstance(interval, dict):
            validation_errors.append("online_interval must be an object")
        else:
            for key in interval_fields:
                value = interval.get(key)
                valid = (
                    isinstance(value, str) and bool(value.strip())
                    if key.endswith("_event")
                    else isinstance(value, list) and bool(value)
                    and all(isinstance(item, str) and item.strip() for item in value)
                )
                if not valid:
                    validation_errors.append(f"online_interval.{key} is missing or invalid")
            for included_key, required_stages in ONLINE_REQUIRED_STAGES.items():
                included = interval.get(included_key)
                if isinstance(included, list):
                    missing = sorted(set(required_stages) - set(included))
                    if missing:
                        validation_errors.append(
                            f"online_interval.{included_key} is missing required stages: {', '.join(missing)}"
                        )

        policy = report.get("rho_policy")
        if not isinstance(policy, dict):
            validation_errors.append("rho_policy must be an object")
        else:
            for key in stage_schema.get("rho_policy_integer_fields", []):
                value = policy.get(key)
                minimum = stage_schema.get("rho_policy_minimums", {}).get(key)
                if type(value) is not int or (minimum is not None and value < minimum):
                    validation_errors.append(f"rho_policy.{key} must be an integer >= {minimum}")
            for key in stage_schema.get("rho_policy_string_fields", []):
                value = policy.get(key)
                if not isinstance(value, str) or not value.strip():
                    validation_errors.append(f"rho_policy.{key} must be a nonempty string")

    ok = not missing_stage and not missing_global and not pairing_errors and not validation_errors
    return {
        "schema_version": ledger.get("schema_version"),
        "stage": stage,
        "fail_closed": True,
        "status": "PASS" if ok else "FAIL",
        "missing_stage_fields": missing_stage,
        "missing_global_provenance": missing_global,
        "pairing_errors": pairing_errors,
        "validation_errors": validation_errors,
        "required_stage_fields": required,
        "required_global_provenance": global_required,
    }


class RunnerLock:
    """Exclusive advisory lock using an atomic create + live-pid check."""

    def __enter__(self) -> "RunnerLock":
        RUNS_DIR.mkdir(parents=True, exist_ok=True)
        while True:
            try:
                fd = os.open(LOCK_PATH, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except FileExistsError:
                try:
                    record = read_json(LOCK_PATH)
                except Exception:
                    record = {}
                pid = int(record.get("pid") or 0)
                if pid and _pid_alive(pid):
                    raise AutolabError(f"autolab lock is held by pid {pid}")
                LOCK_PATH.unlink(missing_ok=True)
                continue
            payload = {"pid": os.getpid(), "created_at": now(), "task_id": TASK_ID}
            os.write(fd, json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n")
            os.close(fd)
            return self

    def __exit__(self, *_: Any) -> None:
        try:
            record = read_json(LOCK_PATH)
        except Exception:
            record = {}
        if record.get("pid") == os.getpid():
            LOCK_PATH.unlink(missing_ok=True)


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def make_run_id(beat_id: str) -> str:
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    digest = hashlib.sha256(f"{TASK_ID}|{beat_id}|{stamp}".encode()).hexdigest()[:10]
    return f"{stamp}-{digest}"


def current_run_id() -> str:
    require(CURRENT_PATH.is_file(), "no current autolab run")
    run_id = read_json(CURRENT_PATH).get("run_id")
    require(bool(run_id), "current.json missing run_id")
    return str(run_id)


def resolve_run(run_id: str | None) -> Path:
    resolved = run_id or current_run_id()
    path = RUNS_DIR / resolved
    require(path.is_dir(), f"run directory missing: {path}")
    return path


def host_record() -> dict[str, Any]:
    return {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "machine": platform.machine(),
        "node": platform.node(),
        "cpu_count": os.cpu_count(),
    }


def preflight(protocol: dict[str, Any], ledger: dict[str, Any]) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []

    def add(name: str, ok: bool, detail: str) -> None:
        checks.append({"name": name, "ok": ok, "detail": detail})

    add(
        "ledger_schema_v2",
        ledger.get("schema_version") == 2 and "measurement_schema" in ledger,
        f"schema_version={ledger.get('schema_version')}",
    )
    cargo = shutil.which("cargo")
    add("cargo", cargo is not None, cargo or "cargo not on PATH")
    rustc = shutil.which("rustc")
    add("rustc", rustc is not None, rustc or "rustc not on PATH")
    cms = Path(os.environ.get("KIC_AUTOLAB_CMS", "/opt/homebrew/bin/cryptominisat5"))
    add(
        "cryptominisat5_optional",
        cms.is_file() or shutil.which("cryptominisat5") is not None,
        str(cms if cms.is_file() else shutil.which("cryptominisat5") or "absent"),
    )
    for key, producer in protocol["producers"].items():
        source = REPO / producer["source"]
        add(f"producer_source_{key}", source.is_file(), str(source))
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    add("git_head", git_head.returncode == 0, git_head.stdout.strip())
    dirty = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    add(
        "git_status_recorded",
        dirty.returncode == 0,
        f"{len(dirty.stdout.splitlines())} dirty paths",
    )
    ok = all(
        check["ok"]
        for check in checks
        if check["name"] != "cryptominisat5_optional"
    )
    return {
        "schema_version": "1.0",
        "task_id": TASK_ID,
        "kind": "local_autolab_preflight",
        "ok": ok,
        "host": host_record(),
        "ledger_sha256": sha256(ledger_path(protocol)),
        "protocol_sha256": sha256(PROTOCOL_PATH),
        "agent_priorities": ledger.get("agent_priorities", []),
        "checks": checks,
        "claim_boundary": (
            "preflight only; not a fixed-arm cost, relation, rank, memory, rho, "
            "or crossover result"
        ),
        "created_at": now(),
    }


def plan(protocol: dict[str, Any], ledger: dict[str, Any]) -> dict[str, Any]:
    contract = protocol.get("primary_speedup_contract", {})
    contract_ready = contract.get("status") == "ready" and contract.get("target_count") == 1
    beats = []
    for beat_id, beat in sorted(
        protocol["beats"].items(),
        key=lambda item: (item[1].get("priority", 99), item[0]),
    ):
        eligible = contract_ready and bool(beat.get("primary_speedup_eligible")) and (
            beat.get("timing_class_goal") == "single_target_online_wall"
        )
        beats.append(
            {
                "beat_id": beat_id,
                "priority": beat.get("priority"),
                "label": beat.get("label"),
                "regime": beat.get("regime"),
                "stage": beat.get("stage"),
                "n": beat.get("n"),
                "timing_class_goal": beat.get("timing_class_goal"),
                "primary_speedup_eligible": eligible,
                "primary_speedup_blocker": None if eligible else (
                    beat.get("primary_speedup_blocker")
                    or protocol.get("primary_speedup_contract", {}).get("blocker")
                    or "beat is diagnostic-only"
                ),
            }
        )
    return {
        "schema_version": "1.0",
        "task_id": TASK_ID,
        "ledger": str(protocol["ledger"]["path"]),
        "ledger_schema_version": ledger.get("schema_version"),
        "primary_speedup_contract": contract,
        "agent_priorities": ledger.get("agent_priorities", []),
        "beats": beats,
        "how_to_beat": ledger.get("how_to_beat", []),
        "commands": {
            "preflight": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py preflight",
            "smoke": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py launch --beat smoke.koblitz.vs_rho.n13",
            "n37_wall": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py launch --beat koblitz.vs_rho.n37_wall",
            "n41_charged": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py launch --beat koblitz.vs_rho.n41_charged",
            "n53_factor_base": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py launch --beat koblitz.factor_base.n53",
            "n37_single_target": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py launch --beat koblitz.vs_rho.n37_wall --fixtures 1",
            "n61_bl_reference": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py launch-equal-advice --beat koblitz.compact_orbit.n61_bl_reference --targets 1024",
            "claim_check": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py claim-check --report PATH --stage vs_rho",
        },
    }


def build_producers() -> dict[str, str]:
    examples = ["koblitz_rank_fixture", "koblitz_rho_fixture"]
    for example in examples:
        source = REPO / "examples" / f"{example}.rs"
        require(source.is_file(), f"missing producer source: {source}")
    command = [
        "cargo",
        "build",
        "--release",
        "--example",
        "koblitz_rank_fixture",
        "--example",
        "koblitz_rho_fixture",
    ]
    completed = subprocess.run(command, cwd=REPO, capture_output=True, text=True)
    require(
        completed.returncode == 0,
        "cargo build failed:\n" + completed.stderr[-4000:],
    )
    binaries = {
        "direct": str((REPO / "target/release/examples/koblitz_rank_fixture").resolve()),
        "rho": str((REPO / "target/release/examples/koblitz_rho_fixture").resolve()),
    }
    for path in binaries.values():
        require(Path(path).is_file(), f"built binary missing: {path}")
    return binaries


def build_example(example: str) -> str:
    source = REPO / "examples" / f"{example}.rs"
    require(source.is_file(), f"missing producer source: {source}")
    completed = subprocess.run(
        ["cargo", "build", "--release", "--example", example],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    require(
        completed.returncode == 0,
        f"cargo build --example {example} failed:\n" + completed.stderr[-4000:],
    )
    binary = (REPO / "target/release/examples" / example).resolve()
    require(binary.is_file(), f"built binary missing: {binary}")
    return str(binary)


RESOURCE_RECEIPT_FIELDS = (
    "exit_code",
    "whole_process_wall_ms",
    "whole_process_wall_samples_ms",
    "cold_start_wall_ms",
    "timed_repeats",
    "stdout_stable",
    "children_cpu_user_ms",
    "children_cpu_system_ms",
    "children_peak_rss_bytes",
    "command",
    "seed",
)

DEFAULT_TIMED_REPEATS = 3
# Past this, further repeats cost more than the spread they resolve, so stop
# after the first timed execution.
REPEAT_BUDGET_MS = 20_000.0


def _run_once(
    command: list[str], *, env: dict[str, str], cwd: Path
) -> tuple[subprocess.CompletedProcess[str], float, float, float]:
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
    )
    elapsed_ms = (time.perf_counter() - started) * 1000.0
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    # Children CPU deltas (seconds -> ms). Wall is the outer process wait.
    cpu_user_ms = (usage_after.ru_utime - usage_before.ru_utime) * 1000.0
    cpu_system_ms = (usage_after.ru_stime - usage_before.ru_stime) * 1000.0
    # Peak RSS after wait is cumulative for this process's children; treat as
    # an observed upper bound for single-producer launches.
    peak_rss = children_rss_bytes(int(usage_after.ru_maxrss))
    return completed, elapsed_ms, cpu_user_ms, cpu_system_ms


def _timing_free(text: str) -> str:
    """Producer stdout with its own timing fields dropped.

    Every row carries measured durations, so raw stdout never repeats exactly.
    What should repeat is the computation: base hash, rank, status, solution.
    """
    rows = []
    for row in parse_json_lines(text):
        rows.append(
            {
                k: v
                for k, v in row.items()
                if not k.endswith("_ms") and k != "timing_breakdown_ms"
            }
        )
    return json.dumps(rows, sort_keys=True)


def run_timed(
    command: list[str],
    *,
    env: dict[str, str],
    cwd: Path,
    repeats: int = DEFAULT_TIMED_REPEATS,
) -> dict[str, Any]:
    """Time a producer without charging it for first-execution cost.

    The first execution of a freshly linked binary pays page-in and, on macOS,
    signature validation: measured here at 183-484 ms against 7 ms warm, which
    is 30x the whole n=13 comparison and larger than either arm at n=37. Timing
    a single run therefore adds a roughly constant per-process term to both
    arms, which drags every ratio toward parity and flatters whichever arm is
    slower.

    The warmup runs the binary with no arguments, which both producers reject
    immediately. That still pays the whole load cost -- after it the real run
    drops from 190 ms to 7.7 ms -- so a rung that takes a quarter of an hour is
    not run twice to save 200 ms.

    The producers take their seed in argv and are deterministic, so repeats are
    the same computation; `stdout_stable` records whether that held, comparing
    the rows with their own timing fields dropped.
    """
    require(repeats >= 1, f"repeats must be >= 1, got {repeats}")
    _, cold_wall_ms, _, _ = _run_once([command[0]], env=env, cwd=cwd)

    walls: list[float] = []
    users: list[float] = []
    systems: list[float] = []
    outputs: list[str] = []
    completed: subprocess.CompletedProcess[str] | None = None
    for _ in range(repeats):
        completed, wall_ms, user_ms, system_ms = _run_once(command, env=env, cwd=cwd)
        walls.append(wall_ms)
        users.append(user_ms)
        systems.append(system_ms)
        outputs.append(_timing_free(completed.stdout))
        if completed.returncode != 0 or wall_ms > REPEAT_BUDGET_MS:
            break

    return {
        "command": command,
        "exit_code": completed.returncode,
        "whole_process_wall_ms": statistics.median(walls),
        "whole_process_wall_samples_ms": walls,
        "cold_start_wall_ms": cold_wall_ms,
        "timed_repeats": len(walls),
        "stdout_stable": len(set(outputs)) <= 1,
        "children_cpu_user_ms": statistics.median(users),
        "children_cpu_system_ms": statistics.median(systems),
        "children_peak_rss_bytes": children_rss_bytes(
            int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
        ),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def parse_json_lines(text: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            rows.append(value)
    return rows


def parse_json_objects(text: str) -> list[dict[str, Any]]:
    """Parse concatenated / pretty-printed JSON objects from producer stdout."""
    rows = parse_json_lines(text)
    if rows:
        return rows
    rows = []
    decoder = json.JSONDecoder()
    index = 0
    while index < len(text):
        while index < len(text) and text[index].isspace():
            index += 1
        if index >= len(text):
            break
        try:
            value, end = decoder.raw_decode(text, index)
        except json.JSONDecodeError:
            break
        if isinstance(value, dict):
            rows.append(value)
        index = end
    return rows


def children_rss_bytes(ru_maxrss: int) -> int:
    # Darwin reports bytes; Linux reports kilobytes.
    if sys.platform == "darwin":
        return int(ru_maxrss)
    return int(ru_maxrss) * 1024


def seed_for(beat_id: str, arm: str, repetition: int) -> int:
    material = f"{TASK_ID}|{beat_id}|arm={arm}|repetition={repetition}".encode()
    return int.from_bytes(hashlib.sha256(material).digest()[:8], "big")


def exclusive_online_timing_ok(row: dict[str, Any] | None) -> bool:
    """Reject a phase sum that re-adds nested LA/replay to collection_ms."""
    if not isinstance(row, dict):
        return False
    phases = row.get("target_online_phase_ms")
    if not isinstance(phases, dict):
        return False
    keys = ("fixture_setup_ms", "collection_ms", "linear_solve_ms",
            "solution_validation_ms", "reference_validation_ms", "target_online_wall_ms")
    values = {key: row.get(key) for key in keys}
    if any(isinstance(value, bool) or not isinstance(value, (int, float))
           or not math.isfinite(value) or value < 0 for value in values.values()):
        return False
    breakdown = row.get("timing_breakdown_ms")
    if not isinstance(breakdown, dict):
        return False
    target_generation_ms = breakdown.get("target_generation")
    packed_verification_ms = breakdown.get("packed_verification")
    if any(isinstance(value, bool) or not isinstance(value, (int, float))
           or not math.isfinite(value) or value < 0
           for value in (target_generation_ms, packed_verification_ms)):
        return False
    expected_pdp = (values["collection_ms"] - values["linear_solve_ms"]
                    - values["solution_validation_ms"] - target_generation_ms
                    - packed_verification_ms)
    tolerance = max(1e-6, values["target_online_wall_ms"] * 1e-8)
    if expected_pdp < -tolerance:
        return False
    expected = {
        "target_query": values["fixture_setup_ms"] + target_generation_ms,
        "target_pdp": expected_pdp,
        "target_relation_check": values["reference_validation_ms"]
                                 + packed_verification_ms,
        "target_descent": 0.0,
        "target_recovery_check": values["linear_solve_ms"] + values["solution_validation_ms"],
    }
    if set(phases) != set(expected):
        return False
    if any(isinstance(phases[key], bool) or not isinstance(phases[key], (int, float))
           or not math.isfinite(phases[key]) or phases[key] < 0
           or abs(phases[key] - value) > tolerance for key, value in expected.items()):
        return False
    charged = (values["fixture_setup_ms"] + values["collection_ms"]
               + values["reference_validation_ms"])
    return abs(charged - values["target_online_wall_ms"]) <= tolerance


def exclusive_precomputation_ms(row: dict[str, Any] | None) -> float | None:
    """Collection already contains final LA and solution validation."""
    if not isinstance(row, dict):
        return None
    keys = ("setup_ms", "fixture_setup_ms", "collection_ms", "reference_validation_ms")
    values = [row.get(key) for key in keys]
    if any(isinstance(value, bool) or not isinstance(value, (int, float))
           or not math.isfinite(value) or value < 0 for value in values):
        return None
    return sum(values)


def positive_cost(value):
    if type(value) not in (int, float) or value <= 0:
        return None
    try:
        result = float(value)
    except (OverflowError, ValueError):
        return None
    return result if math.isfinite(result) else None


def extract_ic_cost(rows: list[dict[str, Any]], timing_class: str) -> float | None:
    if not rows:
        return None
    row = rows[-1]
    if timing_class == "whole_process_wall":
        return None  # filled from outer wait4
    # Online-only and projection-matched slices cannot stand in for a full job.
    for key in ("full_algorithm_charged_total_ms",):
        if key in row:
            return positive_cost(row[key])
        timing = row.get("timing_breakdown_ms")
        if isinstance(timing, dict) and key in timing:
                return positive_cost(timing[key])
    return None


def extract_rho_cost(rows: list[dict[str, Any]]) -> float | None:
    """Return a cost only for exactly one explicit target result row."""
    targets = [row for row in rows if row.get("kind") == "rho_public_fixture"]
    if len(targets) != 1:
        return None
    row = targets[0]
    if "total_ms" in row:
        return positive_cost(row["total_ms"])
    timing = row.get("timing_breakdown_ms")
    if isinstance(timing, dict) and "total_ms" in timing:
        return positive_cost(timing["total_ms"])
    return None


def comparison_integrity(direct_rows, rho_rows, expected_count=None):
    """Bind a comparison to actual complete public targets, never source bytes.

    This checks producer consistency only. Independent group/rank replay is
    still required for a scientific result.
    """
    try:
        require(expected_count == 1, "primary comparison requires exactly one target")
    except AutolabError as exc:
        return {"status": "INVALID_COMPARISON", "reason": str(exc), "fixture_hash": None,
                "independent_validation": False}

    def corpus(rows, kind, verified_key):
        selected = [r for r in rows if r.get('kind') == kind]
        require(bool(selected), 'no complete fixture records')
        headers = [r for r in rows if r.get('kind') == 'point_defined_factor_base']
        if kind == 'relation_rank_summary':
            require(len(headers) == 1, 'missing or ambiguous curve header')
        entries = []
        for row in selected:
            metadata = headers[0] if kind == 'relation_rank_summary' else row
            require(row.get(verified_key) is True, 'unverified fixture')
            require(row.get('recovered_fixture_scalar') is not None, 'missing recovered scalar')
            require(int(row['recovered_fixture_scalar']) == int(row['published_fixture_scalar']), 'wrong scalar')
            entries.append({'index': int(row['fixture_index']), 'n': int(row['n']), 'a': int(row['a']),
                'subgroup_order': int(metadata['subgroup_order']),
                'field_modulus_low_terms': sorted(int(t) for t in metadata['field_modulus_low_terms']),
                'generator': [int(x) for x in row['generator_point_key']],
                'target': [int(x) for x in row['published_q_point_key']]})
        entries.sort(key=lambda r: r['index'])
        count = len(entries) if expected_count is None else expected_count
        require([r['index'] for r in entries] == list(range(count)), 'missing or duplicate fixture')
        return entries
    try:
        direct = corpus(direct_rows, 'relation_rank_summary', 'linear_solution_verified')
        rho = corpus(rho_rows, 'rho_public_fixture', 'verified')
        require(direct == rho, 'IC and rho solved different public targets')
        fingerprint = hashlib.sha256(json.dumps(direct, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        return {'status': 'MATCHED', 'fixture_hash': fingerprint, 'fixtures': len(direct),
                'independent_validation': False}
    except (AutolabError, KeyError, TypeError, ValueError) as exc:
        return {'status': 'INVALID_COMPARISON', 'reason': str(exc), 'fixture_hash': None,
                'independent_validation': False}


def automorphism_discount(n: int) -> dict[str, Any]:
    return {
        "formula": "sqrt(2*n)",
        "A": 2 * n,
        "description": "signed Frobenius automorphism discount family used by Koblitz rho control",
        "n": n,
    }


def public_target_hash(beat: dict[str, Any], point: Any) -> str | None:
    """Bind a producer's public point to the beat's field and curve parameters."""
    if not isinstance(point, list) or len(point) != 2:
        return None
    try:
        coordinates = []
        for value in point:
            if type(value) is int and value >= 0:
                coordinates.append(str(value))
            elif isinstance(value, str) and value.isascii() and value.isdecimal():
                coordinates.append(str(int(value)))
            else:
                return None
        payload = {"n": int(beat["n"]), "a": int(beat["a"]), "q": coordinates}
    except (KeyError, TypeError, ValueError):
        return None
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(canonical).hexdigest()


def draft_vs_rho_claim(
    *,
    beat: dict[str, Any],
    beat_id: str,
    run: Path,
    direct_obs: dict[str, Any],
    rho_obs: dict[str, Any],
    binaries: dict[str, str],
) -> dict[str, Any]:
    timing_class = beat["timing_class_goal"]
    direct_rows = parse_json_objects(direct_obs["stdout"])
    rho_rows = parse_json_objects(rho_obs["stdout"])
    direct_summaries = [
        row for row in direct_rows if row.get("kind") == "relation_rank_summary"
    ]
    precompute_row = next(
        (row for row in direct_summaries if row.get("precomputation_fixture")), None
    )
    online_row = next(
        (row for row in reversed(direct_summaries)
         if row.get("online_target_count") == 1), None
    )
    rho_row = next(
        (row for row in reversed(rho_rows) if row.get("kind") == "rho_public_fixture"),
        None,
    )
    ic_q = online_row.get("published_q") if online_row else None
    rho_q = rho_row.get("published_q") if rho_row else None
    ic_target_hash = public_target_hash(beat, ic_q)
    rho_target_hash = public_target_hash(beat, rho_q)
    same_public_point = ic_target_hash is not None and ic_target_hash == rho_target_hash
    ic_scalar = online_row.get("published_fixture_scalar") if online_row else None
    rho_scalar = rho_row.get("published_fixture_scalar") if rho_row else None
    same_fixture_scalar = ic_scalar is not None and ic_scalar == rho_scalar

    ic_phases = online_row.get("target_online_phase_ms") if online_row else None
    ic_online = online_row.get("target_online_wall_ms") if online_row else None
    ic_phase_exclusive_ok = exclusive_online_timing_ok(online_row)
    ic_phase_sum_ok = False
    if ic_phase_exclusive_ok:
        phase_sum = math.fsum(float(value) for value in ic_phases.values())
        ic_phase_sum_ok = math.isclose(
            phase_sum, float(ic_online), rel_tol=1e-8, abs_tol=0.02
        )
    rho_phases = None
    rho_online = None
    if rho_row is not None:
        walk_ms = rho_row.get("walk_ms")
        validation_ms = rho_row.get("validation_ms")
        if all(
            type(value) in (int, float) and math.isfinite(value) and value >= 0
            for value in (walk_ms, validation_ms)
        ):
            rho_phases = {
                "target_walk": float(walk_ms),
                "scalar_recovery_check": float(validation_ms),
            }
            rho_online = sum(rho_phases.values())

    ic_verified = bool(
        online_row
        and online_row.get("status") == "SHARED_FACTOR_LOG_ONE_RELATION"
        and online_row.get("uses_retained_factor_logs") is True
        and online_row.get("online_target_count") == 1
        and online_row.get("linear_solution_verified") is True
        and online_row.get("recovered_fixture_scalar") == ic_scalar
        and ic_phase_sum_ok
        and ic_phase_exclusive_ok
    )
    rho_verified = bool(
        rho_row
        and rho_row.get("verified") is True
        and rho_row.get("recovered_fixture_scalar") == rho_scalar
        and rho_online is not None
    )
    direct_rows = parse_json_lines(direct_obs["stdout"])
    rho_rows = parse_json_lines(rho_obs["stdout"])
    producers_ok = direct_obs["exit_code"] == 0 and rho_obs["exit_code"] == 0
    integrity = comparison_integrity(direct_rows, rho_rows, direct_obs.get('fixtures'))
    ic_cost = (
        float(direct_obs["whole_process_wall_ms"])
        if timing_class == "whole_process_wall"
        else extract_ic_cost(direct_rows, timing_class)
        and positive_cost(rho_online) is not None
    )
    target_count = int(online_row.get("online_target_count", 0)) if online_row else 0
    pairing_ok = same_public_point and same_fixture_scalar and target_count == 1
    stable_outputs = all(
        observation.get("stdout_stable", True) is True
        for observation in (direct_obs, rho_obs)
    )
    online_speedup = (
        float(rho_online) / float(ic_online)
        if pairing_ok and ic_verified and rho_verified and stable_outputs
        and positive_cost(ic_online) is not None
        and positive_cost(rho_online) is not None
        else None
    )
    producers_ok = (
        direct_obs["exit_code"] == 0 and rho_obs["exit_code"] == 0
        and pairing_ok and ic_verified and rho_verified and stable_outputs
    )
    precompute_ms = exclusive_precomputation_ms(precompute_row)
    if timing_class in ("single_target_online", "single_target_online_wall"):
        ic_cost = positive_cost(ic_online)
        rho_cost = positive_cost(rho_online)
    else:
        ic_cost = (
            float(direct_obs["whole_process_wall_ms"])
            if timing_class == "whole_process_wall"
            else extract_ic_cost(direct_rows, timing_class)
        )
        rho_cost = (
            float(rho_obs["whole_process_wall_ms"])
            if timing_class == "whole_process_wall"
            else extract_rho_cost(rho_rows)
        )
    ic_operations = online_row.get("target_trials") if online_row else None
    rho_operations = rho_row.get("walk_steps") if rho_row else None
    operation_accounting = None
    if all(
        isinstance(value, (int, float)) and not isinstance(value, bool) and value > 0
        for value in (ic_operations, rho_operations)
    ):
        operation_accounting = {
            "schema_version": "1.0",
            "unit_assumption": "IC target trials and rho walk steps are distinct uncalibrated counters",
            "comparison_status": "native_counters_only",
            "operation_units": {"ic": "target relation trials", "rho": "walk steps"},
            "ic_online_operations": ic_operations,
            "ic_online_operations_basis": "target_trials of the one online target",
            "rho_online_operations": rho_operations,
            "rho_online_operations_basis": "measured walk steps of this run",
            "rho_per_ic_native_counter": float(rho_operations) / float(ic_operations),
            "ops_speedup_online": None,
        }
    claim = {
        "schema_version": 2,
        "task_id": TASK_ID,
        "beat_id": beat_id,
        "regime": beat["regime"],
        "stage": "vs_rho",
        "n": beat["n"],
        "n_or_bits": beat["n"],
        "timing_class": timing_class,
        "candidate_id": None,
        "candidate_manifest_sha256": None,
        "workload_id": None,
        "workload_manifest_sha256": None,
        "run_id": None,
        "autolab_run_id": run.name,
        "target_count": target_count,
        "ic_target_hash": ic_target_hash,
        "rho_target_hash": rho_target_hash,
        "ic_online_wall_ms": positive_cost(ic_online),
        "rho_online_wall_ms": positive_cost(rho_online),
        "online_interval": None,
        "same_resource_envelope": None,
        "independent_validation": False,
        "ic_replay_certificate_sha256": None,
        "rho_replay_certificate_sha256": None,
        "ic_resource_envelope": None,
        "rho_resource_envelope": None,
        "ic_scalar_verified": ic_verified,
        "rho_scalar_verified": rho_verified,
        "rho_policy": None,
        "ic_cost": ic_cost,
        "rho_cost": rho_cost,
        "automorphism_discount": automorphism_discount(int(beat["n"])),
        "all_stages_charged_same_series": bool(
            pairing_ok and ic_verified and rho_verified and ic_phase_sum_ok
            and ic_phase_exclusive_ok and stable_outputs
        ),
        "verdict": (
            "DRAFT_PENDING_INDEPENDENT_VALIDATION"
            if producers_ok
            else "PRODUCER_FAILURE"
        ),
        "claim_boundary": (
            "Public synthetic Koblitz, one-target online comparison after reusable "
            "IC setup. No batch throughput, key recovery, asymptotic sub-rho claim, "
            "or ledger promotion before independent validation."
        ),
        "claim_boundary_non_claims": [
            "not key recovery",
            "not asymptotic sub-sqrt",
            "not imported/external points",
            "not a multi-target batch",
            "not ledger promotion until independent validation",
            "not compared against Pollard rho with precomputation at equal "
            "precompute and memory",
        ],
        "paired_target": {
            "ic_public_q": ic_q,
            "rho_public_q": rho_q,
            "same_public_point": same_public_point,
            "ic_fixture_scalar": ic_scalar,
            "rho_fixture_scalar": rho_scalar,
            "same_fixture_scalar": same_fixture_scalar,
            "point_generation_excluded_from_online": bool(
                online_row and rho_row
                and "fixture_generation_ms" in online_row
                and "target_generation_ms_excluded" in rho_row
            ),
        },
        "ic_verified": ic_verified,
        "ic_exclusive_timing_verified": ic_phase_exclusive_ok,
        "rho_verified": rho_verified,
        "ic_online_ms": float(ic_online) if isinstance(ic_online, (int, float)) else None,
        "rho_online_ms": float(rho_online) if isinstance(rho_online, (int, float)) else None,
        "online_speedup": online_speedup,
        "operation_accounting": operation_accounting,
        "ic_online_phase_ms": ic_phases,
        "rho_online_phase_ms": rho_phases,
        "ic_online_interval": online_row.get("target_online_interval") if online_row else None,
        "rho_online_interval": (
            "first target-dependent walk step through recovered scalar and verification; "
            "excludes launch, walk setup, and fixture point generation"
            if rho_row else None
        ),
        "ic_precomputation_fixture_ms": precompute_ms,
        "ic_target_relation_attempts": online_row.get("target_trials") if online_row else None,
        "independent_replay_pointer": independent_replay_pointer(run),
        "fixture_hash": ic_target_hash if same_public_point else None,
        "executable_or_source_hash": {
            "direct": sha256(binaries["direct"]),
            "rho": sha256(binaries["rho"]),
            "rank_fixture_source": sha256(REPO / "examples/koblitz_rank_fixture.rs"),
            "rho_fixture_source": sha256(REPO / "examples/koblitz_rho_fixture.rs"),
        },
        "host_id": host_record(),
        "resource_caps": {"common_cap_bytes": beat.get("resource_cap_bytes")},
        "seeds": {
            "direct": direct_obs.get("seed"),
            "rho": rho_obs.get("seed"),
        },
        "producer_exit_codes": {
            "direct": direct_obs["exit_code"],
            "rho": rho_obs["exit_code"],
        },
        "whole_process_wall_ms": {
            "direct": direct_obs["whole_process_wall_ms"],
            "rho": rho_obs["whole_process_wall_ms"],
        },
        "whole_process_wall_method": {
            "description": "median of timed executions after one discarded warmup; first-execution cost is reported separately",
            "cold_start_wall_ms": {
                "direct": direct_obs.get("cold_start_wall_ms"),
                "rho": rho_obs.get("cold_start_wall_ms"),
            },
            "samples_ms": {
                "direct": direct_obs.get("whole_process_wall_samples_ms"),
                "rho": rho_obs.get("whole_process_wall_samples_ms"),
            },
            "timed_repeats": {
                "direct": direct_obs.get("timed_repeats"),
                "rho": rho_obs.get("timed_repeats"),
            },
            "stdout_stable": {
                "direct": direct_obs.get("stdout_stable"),
                "rho": rho_obs.get("stdout_stable"),
            },
        },
        "direct_rows": len(direct_rows),
        "rho_rows": len(rho_rows),
        "direct_target_summary_rows": len(direct_summaries),
    }
    return claim


def select_factor_base_row(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    for row in rows:
        if row.get("kind") == "point_defined_factor_base":
            return row
        if row.get("evidence_class") == "measured_factor_base_construction":
            return row
    for row in rows:
        if "factor_base_points" in row or "orbit_columns" in row:
            return row
    return rows[0] if rows else None


def parse_paired_oracle_stdout(stdout: str) -> dict[str, Any]:
    """Parse symmetrised_oracle_bench / format_paired text into structured fields."""
    import re

    header = None
    for line in stdout.splitlines():
        m = re.match(
            r"== K_(\d+)/F_2\^(\d+)\s+m = (\d+)\s+dim V = (\d+)\s+"
            r"\|F_x\| = (\d+)\s+\|F_u\| = (\d+)\s+targets = (\d+)",
            line.strip(),
        )
        if m:
            header = {
                "a": int(m.group(1)),
                "n": int(m.group(2)),
                "m": int(m.group(3)),
                "dim": int(m.group(4)),
                "F_x": int(m.group(5)),
                "F_u": int(m.group(6)),
                "targets": int(m.group(7)),
            }
            break
    require(header is not None, "paired oracle stdout missing K_/F_2 header")

    arms: dict[str, Any] = {}
    divisor_indices: list[int] | None = None
    wall_s: float | None = None
    arm_map = {
        "x-chained F4": "x_chained_f4",
        "sym F4": "sym_f4",
        "x-direct F4": "x_direct_f4",
    }
    for line in stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("divisor indices"):
            dm = re.search(r"divisor indices\s+(\[.*?\])\s+\(([0-9.]+)\s*s\)", stripped)
            if dm:
                divisor_indices = json.loads(dm.group(1).replace("'", '"'))
                wall_s = float(dm.group(2))
            continue
        for label, key in arm_map.items():
            if not stripped.startswith(label):
                continue
            rest = stripped[len(label) :].split()
            # vars eqs deg found refuted inconc viaT found_ms refuted_ms effort [built] ffd gate...
            require(len(rest) >= 12, f"paired oracle arm row too short: {stripped}")
            # Current format_paired inserts a 'built' column before ffd.
            if len(rest) >= 13:
                vars_n, eqs, deg, found, refuted, inconc, via_t = map(int, rest[0:7])
                found_ms = float(rest[7]) if rest[7] != "-" else float("nan")
                refuted_ms = float(rest[8]) if rest[8] != "-" else float("nan")
                effort = float(rest[9]) if rest[9] != "-" else float("nan")
                ffd_tok = rest[11]
                gate = rest[12]
            else:
                vars_n, eqs, deg, found, refuted, inconc, via_t = map(int, rest[0:7])
                found_ms = float(rest[7]) if rest[7] != "-" else float("nan")
                refuted_ms = float(rest[8]) if rest[8] != "-" else float("nan")
                effort = float(rest[9]) if rest[9] != "-" else float("nan")
                ffd_tok = rest[10]
                gate = rest[11]
            ffd = None if ffd_tok == "-" else int(ffd_tok)
            arms[key] = {
                "vars": vars_n,
                "eqs": eqs,
                "deg": deg,
                "found": found,
                "refuted": refuted,
                "inconclusive": inconc,
                "via_t": via_t,
                "found_ms_median": found_ms,
                "refuted_ms_median": refuted_ms,
                "effort_splits": effort,
                "ffd": ffd,
                "gate": gate,
            }
            break

    require(arms, "paired oracle stdout missing arm rows")
    ffds = [arm["ffd"] for arm in arms.values() if isinstance(arm.get("ffd"), int)]
    # Honesty gate: max of class medians (found/refuted) across arms.
    class_medians = []
    for arm in arms.values():
        for key in ("found_ms_median", "refuted_ms_median"):
            val = arm.get(key)
            if isinstance(val, (int, float)) and val == val:  # not NaN
                class_medians.append(float(val))
    median_ms = max(class_medians) if class_medians else None
    disagreements = 0
    for arm in arms.values():
        gate = str(arm.get("gate") or "")
        if gate != "ok":
            # FAIL N
            parts = gate.split()
            if len(parts) == 2 and parts[0] == "FAIL" and parts[1].isdigit():
                disagreements += int(parts[1])
            else:
                disagreements += 1
    return {
        "header": header,
        "arms": arms,
        "divisor_indices": divisor_indices,
        "wall_s": wall_s,
        "ffds": ffds,
        "median_ms_per_target": median_ms,
        "disagreements": disagreements,
    }


def parse_semaev_decomp_bench_stdout(stdout: str) -> dict[str, Any]:
    """Parse examples/semaev_decomp_bench text ladder into structured fields."""
    import re

    ladder: list[dict[str, Any]] = []
    for line in stdout.splitlines():
        m = re.match(
            r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+([0-9.]+)s\s+([0-9.]+)s\s+([0-9.]+)s\s+([0-9.]+)x\s*$",
            line,
        )
        if not m:
            continue
        ladder.append(
            {
                "l": int(m.group(1)),
                "n": int(m.group(2)),
                "triples": int(m.group(3)),
                "triples_gen_s": float(m.group(4)),
                "triples_1wd_s": float(m.group(5)),
                "pairs_s": float(m.group(6)),
                "speedup_vs_triples": float(m.group(7)),
            }
        )
    require(ladder, "semaev_decomp_bench stdout missing pairs ladder rows")
    return {"ladder": ladder}


def parse_j0_known_answer_stdout(stdout: str) -> dict[str, Any]:
    """Parse examples/j0_known_answer_bits JSON result object from stdout."""
    rows = parse_json_objects(stdout)
    row = next(
        (
            item
            for item in reversed(rows)
            if isinstance(item.get("bits"), int)
            and "ic_ok" in item
            and "rho_ok" in item
        ),
        None,
    )
    require(row is not None, "j0_known_answer_bits stdout missing known-answer JSON row")
    return row


def draft_j0_end_to_end_dlp_claim(
    *,
    beat: dict[str, Any],
    beat_id: str,
    run: Path,
    direct_obs: dict[str, Any],
    binaries: dict[str, str],
) -> dict[str, Any]:
    parsed = parse_j0_known_answer_stdout(direct_obs.get("stdout") or "")
    bits = int(beat.get("bits") or parsed["bits"])
    require(int(parsed["bits"]) == bits, f"producer bits {parsed['bits']} != beat bits {bits}")
    producers_ok = (
        direct_obs["exit_code"] == 0
        and parsed.get("ic_ok") is True
        and parsed.get("rho_ok") is True
        and parsed.get("ic_agrees_rho") is True
    )
    peak_rss = direct_obs.get("children_peak_rss_bytes")
    cap = beat.get("resource_cap_bytes")
    under_cap = (
        peak_rss is not None and cap is not None and int(peak_rss) <= int(cap)
        if cap is not None
        else None
    )
    ic_ms = float(parsed["ic_ms"])
    rho_ms = float(parsed["rho_ms"])
    source = REPO / "examples" / f"{beat['producer']['example']}.rs"
    binary_path = binaries.get("custom") or binaries.get("direct")
    claim = {
        "schema_version": 2,
        "task_id": TASK_ID,
        "beat_id": beat_id,
        "regime": beat["regime"],
        "stage": "end_to_end_dlp",
        "n": bits,
        "n_or_bits": bits,
        "curve": parsed.get("curve"),
        "recovered_d_verified": bool(parsed.get("ic_ok")),
        "ic_ok": parsed.get("ic_ok"),
        "rho_ok": parsed.get("rho_ok"),
        "ic_matches_truth": parsed.get("ic_matches_truth"),
        "ic_agrees_rho": parsed.get("ic_agrees_rho"),
        "ic_ms": ic_ms,
        "rho_ms": rho_ms,
        "stage_timers": {
            "ic_ms": ic_ms,
            "rho_ms": rho_ms,
            "note": (
                "producer does not yet split FB/relations/LA; "
                "wall is full j0_index_calculus_dlp vs pollard_rho"
            ),
        },
        "target_orbits": parsed.get("target_orbits"),
        "max_trials": parsed.get("max_trials"),
        "rho_steps": parsed.get("rho_steps"),
        "wall_s": float(direct_obs.get("whole_process_wall_ms") or 0.0) / 1000.0,
        "peak_rss_bytes": peak_rss,
        "verdict": (
            "J0_16BIT_KA_IC_AGREES_RHO_PENDING_INDEPENDENT_VALIDATION"
            if producers_ok
            else "PRODUCER_FAILURE"
        ),
        "claim_boundary": (
            "Public synthetic known-answer j=0 IC vs rho only. "
            "Not key recovery, not asymptotic sub-rho, not ledger promotion."
        ),
        "claim_boundary_non_claims": [
            "not key recovery on deployed curves",
            "not asymptotic sub-rho",
            "not ledger promotion until independent validation",
            "IC slower than rho on this instance (expected; not a vs_rho crossover)",
        ],
        "independent_replay_pointer": independent_replay_pointer(run),
        "fixture_hash": sha256(source),
        "executable_or_source_hash": {
            "executable": sha256(binary_path) if binary_path else None,
            "source": sha256(source),
        },
        "host_id": host_record(),
        "resource_caps": {
            "common_cap_bytes": cap,
            "observed_peak_rss_bytes": peak_rss,
            "under_cap": under_cap,
        },
        "seeds": {"random_scalar": "runtime"},
        "evidence_paths": [
            str((run / "logs/direct.stdout.jsonl").resolve()),
            str((run / "receipts/direct.resource.json").resolve()),
        ],
        "created_at": now(),
        "autolab_run_id": run.name,
        "producer_exit_codes": {"direct": direct_obs["exit_code"]},
    }
    return claim


def draft_binary_pairs_decomposition_claim(
    *,
    beat: dict[str, Any],
    beat_id: str,
    run: Path,
    direct_obs: dict[str, Any],
    binaries: dict[str, str],
) -> dict[str, Any]:
    parsed = parse_semaev_decomp_bench_stdout(direct_obs.get("stdout") or "")
    ladder = parsed["ladder"]
    l_goal = int(beat.get("l") or 8)
    row = next((item for item in ladder if item["l"] == l_goal), None)
    require(row is not None, f"semaev_decomp_bench stdout missing l={l_goal} ladder row")
    producers_ok = direct_obs["exit_code"] == 0
    peak_rss = direct_obs.get("children_peak_rss_bytes")
    cap = beat.get("resource_cap_bytes")
    under_cap = (
        peak_rss is not None and cap is not None and int(peak_rss) <= int(cap)
        if cap is not None
        else None
    )
    pairs_s = float(row["pairs_s"])
    median_ms = pairs_s * 1000.0
    source = REPO / "examples" / f"{beat['producer']['example']}.rs"
    binary_path = binaries.get("custom") or binaries.get("direct")
    claim = {
        "schema_version": 2,
        "task_id": TASK_ID,
        "beat_id": beat_id,
        "regime": beat["regime"],
        "stage": "decomposition",
        "n": beat.get("n", row["n"]),
        "n_or_bits": beat.get("n", row["n"]),
        "m_summands": beat.get("m", 3),
        "dimension_l_or_dim": l_goal,
        "unknowns": "2*l_bit_pair_space",
        "system_degree": "quartic_s4_symmetrised",
        "eq_var_ratio": "inapplicable_pairs_and_solve_no_boolean_system",
        "ffd_or_degree_of_regularity": {
            "status": "inapplicable",
            "reason": (
                "oracle_class=pairs_and_solve enumerates pair space; "
                "no Boolean Macaulay/GB system in this arm"
            ),
            "note": (
                "Ledger sub-2^(2l) beat still wants subspace/oracle FFD at l=8; "
                "see runs_manual/binary_s4_ffd_l8_20260920 for related S4 subspace FFD."
            ),
        },
        "oracle_class": "pairs_and_solve",
        "median_ms_per_target": median_ms,
        "largest_solvable": {
            "budget": "pairs_baseline",
            "l": l_goal,
            "n": row["n"],
            "pairs_solve_s": pairs_s,
        },
        "ladder": [
            {
                "l": item["l"],
                "pairs_s": item["pairs_s"],
                "speedup_vs_triples": item["speedup_vs_triples"],
                "triples_gen_s": item["triples_gen_s"],
            }
            for item in ladder
        ],
        "ledger_gate_assessment": {
            "ffd_logged": False,
            "median_le_half_pairs": "n/a_until_faster_oracle_measured",
            "note": (
                "pairs-and-solve IS Theta(2^(2l)); this measures the baseline, "
                "not a sub-2^(2l) win"
            ),
            "sub_2_to_2l_oracle": False,
            "targets_le_64": True,
        },
        "wall_s": float(direct_obs.get("whole_process_wall_ms") or 0.0) / 1000.0,
        "peak_rss_bytes": peak_rss,
        "verdict": (
            "L8_PAIRS_BASELINE_MEASURED_NO_SUB_2_2L_PENDING_INDEPENDENT_VALIDATION"
            if producers_ok
            else "PRODUCER_FAILURE"
        ),
        "claim_boundary": (
            "Public synthetic binary Semaev pairs-and-solve baseline only. "
            "Not key recovery, not a sub-2^(2l) oracle claim, not ledger promotion."
        ),
        "claim_boundary_non_claims": [
            "not key recovery",
            "not sub-2^(2l) oracle claim",
            "not ledger promotion",
            "FFD unmeasured on pairs arm",
        ],
        "independent_replay_pointer": independent_replay_pointer(run),
        "fixture_hash": sha256(source),
        "executable_or_source_hash": {
            "executable": sha256(binary_path) if binary_path else None,
            "source": sha256(source),
        },
        "host_id": host_record(),
        "resource_caps": {
            "common_cap_bytes": cap,
            "observed_peak_rss_bytes": peak_rss,
            "under_cap": under_cap,
        },
        "seeds": {"decomp_bench_internal": "0x9E3779B9"},
        "evidence_paths": [
            str((run / "logs/direct.stdout.jsonl").resolve()),
            str((run / "receipts/direct.resource.json").resolve()),
        ],
        "created_at": now(),
        "autolab_run_id": run.name,
        "producer_exit_codes": {"direct": direct_obs["exit_code"]},
    }
    return claim


def draft_decomposition_claim(
    *,
    beat: dict[str, Any],
    beat_id: str,
    run: Path,
    direct_obs: dict[str, Any],
    binaries: dict[str, str],
    fixtures: int,
) -> dict[str, Any]:
    producer = beat.get("producer") or {}
    if str(producer.get("example") or "") == "semaev_decomp_bench":
        return draft_binary_pairs_decomposition_claim(
            beat=beat,
            beat_id=beat_id,
            run=run,
            direct_obs=direct_obs,
            binaries=binaries,
        )
    parsed = parse_paired_oracle_stdout(direct_obs.get("stdout") or "")
    header = parsed["header"]
    arms = parsed["arms"]
    producers_ok = direct_obs["exit_code"] == 0
    peak_rss = direct_obs.get("children_peak_rss_bytes")
    cap = beat.get("resource_cap_bytes")
    under_cap = peak_rss is not None and cap is not None and int(peak_rss) <= int(cap)
    ffds = parsed["ffds"]
    ffd_block = {
        "draws": 1,
        "min": min(ffds) if ffds else None,
        "max": max(ffds) if ffds else None,
        "mean": (sum(ffds) / len(ffds)) if ffds else None,
        "note": "FFD from first target draw in paired_bench; control-plane decomposition launch",
    }
    for key, arm in arms.items():
        if isinstance(arm.get("ffd"), int):
            ffd_block[key] = arm["ffd"]
    x_arm = arms.get("x_chained_f4") or next(iter(arms.values()))
    unknowns = int(beat.get("unknowns_goal") or x_arm.get("vars") or 0)
    dim = int(beat.get("dim_goal") or header["dim"])
    median_ms = parsed["median_ms_per_target"]
    under_1h = isinstance(median_ms, (int, float)) and median_ms <= 3_600_000.0
    source = REPO / "examples" / f"{beat['producer']['example']}.rs"
    binary_path = binaries.get("custom") or binaries.get("direct")
    claim = {
        "schema_version": 2,
        "task_id": TASK_ID,
        "beat_id": beat_id,
        "regime": beat["regime"],
        "stage": "decomposition",
        "n": beat["n"],
        "a": beat.get("a", header["a"]),
        "n_or_bits": beat["n"],
        "m_summands": beat.get("m", header["m"]),
        "dimension_l_or_dim": dim,
        "unknowns": unknowns,
        "system_degree": beat.get("system_degree", "quadratic"),
        "eq_var_ratio": (x_arm["eqs"] / x_arm["vars"]) if x_arm.get("vars") else None,
        "ffd_or_degree_of_regularity": ffd_block,
        "oracle_class": "groebner",
        "median_ms_per_target": median_ms,
        "median_ms_per_target_detail": {
            "reported_as": "max(found, refuted) class medians for 1h gate honesty",
            **{
                f"{key}_found_median_ms": arm.get("found_ms_median")
                for key, arm in arms.items()
            },
            **{
                f"{key}_refuted_median_ms": arm.get("refuted_ms_median")
                for key, arm in arms.items()
            },
        },
        "largest_solvable": {
            "budget": "1h",
            "dim": dim,
            "m": beat.get("m", header["m"]),
            "n": beat["n"],
            "observed_median_ms": median_ms,
            "under_budget": under_1h,
            "unknowns": unknowns,
            "verdict_mix": {
                "found": x_arm.get("found"),
                "inconclusive": x_arm.get("inconclusive"),
                "refuted": x_arm.get("refuted"),
            },
        },
        "arms": arms,
        "disagreements": parsed["disagreements"],
        "divisor_indices": parsed["divisor_indices"],
        "fb_sizes": {"F_x": header["F_x"], "F_u": header["F_u"]},
        "targets": header["targets"],
        "node_budget": 50000,
        "sat": False,
        "wall_s": parsed["wall_s"]
        if parsed["wall_s"] is not None
        else (float(direct_obs.get("whole_process_wall_ms") or 0.0) / 1000.0),
        "peak_rss_bytes": peak_rss,
        "acceptance_vs_ledger": {
            "dim_16_m_2_quadratic": dim == 16 and int(beat.get("m", 2)) == 2,
            "disagreements_zero": parsed["disagreements"] == 0,
            "ffd_logged": bool(ffds),
            "median_under_1h": under_1h,
        },
        "verdict": (
            "N31_M2_DIM16_PANEL_PASS_UNDER_1H_PENDING_INDEPENDENT_VALIDATION"
            if producers_ok
            and parsed["disagreements"] == 0
            and under_1h
            and ffds
            else (
                "PRODUCER_FAILURE"
                if not producers_ok
                else "DRAFT_INCOMPLETE_OR_GATE_FAIL"
            )
        ),
        "claim_boundary": (
            "Public synthetic Koblitz Semaev decomposition oracle measure only. "
            "Not key recovery, not asymptotic sub-rho, not imported points, not "
            "ledger promotion until independent validation."
        ),
        "claim_boundary_non_claims": [
            "not key recovery",
            "not asymptotic sub-sqrt",
            "not imported/external points",
            "not ledger promotion until independent validation",
        ],
        "independent_replay_pointer": independent_replay_pointer(run),
        "fixture_hash": sha256(source),
        "executable_or_source_hash": {
            "executable": sha256(binary_path) if binary_path else None,
            "source": sha256(source),
        },
        "host_id": host_record(),
        "resource_caps": {
            "common_cap_bytes": cap,
            "node_budget": 50000,
            "observed_peak_rss_bytes": peak_rss,
            "under_cap": under_cap,
        },
        "seeds": {"paired_bench_default": hex(int((beat.get("env") or {}).get("seed", 0x5EED)))},
        "evidence_paths": [
            str((run / "logs/direct.stdout.jsonl").resolve()),
            str((run / "receipts/direct.resource.json").resolve()),
        ],
        "created_at": now(),
        "autolab_run_id": run.name,
        "producer_exit_codes": {"direct": direct_obs["exit_code"]},
    }
    return claim


def draft_factor_base_claim(
    *,
    beat: dict[str, Any],
    beat_id: str,
    run: Path,
    direct_obs: dict[str, Any],
    binaries: dict[str, str],
) -> dict[str, Any]:
    rows = parse_json_objects(direct_obs["stdout"])
    base = select_factor_base_row(rows) or {}
    producers_ok = direct_obs["exit_code"] == 0
    factor_base_size = base.get("factor_base_points")
    if isinstance(factor_base_size, list):
        factor_base_size = len(factor_base_size)
    orbit_count = base.get("orbit_columns")
    retained = base.get("support_payload_lower_bound_bytes")
    construction_wall = base.get("total_setup_ms")
    if construction_wall is None:
        construction_wall = direct_obs.get("whole_process_wall_ms")
    peak_rss = direct_obs.get("children_peak_rss_bytes")
    cap = beat.get("resource_cap_bytes")
    under_cap = (
        peak_rss is not None and cap is not None and int(peak_rss) <= int(cap)
    )
    pair_mode = beat.get("direct", {}).get("pair_mode", "unknown")
    eta = beat.get("eta")
    eta_label = (
        f"eta_{eta[0]}_{eta[1]}"
        if isinstance(eta, list) and len(eta) == 2
        else "eta_unknown"
    )
    claim = {
        "schema_version": 2,
        "task_id": TASK_ID,
        "beat_id": beat_id,
        "regime": beat["regime"],
        "stage": "factor_base",
        "n": beat["n"],
        "n_or_bits": beat["n"],
        "factor_base_size_F": factor_base_size,
        "orbit_count_K": orbit_count,
        "dimension_l_or_dim": orbit_count,
        "construction_method": (
            f"point_defined_{pair_mode}_{eta_label}"
            + ("_summary_only" if (beat.get("env") or {}).get("KIC_SUMMARY_ONLY") == "1" else "")
        ),
        "materialized": True,
        "construction_wall_ms": construction_wall,
        "retained_bytes": retained,
        "frobenius_closed": base.get("frobenius_closed"),
        "negation_closed": base.get("negation_closed"),
        "subgroup_membership_verified": base.get("subgroup_membership_verified"),
        "pair_index_mode": base.get("pair_index_mode", pair_mode),
        "eta": base.get("eta")
        or (
            {"numerator": eta[0], "denominator": eta[1]}
            if isinstance(eta, list) and len(eta) == 2
            else eta
        ),
        "base_hash": base.get("base_hash"),
        "support_table_allocated_bytes": base.get("support_table_allocated_bytes"),
        "producer_timings_ms": {
            "curve_setup_ms": base.get("curve_setup_ms"),
            "base_construction_ms": base.get("base_construction_ms"),
            "support_index_ms": base.get("support_index_ms"),
            "independent_base_validation_ms": base.get(
                "independent_base_validation_ms"
            ),
            "total_setup_ms": base.get("total_setup_ms"),
            "whole_process_wall_ms": direct_obs.get("whole_process_wall_ms"),
        },
        "verdict": (
            "DRAFT_PENDING_INDEPENDENT_VALIDATION"
            if producers_ok and under_cap
            else (
                "PRODUCER_FAILURE"
                if not producers_ok
                else "DRAFT_RESOURCE_CAP_EXCEEDED_OR_INCOMPLETE"
            )
        ),
        "claim_boundary": (
            "Public synthetic Koblitz factor-base construction measure only. "
            "Not key recovery, not asymptotic sub-rho, not an imported-point "
            "attack, and not a ledger promotion until independent validation "
            "and schema PASS."
        ),
        "claim_boundary_non_claims": [
            "not key recovery",
            "not asymptotic sub-sqrt",
            "not imported/external points",
            "not ledger promotion until independent validation",
        ],
        "independent_replay_pointer": independent_replay_pointer(run),
        "fixture_hash": sha256(REPO / "examples/koblitz_rank_fixture.rs"),
        "executable_or_source_hash": {
            "direct": sha256(binaries["direct"]),
            "rank_fixture_source": sha256(REPO / "examples/koblitz_rank_fixture.rs"),
        },
        "host_id": host_record(),
        "resource_caps": {
            "common_cap_bytes": cap,
            "observed_peak_rss_bytes": peak_rss,
            "under_cap": under_cap,
        },
        "seeds": {"direct": direct_obs.get("seed")},
        "producer_exit_codes": {"direct": direct_obs["exit_code"]},
        "whole_process_wall_ms": direct_obs["whole_process_wall_ms"],
        "direct_rows": len(rows),
    }
    return claim


def launch(arguments: argparse.Namespace) -> dict[str, Any]:
    protocol = load_protocol()
    ledger = load_ledger(protocol)
    beat_id = arguments.beat
    require(beat_id in protocol["beats"], f"unknown beat id: {beat_id}")
    beat = protocol["beats"][beat_id]
    require(
        beat.get("launch_mode") != "single_target_panel",
        f"{beat_id} is a single-target panel; use launch-single",
    )
    require(
        beat.get("launch_mode") != "equal_advice_bl",
        f"{beat_id} is an equal-advice BL reference beat; use launch-equal-advice",
    )
    require(
        beat.get("launch_mode", "single_target") == "single_target",
        f"{beat_id} is a multi-target batch panel; use launch-panel",
    )
    fixtures = arguments.fixtures
    if fixtures is None:
        fixtures = int(beat.get("fixtures", beat.get("fixtures_default", 1)))
    require(fixtures > 0, "fixtures must be positive")
    if beat.get("stage") == "vs_rho":
        require(
            fixtures == 1,
            "vs_rho primary launches require exactly one target; define a separate secondary workload before using multiple targets",
        )

    with RunnerLock():
        preflight_receipt = preflight(protocol, ledger)
        require(preflight_receipt["ok"], "preflight failed; see artifacts after launch dir create")
        run_id = arguments.run_id or make_run_id(beat_id)
        run = RUNS_DIR / run_id
        require(not run.exists(), f"run already exists: {run}")
        for name in ("artifacts", "inputs", "logs", "receipts"):
            (run / name).mkdir(parents=True, exist_ok=True)
        write_json(CURRENT_PATH, {"run_id": run_id, "beat_id": beat_id, "updated_at": now()})
        write_json(run / "artifacts/preflight.json", preflight_receipt)
        write_json(run / "inputs/protocol.json", protocol)
        write_json(run / "inputs/boundary_targets.json", ledger)
        write_json(
            run / "inputs/ledger_pin.json",
            {
                "path": protocol["ledger"]["path"],
                "sha256": sha256(ledger_path(protocol)),
                "schema_version": ledger["schema_version"],
            },
        )

        state: dict[str, Any] = {
            "schema_version": "1.0",
            "task_id": TASK_ID,
            "run_id": run_id,
            "beat_id": beat_id,
            "status": "ACTIVE",
            "phase": "build",
            "fixtures": fixtures,
            "online_target_count": fixtures,
            "precompute_fixture_count": int(beat.get("precompute_fixtures", 0)),
            "created_at": now(),
            "updated_at": now(),
        }
        write_json(run / "state.json", state)

        if arguments.prepare_only:
            state.update(status="PREPARED", phase="prepared", updated_at=now())
            write_json(run / "state.json", state)
            commands = producer_commands(beat, beat_id, fixtures, binaries=None)
            write_json(run / "artifacts/prepared_commands.json", commands)
            return state

        stage = str(beat.get("stage") or "vs_rho")
        if stage in ("decomposition", "end_to_end_dlp") and isinstance(
            beat.get("producer"), dict
        ) and beat["producer"].get("example"):
            require(
                beat["producer"].get("example"),
                f"beat {beat_id} stage {stage} requires producer.example",
            )
            custom = build_example(str(beat["producer"]["example"]))
            binaries = {"custom": custom, "direct": custom}
        elif stage == "decomposition":
            raise AutolabError(
                f"beat {beat_id} stage decomposition requires producer.example"
            )
        else:
            binaries = build_producers()
        write_json(
            run / "artifacts/binaries.json",
            {key: {"path": path, "sha256": sha256(path)} for key, path in binaries.items()},
        )
        state.update(phase="measurement", updated_at=now())
        write_json(run / "state.json", state)

        env = os.environ.copy()
        # Default to incremental crosscheck for small rungs; beats may override
        # (n=53 dense recompute-after-every-relation is hour-class).
        env.setdefault("KIC_INCREMENTAL_RANK_CROSSCHECK", "1")
        for key, value in (beat.get("env") or {}).items():
            env[str(key)] = str(value)

        commands = producer_commands(beat, beat_id, fixtures, binaries=binaries)
        write_json(run / "artifacts/commands.json", commands)

        direct_seed = seed_for(beat_id, "direct", 0)
        direct_cmd = commands["direct_argv"]
        direct_obs = run_timed(direct_cmd, env=env, cwd=REPO)
        direct_obs["seed"] = direct_seed
        (run / "logs/direct.stdout.jsonl").write_text(direct_obs["stdout"])
        (run / "logs/direct.stderr.txt").write_text(direct_obs["stderr"])
        write_json(
            run / "receipts/direct.resource.json",
            {key: direct_obs.get(key) for key in RESOURCE_RECEIPT_FIELDS},
        )

        if stage == "factor_base":
            claim = draft_factor_base_claim(
                beat=beat,
                beat_id=beat_id,
                run=run,
                direct_obs=direct_obs,
                binaries=binaries,
            )
            write_json(run / "artifacts/claim_draft.json", claim)
            validation = validate_claim(claim, stage="factor_base", ledger=ledger)
            write_json(run / "artifacts/claim_check.json", validation)
            producers_ok = direct_obs["exit_code"] == 0
            status = (
                "PENDING_INDEPENDENT_VALIDATION" if producers_ok else "PRODUCER_FAILURE"
            )
            if producers_ok and validation["status"] != "PASS":
                status = "SCHEMA_INCOMPLETE"
            state.update(
                status=status,
                phase="analysis",
                updated_at=now(),
                direct_exit_code=direct_obs["exit_code"],
                claim_check=validation["status"],
            )
        elif stage == "vs_rho":
            require(
                isinstance(beat.get("rho"), dict),
                f"beat {beat_id} stage vs_rho requires a rho config object",
            )
            rho_seed = seed_for(beat_id, "rho", 0)
            rho_cmd = commands["rho_argv"]
            rho_env = env.copy()
            rho_fixed_target_scalar = None
            if beat.get("pair_rho_to_direct_target", True):
                direct_rows = parse_json_objects(direct_obs["stdout"])
                direct_online = next(
                    (row for row in reversed(direct_rows)
                     if row.get("kind") == "relation_rank_summary"
                     and row.get("online_target_count") == 1),
                    None,
                )
                if direct_online is not None:
                    rho_fixed_target_scalar = direct_online.get("published_fixture_scalar")
                if isinstance(rho_fixed_target_scalar, int):
                    rho_env["KIC_RHO_FIXED_TARGET_SCALAR"] = str(rho_fixed_target_scalar)
                    commands["rho_environment_overrides"] = {
                        "KIC_RHO_FIXED_TARGET_SCALAR": str(rho_fixed_target_scalar),
                        "source": "the one IC online target's published synthetic fixture scalar; target generation excluded from both online intervals",
                    }
                else:
                    commands["rho_environment_overrides"] = None
            write_json(run / "artifacts/commands.json", commands)
            if rho_fixed_target_scalar is None and beat.get("pair_rho_to_direct_target", True):
                rho_obs = {
                    "exit_code": 1,
                    "whole_process_wall_ms": 0.0,
                    "children_cpu_user_ms": 0.0,
                    "children_cpu_system_ms": 0.0,
                    "children_peak_rss_bytes": 0,
                    "command": rho_cmd,
                    "stdout": "",
                    "stderr": "IC producer did not emit the single online target fixture scalar; rho was not launched",
                }
            else:
                rho_obs = run_timed(rho_cmd, env=rho_env, cwd=REPO)
            rho_obs["seed"] = rho_seed
            rho_obs["fixed_target_scalar"] = rho_fixed_target_scalar
            (run / "logs/rho.stdout.jsonl").write_text(rho_obs["stdout"])
            (run / "logs/rho.stderr.txt").write_text(rho_obs["stderr"])
            write_json(
                run / "receipts/rho.resource.json",
                {
                    **{key: rho_obs.get(key) for key in RESOURCE_RECEIPT_FIELDS},
                    "fixed_target_scalar": rho_obs.get("fixed_target_scalar"),
                },
            )

            claim = draft_vs_rho_claim(
                beat=beat,
                beat_id=beat_id,
                run=run,
                direct_obs=direct_obs,
                rho_obs=rho_obs,
                binaries=binaries,
            )
            write_json(run / "artifacts/claim_draft.json", claim)
            validation = validate_claim(claim, stage="vs_rho", ledger=ledger)
            write_json(run / "artifacts/claim_check.json", validation)

            producers_ok = direct_obs["exit_code"] == 0 and rho_obs["exit_code"] == 0
            status = (
                "PENDING_INDEPENDENT_VALIDATION" if producers_ok else "PRODUCER_FAILURE"
            )
            pairing_errors = []
            if producers_ok and claim["paired_target"]["same_public_point"] is not True:
                pairing_errors.append("IC and rho public targets differ or are missing")
            if producers_ok and claim["paired_target"]["same_fixture_scalar"] is not True:
                pairing_errors.append("IC and rho fixture scalars differ or are missing")
            if producers_ok and pairing_errors:
                status = "PAIRING_FAILURE"
            elif producers_ok and not claim["all_stages_charged_same_series"]:
                status = "ACCOUNTING_INCOMPLETE"
            elif producers_ok and validation["status"] != "PASS":
                status = "SCHEMA_INCOMPLETE"
            state.update(
                status=status,
                phase="analysis",
                updated_at=now(),
                direct_exit_code=direct_obs["exit_code"],
                rho_exit_code=rho_obs["exit_code"],
                claim_check=validation["status"],
                pairing_errors=pairing_errors,
            )
        elif stage == "decomposition":
            claim = draft_decomposition_claim(
                beat=beat,
                beat_id=beat_id,
                run=run,
                direct_obs=direct_obs,
                binaries=binaries,
                fixtures=fixtures,
            )
            write_json(run / "artifacts/claim_draft.json", claim)
            validation = validate_claim(claim, stage="decomposition", ledger=ledger)
            write_json(run / "artifacts/claim_check.json", validation)
            producers_ok = direct_obs["exit_code"] == 0
            status = (
                "PENDING_INDEPENDENT_VALIDATION" if producers_ok else "PRODUCER_FAILURE"
            )
            if producers_ok and validation["status"] != "PASS":
                status = "SCHEMA_INCOMPLETE"
            state.update(
                status=status,
                phase="analysis",
                updated_at=now(),
                direct_exit_code=direct_obs["exit_code"],
                claim_check=validation["status"],
            )
        elif stage == "end_to_end_dlp":
            producer = beat.get("producer") or {}
            require(
                str(producer.get("example") or "") == "j0_known_answer_bits",
                f"beat {beat_id} end_to_end_dlp launch supports j0_known_answer_bits only",
            )
            claim = draft_j0_end_to_end_dlp_claim(
                beat=beat,
                beat_id=beat_id,
                run=run,
                direct_obs=direct_obs,
                binaries=binaries,
            )
            write_json(run / "artifacts/claim_draft.json", claim)
            validation = validate_claim(claim, stage="end_to_end_dlp", ledger=ledger)
            write_json(run / "artifacts/claim_check.json", validation)
            producers_ok = (
                direct_obs["exit_code"] == 0
                and claim.get("ic_ok") is True
                and claim.get("rho_ok") is True
                and claim.get("ic_agrees_rho") is True
            )
            status = (
                "PENDING_INDEPENDENT_VALIDATION" if producers_ok else "PRODUCER_FAILURE"
            )
            if producers_ok and validation["status"] != "PASS":
                status = "SCHEMA_INCOMPLETE"
            state.update(
                status=status,
                phase="analysis",
                updated_at=now(),
                direct_exit_code=direct_obs["exit_code"],
                claim_check=validation["status"],
            )
        else:
            raise AutolabError(f"unsupported launch stage: {stage}")

        write_json(run / "state.json", state)
        write_json(
            run / "artifacts/candidate.json",
            {
                "schema_version": "1.0",
                "task_id": TASK_ID,
                "run_id": run_id,
                "beat_id": beat_id,
                "status": status,
                "claim_draft_sha256": sha256(run / "artifacts/claim_draft.json"),
                "claim_check": validation,
                "ledger_sha256": sha256(ledger_path(protocol)),
                "created_at": now(),
                "note": (
                    "Draft only. Promote the ledger only after independent validation "
                    "and a claim-check PASS with every required measurement field."
                ),
            },
        )
        files = {
            str(path.relative_to(run)): sha256(path)
            for path in sorted(run.rglob("*"))
            if path.is_file()
        }
        write_json(
            run / "artifacts/review_manifest.json",
            {"schema_version": "1.0", "task_id": TASK_ID, "files": files},
        )
        return state


def producer_commands(
    beat: dict[str, Any],
    beat_id: str,
    fixtures: int,
    binaries: dict[str, str] | None,
) -> dict[str, Any]:
    producer = beat.get("producer")
    if isinstance(producer, dict) and producer.get("args_template"):
        example = str(producer["example"])
        direct_bin = (
            binaries.get("custom")
            if binaries and binaries.get("custom")
            else f"target/release/examples/{example}"
        )
        rendered = str(producer["args_template"]).format(
            a=beat.get("a", 0),
            n=beat.get("n"),
            m=beat.get("m"),
            l=beat.get("l"),
            bits=beat.get("bits"),
            fixtures=fixtures,
        )
        direct_argv = [direct_bin, *shlex.split(rendered)]
        return {
            "build": f"cargo build --release --example {example}",
            "stage": beat.get("stage"),
            "direct": " ".join(direct_argv),
            "direct_argv": direct_argv,
            "fixtures": fixtures,
            "online_target_count": fixtures,
            "precompute_fixtures": 0,
            "direct_producer_fixture_count": fixtures,
            "rho_producer_fixture_count": 0,
            "env": beat.get("env") or {},
            "rho": None,
            "rho_argv": None,
            "custom_example": example,
        }

    direct_bin = (
        binaries["direct"]
        if binaries
        else "target/release/examples/koblitz_rank_fixture"
    )
    rho_bin = (
        binaries["rho"] if binaries else "target/release/examples/koblitz_rho_fixture"
    )
    direct_seed = seed_for(beat_id, "direct", 0)
    precompute_fixtures = int(beat.get("precompute_fixtures", 0))
    direct_fixtures = fixtures + precompute_fixtures
    direct_argv = [
        direct_bin,
        str(beat["n"]),
        str(beat["a"]),
        str(beat["eta"][0]),
        str(beat["eta"][1]),
        str(direct_seed),
        beat["direct"]["pair_mode"],
        beat["direct"]["target_mode"],
        beat["direct"]["query_mode"],
        str(direct_fixtures),
    ]
    result: dict[str, Any] = {
        "build": (
            "cargo build --release --example koblitz_rank_fixture "
            "--example koblitz_rho_fixture"
        ),
        "stage": beat.get("stage"),
        "direct": " ".join(direct_argv),
        "direct_argv": direct_argv,
        "fixtures": fixtures,
        "online_target_count": fixtures,
        "precompute_fixtures": precompute_fixtures,
        "direct_producer_fixture_count": direct_fixtures,
        "rho_producer_fixture_count": fixtures,
        "env": beat.get("env") or {},
    }
    rho_cfg = beat.get("rho")
    if isinstance(rho_cfg, dict):
        rho_seed = seed_for(beat_id, "rho", 0)
        rho_argv = [
            rho_bin,
            str(beat["n"]),
            str(beat["a"]),
            rho_cfg["quotient_mode"],
            str(fixtures),
            rho_cfg["backend"],
            str(rho_seed),
        ]
        result["rho"] = " ".join(rho_argv)
        result["rho_argv"] = rho_argv
    else:
        result["rho"] = None
        result["rho_argv"] = None
    return result


TIME_L_FIELDS = {
    "maximum resident set size": "max_rss_bytes",
    "instructions retired": "instructions_retired",
    "cycles elapsed": "cycles_elapsed",
    "peak memory footprint": "peak_memory_footprint_bytes",
    "page reclaims": "page_reclaims",
    "page faults": "page_faults",
    "swaps": "swaps",
    "voluntary context switches": "voluntary_context_switches",
    "involuntary context switches": "involuntary_context_switches",
}


def parse_time_l(text: str) -> dict[str, int]:
    """Counters from macOS `/usr/bin/time -l` (stderr may hold producer lines too)."""
    record: dict[str, int] = {}
    for line in text.splitlines():
        stripped = line.strip()
        for label, key in TIME_L_FIELDS.items():
            if stripped.endswith(label) and stripped.split()[0].isdigit():
                record[key] = int(stripped.split()[0])
    return record


def host_conditions() -> dict[str, Any]:
    conditions: dict[str, Any] = {"at": now(), "loadavg": list(os.getloadavg())}
    if sys.platform == "darwin":
        vm = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
        pages: dict[str, int] = {}
        page_size = 4096
        for line in vm.splitlines():
            if "page size of" in line:
                page_size = int(line.split("page size of")[1].split()[0])
            elif ":" in line:
                key, _, value = line.partition(":")
                value = value.strip().rstrip(".")
                if value.isdigit():
                    pages[key.strip()] = int(value)
        reclaimable = sum(
            pages.get(key, 0)
            for key in ("Pages free", "Pages inactive", "Pages purgeable", "Pages speculative")
        )
        conditions["available_memory_bytes"] = reclaimable * page_size
        conditions["swap"] = subprocess.run(
            ["sysctl", "-n", "vm.swapusage"], capture_output=True, text=True
        ).stdout.strip()
    elif Path("/proc/meminfo").is_file():
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                conditions["available_memory_bytes"] = int(line.split()[1]) * 1024
    return conditions


def run_panel_process(
    command: list[str],
    *,
    env: dict[str, str],
    stdout_path: Path,
    stderr_path: Path,
    cpu: int | None,
) -> dict[str, Any]:
    """One whole-process run: wall from the parent, CPU and RSS from wait4.

    On macOS the command runs under `/usr/bin/time -l` for the retired
    instruction count; `cpu` pins it with taskset (Linux only).
    """
    wrapped = list(command)
    if cpu is not None:
        require(shutil.which("taskset") is not None, "--cpu needs taskset (Linux)")
        wrapped = ["taskset", "-c", str(cpu)] + wrapped
    if sys.platform == "darwin" and Path("/usr/bin/time").is_file():
        wrapped = ["/usr/bin/time", "-l"] + wrapped
    before = host_conditions()
    with stdout_path.open("wb") as out, stderr_path.open("wb") as err:
        started = time.perf_counter()
        process = subprocess.Popen(wrapped, cwd=REPO, env=env, stdout=out, stderr=err)
        _, status_code, usage = os.wait4(process.pid, 0)
        wall_s = time.perf_counter() - started
    process.returncode = os.waitstatus_to_exitcode(status_code)
    rss_scale = 1 if sys.platform == "darwin" else 1024
    record: dict[str, Any] = {
        "command": command,
        "pinned_cpu": cpu,
        "exit_code": process.returncode,
        "wall_s": wall_s,
        "user_s": usage.ru_utime,
        "sys_s": usage.ru_stime,
        "max_rss_bytes": usage.ru_maxrss * rss_scale,
        "voluntary_context_switches": usage.ru_nvcsw,
        "involuntary_context_switches": usage.ru_nivcsw,
        "conditions_before": before,
        "conditions_after": host_conditions(),
    }
    if sys.platform == "darwin":
        record.update(
            {k: v for k, v in parse_time_l(stderr_path.read_text(errors="replace")).items()
             if k not in record or k == "max_rss_bytes"}
        )
    return record


def panel_corpus_name(template: str, targets: int) -> str:
    return template.format(L=targets)


def estimate_ic_rss_bytes(beat: dict[str, Any], k: int) -> int:
    """Peak RSS of the compact-orbit IC at K orbit columns.

    `compact_orbit_rss` follows the producer's allocations: K^2*n regular states
    of `state_bytes` each (the vector's touched part) plus a root table of
    `table_slot_bytes` per slot, sized to the next power of two of
    `table_slots_per_state` * K^2*n and fully initialised. The table doubles at
    power-of-two boundaries, which a flat bytes-per-state constant misses.
    """
    states = int(beat["n"]) * k * k
    model = beat.get("ic_rss_model")
    if model is None:
        return int(beat["ic_bytes_per_regular_state"]) * states
    require(model["kind"] == "compact_orbit_rss", f"unknown IC RSS model {model['kind']}")
    slots = 1 << max(4, (int(model["table_slots_per_state"]) * states - 1).bit_length())
    return (int(model["state_bytes"]) * states + int(model["table_slot_bytes"]) * slots
            + int(model["fixed_bytes"]))


def k_fits(beat: dict[str, Any], k: int) -> tuple[bool, int, int | None]:
    need = estimate_ic_rss_bytes(beat, k)
    headroom = int((beat.get("ic_rss_model") or {}).get("headroom_bytes", 0))
    available = host_conditions().get("available_memory_bytes")
    return (available is None or need + headroom < available), need, available


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return parse_json_lines(path.read_text()) if path.is_file() else []


def untimed_digest(records: list[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for record in records:
        digest.update(
            json.dumps({k: v for k, v in record.items() if "_ms" not in k}, sort_keys=True).encode()
        )
    return digest.hexdigest()


def check_panel_block(
    ic_records: list[dict[str, Any]],
    rho_records: list[dict[str, Any]],
    corpus: list[int],
) -> dict[str, Any]:
    ic = [r for r in ic_records if r.get("kind") == "compact_orbit_dlp_target"]
    rho = [r for r in rho_records if r.get("kind") == "rho_ks_batch_fixture"]
    ic_scalars = [r.get("published_fixture_scalar") for r in ic]
    rho_scalars = [r.get("published_fixture_scalar") for r in sorted(rho, key=lambda r: r["fixture_index"])]
    ic_points = [r.get("target") for r in ic]
    rho_points = [r.get("published_q") for r in sorted(rho, key=lambda r: r["fixture_index"])]
    return {
        "ic_targets": len(ic),
        "rho_targets": len(rho),
        "ic_all_verified": bool(ic) and all(
            r.get("recovered_matches_published") is True and r.get("group_verified") is True for r in ic
        ),
        "rho_all_verified": bool(rho) and all(
            r.get("verified") is True
            and r.get("recovered_fixture_scalar") == r.get("published_fixture_scalar")
            for r in rho
        ),
        "ic_matches_corpus": ic_scalars == corpus,
        "rho_matches_corpus": rho_scalars == corpus,
        "same_target_points": bool(ic_points) and ic_points == rho_points,
        "ic_untimed_sha256": untimed_digest(ic),
        "rho_untimed_sha256": untimed_digest(rho),
        "ic_target_probes_total": sum(r.get("probes") or 0 for r in ic),
    }


def panel_ratio(ic_run: dict[str, Any], rho_run: dict[str, Any], key: str) -> float | None:
    ic_value, rho_value = ic_run.get(key), rho_run.get(key)
    if isinstance(ic_value, (int, float)) and isinstance(rho_value, (int, float)) and rho_value > 0:
        return ic_value / rho_value
    return None


def median_range(values: list[float]) -> dict[str, Any] | None:
    if not values:
        return None
    return {"median": statistics.median(values), "range": [min(values), max(values)], "values": values}


def draft_panel_claims(
    *,
    beat_id: str,
    beat: dict[str, Any],
    run: Path,
    summary: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """end_to_end_dlp (primary) and vs_rho (supplementary, fails closed) drafts.

    Only measured values are written.  The vs_rho schema's single-target
    online fields are absent because the batch producers do not emit them.
    """
    common = {
        "schema_version": 2,
        "task_id": TASK_ID,
        "beat_id": beat_id,
        "autolab_run_id": run.name,
        "regime": beat["regime"],
        "result_class": beat["result_class"],
        "n_or_bits": {"n": beat["n"], "a": beat["a"], "subgroup_order_bits": summary.get("subgroup_order_bits")},
        "fixture_hash": summary["corpora"],
        "executable_or_source_hash": summary["executables"],
        "host_id": summary["host"],
        "resource_caps": {
            "threads": 1,
            "memory": "no cap; K candidates skipped when the estimated IC RSS exceeded available memory",
            "isolation": summary["isolation"],
        },
        "seeds": {
            "batch_seed": beat["batch_seed"],
            "ic_rank_seed": beat["ic_rank_seed"],
            "rho_env": beat["panel_producers"]["rho"].get("env", {}),
        },
        "claim_boundary_non_claims": [
            "public synthetic known-answer fixtures only",
            "no key recovery",
            "no asymptotic sub-rho claim; both arms scale as sqrt(L*r/n) at optimal K",
            f"multi-target L={summary['targets']} batch: diagnostic, not the one-target primary comparison",
            beat["comparator_status"],
            "same-host replay; independent-host validation still required",
        ],
    }
    replay_pointer = [
        str(p.relative_to(REPO)) for p in sorted((run / "artifacts").glob("replay_*.json"))
    ]
    verified = summary["verification"]
    end_to_end = {
        **common,
        "stage": "end_to_end_dlp",
        "status": "PENDING_INDEPENDENT_VALIDATION",
        "targets_per_block": summary["targets"],
        "K": summary["K"],
        "recovered_d_verified": verified,
        "stage_timers": {
            f"b{row['block']}": row.get("ic_timing_ms") for row in summary["blocks"]
        },
        "claim_boundary": (
            f"Known-answer shared-log DLP on public synthetic a={beat['a']} n={beat['n']} fixtures, "
            f"L={summary['targets']} per batch; every recovered log checked as [d]G = Q in-process "
            "and by the pure-Python replay."
        ),
        "independent_replay_pointer": replay_pointer,
        "verdict": "DRAFT_PENDING_INDEPENDENT_VALIDATION",
    }
    vs_rho = {
        **common,
        "stage": "vs_rho",
        "status": "PENDING_INDEPENDENT_VALIDATION",
        "comparator": beat["panel_producers"]["rho"]["example"],
        "target_count": summary["targets"],
        "timing_class": "whole_process_wall",
        "K": summary["K"],
        "wall_ratio_compact_over_rho": summary["wall_ratio"],
        "user_cpu_ratio_compact_over_rho": summary["user_ratio"],
        "instructions_retired_ratio_compact_over_rho": summary["instructions_ratio"],
        "ic_scalar_verified": verified["ic_all_verified"],
        "rho_scalar_verified": verified["rho_all_verified"],
        "independent_validation": False,
        "verdict": "BATCH_DIAGNOSTIC_ONLY_NOT_SINGLE_TARGET",
        "claim_boundary": (
            f"L={summary['targets']} batched whole-process diagnostic; ineligible for the "
            "single-target online vs_rho schema."
        ),
        "independent_replay_pointer": replay_pointer,
    }
    return end_to_end, vs_rho


def launch_panel(arguments: argparse.Namespace) -> dict[str, Any]:
    protocol = load_protocol()
    ledger = load_ledger(protocol)
    beat_id = arguments.beat
    require(beat_id in protocol["beats"], f"unknown beat id: {beat_id}")
    beat = protocol["beats"][beat_id]
    require(beat.get("launch_mode") == "batch_panel", f"{beat_id} is not a batch panel beat")
    resume = getattr(arguments, "resume", None)
    if resume:
        run = RUNS_DIR / resume
        require((run / "state.json").is_file(), f"no run to resume: {run}")
        state: dict[str, Any] = read_json(run / "state.json")
        require(state.get("beat_id") == beat_id, "resume beat differs from the run's beat")
        require(state.get("phase") != "done", "run already finished")
        targets, blocks = int(state["targets"]), int(state["blocks"])
        tune_targets = int(state.get("tune_targets", targets))
        candidates = [int(k) for k in state["k_candidates"]]
        state.setdefault("resumed_at", []).append(now())
        if state.get("status") == "TUNE_ONLY":
            state["status"] = "ACTIVE"
    else:
        targets = arguments.targets or int(beat["targets_default"])
        blocks = arguments.blocks or int(beat["blocks_default"])
        tune_targets = arguments.tune_targets or targets
        if arguments.k is not None:
            candidates = [arguments.k]
        elif arguments.k_candidates:
            candidates = [int(v) for v in arguments.k_candidates.split(",") if v.strip()]
        elif tune_targets != targets:
            by_tune = beat.get("k_candidates_by_tune_targets", {}).get(str(tune_targets), {})
            candidates = [int(v) for v in by_tune.get(str(targets), [])]
        else:
            candidates = [int(v) for v in beat["k_candidates"].get(str(targets), [])]
    minimum = int(beat["targets_minimum"])
    require(targets >= minimum, f"panel targets must be >= {minimum} (no small-L panels)")
    require(tune_targets >= minimum, f"tune targets must be >= {minimum} (no small-L tunes)")
    require(blocks >= 1, "blocks must be positive")
    require(bool(candidates) and all(k > 0 for k in candidates), "no K candidates for this L")
    n, a = int(beat["n"]), int(beat["a"])
    producers = beat["panel_producers"]

    with RunnerLock():
        preflight_receipt = preflight(protocol, ledger)
        for arm, producer in producers.items():
            source = REPO / producer["source"]
            preflight_receipt["checks"].append(
                {"name": f"panel_producer_source_{arm}", "ok": source.is_file(), "detail": str(source)}
            )
        preflight_receipt["ok"] = all(
            c["ok"] for c in preflight_receipt["checks"] if c["name"] != "cryptominisat5_optional"
        )
        require(preflight_receipt["ok"], "preflight failed")
        git_head = next(c["detail"] for c in preflight_receipt["checks"] if c["name"] == "git_head")
        if resume:
            run_id = resume
            require(state.get("git_head") in (None, git_head),
                    f"resume needs the run's commit {state.get('git_head')}, not {git_head}")
            write_json(run / f"artifacts/preflight_resume{len(state['resumed_at'])}.json", preflight_receipt)
        else:
            run_id = arguments.run_id or make_run_id(beat_id)
            run = RUNS_DIR / run_id
            require(not run.exists(), f"run already exists: {run}")
        for name in ("artifacts", "inputs", "logs", "receipts"):
            (run / name).mkdir(parents=True, exist_ok=True)
        write_json(CURRENT_PATH, {"run_id": run_id, "beat_id": beat_id, "updated_at": now()})
        if not resume:
            write_json(run / "artifacts/preflight.json", preflight_receipt)
            write_json(run / "inputs/protocol.json", protocol)
            write_json(run / "inputs/boundary_targets.json", ledger)
            write_json(
                run / "inputs/ledger_pin.json",
                {"path": protocol["ledger"]["path"], "sha256": sha256(ledger_path(protocol)),
                 "schema_version": ledger["schema_version"]},
            )
            state = {
                "schema_version": "1.0", "task_id": TASK_ID, "run_id": run_id, "beat_id": beat_id,
                "launch_mode": "batch_panel", "status": "ACTIVE", "phase": "build",
                "targets": targets, "tune_targets": tune_targets, "blocks": blocks,
                "k_candidates": candidates,
                "pinned_cpu": arguments.cpu, "git_head": git_head, "created_at": now(), "updated_at": now(),
            }

        def advance(phase: str, **extra: Any) -> None:
            state.update(phase=phase, updated_at=now(), **extra)
            write_json(run / "state.json", state)

        advance("build")
        command = ["cargo", "build", "--release"]
        for producer in producers.values():
            command += ["--example", producer["example"]]
        completed = subprocess.run(command, cwd=REPO, capture_output=True, text=True)
        require(completed.returncode == 0, "cargo build failed:\n" + completed.stderr[-4000:])
        binaries = {
            arm: str((REPO / "target/release/examples" / producer["example"]).resolve())
            for arm, producer in producers.items()
        }
        executables = {
            f"{arm}_binary_sha256": sha256(path) for arm, path in binaries.items()
        } | {
            f"{Path(p['source']).name}_sha256": sha256(REPO / p["source"]) for p in producers.values()
        }
        executables["git_head"] = git_head
        executables["rustc"] = subprocess.run(
            ["rustc", "--version"], capture_output=True, text=True
        ).stdout.strip()
        if resume and (run / "artifacts/binaries.json").is_file():
            previous = read_json(run / "artifacts/binaries.json")["hashes"]
            require(all(previous.get(k) == v for k, v in executables.items() if k.endswith("_sha256")),
                    "resume rebuilt different binaries or sources")
        write_json(run / "artifacts/binaries.json", {"paths": binaries, "hashes": executables})

        env = os.environ.copy()
        env["RAYON_NUM_THREADS"] = "1"
        rho_env = dict(env, **{k: str(v) for k, v in producers["rho"].get("env", {}).items()})

        advance("corpora")
        corpora: dict[str, Any] = {}
        corpus_values: dict[str, list[int]] = {}
        for role in ("tune", "eval"):
            size = tune_targets if role == "tune" else targets
            name = panel_corpus_name(beat["corpora"][role], size)
            generated = subprocess.run(
                [binaries["rho"], str(n), str(a), "signed_frobenius", str(size), str(beat["batch_seed"])],
                cwd=REPO, capture_output=True, text=True,
                env=dict(env, KIC_RHO_GENERATE_ONLY="1", KIC_RHO_BATCH_CORPUS=name),
            )
            require(generated.returncode == 0, f"corpus generation failed for {name}")
            scalars = [
                int(row["published_fixture_scalar"])
                for row in parse_json_lines(generated.stdout)
                if row.get("kind") == "rho_ks_public_fixture"
            ]
            require(len(scalars) == size, f"corpus {name} has {len(scalars)} != {size} targets")
            path = run / f"inputs/scalars_{name}.txt"
            text = "".join(f"{s}\n" for s in scalars)
            if path.is_file():
                require(path.read_text() == text, f"regenerated corpus {name} differs from the run's copy")
            path.write_text(text)
            corpora[role] = {"name": name, "scalars_sha256": sha256(path), "targets": len(scalars)}
            corpus_values[role] = scalars
        overlap = set(corpus_values["tune"]) & set(corpus_values["eval"])
        corpora["disjoint"] = not overlap
        if tune_targets != targets:
            corpora["tune_eval_shared_targets"] = len(overlap)
        require(corpora["disjoint"], "tune and eval corpora share targets")
        tune_path = run / f"inputs/scalars_{corpora['tune']['name']}.txt"
        eval_path = run / f"inputs/scalars_{corpora['eval']['name']}.txt"

        def ic_command(k: int, scalars_path: Path, out_path: str) -> list[str]:
            return [binaries["ic"], f"construct:{n}:{a}:{k}", str(scalars_path),
                    str(beat["ic_rank_seed"]), out_path]

        advance("k_tune")
        tune_file = run / "artifacts/k_tune.json"
        tune_rows = read_json(tune_file)["rows"] if resume and tune_file.is_file() else []
        finished_k = {int(row["K"]) for row in tune_rows}
        for k in candidates:
            if k in finished_k:
                continue
            fits, need, available = k_fits(beat, k)
            row: dict[str, Any] = {"K": k, "estimated_ic_rss_bytes": need, "available_memory_bytes": available}
            if not fits:
                row["state"] = "skipped_memory"
                tune_rows.append(row)
                write_json(run / "artifacts/k_tune.json", {"rows": tune_rows})
                continue
            if len(candidates) == 1:
                row["state"] = "fixed"
                tune_rows.append(row)
                break
            if resume:
                for partial in (run / "logs").glob(f"tune_K{k}.*"):
                    partial.rename(partial.with_name(f"{partial.name}.aborted{len(state['resumed_at'])}"))
            measured = run_panel_process(
                ic_command(k, tune_path, os.devnull), env=env,
                stdout_path=run / f"logs/tune_K{k}.summary.json",
                stderr_path=run / f"logs/tune_K{k}.stderr.txt", cpu=arguments.cpu,
            )
            summary_rows = read_jsonl(run / f"logs/tune_K{k}.summary.json")
            ic_summary = summary_rows[-1] if summary_rows else {}
            ok = measured["exit_code"] == 0 and ic_summary.get("targets_failed") == 0 and (
                ic_summary.get("targets_solved") == tune_targets
            )
            row.update(state="ok" if ok else "failed", run=measured, ic_summary=ic_summary)
            tune_rows.append(row)
            write_json(run / "artifacts/k_tune.json", {"rows": tune_rows})
        completed_rows = [r for r in tune_rows if r["state"] in ("ok", "fixed")]
        require(bool(completed_rows), "no K candidate completed")
        chosen = (
            completed_rows[0]
            if completed_rows[0]["state"] == "fixed"
            else min(completed_rows, key=lambda r: r["run"]["wall_s"])
        )
        k_choice = int(chosen["K"])
        write_json(
            run / "artifacts/k_tune.json",
            {"rows": tune_rows, "tune_corpus": corpora["tune"], "selection": "lowest wall_s",
             "chosen_K": k_choice, "k_source": "fixed" if chosen["state"] == "fixed" else "tuned"}
            | ({"tune_targets": tune_targets, "panel_targets": targets,
                "tune_eval_disjoint": corpora["disjoint"]} if tune_targets != targets else {}),
        )
        if getattr(arguments, "tune_only", False):
            state.update(status="TUNE_ONLY", phase="tuned", updated_at=now(), chosen_K=k_choice)
            write_json(run / "state.json", state)
            write_json(
                run / "artifacts/candidate.json",
                {"schema_version": "1.0", "task_id": TASK_ID, "run_id": run_id, "beat_id": beat_id,
                 "status": "TUNE_ONLY", "chosen_K": k_choice, "tune_targets": tune_targets,
                 "panel_targets": targets, "tune_corpus": corpora["tune"],
                 "tune_eval_disjoint": corpora["disjoint"], "created_at": now(),
                 "note": "K tune only; no panel blocks. Resume without --tune-only to run the panel."},
            )
            write_panel_manifest(run)
            return state

        tune_manifest = run / "artifacts/review_manifest.json"
        if resume and tune_manifest.is_file():
            tune_manifest.rename(run / "artifacts/review_manifest_tune_only.json")

        advance("blocks", chosen_K=k_choice)
        base_path = run / f"logs/base_n{n}_K{k_choice}.jsonl"
        blocks_file = run / "artifacts/blocks.json"
        block_rows: list[dict[str, Any]] = read_json(blocks_file) if resume and blocks_file.is_file() else []
        exit_codes: list[int] = [row[arm]["exit_code"] for row in block_rows for arm in ("ic", "rho")]
        for b in range(blocks):
            if any(row["block"] == b for row in block_rows):
                continue
            if resume:
                # A block interrupted between its two arms is rerun whole; its
                # partial logs are kept beside the new ones, never overwritten.
                attempt = len(state["resumed_at"])
                for partial in list((run / "logs").glob(f"*_n{n}_b{b}.*")) + list(
                        (run / "receipts").glob(f"*_b{b}.resource.json")):
                    partial.rename(partial.with_name(f"{partial.name}.aborted{attempt}"))
            order = ("ic", "rho") if b % 2 == 0 else ("rho", "ic")
            row = {"block": b, "order": "_then_".join(order)}
            for arm in order:
                if arm == "ic":
                    measured = run_panel_process(
                        ic_command(k_choice, eval_path, str(run / f"logs/ic_n{n}_b{b}.jsonl")),
                        env=dict(env, KIC_DUMP_BASE=str(base_path)),
                        stdout_path=run / f"logs/ic_n{n}_b{b}.summary.json",
                        stderr_path=run / f"logs/ic_n{n}_b{b}.stderr.txt", cpu=arguments.cpu,
                    )
                else:
                    measured = run_panel_process(
                        [binaries["rho"], str(n), str(a), "signed_frobenius", str(targets),
                         str(beat["batch_seed"])],
                        env=dict(rho_env, KIC_RHO_BATCH_CORPUS=corpora["eval"]["name"]),
                        stdout_path=run / f"logs/ks_n{n}_b{b}.jsonl",
                        stderr_path=run / f"logs/ks_n{n}_b{b}.stderr.txt", cpu=arguments.cpu,
                    )
                exit_codes.append(measured["exit_code"])
                row[arm] = measured
                write_json(run / f"receipts/{arm}_b{b}.resource.json", measured)
            ic_rows = read_jsonl(run / f"logs/ic_n{n}_b{b}.jsonl")
            rho_rows = read_jsonl(run / f"logs/ks_n{n}_b{b}.jsonl")
            ic_summary_rows = read_jsonl(run / f"logs/ic_n{n}_b{b}.summary.json")
            rho_summary = [r for r in rho_rows if r.get("kind") == "rho_ks_batch_summary"]
            row["ic_summary"] = ic_summary_rows[-1] if ic_summary_rows else None
            row["ic_timing_ms"] = (row["ic_summary"] or {}).get("timing_ms")
            row["rho_summary"] = rho_summary[-1] if rho_summary else None
            row["checks"] = check_panel_block(ic_rows, rho_rows, corpus_values["eval"])
            for key in ("wall_s", "user_s", "instructions_retired"):
                row[f"{key}_ratio"] = panel_ratio(row["ic"], row["rho"], key)
            block_rows.append(row)
            write_json(run / "artifacts/blocks.json", block_rows)
        producers_ok = all(code == 0 for code in exit_codes)

        advance("replay")
        replay_dir = REPO / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924"
        replays: dict[str, Any] = {}
        if producers_ok:
            ic_replay = subprocess.run(
                [sys.executable, str(replay_dir / "independent_replay.py"), "--dlp",
                 str(run / "logs"), str(run / "artifacts/replay_ic.json"), str(n)],
                cwd=REPO, capture_output=True, text=True, env=dict(env, REPLAY_BASE=str(base_path)),
            )
            (run / "logs/replay_ic.log").write_text(ic_replay.stdout + ic_replay.stderr)
            replays["ic"] = read_json(run / "artifacts/replay_ic.json")
            for b in range(blocks):
                report = run / f"artifacts/replay_rho_b{b}.json"
                rho_replay = subprocess.run(
                    [sys.executable, str(replay_dir / "growing_n_n61_L65536_20260930_rho_replay.py"),
                     str(base_path), str(eval_path), str(run / f"logs/ic_n{n}_b{b}.jsonl"),
                     str(run / f"logs/ks_n{n}_b{b}.jsonl"), str(report)],
                    cwd=REPO, capture_output=True, text=True, env=env,
                )
                (run / f"logs/replay_rho_b{b}.log").write_text(rho_replay.stdout + rho_replay.stderr)
                replays[f"rho_b{b}"] = read_json(report) if report.is_file() else {"all_pass": False}
        ic_field = (replays.get("ic") or {}).get("fields", {}).get(str(n), {})
        rho_reports = [v for k, v in replays.items() if k.startswith("rho_")]
        verification = {
            "targets_per_block": targets,
            "blocks": blocks,
            "ic_all_verified": all(r["checks"]["ic_all_verified"] for r in block_rows),
            "rho_all_verified": all(r["checks"]["rho_all_verified"] for r in block_rows),
            "scalars_match_corpus": all(
                r["checks"]["ic_matches_corpus"] and r["checks"]["rho_matches_corpus"] for r in block_rows
            ),
            "same_target_points_all_blocks": all(r["checks"]["same_target_points"] for r in block_rows),
            "untimed_records_identical_across_blocks": {
                arm: len({r["checks"][f"{arm}_untimed_sha256"] for r in block_rows}) == 1
                for arm in ("ic", "rho")
            },
            "ic_replay": {k: ic_field.get(k) for k in ("records", "pass", "fail")},
            "ic_replay_all_pass": bool((replays.get("ic") or {}).get("all_pass"))
            and ic_field.get("records") == targets * blocks,
            "rho_replay": {
                "records": sum(r.get("records", 0) for r in rho_reports),
                "pass": sum(r.get("pass", 0) for r in rho_reports),
                "fail": sum(r.get("fail", 0) for r in rho_reports),
            },
            "rho_replay_all_pass": len(rho_reports) == blocks and all(r.get("all_pass") for r in rho_reports),
        }
        verified_ok = producers_ok and all(
            verification[key] for key in (
                "ic_all_verified", "rho_all_verified", "scalars_match_corpus",
                "same_target_points_all_blocks", "ic_replay_all_pass", "rho_replay_all_pass",
            )
        )

        advance("analysis")
        ic_summary0 = block_rows[0].get("ic_summary") or {}
        summary = {
            "schema_version": "1.0",
            "beat_id": beat_id,
            "run_id": run_id,
            "n": n, "a": a, "targets": targets, "blocks": blocks, "K": k_choice,
            "subgroup_order_bits": None,
            "comparator": beat["comparator_status"],
            "corpora": corpora | {"base_hash": ic_summary0.get("base_hash")},
            "executables": executables,
            "host": host_record() | {
                "cpu": subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                      capture_output=True, text=True).stdout.strip()
                if sys.platform == "darwin" else platform.processor(),
            },
            "isolation": (
                f"children pinned to CPU {arguments.cpu} with taskset; wrap this command in "
                "tools/isolated_bench.py reserve for section-10 evidence"
                if arguments.cpu is not None
                else "none: unpinned shared host; per-run load, memory and swap recorded"
            ),
            "k_tune": read_json(run / "artifacts/k_tune.json"),
            "blocks": [
                {k: v for k, v in row.items() if k not in ("ic_summary", "rho_summary")}
                | {"ic_peak_rss_bytes": (row.get("ic_summary") or {}).get("peak_rss_bytes")}
                for row in block_rows
            ],
            "wall_ratio": median_range([r["wall_s_ratio"] for r in block_rows if r["wall_s_ratio"]]),
            "user_ratio": median_range([r["user_s_ratio"] for r in block_rows if r["user_s_ratio"]]),
            "instructions_ratio": median_range(
                [r["instructions_retired_ratio"] for r in block_rows if r["instructions_retired_ratio"]]
            ),
            "verification": verification,
        }
        base_rows = read_jsonl(base_path)[:1]
        if base_rows and base_rows[0].get("subgroup_order"):
            summary["subgroup_order_bits"] = round(math.log2(int(base_rows[0]["subgroup_order"])), 2)
        write_json(run / "artifacts/panel_summary.json", summary)

        end_to_end, vs_rho = draft_panel_claims(beat_id=beat_id, beat=beat, run=run, summary=summary)
        write_json(run / "artifacts/claim_draft.json", end_to_end)
        validation = validate_claim(end_to_end, stage="end_to_end_dlp", ledger=ledger)
        write_json(run / "artifacts/claim_check.json", validation)
        write_json(run / "artifacts/claim_draft_vs_rho.json", vs_rho)
        vs_rho_validation = validate_claim(vs_rho, stage="vs_rho", ledger=ledger)
        write_json(run / "artifacts/claim_check_vs_rho.json", vs_rho_validation)

        if not producers_ok:
            status_value = "PRODUCER_FAILURE"
        elif not verified_ok:
            status_value = "VERIFICATION_FAILURE"
        elif validation["status"] != "PASS":
            status_value = "SCHEMA_INCOMPLETE"
        else:
            status_value = "PENDING_INDEPENDENT_VALIDATION"
        state.update(
            status=status_value, phase="done", updated_at=now(), chosen_K=k_choice,
            claim_check=validation["status"], claim_check_vs_rho=vs_rho_validation["status"],
            wall_ratio=summary["wall_ratio"], producer_exit_codes=exit_codes,
        )
        write_json(run / "state.json", state)
        write_json(
            run / "artifacts/candidate.json",
            {"schema_version": "1.0", "task_id": TASK_ID, "run_id": run_id, "beat_id": beat_id,
             "status": status_value, "claim_draft_sha256": sha256(run / "artifacts/claim_draft.json"),
             "claim_check": validation, "claim_check_vs_rho": vs_rho_validation,
             "chosen_K": k_choice, "tune_targets": tune_targets, "tune_corpus": corpora["tune"],
             "tune_eval_disjoint": corpora["disjoint"],
             "ledger_sha256": sha256(ledger_path(protocol)), "created_at": now(),
             "note": "Multi-target batch diagnostic; never a vs_rho ledger promotion."},
        )
        write_panel_manifest(run)
        return state


def write_panel_manifest(run: Path) -> None:
    manifest = run / "artifacts/review_manifest.json"
    files = {
        str(path.relative_to(run)): sha256(path)
        for path in sorted(run.rglob("*")) if path.is_file() and path != manifest
    }
    write_json(manifest, {"schema_version": "1.0", "task_id": TASK_ID, "files": files})


def status(arguments: argparse.Namespace) -> None:
    run = resolve_run(arguments.run_id)
    state = read_json(run / "state.json")
    pairing_audit_path = run / "artifacts/paired_target_audit.json"
    if pairing_audit_path.is_file():
        audit = read_json(pairing_audit_path)
        if audit.get("status") == "PAIRING_REJECTED" or audit.get("speedup_claim_valid") is False:
            state["recorded_status"] = state.get("status")
            state["recorded_claim_check"] = state.get("claim_check")
            state["status"] = "PAIRING_FAILURE"
            state["claim_check"] = "FAIL"
            state["pairing_audit_status"] = audit.get("status")
            state["pairing_errors"] = [audit.get("reason") or "paired-target audit rejected the run"]
            state["status_note"] = (
                "Effective status includes the post-run pairing audit; "
                "the original state.json remains unchanged."
            )
    print(json.dumps(state, indent=2, sort_keys=True) + "\n", end="")


def verify(arguments: argparse.Namespace) -> dict[str, Any]:
    run = resolve_run(arguments.run_id)
    manifest_path = run / "artifacts/review_manifest.json"
    require(manifest_path.is_file(), f"manifest missing in {run}")
    manifest = read_json(manifest_path)
    mismatches = [
        relative
        for relative, digest in manifest["files"].items()
        if not (run / relative).is_file() or sha256(run / relative) != digest
    ]
    result = {
        "schema_version": "1.0",
        "task_id": TASK_ID,
        "run_id": run.name,
        "status": "PASS" if not mismatches else "FAIL",
        "files": len(manifest["files"]),
        "mismatches": mismatches,
    }
    write_json(run / "artifacts/verification.json", result)
    return result


def claim_check(arguments: argparse.Namespace) -> dict[str, Any]:
    protocol = load_protocol()
    ledger = load_ledger(protocol)
    report = read_json(Path(arguments.report))
    stage = arguments.stage or report.get("stage")
    require(bool(stage), "stage required via --stage or report.stage")
    result = validate_claim(report, stage=str(stage), ledger=ledger)
    if arguments.out:
        write_json(Path(arguments.out), result)
    return result


def launch_equal_advice(arguments: argparse.Namespace) -> dict[str, Any]:
    """Diagnostic Bernstein–Lange equal-advice reference arm (not KS batch rho).

    Drives ``ic_next10_09_equal_advice`` via its ``targets`` / ``bl`` subcommands.
    Optional ``--reuse-bl-archive`` + ``--pair-ic`` reuses a hash-verified BL
    corpus and measures compact-orbit IC online on the same public points.
    Fail-closed for vs_rho ledger promotion: claim draft stays SCHEMA_INCOMPLETE
    until a schema-complete equal-advice paired report exists.
    """
    import equal_advice_pair

    protocol = load_protocol()
    ledger = load_ledger(protocol)
    beat_id = arguments.beat
    require(beat_id in protocol["beats"], f"unknown beat id: {beat_id}")
    beat = protocol["beats"][beat_id]
    require(
        beat.get("launch_mode") == "equal_advice_bl",
        f"{beat_id} is not an equal_advice_bl beat",
    )
    targets = int(arguments.targets or beat.get("targets_default", 1024))
    minimum = int(beat.get("targets_minimum", 1024))
    require(targets >= minimum, f"equal-advice targets must be >= {minimum} (no small-L panels)")
    table_bytes = int(arguments.table_bytes or beat.get("table_bytes_default", 5711008))
    label = str(arguments.label or beat.get("target_label_default", "autolab-n61-bl-ref-eval-v1"))
    pair_ic = bool(getattr(arguments, "pair_ic", False))
    reuse_archive = getattr(arguments, "reuse_bl_archive", None)
    ic_k = int(getattr(arguments, "ic_k", None) or beat.get("ic_k_default", 256))
    ic_rank_seed = int(getattr(arguments, "ic_rank_seed", None) or beat.get("ic_rank_seed", 7))
    ic_target_count = getattr(arguments, "ic_target_count", None)
    n, a = int(beat["n"]), int(beat["a"])
    producers = beat["panel_producers"]
    bl_prod = producers["bl"]
    ic_prod = producers["ic"]
    if pair_ic:
        require("ic" in producers, f"{beat_id} has no IC panel producer")
    if reuse_archive:
        require(pair_ic or arguments.prepare_only, "--reuse-bl-archive requires --pair-ic (or --prepare-only)")

    with RunnerLock():
        # prepare-only must work in the Python-only harness (no cargo/rustc).
        # Full execution still requires a complete preflight.
        if arguments.prepare_only:
            checks: list[dict[str, Any]] = []
            for arm, producer in producers.items():
                source = REPO / producer["source"]
                checks.append(
                    {
                        "name": f"equal_advice_producer_source_{arm}",
                        "ok": source.is_file(),
                        "detail": str(source),
                    }
                )
            git_head_proc = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=REPO,
                capture_output=True,
                text=True,
                check=False,
            )
            checks.append(
                {
                    "name": "git_head",
                    "ok": git_head_proc.returncode == 0,
                    "detail": git_head_proc.stdout.strip(),
                }
            )
            if reuse_archive:
                try:
                    archive = equal_advice_pair.resolve_bl_archive(Path(reuse_archive))
                    checks.append(
                        {
                            "name": "reuse_bl_archive",
                            "ok": True,
                            "detail": str(archive),
                        }
                    )
                except (OSError, ValueError) as exc:
                    checks.append(
                        {
                            "name": "reuse_bl_archive",
                            "ok": False,
                            "detail": str(exc),
                        }
                    )
            preflight_receipt = {
                "schema_version": "1.0",
                "task_id": TASK_ID,
                "kind": "equal_advice_prepare_only_preflight",
                "ok": all(c["ok"] for c in checks),
                "checks": checks,
                "claim_boundary": (
                    "prepare-only scaffolding; cargo/rustc not required; "
                    "not a fixed-arm measurement"
                ),
            }
        else:
            preflight_receipt = preflight(protocol, ledger)
            for arm, producer in producers.items():
                source = REPO / producer["source"]
                preflight_receipt["checks"].append(
                    {
                        "name": f"equal_advice_producer_source_{arm}",
                        "ok": source.is_file(),
                        "detail": str(source),
                    }
                )
            preflight_receipt["ok"] = all(
                c["ok"]
                for c in preflight_receipt["checks"]
                if c["name"] != "cryptominisat5_optional"
            )
        if not preflight_receipt["ok"]:
            failed = [c for c in preflight_receipt["checks"] if not c["ok"]]
            require(False, "preflight failed: " + "; ".join(f"{c['name']}={c['detail']}" for c in failed))
        git_head = next(c["detail"] for c in preflight_receipt["checks"] if c["name"] == "git_head")
        run_id = arguments.run_id or make_run_id(beat_id)
        run = RUNS_DIR / run_id
        require(not run.exists(), f"run already exists: {run}")
        for name in ("artifacts", "inputs", "logs", "receipts"):
            (run / name).mkdir(parents=True, exist_ok=True)
        write_json(CURRENT_PATH, {"run_id": run_id, "beat_id": beat_id, "updated_at": now()})
        write_json(run / "artifacts/preflight.json", preflight_receipt)
        write_json(run / "inputs/protocol.json", protocol)
        write_json(run / "inputs/boundary_targets.json", ledger)
        write_json(
            run / "inputs/ledger_pin.json",
            {
                "path": protocol["ledger"]["path"],
                "sha256": sha256(ledger_path(protocol)),
                "schema_version": ledger["schema_version"],
            },
        )

        points_path = run / "inputs/points.txt"
        scalars_path = run / "inputs/scalars.txt"
        bl_out = run / "logs/bl.jsonl"
        # Prefer CARGO_TARGET_DIR when set (autolab worktrees often redirect builds).
        cargo_target = Path(os.environ.get("CARGO_TARGET_DIR", str(REPO / "target")))
        bl_bin = str((cargo_target / "release/examples" / bl_prod["example"]).resolve())
        ic_bin = str((cargo_target / "release/examples" / ic_prod["example"]).resolve())
        targets_cmd = [
            bl_bin,
            "targets",
            str(n),
            str(a),
            str(targets),
            label,
            str(points_path),
            str(scalars_path),
        ]
        bl_cmd = [
            bl_bin,
            "bl",
            str(n),
            str(a),
            str(table_bytes),
            str(points_path),
            str(bl_out),
        ]
        build_examples = [bl_prod["example"]]
        if pair_ic:
            build_examples.append(ic_prod["example"])
        commands: dict[str, Any] = {
            "build": ["cargo", "build", "--release"]
            + [item for example in build_examples for item in ("--example", example)],
            "targets": targets_cmd,
            "bl": bl_cmd,
        }
        if pair_ic:
            commands["ic"] = [
                ic_bin,
                f"construct:{n}:{a}:{ic_k}",
                str(run / "inputs/ic_points.txt"),
                str(ic_rank_seed),
                str(run / "logs/ic.jsonl"),
            ]
        if reuse_archive:
            commands["reuse_bl_archive"] = str(Path(reuse_archive).resolve())
        write_json(run / "artifacts/commands.json", commands)

        state: dict[str, Any] = {
            "schema_version": "1.0",
            "task_id": TASK_ID,
            "run_id": run_id,
            "beat_id": beat_id,
            "launch_mode": "equal_advice_bl",
            "status": "PREPARE_ONLY" if arguments.prepare_only else "ACTIVE",
            "phase": "prepare",
            "n": n,
            "a": a,
            "targets": targets,
            "table_bytes": table_bytes,
            "target_label": label,
            "pair_ic": pair_ic,
            "reuse_bl_archive": bool(reuse_archive),
            "ic_k": ic_k if pair_ic else None,
            "ic_rank_seed": ic_rank_seed if pair_ic else None,
            "diagnostic_only": True,
            "reference_arm": "bernstein_lange_equal_advice",
            "existing_study": beat.get("existing_study"),
            "git_head": git_head,
            "created_at": now(),
            "updated_at": now(),
            "claim_boundary": (
                "Diagnostic Bernstein–Lange equal-advice reference arm; public synthetic "
                "only; fail-closed for vs_rho ledger promotion; PENDING_INDEPENDENT_VALIDATION "
                "at most; not a KS v2 batch comparator and not an asymptotic sub-rho claim."
            ),
        }
        write_json(run / "state.json", state)
        if arguments.prepare_only:
            write_json(
                run / "artifacts/claim_draft.json",
                {
                    "beat_id": beat_id,
                    "stage": "vs_rho",
                    "status": "SCHEMA_INCOMPLETE",
                    "diagnostic_only": True,
                    "reason": "prepare-only; producers not executed",
                    "n": n,
                    "targets": targets,
                    "table_bytes": table_bytes,
                    "pair_ic": pair_ic,
                    "reuse_bl_archive": bool(reuse_archive),
                    "claim_boundary": state["claim_boundary"],
                },
            )
            write_json(
                run / "artifacts/claim_check.json",
                {
                    "status": "FAIL",
                    "stage": "vs_rho",
                    "reason": "prepare-only scaffolding; vs_rho remains fail-closed",
                },
            )
            return state

        pin: dict[str, Any] | None = None
        if reuse_archive:
            state.update(phase="reuse_bl", updated_at=now())
            write_json(run / "state.json", state)
            try:
                pin = equal_advice_pair.materialize_bl_archive(
                    archive=Path(reuse_archive),
                    run=run,
                    lab=sys.modules[__name__],
                    targets=targets,
                    table_bytes=table_bytes,
                )
            except (OSError, ValueError, FileNotFoundError) as exc:
                require(False, f"reuse-bl-archive failed: {exc}")
        else:
            state.update(phase="build", updated_at=now())
            write_json(run / "state.json", state)
            built = subprocess.run(commands["build"], cwd=REPO, capture_output=True, text=True)
            require(built.returncode == 0, "cargo build failed:\n" + built.stderr[-4000:])
            write_json(
                run / "artifacts/binaries.json",
                {
                    "paths": {"bl": bl_bin, **({"ic": ic_bin} if pair_ic else {})},
                    "hashes": {
                        "bl_binary_sha256": sha256(bl_bin),
                        "bl_source_sha256": sha256(REPO / bl_prod["source"]),
                        **(
                            {
                                "ic_binary_sha256": sha256(ic_bin),
                                "ic_source_sha256": sha256(REPO / ic_prod["source"]),
                            }
                            if pair_ic
                            else {}
                        ),
                        "git_head": git_head,
                    },
                },
            )

            state.update(phase="targets", updated_at=now())
            write_json(run / "state.json", state)
            gen = subprocess.run(targets_cmd, cwd=REPO, capture_output=True, text=True)
            (run / "logs/targets.stdout").write_text(gen.stdout)
            (run / "logs/targets.stderr").write_text(gen.stderr)
            require(gen.returncode == 0, "targets generation failed:\n" + gen.stderr[-2000:])
            require(points_path.is_file() and scalars_path.is_file(), "targets mode did not write point/scalar files")
            point_lines = [ln for ln in points_path.read_text().splitlines() if ln.strip()]
            require(len(point_lines) == targets, f"expected {targets} points, got {len(point_lines)}")

            state.update(phase="bl", updated_at=now())
            write_json(run / "state.json", state)
            bl_run = subprocess.run(bl_cmd, cwd=REPO, capture_output=True, text=True)
            (run / "logs/bl.stdout").write_text(bl_run.stdout)
            (run / "logs/bl.stderr").write_text(bl_run.stderr)
            require(bl_run.returncode == 0, "BL producer failed:\n" + bl_run.stderr[-4000:])
            require(bl_out.is_file(), "BL producer did not write jsonl")

        records = parse_json_lines(bl_out.read_text())
        summary = next((r for r in records if r.get("kind") == "bl_summary" or r.get("kind") == "bl_level_summary"), None)
        precompute = next((r for r in records if r.get("kind") == "bl_precompute"), None)
        online_rows = [
            r
            for r in records
            if r.get("kind") in ("bl_online", "bl_online_target", "bl_target", "bl_level_target")
        ]
        solved = sum(
            1
            for r in online_rows
            if r.get("solved") or r.get("verified") or r.get("verified_vs_truth")
        )

        claim_draft: dict[str, Any]
        check: dict[str, Any]
        if pair_ic:
            state.update(phase="build_ic", updated_at=now())
            write_json(run / "state.json", state)
            built = subprocess.run(commands["build"], cwd=REPO, capture_output=True, text=True)
            require(built.returncode == 0, "cargo build failed:\n" + built.stderr[-4000:])
            binaries = {"bl": bl_bin, "ic": ic_bin}
            write_json(
                run / "artifacts/binaries.json",
                {
                    "paths": binaries,
                    "hashes": {
                        "bl_binary_sha256": sha256(bl_bin) if Path(bl_bin).is_file() else None,
                        "ic_binary_sha256": sha256(ic_bin),
                        "bl_source_sha256": sha256(REPO / bl_prod["source"]),
                        "ic_source_sha256": sha256(REPO / ic_prod["source"]),
                        "git_head": git_head,
                    },
                },
            )
            points = equal_advice_pair.parse_points_file(points_path)
            scalars = equal_advice_pair.parse_scalars_file(scalars_path)
            count = int(ic_target_count) if ic_target_count is not None else len(points)
            require(1 <= count <= len(points), f"ic_target_count must be in 1..{len(points)}")
            state.update(phase="ic", updated_at=now(), ic_target_count=count)
            write_json(run / "state.json", state)
            try:
                ic_obs = equal_advice_pair.run_ic_arm(
                    lab=sys.modules[__name__],
                    run=run,
                    binaries=binaries,
                    n=n,
                    a=a,
                    k=ic_k,
                    rank_seed=ic_rank_seed,
                    points=points,
                    ic_target_count=count,
                )
            except (OSError, RuntimeError, ValueError) as exc:
                require(False, f"paired IC arm failed: {exc}")
            write_json(
                run / "artifacts/ic_arm.json",
                {
                    "exit_code": ic_obs["exit_code"],
                    "ic_target_count": ic_obs["ic_target_count"],
                    "k": ic_obs["k"],
                    "rank_seed": ic_obs["rank_seed"],
                    "points_sha256": ic_obs["points_sha256"],
                    "ic_jsonl_sha256": ic_obs["ic_jsonl_sha256"],
                    "summary": ic_obs["summary"],
                },
            )
            if pin is None:
                pin = {
                    "kind": "fresh_bl_measurement",
                    "archive_run_id": None,
                    "points_sha256": sha256(points_path),
                    "scalars_sha256": sha256(scalars_path),
                    "bl_jsonl_sha256": sha256(bl_out),
                    "remeasured": True,
                }
            claim_draft = equal_advice_pair.draft_paired_vs_rho_slice(
                lab=sys.modules[__name__],
                beat=beat,
                beat_id=beat_id,
                run=run,
                pin=pin,
                ic_obs=ic_obs,
                bl_records=records,
                points=points,
                scalars=scalars,
                table_bytes=table_bytes,
                binaries=binaries,
                git_head=git_head,
            )
            claim_draft["bl_records"] = len(records)
            claim_draft["bl_online_rows"] = len(online_rows)
            claim_draft["bl_solved_or_verified_rows"] = solved
            claim_draft["bl_precompute"] = precompute
            claim_draft["bl_summary"] = summary
            claim_draft["target_label"] = label
            panel = equal_advice_pair.panel_pairing_summary(
                ic_obs=ic_obs,
                bl_records=records,
                points=points,
                scalars=scalars,
                table_bytes=table_bytes,
            )
            write_json(run / "artifacts/panel_pairing.json", panel)
            claim_draft["panel_pairing"] = {
                "paired_targets": panel["paired_targets"],
                "verified_pairs": panel["verified_pairs"],
                "online_speedup_bl_over_ic": panel["online_speedup_bl_over_ic"],
            }
            write_json(run / "artifacts/claim_draft.json", claim_draft)
            check = validate_claim(claim_draft, stage="vs_rho", ledger=ledger)
            write_json(run / "artifacts/claim_check.json", check)
            write_json(run / "artifacts/claim_check_tool.json", check)
        else:
            claim_draft = {
                "beat_id": beat_id,
                "stage": "vs_rho",
                "status": "SCHEMA_INCOMPLETE",
                "diagnostic_only": True,
                "reference_arm": "bernstein_lange_equal_advice",
                "n": n,
                "a": a,
                "targets": targets,
                "table_bytes": table_bytes,
                "target_label": label,
                "bl_records": len(records),
                "bl_online_rows": len(online_rows),
                "bl_solved_or_verified_rows": solved,
                "bl_precompute": precompute,
                "bl_summary": summary,
                "existing_study": beat.get("existing_study"),
                "timing_class": beat.get("timing_class_goal"),
                "claim_boundary": state["claim_boundary"],
                "missing_vs_rho_fields_note": (
                    "Equal-advice BL reference smoke does not emit the full vs_rho "
                    "paired IC/rho schema; fail-closed until a schema-complete report."
                ),
                "independent_replay_pointer": None,
                "verdict": "DIAGNOSTIC_REFERENCE_ARM_ONLY",
            }
            write_json(run / "artifacts/claim_draft.json", claim_draft)
            check = {
                "status": "FAIL",
                "stage": "vs_rho",
                "reason": "diagnostic equal-advice BL reference; vs_rho schema incomplete by design",
                "claim_draft_status": claim_draft["status"],
            }
            write_json(run / "artifacts/claim_check.json", check)

        state.update(
            status="PENDING_INDEPENDENT_VALIDATION",
            phase="done",
            updated_at=now(),
            claim_check_status=check["status"],
            claim_draft_status=claim_draft.get("status"),
            remaining_schema_gaps=claim_draft.get("remaining_schema_gaps"),
        )
        write_json(run / "state.json", state)
        write_json(CURRENT_PATH, {"run_id": run_id, "beat_id": beat_id, "updated_at": now()})
        return state


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    sub = result.add_subparsers(dest="command", required=True)

    sub.add_parser("plan", help="Show ledger priorities and beat launch commands")
    sub.add_parser("preflight", help="Check ledger, deps, and producer sources")

    launch_parser = sub.add_parser("launch", help="Create a run and execute producers")
    launch_parser.add_argument(
        "--beat",
        required=True,
        help="Beat id from protocol.json (e.g. koblitz.vs_rho.n37_wall)",
    )
    launch_parser.add_argument("--run-id")
    launch_parser.add_argument(
        "--fixtures",
        type=int,
        help="Compatibility option; the only accepted value is 1 target per workload",
    )
    launch_parser.add_argument(
        "--repeats",
        type=int,
        help=(
            "Timed executions per arm after a discarded warmup exec "
            f"(default {DEFAULT_TIMED_REPEATS}; stops early once a run "
            f"exceeds {int(REPEAT_BUDGET_MS / 1000)}s)"
        ),
    )
    launch_parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Create run scaffolding and print commands without executing producers",
    )

    panel = sub.add_parser(
        "launch-panel",
        help="Run a multi-target batch panel beat (diagnostic only; never a vs_rho promotion)",
    )
    panel.add_argument("--beat", required=True, help="Batch panel beat id from protocol.json")
    panel.add_argument("--run-id")
    panel.add_argument("--targets", type=int, help="Targets per batch (default: beat targets_default)")
    panel.add_argument("--blocks", type=int, help="Paired blocks (default: beat blocks_default)")
    panel.add_argument("--k", type=int, help="Use this K and skip the tune")
    panel.add_argument("--k-candidates", help="Comma-separated K tune candidates (default: beat list for this L)")
    panel.add_argument("--tune-targets", type=int,
                       help="Targets in the disjoint tune corpus (default: the panel L; e.g. 1024)")
    panel.add_argument("--tune-only", action="store_true",
                       help="Stop after the K tune (status TUNE_ONLY); --resume the run to continue")
    panel.add_argument("--cpu", type=int, help="Pin every producer to this CPU with taskset (Linux)")
    panel.add_argument("--resume", metavar="RUN_ID",
                       help="Continue an interrupted panel run; finished K rows and blocks are kept")

    single = sub.add_parser(
        "launch-single",
        help="Run W one-target workloads (fresh processes per arm) for a single-target panel beat",
    )
    single.add_argument("--beat", required=True, help="Single-target panel beat id from protocol.json")
    single.add_argument("--run-id")
    single.add_argument("--workloads", type=int, help="One-target eval workloads (default: beat workloads_default)")
    single.add_argument("--tune-workloads", type=int, help="Disjoint one-target tune workloads per K")
    single.add_argument("--k", type=int, help="Use this K and skip the tune")
    single.add_argument("--k-candidates", help="Comma-separated K tune candidates (default: beat list)")
    single.add_argument("--cpu", type=int, help="Pin every producer to this CPU with taskset (Linux)")
    single.add_argument("--resume", metavar="RUN_ID", help="Continue an interrupted single-target run")

    equal = sub.add_parser(
        "launch-equal-advice",
        help=(
            "Diagnostic Bernstein–Lange equal-advice reference arm "
            "(ic_next10_09_equal_advice); not a KS batch panel; fail-closed vs_rho"
        ),
    )
    equal.add_argument("--beat", required=True, help="equal_advice_bl beat id from protocol.json")
    equal.add_argument("--run-id")
    equal.add_argument("--targets", type=int, help="Public synthetic targets (default: beat targets_default, >=1024)")
    equal.add_argument("--table-bytes", type=int, help="BL table budget bytes (default: beat table_bytes_default)")
    equal.add_argument("--label", help="Target corpus label for targets mode")
    equal.add_argument(
        "--prepare-only",
        action="store_true",
        help="Scaffold run dir + commands.json without building or executing producers",
    )
    equal.add_argument(
        "--reuse-bl-archive",
        metavar="PATH",
        help=(
            "Hash-verify and reuse an archived BL run (inputs/points.txt + logs/bl.jsonl) "
            "instead of remeasuring BL; requires --pair-ic"
        ),
    )
    equal.add_argument(
        "--pair-ic",
        action="store_true",
        help="Measure koblitz_orbit_dlp_fast_online on the same public points as the BL arm",
    )
    equal.add_argument(
        "--ic-k",
        type=int,
        help="Orbit columns K for construct:n:a:K (default: beat ic_k_default / 256 for n61 equal-advice)",
    )
    equal.add_argument(
        "--ic-rank-seed",
        type=int,
        help="Rank seed for the IC online producer (default: beat ic_rank_seed / 7)",
    )
    equal.add_argument(
        "--ic-target-count",
        type=int,
        help="How many leading archived points to solve with IC (default: all)",
    )

    for name in ("status", "verify"):
        command = sub.add_parser(name)
        command.add_argument("--run-id")

    claim = sub.add_parser(
        "claim-check",
        help="Fail-closed validate a JSON report against measurement_schema",
    )
    claim.add_argument("--report", required=True)
    claim.add_argument("--stage", help="Stage id (default: report.stage)")
    claim.add_argument("--out", help="Optional path to write the check receipt")
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        if arguments.command == "plan":
            protocol = load_protocol()
            ledger = load_ledger(protocol)
            print(json.dumps(plan(protocol, ledger), indent=2, sort_keys=True))
            return 0
        if arguments.command == "preflight":
            protocol = load_protocol()
            ledger = load_ledger(protocol)
            receipt = preflight(protocol, ledger)
            print(json.dumps(receipt, indent=2, sort_keys=True))
            return 0 if receipt["ok"] else 2
        if arguments.command == "launch":
            state = launch(arguments)
            print(json.dumps(state, indent=2, sort_keys=True))
            return 0 if state.get("status") != "PRODUCER_FAILURE" else 3
        if arguments.command == "launch-panel":
            state = launch_panel(arguments)
            print(json.dumps(state, indent=2, sort_keys=True))
            return 0 if state.get("status") in ("PENDING_INDEPENDENT_VALIDATION", "TUNE_ONLY") else 3
        if arguments.command == "launch-single":
            import single_target_panel

            state = single_target_panel.launch_single(arguments, sys.modules[__name__])
            print(json.dumps(state, indent=2, sort_keys=True))
            return 0 if state.get("status") == "PENDING_INDEPENDENT_VALIDATION" else 3
        if arguments.command == "launch-equal-advice":
            state = launch_equal_advice(arguments)
            print(json.dumps(state, indent=2, sort_keys=True))
            return 0 if state.get("status") in (
                "PENDING_INDEPENDENT_VALIDATION",
                "PREPARE_ONLY",
            ) else 3
        if arguments.command == "status":
            status(arguments)
            return 0
        if arguments.command == "verify":
            print(json.dumps(verify(arguments), indent=2, sort_keys=True))
            return 0
        if arguments.command == "claim-check":
            result = claim_check(arguments)
            print(json.dumps(result, indent=2, sort_keys=True))
            return 0 if result["status"] == "PASS" else 4
        raise AutolabError(f"unknown command {arguments.command}")
    except AutolabError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
