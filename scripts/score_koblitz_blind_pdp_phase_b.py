#!/usr/bin/env python3
"""Open the sealed Phase-A truth only after a Phase-B solver run is frozen."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import stat
import statistics
from typing import Any

import run_koblitz_blind_pdp_phase_b as phase_b


SCORE_SCHEMA = "koblitz_pdp_phase_b_score.v1"
SCORE_SEAL_SCHEMA = "koblitz_pdp_phase_b_score_seal.v1"
ORACLE_SCHEMA = "koblitz_pdp_sealed_oracle_ledger.v1"
TERMINAL_STATUSES = {"sat", "unsat"}
KNOWN_SOLVER_STATUSES = {
    "sat",
    "unsat",
    "unknown_inconclusive",
    "timeout_inconclusive",
    "model_cap_inconclusive",
    "sat_nonlifting_model_inconclusive",
    "sat_invalid_model",
    "sat_invalid_terminal_status",
    "unsat_invalid_terminal_status",
    "sat_point_validation_error",
    "sat_point_validation_timeout_inconclusive",
    "solver_error",
    "backend_error",
    "backend_contract_error",
}


def validate_build_provenance(plan: dict, protocol: dict, *, allow_smoke: bool) -> dict:
    import build_koblitz_phase_b_tools as tool_builder

    tool_builds = plan.get("additional_tool_builds", {})
    if not isinstance(tool_builds, dict) or set(tool_builds) - {"rust", "cryptominisat"}:
        raise phase_b.PhaseBError("Phase-B plan has unexpected tool build receipts")
    if not allow_smoke and set(tool_builds) != {"rust", "cryptominisat"}:
        raise phase_b.PhaseBError("production scoring requires bound Rust and CryptoMiniSat build receipts")
    source_state = plan.get("source_revision")
    if not isinstance(source_state, dict):
        raise phase_b.PhaseBError("Phase-B plan lacks implementation provenance")
    phase_b.require_hex40(source_state.get("commit"), "Phase-B implementation revision")
    evidence_class = plan.get("evidence_class")
    if evidence_class not in {"operational_smoke", phase_b.PRODUCTION_EVIDENCE_CLASS}:
        raise phase_b.PhaseBError("Phase-B plan has an unknown evidence class")
    if evidence_class != "operational_smoke" and (
        source_state.get("dirty") is not False or source_state.get("porcelain") != []
        or plan.get("allow_dirty_requested") is not False
        or set(tool_builds) != {"rust", "cryptominisat"}
    ):
        raise phase_b.PhaseBError("Phase-B measurement requires a clean implementation and bound builds")
    if "rust" in tool_builds:
        expected_objects = tool_builder.rust_source_objects(phase_b.REPO, source_state["commit"])
        if expected_objects != plan.get("rust_build_source_objects"):
            raise phase_b.PhaseBError("Rust build provenance differs from the recorded implementation revision")
        tool_builder.git(phase_b.REPO, "merge-base", "--is-ancestor", protocol["phase_a_binding"]["source_revision"], source_state["commit"])
    return tool_builds


def validate_process_summary(process: Any, context: str) -> None:
    if not isinstance(process, dict):
        raise phase_b.PhaseBError(f"{context} lacks its process receipt")
    if (
        not isinstance(process.get("command"), list)
        or not all(isinstance(value, str) for value in process["command"])
        or not isinstance(process.get("returncode"), int)
        or not isinstance(process.get("timed_out"), bool)
        or not isinstance(process.get("orphan_group_terminated"), bool)
    ):
        raise phase_b.PhaseBError(f"{context} process terminal fields are invalid")
    if process["orphan_group_terminated"]:
        raise phase_b.PhaseBError(f"{context} left descendant processes after its leader exited")
    metrics = process.get("metrics")
    required = {
        "wall_seconds",
        "user_seconds",
        "system_seconds",
        "total_core_seconds",
        "single_core_seconds",
        "peak_rss_bytes",
        "meter",
    }
    if not isinstance(metrics, dict) or set(metrics) != required:
        raise phase_b.PhaseBError(f"{context} process resources are incomplete")
    for name in required - {"meter", "peak_rss_bytes"}:
        value = metrics[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise phase_b.PhaseBError(f"{context} has invalid numeric process resources")
    if (
        not isinstance(metrics["peak_rss_bytes"], int)
        or isinstance(metrics["peak_rss_bytes"], bool)
        or metrics["peak_rss_bytes"] < 0
        or metrics["meter"] != "fresh-process getrusage(RUSAGE_CHILDREN)"
        or not math.isclose(
            metrics["total_core_seconds"],
            metrics["user_seconds"] + metrics["system_seconds"],
            rel_tol=0,
            abs_tol=1e-9,
        )
        or metrics["single_core_seconds"] != metrics["total_core_seconds"]
    ):
        raise phase_b.PhaseBError(f"{context} process resources are inconsistent")


def sealed_stdout(row: dict[str, Any], task_root: Path, context: str) -> str:
    process = row["process"]
    relative = Path(process.get("stdout_path", ""))
    if relative.is_absolute() or len(relative.parts) != 1 or relative.parts[0] in {"", ".", ".."}:
        raise phase_b.PhaseBError(f"{context} has an unsafe stdout path")
    data = phase_b.regular_file_bytes(task_root / relative, f"{context} stdout")
    if (
        len(data) != process.get("stdout_bytes")
        or phase_b.sha256_bytes(data) != process.get("stdout_sha256")
    ):
        raise phase_b.PhaseBError(f"{context} stdout differs from its process receipt")
    return data.decode(errors="replace")


def validate_backend_admission(
    row: dict[str, Any],
    blind_id: str,
    task: dict[str, Any] | None = None,
    task_root: Path | None = None,
) -> None:
    backend = row["solver"]
    status = row.get("status")
    validate_process_summary(row.get("process"), f"{blind_id} {backend}")
    process = row["process"]
    if process["timed_out"] and status != "timeout_inconclusive":
        raise phase_b.PhaseBError(f"{blind_id} {backend} relabels a timeout")
    if status == "sat":
        if row.get("source_model_valid") is not True or row.get("source_witness_valid") is not True:
            raise phase_b.PhaseBError(f"{blind_id} {backend} SAT lacks exact model and witness validation")
        if backend == "native-xor":
            report = row.get("backend_report")
            if (
                not isinstance(report, dict)
                or report.get("schema") != "koblitz_pdp_isolated_backend.v1"
                or report.get("backend") != "native-sat"
                or report.get("status") != "sat"
                or report.get("source_instance_verified") is not True
                or report.get("regenerated_source_exact") is not True
            ):
                raise phase_b.PhaseBError(f"{blind_id} native SAT report is not independently admissible")
            if process["returncode"] != 0:
                raise phase_b.PhaseBError(f"{blind_id} native SAT has an invalid return code")
        else:
            validation = row.get("point_witness_validation")
            if not isinstance(validation, dict) or validation.get("status") != "valid_point_witness":
                raise phase_b.PhaseBError(f"{blind_id} {backend} SAT lacks its validator receipt")
            validate_process_summary(
                validation.get("process"), f"{blind_id} {backend} point validation"
            )
            if validation["process"]["returncode"] != 0:
                raise phase_b.PhaseBError(f"{blind_id} {backend} validator did not terminate cleanly")
            expected_returncode = 0 if backend == "wdsat" else 10
            if process["returncode"] != expected_returncode:
                raise phase_b.PhaseBError(f"{blind_id} {backend} SAT has an invalid return code")
    if status == "unsat":
        if backend == "native-xor":
            report = row.get("backend_report")
            if (
                process["returncode"] != 0
                or not isinstance(report, dict)
                or report.get("schema") != "koblitz_pdp_isolated_backend.v1"
                or report.get("backend") != "native-sat"
                or report.get("status") != "unsat"
                or report.get("source_instance_verified") is not True
                or report.get("regenerated_source_exact") is not True
            ):
                raise phase_b.PhaseBError(f"{blind_id} native UNSAT lacks its authenticated backend report")
        elif task is None or task_root is None:
            raise phase_b.PhaseBError(f"{blind_id} {backend} UNSAT lacks sealed terminal context")
        else:
            stdout = sealed_stdout(row, task_root, f"{blind_id} {backend}")
            lines = [line.strip() for line in stdout.splitlines() if line.strip()]
            if backend == "wdsat":
                source_variables = task.get("wdsat_requirements", {}).get("max_anf_id", 0) - 1
                models = [
                    line
                    for line in lines
                    if source_variables > 0
                    and len(line) >= source_variables
                    and set(line) <= {"0", "1"}
                ]
                if (
                    process["returncode"] != 0
                    or lines.count("UNSAT") != 1
                    or any(line in {"UNKNOWN", "s UNKNOWN"} for line in lines)
                    or models
                ):
                    raise phase_b.PhaseBError(f"{blind_id} WDSat UNSAT lacks one exact terminal")
            elif backend == "cryptominisat":
                if (
                    process["returncode"] != 20
                    or lines.count("s UNSATISFIABLE") != 1
                    or any(line in {"s SATISFIABLE", "UNKNOWN", "s UNKNOWN"} for line in lines)
                    or any(line.startswith("v ") for line in lines)
                ):
                    raise phase_b.PhaseBError(
                        f"{blind_id} CryptoMiniSat UNSAT lacks one exact terminal"
                    )
    if backend == "cryptominisat":
        command = process["command"]
        thread_index = command.index("--threads") if "--threads" in command else -1
        if thread_index < 0 or thread_index + 1 >= len(command) or command[thread_index + 1] != "1":
            raise phase_b.PhaseBError(f"{blind_id} CryptoMiniSat did not bind one thread")


def verify_run_tree(
    run_root: Path,
    solver_root: Path,
    protocol: dict[str, Any],
    protocol_bytes: bytes,
    *,
    allow_smoke: bool,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    metadata = run_root.lstat()
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise phase_b.PhaseBError("run root must be a real directory")
    seal, _ = phase_b.read_json(run_root / "run-seal.json", "Phase-B run seal")
    if not isinstance(seal, dict) or seal.get("schema") != phase_b.RUN_SEAL_SCHEMA:
        raise phase_b.PhaseBError("Phase-B run seal schema changed")
    if seal.get("status") != "solver_outputs_frozen":
        raise phase_b.PhaseBError("Phase-B solver outputs are not frozen")
    seal_payload = {key: value for key, value in seal.items() if key != "seal_payload_sha256"}
    phase_b.require_hex64(seal.get("seal_payload_sha256"), "Phase-B run seal self-hash")
    if phase_b.canonical_sha256(seal_payload) != seal["seal_payload_sha256"]:
        raise phase_b.PhaseBError("Phase-B run seal self-hash is invalid")
    inventory = phase_b.all_regular_inventory(run_root, {"run-seal.json"})
    if inventory != seal.get("inventory"):
        raise phase_b.PhaseBError("Phase-B result inventory differs from its terminal seal")
    if phase_b.canonical_sha256(inventory) != seal.get("inventory_sha256"):
        raise phase_b.PhaseBError("Phase-B result inventory self-hash is invalid")
    protocol_hash = phase_b.sha256_bytes(protocol_bytes)
    if seal.get("protocol_sha256") != protocol_hash:
        raise phase_b.PhaseBError("Phase-B run seal uses a different protocol")
    binding, bundle = phase_b.load_solver_root(solver_root, protocol)
    if seal.get("source_phase_a_seal_sha256") != binding["source_phase_a_seal_sha256"]:
        raise phase_b.PhaseBError("Phase-B run seal uses a different Phase-A seal commitment")
    if seal.get("blind_bundle_sha256") != binding["blind_bundle_sha256"]:
        raise phase_b.PhaseBError("Phase-B run seal uses a different blind bundle")
    plan, _ = phase_b.read_json(run_root / "execution-plan.json", "Phase-B execution plan")
    if not isinstance(plan, dict) or plan.get("schema") != phase_b.RUN_PLAN_SCHEMA:
        raise phase_b.PhaseBError("Phase-B execution plan schema changed")
    if plan.get("protocol_sha256") != protocol_hash:
        raise phase_b.PhaseBError("Phase-B execution plan uses a different protocol")
    if plan.get("solver_input_binding") != binding:
        raise phase_b.PhaseBError("Phase-B execution plan uses a different solver input binding")
    if plan.get("tool_build_accounting") != protocol["tool_build_accounting"]:
        raise phase_b.PhaseBError("Phase-B execution plan changed the incomplete tool-build accounting")
    import build_koblitz_phase_b_tools as tool_builder
    import build_koblitz_phase_b_wdsat as wdsat_builder

    tool_builds = validate_build_provenance(plan, protocol, allow_smoke=allow_smoke)
    evidence_class = plan["evidence_class"]
    if plan.get("wdsat_build_capsule_bound") is not True:
        raise phase_b.PhaseBError("Phase-B execution plan lacks its WDSat build capsule binding")
    checked_wdsat = wdsat_builder.validate_capsule(
        run_root / "tool-builds" / "wdsat" / "build-seal.json",
        protocol,
        plan["tool_identities"]["wdsat"],
        plan["source_revision"],
        require_clean_implementation=evidence_class != "operational_smoke",
    )
    if checked_wdsat != plan.get("wdsat_build"):
        raise phase_b.PhaseBError("Phase-B WDSat build capsule differs from the execution binding")
    for name, bound in tool_builds.items():
        checked = tool_builder.validate_receipt(
            run_root / "tool-builds" / name / "receipt.json", name, plan["tool_identities"],
            plan["rust_build_source_objects"] if name == "rust" else None,
            plan["source_revision"],
            require_clean_implementation=evidence_class != "operational_smoke",
        )
        if checked != bound:
            raise phase_b.PhaseBError("Phase-B tool build receipt differs from the execution binding")
    selected_ids = plan.get("selected_blind_instance_ids")
    if not isinstance(selected_ids, list) or not all(isinstance(value, str) for value in selected_ids):
        raise phase_b.PhaseBError("Phase-B execution plan has an invalid selected-id list")
    if len(selected_ids) != len(set(selected_ids)) or len(selected_ids) != plan.get(
        "selected_instance_count"
    ):
        raise phase_b.PhaseBError("Phase-B execution plan selected ids are inconsistent")
    bundle_ids = [item["blind_instance_id"] for item in bundle["instances"]]
    if selected_ids != bundle_ids[: len(selected_ids)]:
        raise phase_b.PhaseBError("Phase-B task selection is not the frozen bundle prefix")
    task_results = []
    for index, blind_id in enumerate(selected_ids):
        task_path = phase_b.task_directory(run_root, index, blind_id) / "task-result.json"
        task, _ = phase_b.read_json(task_path, "Phase-B task result")
        if not isinstance(task, dict) or task.get("schema") != phase_b.TASK_SCHEMA:
            raise phase_b.PhaseBError(f"task result schema changed for {blind_id}")
        if task.get("ordinal") != index or task.get("blind_instance_id") != blind_id:
            raise phase_b.PhaseBError(f"task result identity changed for {blind_id}")
        rows = task.get("backends")
        if not isinstance(rows, list) or [row.get("solver") for row in rows] != list(
            phase_b.BACKENDS
        ):
            raise phase_b.PhaseBError(f"task result backend inventory changed for {blind_id}")
        if (
            task.get("status") != "terminal_outputs_recorded"
            or task.get("source_artifacts_unchanged") is not True
            or task.get("source_artifacts_before") != task.get("source_artifacts_after")
            or task.get("export", {}).get("returncode") != 0
            or task.get("export", {}).get("timed_out") is not False
        ):
            raise phase_b.PhaseBError(f"task result is not a complete frozen export for {blind_id}")
        validate_process_summary(task.get("export"), f"{blind_id} explicit export")
        source_verification = task.get("source_verification")
        source_report = (
            source_verification.get("report") if isinstance(source_verification, dict) else None
        )
        if (
            not isinstance(source_verification, dict)
            or source_verification.get("status") != "verified"
            or source_verification.get("source_instance_id") != task.get("source_instance_id")
            or source_verification.get("source_instance_verified") is not True
            or source_verification.get("regenerated_source_exact") is not True
            or not isinstance(source_report, dict)
            or source_report.get("schema") != "koblitz_pdp_source_verification.v1"
            or source_report.get("source_instance_id") != task.get("source_instance_id")
        ):
            raise phase_b.PhaseBError(f"task lacks isolated algebraic source verification for {blind_id}")
        validate_process_summary(
            source_verification.get("process"), f"{blind_id} source verification"
        )
        for row in rows:
            if row.get("status") not in KNOWN_SOLVER_STATUSES:
                raise phase_b.PhaseBError(f"task result has an unknown solver status for {blind_id}")
            validate_backend_admission(
                row,
                blind_id,
                task,
                phase_b.task_directory(run_root, index, blind_id),
            )
        task_results.append(task)
    if seal.get("selected_instance_count") != len(task_results):
        raise phase_b.PhaseBError("Phase-B run seal selected count changed")
    if seal.get("backend_outcomes") != len(task_results) * len(phase_b.BACKENDS):
        raise phase_b.PhaseBError("Phase-B run seal backend count changed")
    summary, _ = phase_b.read_json(run_root / "run-summary.json", "Phase-B run summary")
    if not isinstance(summary, dict) or summary.get("schema") != phase_b.RUN_SUMMARY_SCHEMA:
        raise phase_b.PhaseBError("Phase-B run summary schema changed")
    independently_summarized = phase_b.summarize_run(
        run_root,
        protocol,
        len(selected_ids),
        len(bundle_ids),
        task_results,
        plan["wdsat_build"],
        tool_builds,
    )
    if summary != independently_summarized:
        raise phase_b.PhaseBError("Phase-B run summary is not independently reproducible")
    independently_full = (
        selected_ids == bundle_ids
        and len(task_results) == len(bundle_ids) == protocol["phase_a_binding"]["instance_count"]
        and plan.get("selected_instance_count") == len(bundle_ids)
        and plan.get("full_instance_count") == len(bundle_ids)
        and plan.get("planned_backend_outcomes") == protocol["execution"]["expected_backend_runs"]
        and seal.get("selected_instance_count") == len(bundle_ids)
        and seal.get("full_instance_count") == len(bundle_ids)
        and seal.get("backend_outcomes") == protocol["execution"]["expected_backend_runs"]
        and summary.get("selected_instances") == len(bundle_ids)
        and summary.get("frozen_instances") == len(bundle_ids)
        and summary.get("task_records") == len(bundle_ids)
        and summary.get("backend_outcomes") == protocol["execution"]["expected_backend_runs"]
        and summary.get("selection_complete") is True
        and summary.get("full_panel_complete") is True
    )
    if seal.get("full_panel_complete") is not independently_full:
        raise phase_b.PhaseBError("Phase-B run seal full-panel label is not independently derived")
    if not allow_smoke and not independently_full:
        raise phase_b.PhaseBError("scoring requires the independently complete full panel")
    if evidence_class != "operational_smoke" and not independently_full:
        raise phase_b.PhaseBError("a partial Phase-B panel cannot carry the full measurement label")
    exports, _ = phase_b.read_json(
        run_root / "export-inventory.json", "Phase-B frozen export inventory"
    )
    if (
        not isinstance(exports, dict)
        or exports.get("schema") != "koblitz_pdp_phase_b_frozen_export_inventory.v1"
        or exports.get("status")
        != "all_selected_explicit_targets_exported_before_solver_execution"
        or exports.get("selected_instances") != len(task_results)
        or exports.get("source_instance_ids")
        != [task["source_instance_id"] for task in task_results]
        or exports.get("wdsat_capacity", {}).get("capacity_verified") is not True
    ):
        raise phase_b.PhaseBError("Phase-B frozen export inventory is inconsistent")
    return seal, plan, bundle, task_results, summary


def validate_outer_metrics(
    path: Path,
    run_root: Path,
    solver_root: Path,
    plan: dict[str, Any],
    run_summary: dict[str, Any],
) -> dict[str, Any]:
    record, data = phase_b.read_json(path, "outer driver metrics")
    if not isinstance(record, dict):
        raise phase_b.PhaseBError("outer driver metrics must be an object")
    command = record.get("command")
    if not isinstance(command, list) or not all(isinstance(item, str) for item in command):
        raise phase_b.PhaseBError("outer driver metrics do not bind a command")
    normalized_command = list(command)
    if normalized_command:
        normalized_command[0] = str(Path(normalized_command[0]).resolve())
    if len(normalized_command) > 1 and normalized_command[1].endswith(".py"):
        normalized_command[1] = str(Path(normalized_command[1]).resolve())
    if normalized_command != plan.get("outer_expected_command"):
        raise phase_b.PhaseBError("outer driver command differs from the frozen execution command")
    if command.count("--output") != 1 or command.count("--solver-root") != 1:
        raise phase_b.PhaseBError("outer driver command has ambiguous output or solver-root binding")
    try:
        output_index = command.index("--output")
        solver_index = command.index("--solver-root")
    except ValueError as error:
        raise phase_b.PhaseBError("outer driver command lacks output or solver-root binding") from error
    if output_index + 1 >= len(command) or Path(command[output_index + 1]).resolve() != run_root.resolve():
        raise phase_b.PhaseBError("outer driver command targets a different run root")
    if solver_index + 1 >= len(command):
        raise phase_b.PhaseBError("outer driver command lacks its solver-root value")
    if Path(command[solver_index + 1]).resolve() != solver_root.resolve():
        raise phase_b.PhaseBError("outer driver command uses a different solver root")
    if "run" not in command or "--max-instances" in command:
        raise phase_b.PhaseBError("outer driver metrics do not cover the full Phase-B run command")
    if (
        record.get("returncode") != 0
        or record.get("timed_out") is not False
        or record.get("orphan_group_terminated") is not False
    ):
        raise phase_b.PhaseBError("outer driver did not terminate cleanly")
    watchdog = record.get("watchdog_seconds")
    if isinstance(watchdog, bool) or not isinstance(watchdog, (int, float)) or not math.isfinite(watchdog) or watchdog <= 0:
        raise phase_b.PhaseBError("outer driver metrics use an invalid watchdog")
    metrics = record.get("metrics")
    required = {
        "wall_seconds",
        "user_seconds",
        "system_seconds",
        "total_core_seconds",
        "single_core_seconds",
        "peak_rss_bytes",
        "meter",
    }
    if not isinstance(metrics, dict) or set(metrics) != required:
        raise phase_b.PhaseBError("outer driver metrics lack required resources")
    for name in required - {"meter", "peak_rss_bytes"}:
        value = metrics[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise phase_b.PhaseBError("outer driver metrics contain invalid numeric resources")
    if (
        isinstance(metrics["peak_rss_bytes"], bool)
        or not isinstance(metrics["peak_rss_bytes"], int)
        or metrics["peak_rss_bytes"] < 0
        or metrics["meter"] != "fresh-process getrusage(RUSAGE_CHILDREN)"
        or not math.isclose(
            metrics["total_core_seconds"],
            metrics["user_seconds"] + metrics["system_seconds"],
            rel_tol=0,
            abs_tol=1e-9,
        )
        or metrics["single_core_seconds"] != metrics["total_core_seconds"]
    ):
        raise phase_b.PhaseBError("outer driver resource accounting is inconsistent")
    inner = run_summary.get("charged_process_resources")
    if not isinstance(inner, dict):
        raise phase_b.PhaseBError("run summary lacks charged process resources")
    if (
        metrics["wall_seconds"] + 1e-9 < inner.get("summed_process_wall_seconds", math.inf)
        or metrics["total_core_seconds"] + 1e-9 < inner.get("total_core_seconds", math.inf)
    ):
        raise phase_b.PhaseBError("outer driver resources do not enclose the charged child processes")
    return {
        "path": str(path.resolve()),
        "bytes": len(data),
        "sha256": phase_b.sha256_bytes(data),
        "command": command,
        "watchdog_seconds": watchdog,
        "metrics": metrics,
        "evidence_class": plan.get("evidence_class"),
    }


def load_authenticated_oracle(
    raw_seal_path: Path,
    oracle_path: Path,
    protocol: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    raw_seal, raw_seal_bytes = phase_b.read_json(raw_seal_path, "raw Phase-A seal")
    phase_b.validate_phase_a_seal(raw_seal, raw_seal_bytes, protocol)
    oracle, oracle_bytes = phase_b.read_json(oracle_path, "Phase-A oracle ledger")
    if not isinstance(oracle, dict) or oracle.get("schema") != ORACLE_SCHEMA:
        raise phase_b.PhaseBError("Phase-A oracle schema changed")
    expected_hash = raw_seal.get("oracle_ledger_sha256")
    phase_b.require_hex64(expected_hash, "raw Phase-A oracle commitment")
    if phase_b.sha256_bytes(oracle_bytes) != expected_hash:
        raise phase_b.PhaseBError("Phase-A oracle bytes do not match the raw seal")
    return raw_seal, oracle


def oracle_index(
    oracle: dict[str, Any], bundle: dict[str, Any], protocol: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    truth: dict[str, dict[str, Any]] = {}
    class_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for cell_record in oracle.get("cells", []):
        if not isinstance(cell_record, dict):
            raise phase_b.PhaseBError("Phase-A oracle cell must be an object")
        cell = cell_record.get("cell")
        entries = cell_record.get("entries")
        if not isinstance(cell, dict) or not isinstance(entries, list):
            raise phase_b.PhaseBError("Phase-A oracle cell lacks its cell or entry list")
        for entry in entries:
            if not isinstance(entry, dict):
                raise phase_b.PhaseBError("Phase-A oracle entry must be an object")
            blind_id = entry.get("blind_instance_id")
            if blind_id in truth:
                raise phase_b.PhaseBError(f"duplicate Phase-A oracle id {blind_id}")
            if entry.get("target_class") not in {"decomposable", "nondecomposable"}:
                raise phase_b.PhaseBError(f"invalid Phase-A class for {blind_id}")
            if entry.get("cell_id") != cell.get("id"):
                raise phase_b.PhaseBError(f"Phase-A oracle cell mismatch for {blind_id}")
            truth[blind_id] = entry
            class_counts[entry["cell_id"]][entry["target_class"]] += 1
    blind = {item["blind_instance_id"]: item for item in bundle["instances"]}
    if set(truth) != set(blind):
        raise phase_b.PhaseBError("Phase-A oracle and blind bundle do not have an exact id bijection")
    for blind_id, instance in blind.items():
        entry = truth[blind_id]
        if entry.get("cell_id") != instance["cell_id"] or entry.get("target") != instance["target"]:
            raise phase_b.PhaseBError(f"Phase-A oracle target mapping changed for {blind_id}")
    for cell in protocol["cells"]:
        expected = cell["expected_instances"]
        if expected % 2 != 0 or class_counts[cell["id"]] != Counter(
            {"decomposable": expected // 2, "nondecomposable": expected // 2}
        ):
            raise phase_b.PhaseBError(f"Phase-A oracle class quota changed for {cell['id']}")
    return truth


def classification(expected: str, status: str) -> str:
    if status == "sat":
        return "true_positive" if expected == "decomposable" else "false_positive"
    if status == "unsat":
        return "true_negative" if expected == "nondecomposable" else "false_negative"
    return "inconclusive"


def metric_values(row: dict[str, Any]) -> list[dict[str, Any]]:
    values = []
    process = row.get("process")
    if isinstance(process, dict) and isinstance(process.get("metrics"), dict):
        values.append(process["metrics"])
    validation = row.get("point_witness_validation")
    if isinstance(validation, dict):
        validation_process = validation.get("process")
        if isinstance(validation_process, dict) and isinstance(validation_process.get("metrics"), dict):
            values.append(validation_process["metrics"])
    return values


def score_run(
    protocol_path: Path,
    solver_root: Path,
    run_root: Path,
    raw_seal_path: Path,
    oracle_path: Path,
    output: Path,
    *,
    allow_smoke: bool,
    outer_metrics_path: Path | None,
) -> dict[str, Any]:
    if output.exists() or output.is_symlink():
        raise phase_b.PhaseBError(f"score output must be new: {output}")
    protocol, protocol_bytes = phase_b.read_json(protocol_path, "Phase-B protocol")
    phase_b.validate_protocol(protocol)
    run_seal, plan, bundle, tasks, run_summary = verify_run_tree(
        run_root, solver_root, protocol, protocol_bytes, allow_smoke=allow_smoke
    )
    outer_metrics = (
        validate_outer_metrics(outer_metrics_path, run_root, solver_root, plan, run_summary)
        if outer_metrics_path is not None
        else None
    )
    if (
        plan.get("evidence_class")
        == phase_b.PRODUCTION_EVIDENCE_CLASS
        and outer_metrics is None
    ):
        raise phase_b.PhaseBError("full scientific scoring requires the bound outer driver receipt")
    if plan.get("evidence_class") == "operational_smoke" and not allow_smoke:
        raise phase_b.PhaseBError("operational-smoke scoring requires --allow-smoke")
    raw_seal, oracle = load_authenticated_oracle(raw_seal_path, oracle_path, protocol)
    truth = oracle_index(oracle, bundle, protocol)
    bundle_by_id = {item["blind_instance_id"]: item for item in bundle["instances"]}
    selected_ids = plan["selected_blind_instance_ids"]
    task_by_id = {task["blind_instance_id"]: task for task in tasks}
    if set(task_by_id) != set(selected_ids):
        raise phase_b.PhaseBError("scoring task ids differ from the frozen selected ids")
    rows = []
    grouped: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    resources: dict[str, dict[str, Any]] = {
        backend: {
            "total_core_seconds": 0.0,
            "single_core_seconds": 0.0,
            "single_core_seconds_alias_of": "total_core_seconds",
            "single_core_elapsed_seconds": None,
            "summed_process_wall_seconds": 0.0,
            "peak_rss_bytes": 0,
            "peak_rss_scope": "maximum fresh-process high-water mark, not aggregate parallel memory",
            "conflicts_sum": 0,
            "conflicts_reported": 0,
            "conflict_values": [],
        }
        for backend in phase_b.BACKENDS
    }
    for blind_id in selected_ids:
        instance = bundle_by_id[blind_id]
        task = task_by_id[blind_id]
        if (
            task.get("cell_id") != instance["cell_id"]
            or task.get("target") != instance["target"]
            or task.get("source_system_id") != instance["source_system_id"]
        ):
            raise phase_b.PhaseBError(f"solver task no longer matches blind instance {blind_id}")
        expected = truth[blind_id]["target_class"]
        for backend_row in task["backends"]:
            backend = backend_row["solver"]
            status = backend_row.get("status")
            outcome = classification(expected, status)
            grouped[(instance["cell_id"], backend, expected)][outcome] += 1
            grouped[(instance["cell_id"], backend, expected)][f"status:{status}"] += 1
            conflicts = backend_row.get("conflicts")
            if isinstance(conflicts, int) and conflicts >= 0:
                resources[backend]["conflicts_sum"] += conflicts
                resources[backend]["conflicts_reported"] += 1
                resources[backend]["conflict_values"].append(conflicts)
            for metrics in metric_values(backend_row):
                resources[backend]["total_core_seconds"] += metrics["total_core_seconds"]
                resources[backend]["single_core_seconds"] += metrics["single_core_seconds"]
                resources[backend]["summed_process_wall_seconds"] += metrics["wall_seconds"]
                resources[backend]["peak_rss_bytes"] = max(
                    resources[backend]["peak_rss_bytes"], metrics["peak_rss_bytes"]
                )
            rows.append(
                {
                    "blind_instance_id": blind_id,
                    "source_system_id": instance["source_system_id"],
                    "cell_id": instance["cell_id"],
                    "target_class": expected,
                    "backend": backend,
                    "solver_status": status,
                    "classification": outcome,
                    "conflicts": conflicts,
                }
            )
    group_rows = []
    for (cell_id, backend, expected), counts in sorted(grouped.items()):
        group_rows.append(
            {
                "cell_id": cell_id,
                "backend": backend,
                "target_class": expected,
                "counts": dict(sorted(counts.items())),
            }
        )
    cluster_sizes = Counter(
        bundle_by_id[blind_id]["source_system_id"] for blind_id in selected_ids
    )
    for backend in phase_b.BACKENDS:
        values = sorted(resources[backend].pop("conflict_values"))
        resources[backend].update(
            {
                "conflict_values": values,
                "conflict_values_sha256": phase_b.canonical_sha256(values),
                "conflicts_min": min(values) if values else None,
                "conflicts_median": statistics.median(values) if values else None,
                "conflicts_max": max(values) if values else None,
            }
        )
    score = {
        "schema": SCORE_SCHEMA,
        "created_at": phase_b.now(),
        "claim_boundary": "Balanced public toy PDP scoring only; class balance is fixed by construction and does not estimate natural prevalence or establish index-calculus SOTA",
        "protocol_sha256": phase_b.sha256_bytes(protocol_bytes),
        "phase_a_seal_sha256": phase_b.sha256_file(raw_seal_path, "raw Phase-A seal"),
        "phase_a_oracle_sha256": raw_seal["oracle_ledger_sha256"],
        "phase_b_run_inventory_sha256": run_seal["inventory_sha256"],
        "outer_driver_accounting": outer_metrics,
        "selected_instances": len(selected_ids),
        "backend_rows": len(rows),
        "full_panel_scored": run_seal["full_panel_complete"],
        "source_system_clusters": {
            "clusters": len(cluster_sizes),
            "largest_cluster": max(cluster_sizes.values(), default=0),
            "multi_target_clusters": sum(size > 1 for size in cluster_sizes.values()),
            "bootstrap_unit": "source_system_id",
        },
        "per_cell_backend_class": group_rows,
        "backend_resources": resources,
        "resource_field_semantics": {
            "single_core_seconds": "legacy alias of total_core_seconds (user plus system CPU), not measured single-core elapsed time",
            "single_core_elapsed_seconds": None,
            "peak_rss_bytes": "maximum fresh-process high-water mark within each backend, not aggregate parallel memory",
        },
        "native_xor_validation_accounting": "Native process resources include source regeneration, solve, source-model validation, and exact lifted point-witness validation",
        "rows": rows,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "full_cost_blockers": phase_b.tool_build_cost_blockers(plan.get("additional_tool_builds", {})),
        "separately_charged_tool_builds": {
            name: {"receipt_sha256": value["receipt_sha256"], **value["receipt"]["resources"]}
            for name, value in sorted(plan.get("additional_tool_builds", {}).items())
        },
    }
    output.mkdir(parents=False)
    phase_b.write_json_new(output / "score.json", score)
    score_bytes = phase_b.regular_file_bytes(output / "score.json", "score")
    score_seal = {
        "schema": SCORE_SEAL_SCHEMA,
        "status": "post_run_truth_scoring_complete",
        "score_path": "score.json",
        "score_bytes": len(score_bytes),
        "score_sha256": phase_b.sha256_bytes(score_bytes),
        "phase_b_run_inventory_sha256": run_seal["inventory_sha256"],
        "phase_a_seal_sha256": score["phase_a_seal_sha256"],
        "phase_a_oracle_sha256": score["phase_a_oracle_sha256"],
    }
    phase_b.write_json_new(output / "score-seal.json", score_seal)
    return score_seal


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=phase_b.DEFAULT_PROTOCOL)
    parser.add_argument("--solver-root", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--phase-a-seal", type=Path, required=True)
    parser.add_argument("--oracle-ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--outer-metrics", type=Path)
    parser.add_argument("--allow-smoke", action="store_true")
    args = parser.parse_args()
    try:
        result = score_run(
            args.protocol.resolve(),
            args.solver_root.resolve(),
            args.run_root.resolve(),
            args.phase_a_seal.resolve(),
            args.oracle_ledger.resolve(),
            args.output.resolve(),
            allow_smoke=args.allow_smoke,
            outer_metrics_path=args.outer_metrics.resolve() if args.outer_metrics else None,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, phase_b.PhaseBError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
