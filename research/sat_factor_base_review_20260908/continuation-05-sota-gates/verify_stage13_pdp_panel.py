#!/usr/bin/env python3
"""Verify and summarize the frozen stage-13 target-matched PDP panel."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics
import tempfile
from typing import Any


HERE = Path(__file__).resolve().parent
DEFAULT_PROTOCOL = HERE / "stage-13-pdp-panel-protocol.json"
EXPORT_KEYS = {
    "wdsat_anf": "instance.anf",
    "cryptominisat_xor_dimacs": "instance.xor.cnf",
    "magma_boolean_f4": "instance.magma",
}
BACKENDS = ("native-xor", "wdsat", "cryptominisat", "magma-f4", "direct-mitm")
ALIASES = {
    "native-sat": "native-xor",
    "wdsat": "wdsat",
    "cryptominisat": "cryptominisat",
    "magma-f4": "magma-f4",
    "direct-mitm": "direct-mitm",
}
ALLOWED_STATUSES = {
    "native-xor": {
        "sat",
        "unsat",
        "unknown_inconclusive",
        "model_cap_inconclusive",
        "timeout_inconclusive",
    },
    "wdsat": {
        "sat",
        "sat_nonlifting_model_inconclusive",
        "unsat",
        "timeout_inconclusive",
        "solver_error",
        "unavailable_operational",
    },
    "cryptominisat": {
        "sat",
        "sat_nonlifting_model_inconclusive",
        "unsat",
        "unknown_inconclusive",
        "timeout_inconclusive",
        "solver_error",
        "unavailable_operational",
    },
    "magma-f4": {
        "unsat",
        "sat_basis_certificate_unverified_model",
        "timeout_inconclusive",
        "solver_error",
        "unavailable_operational",
    },
    "direct-mitm": {"sat", "timeout_inconclusive"},
}
EXTERNAL_EXPORT = {
    "wdsat": "wdsat_anf",
    "cryptominisat": "cryptominisat_xor_dimacs",
    "magma-f4": "magma_boolean_f4",
}


class VerificationError(RuntimeError):
    """A panel artifact violates its frozen protocol or source binding."""


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise VerificationError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise VerificationError(f"expected JSON object in {path}")
    return value


def validate_protocol(protocol: dict) -> None:
    if protocol.get("schema") != "koblitz_target_matched_pdp_panel_protocol.v1":
        raise VerificationError("unexpected protocol schema")
    if protocol.get("status") != "frozen_before_execution":
        raise VerificationError("protocol is not frozen before execution")
    seeds = protocol.get("replicate_seeds")
    cells = protocol.get("cells")
    if not isinstance(seeds, list) or len(seeds) < 2 or len(set(seeds)) != len(seeds):
        raise VerificationError("replicate_seeds must contain at least two unique values")
    if not isinstance(cells, list) or not cells:
        raise VerificationError("cells must be a nonempty list")
    cell_ids = [cell.get("id") for cell in cells]
    if any(not isinstance(cell_id, str) or not cell_id for cell_id in cell_ids):
        raise VerificationError("every cell requires a nonempty id")
    if len(set(cell_ids)) != len(cell_ids):
        raise VerificationError("cell ids are not unique")
    expected = len(seeds) * len(cells)
    if protocol.get("acceptance", {}).get("expected_tasks") != expected:
        raise VerificationError("acceptance.expected_tasks disagrees with seeds times cells")
    for cell in cells:
        if cell.get("m") != 3:
            raise VerificationError(f"{cell['id']}: current matrix runner supports only m=3")
        if cell.get("basis") not in {"standard", "ggmp"}:
            raise VerificationError(f"{cell['id']}: unsupported basis")
        if cell.get("curve_a") not in {0, 1}:
            raise VerificationError(f"{cell['id']}: curve_a must be zero or one")
    required = protocol.get("solver_policy", {}).get("required_backend_bindings")
    if required != list(BACKENDS):
        raise VerificationError("required backend list is not the verifier's exact frozen list")


def tasks(protocol: dict) -> list[tuple[int, dict]]:
    return [(seed, cell) for seed in protocol["replicate_seeds"] for cell in protocol["cells"]]


def task_relpath(seed: int, cell: dict) -> Path:
    return Path("tasks") / f"seed-{seed}" / cell["id"]


def semantic_descriptor(manifest: dict, protocol: dict) -> dict:
    fields = protocol["source_identity"]["semantic_fields"]
    missing = [field for field in fields if field not in manifest]
    if missing:
        raise VerificationError(f"manifest lacks semantic fields: {missing}")
    return {field: manifest[field] for field in fields}


def finite_nonnegative(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and value >= 0
    )


def validated_metrics(record: dict, label: str) -> dict:
    metrics = record.get("metrics")
    if not isinstance(metrics, dict):
        raise VerificationError(f"{label}: missing process metrics")
    numeric = ("wall_seconds", "user_seconds", "system_seconds", "total_core_seconds", "single_core_seconds")
    for name in numeric:
        if not finite_nonnegative(metrics.get(name)):
            raise VerificationError(f"{label}: invalid metric {name}")
    peak = metrics.get("peak_rss_bytes")
    if not isinstance(peak, int) or isinstance(peak, bool) or peak < 0:
        raise VerificationError(f"{label}: invalid peak_rss_bytes")
    if metrics.get("meter") != "fresh-process getrusage(RUSAGE_CHILDREN)":
        raise VerificationError(f"{label}: unexpected process meter")
    total = float(metrics["total_core_seconds"])
    if not math.isclose(total, float(metrics["user_seconds"]) + float(metrics["system_seconds"]), rel_tol=1e-9, abs_tol=1e-9):
        raise VerificationError(f"{label}: total core-seconds does not equal user plus system")
    if not math.isclose(float(metrics["single_core_seconds"]), total, rel_tol=1e-9, abs_tol=1e-9):
        raise VerificationError(f"{label}: single-core and total core-seconds disagree")
    return {name: metrics[name] for name in (*numeric, "peak_rss_bytes", "meter")}


def backend_records(row: dict) -> dict[str, dict]:
    candidates: list[dict] = []
    for key in ("solvers", "isolated_backends", "backend_processes"):
        value = row.get(key, [])
        if isinstance(value, list):
            candidates.extend(item for item in value if isinstance(item, dict))
        elif isinstance(value, dict):
            candidates.extend(
                dict(item, solver=name) if isinstance(item, dict) else {"solver": name}
                for name, item in value.items()
            )
    records: dict[str, dict] = {}
    for item in candidates:
        raw_name = str(item.get("solver", item.get("backend", ""))).lower()
        canonical = ALIASES.get(raw_name)
        if canonical is None:
            raise VerificationError(f"unknown backend record {raw_name!r}")
        if canonical in records:
            raise VerificationError(f"duplicate backend record {canonical}")
        records[canonical] = item
    missing = set(BACKENDS) - set(records)
    if missing:
        raise VerificationError(f"missing backend records: {sorted(missing)}")
    return records


def command_option(command: list, option: str) -> str | None:
    positions = [index for index, value in enumerate(command) if value == option]
    if len(positions) != 1 or positions[0] + 1 >= len(command):
        return None
    return str(command[positions[0] + 1])


def cell_config(cell: dict) -> str:
    return ":".join(
        str(cell[field]) for field in ("n", "ell", "basis", "curve_a", "factor_index")
    )


def archive_relative_path(artifact: Path, panel: Path, label: str) -> Path:
    try:
        panel_root = panel.resolve(strict=True)
        resolved = artifact.resolve(strict=True)
    except OSError as error:
        raise VerificationError(f"{label}: missing archive path: {error}") from error
    if not panel_root.is_dir():
        raise VerificationError("panel archive root is not a directory")
    try:
        return resolved.relative_to(panel_root)
    except ValueError as error:
        raise VerificationError(f"{label}: archive path escapes the panel root") from error


def recorded_archive_path(
    artifact: Path, panel: Path, recorded_panel: Path, label: str
) -> Path:
    return recorded_panel / archive_relative_path(artifact, panel, label)


def recorded_panel_root(
    protocol: dict,
    panel: Path,
    selected_task: tuple[int, dict] | None = None,
) -> Path:
    """Recover the execution-time panel root without requiring it to still exist.

    Archived command paths are custody data.  Verification reads artifacts from
    ``panel`` but checks every path-bearing command against the common absolute
    root recorded by the first frozen task.
    """
    candidates = [selected_task] if selected_task is not None else tasks(protocol)
    located = next(
        (
            (seed, cell, path)
            for seed, cell in candidates
            if (path := panel / task_relpath(seed, cell) / "invocation.json").is_file()
        ),
        None,
    )
    if located is None:
        raise VerificationError("panel has no frozen task invocation from which to recover its root")
    seed, cell, invocation_path = located
    archive_relative_path(invocation_path, panel, "recorded-root invocation")
    invocation = read_json(invocation_path)
    command = invocation.get("runner_command")
    if not isinstance(command, list):
        raise VerificationError("first frozen task lacks a runner command")
    output = command_option(command, "--output")
    if output is None:
        raise VerificationError("first frozen task lacks a unique --output path")
    root = Path(output)
    if not root.is_absolute() or ".." in root.parts:
        raise VerificationError("recorded panel output path is not canonical and absolute")
    suffix = task_relpath(seed, cell) / "matrix"
    for component in reversed(suffix.parts):
        if root.name != component:
            raise VerificationError("recorded panel output path has an unexpected suffix")
        root = root.parent
    if output != str(root / suffix):
        raise VerificationError("recorded panel output path has a noncanonical spelling")
    return root


def repository_root(path: Path) -> Path:
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise VerificationError("cannot locate the repository root for the verifier")


def validate_panel_run(protocol: dict, panel: Path) -> tuple[dict, str]:
    run_path = panel / "panel-run.json"
    archive_relative_path(run_path, panel, "panel execution plan")
    run = read_json(run_path)
    if run.get("schema") != "koblitz_target_matched_pdp_execution_plan.v1":
        raise VerificationError("panel-run.json has the wrong schema")
    if run.get("protocol_sha256") != canonical_sha256(protocol):
        raise VerificationError("panel-run.json is bound to a different protocol")
    if not isinstance(run.get("repo"), str) or not Path(run["repo"]).is_absolute():
        raise VerificationError("panel-run.json lacks an absolute execution repository path")
    expected_tasks = [
        {
            "seed": seed,
            "cell": cell,
            "config": cell_config(cell),
            "relative_output": str(task_relpath(seed, cell) / "matrix"),
        }
        for seed, cell in tasks(protocol)
    ]
    if run.get("frozen_tasks") != expected_tasks:
        raise VerificationError("panel-run.json does not contain the exact frozen task plan")
    if run.get("full_frozen_task_count") != len(expected_tasks):
        raise VerificationError("panel-run.json frozen task count is wrong")
    evidence_class = run.get("evidence_class")
    if evidence_class not in {"scientific_candidate", "operational_smoke"}:
        raise VerificationError("panel-run.json has an invalid evidence class")
    source_revision = run.get("source_revision")
    if not isinstance(source_revision, dict) or not isinstance(source_revision.get("commit"), str):
        raise VerificationError("panel-run.json lacks a source revision")
    if evidence_class == "scientific_candidate" and source_revision.get("dirty") is not False:
        raise VerificationError("scientific-candidate panel started from a dirty checkout")
    implementation = run.get("implementation")
    required_implementation = {"python", "driver", "verifier", "exporter", "backend", "runner", "meter"}
    if not isinstance(implementation, dict) or set(implementation) != required_implementation:
        raise VerificationError("panel-run.json implementation identity is incomplete")
    for name, identity in implementation.items():
        if not isinstance(identity, dict) or not isinstance(identity.get("path"), str):
            raise VerificationError(f"panel-run.json has malformed {name} identity")
        digest = identity.get("sha256")
        if not isinstance(digest, str) or len(digest) != 64:
            raise VerificationError(f"panel-run.json has malformed {name} SHA-256")
    requested_tools = run.get("requested_tools")
    required_tools = {"wdsat_binary", "wdsat_source", "cryptominisat", "magma"}
    if not isinstance(requested_tools, dict) or set(requested_tools) != required_tools:
        raise VerificationError("panel-run.json requested-tool identity is incomplete")
    if evidence_class == "scientific_candidate":
        for name, identity in requested_tools.items():
            git = identity.get("git") if isinstance(identity, dict) else None
            if isinstance(git, dict) and git.get("dirty") is not False:
                raise VerificationError(f"scientific-candidate panel used dirty {name} source")
    if not isinstance(run.get("invocation_history"), list) or not run["invocation_history"]:
        raise VerificationError("panel-run.json lacks invocation history")
    frozen_identity = {
        key: run[key]
        for key in (
            "schema",
            "protocol_sha256",
            "repo",
            "source_revision",
            "implementation",
            "requested_tools",
            "evidence_class",
            "full_frozen_task_count",
            "frozen_tasks",
        )
    }
    return run, canonical_sha256(frozen_identity)


def validate_relocated_custody(panel: Path, recorded_panel: Path, run: dict) -> bytes | None:
    if panel.resolve(strict=True) == recorded_panel:
        return None
    custody_path = panel.parent / "stage-13-custody.json"
    if not custody_path.is_file() or custody_path.is_symlink():
        raise VerificationError("relocated panel requires its adjacent custody record")
    custody = read_json(custody_path)
    source_revision = run["source_revision"]["commit"]
    if (
        custody.get("schema") != "koblitz_stage13_panel_custody.v1"
        or custody.get("panel_archive") != panel.name
        or custody.get("panel_original_root") != str(recorded_panel)
        or custody.get("panel_source_revision") != source_revision
    ):
        raise VerificationError("relocated panel disagrees with its custody identity")
    summary_path = panel / "panel-summary.json"
    archive_relative_path(summary_path, panel, "custody-bound panel summary")
    summary_bytes = summary_path.read_bytes()
    if hashlib.sha256(summary_bytes).hexdigest() != custody.get("panel_summary_sha256"):
        raise VerificationError("relocated panel summary disagrees with its custody hash")
    return summary_bytes


def validate_task_runtime(
    protocol: dict,
    panel: Path,
    recorded_panel: Path,
    run: dict,
    seed: int,
    cell: dict,
    report: dict,
    row: dict,
) -> tuple[dict, dict]:
    task_dir = panel / task_relpath(seed, cell)
    matrix_dir = task_dir / "matrix"
    invocation_path = task_dir / "invocation.json"
    recorded_matrix_dir = recorded_archive_path(
        matrix_dir, panel, recorded_panel, f"{task_dir} matrix directory"
    )
    invocation = read_json(
        panel / archive_relative_path(invocation_path, panel, f"{task_dir} invocation")
    )
    if invocation.get("seed") != seed or invocation.get("cell") != cell:
        raise VerificationError(f"{task_dir}: invocation does not match frozen task")
    key = f"seed-{seed}/{cell['id']}"
    if not any(
        isinstance(item, dict) and key in item.get("selection", [])
        for item in run["invocation_history"]
    ):
        raise VerificationError(f"{task_dir}: task is absent from invocation history")
    command = invocation.get("runner_command")
    if not isinstance(command, list) or len(command) < 2:
        raise VerificationError(f"{task_dir}: malformed runner command")
    implementation = run["implementation"]
    expected_options = {
        "--output": str(recorded_matrix_dir),
        "--exporter": implementation["exporter"]["path"],
        "--backend": implementation["backend"]["path"],
        "--timeout": str(protocol["solver_policy"]["per_process_watchdog_seconds"]),
        "--conflicts": str(protocol["solver_policy"]["native_conflict_budget"]),
        "--configs": cell_config(cell),
        "--seed": str(seed),
    }
    if command[0] != implementation["python"]["path"] or command[1] != implementation["runner"]["path"]:
        raise VerificationError(f"{task_dir}: runner executable identity changed")
    for option, expected in expected_options.items():
        if command_option(command, option) != expected:
            raise VerificationError(f"{task_dir}: runner option {option} differs from protocol")
    tool_options = {
        "wdsat_binary": "--wdsat",
        "wdsat_source": "--wdsat-source",
        "cryptominisat": "--cryptominisat",
        "magma": "--magma",
    }
    for name, option in tool_options.items():
        identity = run["requested_tools"][name]
        actual = command_option(command, option)
        expected = identity.get("path") if identity.get("requested") else None
        if actual != expected:
            raise VerificationError(f"{task_dir}: requested tool option {option} changed")
    expected_meter_command = [
        implementation["python"]["path"],
        implementation["meter"]["path"],
        "--cwd",
        run["repo"],
        "--timeout",
        str(protocol["solver_policy"]["whole_cell_watchdog_seconds"]),
        "--stdout",
        str(
            recorded_archive_path(
                task_dir / "runner.stdout", panel, recorded_panel, f"{task_dir} stdout"
            )
        ),
        "--stderr",
        str(
            recorded_archive_path(
                task_dir / "runner.stderr", panel, recorded_panel, f"{task_dir} stderr"
            )
        ),
        "--metrics",
        str(
            recorded_archive_path(
                task_dir / "outer-metrics.json",
                panel,
                recorded_panel,
                f"{task_dir} outer metrics",
            )
        ),
        "--",
        *command,
    ]
    if invocation.get("meter_command") != expected_meter_command:
        raise VerificationError(f"{task_dir}: outer meter invocation differs from frozen command")
    outer = read_json(task_dir / "outer-metrics.json")
    if outer.get("command") != command:
        raise VerificationError(f"{task_dir}: outer meter command differs from invocation")
    watchdog = protocol["solver_policy"]["whole_cell_watchdog_seconds"]
    if outer.get("watchdog_seconds") != watchdog or invocation.get("outer_watchdog_seconds") != watchdog:
        raise VerificationError(f"{task_dir}: outer watchdog differs from protocol")
    if invocation.get("meter_launcher_returncode") != 0:
        raise VerificationError(f"{task_dir}: process-meter launcher failed")
    if outer.get("returncode") != 0 or outer.get("timed_out") is not False:
        raise VerificationError(f"{task_dir}: outer matrix process did not complete normally")
    if report.get("schema") != "koblitz_pdp_matched_matrix.v1":
        raise VerificationError(f"{task_dir}: matrix report schema differs")
    expected_policy = {
        "unknown_is_unsat": False,
        "solver_subprocess_resources_charged": True,
        "python_orchestration_included_in_child_sum": False,
        "whole_panel_requires_outer_process_meter": True,
        "single_thread_requested": True,
        "missing_tool_is_negative_evidence": False,
        "producer_runs_solver_backends": False,
        "native_sat_and_mitm_isolated": True,
    }
    if report.get("policy") != expected_policy:
        raise VerificationError(f"{task_dir}: matrix policy differs from frozen runtime semantics")
    generator = row.get("generator")
    if not isinstance(generator, dict) or generator.get("status") != "completed" or generator.get("returncode") != 0:
        raise VerificationError(f"{task_dir}: source exporter did not complete")
    generator_command = generator.get("command")
    if (
        not isinstance(generator_command, list)
        or generator_command[0] != implementation["exporter"]["path"]
        or generator_command[-1] != "--export-only"
        or str(seed) not in generator_command
        or str(protocol["solver_policy"]["native_conflict_budget"]) not in generator_command
    ):
        raise VerificationError(f"{task_dir}: source exporter command differs from frozen semantics")
    tools = report.get("tools")
    if not isinstance(tools, dict):
        raise VerificationError(f"{task_dir}: matrix report lacks tool identities")
    isolated = tools.get("isolated_backend", {})
    if (
        isolated.get("available") is not True
        or isolated.get("path") != implementation["backend"]["path"]
        or isolated.get("sha256") != implementation["backend"]["sha256"]
    ):
        raise VerificationError(f"{task_dir}: isolated backend tool identity differs")
    for requested_name, report_name in (
        ("cryptominisat", "cryptominisat"),
        ("magma", "magma"),
    ):
        expected = run["requested_tools"][requested_name]
        observed = tools.get(report_name, {})
        if observed.get("available") is not expected.get("available"):
            raise VerificationError(f"{task_dir}: {report_name} availability differs from frozen plan")
        if expected.get("available") and (
            observed.get("path") != expected.get("path")
            or observed.get("sha256") != expected.get("sha256")
        ):
            raise VerificationError(f"{task_dir}: {report_name} binary identity differs")
    expected_source = run["requested_tools"]["wdsat_source"]
    expected_binary = run["requested_tools"]["wdsat_binary"]
    observed_wdsat = tools.get("wdsat", {})
    if expected_source.get("requested"):
        if (
            observed_wdsat.get("available") is not expected_source.get("available")
            or observed_wdsat.get("path") != expected_source.get("path")
            or observed_wdsat.get("source_commit") != (expected_source.get("git") or {}).get("commit")
        ):
            raise VerificationError(f"{task_dir}: WDSat source identity differs")
    elif expected_binary.get("requested"):
        if (
            observed_wdsat.get("available") is not expected_binary.get("available")
            or observed_wdsat.get("path") != expected_binary.get("path")
            or observed_wdsat.get("sha256") != expected_binary.get("sha256")
        ):
            raise VerificationError(f"{task_dir}: WDSat binary identity differs")
    elif observed_wdsat.get("available") is not False:
        raise VerificationError(f"{task_dir}: unrequested WDSat unexpectedly became available")
    return invocation, outer


def validate_external_point_witness(
    record: dict,
    backend: str,
    manifest_path: Path,
    recorded_manifest_path: Path,
    panel: Path,
    recorded_panel: Path,
    producer_id: str,
    backend_path: str,
    expected_valid: bool,
    label: str,
) -> dict:
    validation = record.get("point_witness_validation")
    if not isinstance(validation, dict):
        raise VerificationError(f"{label}: missing separately metered point-witness validation")
    expected_status = "valid_point_witness" if expected_valid else "nonlifting_source_model"
    expected_returncode = 0 if expected_valid else 2
    if (
        validation.get("status") != expected_status
        or validation.get("returncode") != expected_returncode
        or validation.get("timed_out") is not False
    ):
        raise VerificationError(f"{label}: contradictory point-witness validation status")
    assignment_name = validation.get("assignment_path")
    if (
        not isinstance(assignment_name, str)
        or Path(assignment_name).is_absolute()
        or len(Path(assignment_name).parts) != 1
        or assignment_name != f"{backend}.source-model.json"
    ):
        raise VerificationError(f"{label}: unsafe or unexpected source-model assignment path")
    assignment_path = manifest_path.parent / assignment_name
    if not assignment_path.is_file() or assignment_path.is_symlink():
        raise VerificationError(f"{label}: source-model assignment artifact is missing")
    recorded_assignment_path = recorded_archive_path(
        assignment_path, panel, recorded_panel, f"{label} source-model assignment"
    )
    assignment_bytes = assignment_path.read_bytes()
    if hashlib.sha256(assignment_bytes).hexdigest() != validation.get("assignment_sha256"):
        raise VerificationError(f"{label}: source-model assignment SHA-256 changed")
    try:
        assignment = json.loads(assignment_bytes)
    except json.JSONDecodeError as error:
        raise VerificationError(f"{label}: source-model assignment is not JSON") from error
    report = validation.get("report")
    if (
        not isinstance(report, dict)
        or report.get("schema") != "koblitz_pdp_assignment_validation.v1"
        or report.get("source_instance_id") != producer_id
        or report.get("source_instance_verified") is not True
        or report.get("regenerated_source_exact") is not True
        or report.get("source_assignment") != assignment
        or report.get("assignment_values") != len(assignment)
        or report.get("source_model_valid") is not True
        or report.get("source_witness_valid") is not expected_valid
        or report.get("status") != expected_status
        or report.get("assignment_blake3") != validation.get("assignment_blake3")
    ):
        raise VerificationError(f"{label}: invalid point-witness validation contract")
    expected_command = [
        backend_path,
        "validate-model",
        str(recorded_manifest_path),
        str(recorded_assignment_path),
    ]
    if validation.get("command") != expected_command:
        raise VerificationError(f"{label}: point-witness validation command changed")
    return validated_metrics(validation, f"{label} point-witness validation")


def verify_cell(
    protocol: dict,
    panel: Path,
    seed: int,
    cell: dict,
    write_receipt: bool,
    recorded_panel: Path | None = None,
) -> dict:
    if recorded_panel is None:
        recorded_panel = recorded_panel_root(protocol, panel, (seed, cell))
    run, execution_identity_sha256 = validate_panel_run(protocol, panel)
    validate_relocated_custody(panel, recorded_panel, run)
    task_dir = panel / task_relpath(seed, cell)
    matrix_dir = task_dir / "matrix"
    result_path = matrix_dir / "result.json"
    outer_path = task_dir / "outer-metrics.json"
    if not result_path.is_file():
        raise VerificationError(f"{task_dir}: missing matrix/result.json")
    if not outer_path.is_file():
        raise VerificationError(f"{task_dir}: missing outer-metrics.json")
    archive_relative_path(result_path, panel, f"{task_dir} matrix result")
    archive_relative_path(outer_path, panel, f"{task_dir} outer metrics")
    report = read_json(result_path)
    rows = report.get("instances")
    if not isinstance(rows, list) or len(rows) != 1:
        raise VerificationError(f"{result_path}: expected exactly one matrix instance")
    row = rows[0]
    _, outer = validate_task_runtime(
        protocol, panel, recorded_panel, run, seed, cell, report, row
    )
    actual_cell = row.get("cell", {})
    for field in ("n", "ell", "m", "basis", "curve_a", "factor_index"):
        if actual_cell.get(field) != cell.get(field):
            raise VerificationError(
                f"{cell['id']} seed {seed}: row {field}={actual_cell.get(field)!r}, "
                f"expected {cell.get(field)!r}"
            )
    manifests = list(matrix_dir.glob("*/manifest.json"))
    if len(manifests) != 1:
        raise VerificationError(f"{matrix_dir}: expected one instance manifest, found {len(manifests)}")
    manifest_path = manifests[0]
    recorded_manifest_path = recorded_archive_path(
        manifest_path, panel, recorded_panel, f"{task_dir} manifest"
    )
    manifest = read_json(manifest_path)
    expected_generator_command = [
        run["implementation"]["exporter"]["path"],
        str(cell["n"]),
        str(cell["ell"]),
        cell["basis"],
        str(seed),
        str(protocol["solver_policy"]["native_conflict_budget"]),
        str(recorded_manifest_path.parent),
        str(cell["curve_a"]),
        str(cell["factor_index"]),
        "--export-only",
    ]
    if row.get("generator", {}).get("command") != expected_generator_command:
        raise VerificationError(f"{manifest_path}: exporter command is not the frozen cell command")
    for field in ("n", "ell", "m", "curve_a"):
        if manifest.get(field) != cell.get(field):
            raise VerificationError(f"{manifest_path}: {field} disagrees with protocol")
    if manifest.get("seed") != seed:
        raise VerificationError(
            f"{manifest_path}: seed {manifest.get('seed')!r} does not equal frozen seed {seed}"
        )
    if (
        manifest.get("native_sat", {}).get("status") != "not_run_in_export_process"
        or manifest.get("native_sat", {}).get("conflict_budget")
        != protocol["solver_policy"]["native_conflict_budget"]
        or manifest.get("direct_meet_in_the_middle", {}).get("status")
        != "not_run_in_export_process"
    ):
        raise VerificationError(f"{manifest_path}: producer executed or mis-budgeted a solver backend")
    producer_metrics = validated_metrics(
        row["generator"], f"{cell['id']} seed {seed} source producer"
    )
    timing = manifest.get("timing_ns")
    required_timers = (
        "factor_base_predicate_construction",
        "planted_target_construction",
        "source_system_construction",
        "whole_process_internal",
    )
    if not isinstance(timing, dict) or any(
        not isinstance(timing.get(name), int) or timing[name] < 0 for name in required_timers
    ):
        raise VerificationError(f"{manifest_path}: missing or invalid producer stage timers")
    if timing.get("native_encoding_and_solve") is not None:
        raise VerificationError(f"{manifest_path}: export-only producer contains native solver timing")
    producer_timing_ns = {name: timing[name] for name in required_timers}
    predicate_kind = manifest.get("factor_base_predicate", {}).get("kind")
    expected_kind = "polynomial_subspace" if cell["basis"] == "standard" else "ggmp_linearised_kernel"
    if predicate_kind != expected_kind:
        raise VerificationError(f"{manifest_path}: unexpected predicate kind {predicate_kind!r}")
    if cell["basis"] == "ggmp":
        if manifest["factor_base_predicate"].get("factor_index") != cell["factor_index"]:
            raise VerificationError(f"{manifest_path}: GGMP factor index disagrees with protocol")

    semantic = semantic_descriptor(manifest, protocol)
    semantic_hash = canonical_sha256(semantic)
    exports: dict[str, dict] = {}
    recorded_export_paths: dict[str, Path] = {}
    for export_key, expected_name in EXPORT_KEYS.items():
        declared = manifest.get("exports", {}).get(export_key)
        if not isinstance(declared, dict):
            raise VerificationError(f"{manifest_path}: missing export declaration {export_key}")
        if declared.get("path") != expected_name:
            raise VerificationError(f"{manifest_path}: unexpected path for {export_key}")
        artifact = manifest_path.parent / expected_name
        if not artifact.is_file():
            raise VerificationError(f"missing export artifact {artifact}")
        recorded_export_paths[export_key] = recorded_archive_path(
            artifact, panel, recorded_panel, f"{cell['id']} seed {seed} {export_key} export"
        )
        size = artifact.stat().st_size
        if declared.get("bytes") != size:
            raise VerificationError(f"{artifact}: size disagrees with manifest")
        exports[export_key] = {
            "path": expected_name,
            "bytes": size,
            "sha256": sha256_file(artifact),
            "manifest_blake3": declared.get("blake3"),
        }
    before = row.get("source_artifacts_before")
    after = row.get("source_artifacts_after")
    if row.get("source_artifacts_unchanged") is not True:
        raise VerificationError(f"{cell['id']} seed {seed}: runner reports changed source artifacts")
    if before != exports or after != exports:
        raise VerificationError(
            f"{cell['id']} seed {seed}: pre/post backend SHA-256 custody snapshots "
            "do not equal the currently recomputed exports"
        )
    source_id = canonical_sha256({"semantic_sha256": semantic_hash, "exports": exports})
    producer_source = manifest.get("source_instance")
    if not isinstance(producer_source, dict):
        raise VerificationError(f"{manifest_path}: missing producer source_instance")
    if producer_source.get("schema") != "koblitz_pdp_source_instance.v1":
        raise VerificationError(f"{manifest_path}: unexpected producer source-instance schema")
    producer_id = producer_source.get("id_blake3")
    if (
        not isinstance(producer_id, str)
        or len(producer_id) != 64
        or any(character not in "0123456789abcdef" for character in producer_id)
    ):
        raise VerificationError(f"{manifest_path}: malformed producer BLAKE3 identifier")
    expected_producer_identity = {
        "schema": "koblitz_pdp_source_identity.v1",
        **{
            field: manifest[field]
            for field in (
                "n",
                "ell",
                "m",
                "seed",
                "curve_a",
                "irreducible_low_terms",
                "factor_base_predicate",
                "factor_base_basis_bitmasks",
                "target",
                "representation",
                "source_variables",
                "source_equations",
                "exports",
            )
        },
    }
    if producer_source.get("identity") != expected_producer_identity:
        raise VerificationError(f"{manifest_path}: producer source identity disagrees with manifest")

    records = backend_records(row)
    bindings = []
    for backend in BACKENDS:
        record = records[backend]
        status = record.get("status", "missing")
        if status not in ALLOWED_STATUSES[backend]:
            raise VerificationError(
                f"{cell['id']} seed {seed}: invalid or fail-open {backend} status {status!r}"
            )
        model_valid = record.get("source_model_valid")
        witness_valid = record.get("source_witness_valid")
        if status == "sat" and backend == "direct-mitm" and witness_valid is not True:
            raise VerificationError(
                f"{cell['id']} seed {seed}: direct MITM SAT lacks a valid point witness"
            )
        if status == "sat" and backend != "direct-mitm" and model_valid is not True:
            raise VerificationError(
                f"{cell['id']} seed {seed}: {backend} SAT lacks a valid source-model certificate"
            )
        if status == "sat" and backend in {"native-xor", "wdsat", "cryptominisat"} and witness_valid is not True:
            raise VerificationError(
                f"{cell['id']} seed {seed}: {backend} SAT lacks a valid rational point witness"
            )
        if status == "sat_nonlifting_model_inconclusive" and (
            model_valid is not True or witness_valid is not False
        ):
            raise VerificationError(
                f"{cell['id']} seed {seed}: {backend} nonlifting status lacks its exact certificates"
            )
        timed_out = record.get("timed_out")
        if status == "timeout_inconclusive":
            if timed_out is not True:
                raise VerificationError(f"{cell['id']} seed {seed}: timeout status lacks watchdog receipt")
        elif status != "unavailable_operational" and timed_out is not False:
            raise VerificationError(f"{cell['id']} seed {seed}: non-timeout {backend} has inconsistent watchdog state")
        if backend in {"native-xor", "direct-mitm"} and status != "timeout_inconclusive":
            backend_report = record.get("backend_report")
            if not isinstance(backend_report, dict):
                raise VerificationError(f"{cell['id']} seed {seed}: missing isolated backend report")
            authenticated = backend_report.get("source_artifacts")
            authenticated_exact = (
                isinstance(authenticated, dict)
                and set(authenticated) == set(EXPORT_KEYS)
                and all(
                    isinstance(authenticated[key], dict)
                    and authenticated[key].get("valid") is True
                    and authenticated[key].get("path") == exports[key]["path"]
                    and authenticated[key].get("bytes") == exports[key]["bytes"]
                    and authenticated[key].get("blake3") == exports[key]["manifest_blake3"]
                    for key in EXPORT_KEYS
                )
            )
            isolated_contract = (
                backend_report.get("schema") == "koblitz_pdp_isolated_backend.v1"
                and backend_report.get("backend")
                == ("native-sat" if backend == "native-xor" else "direct-mitm")
                and backend_report.get("source_instance_id") == producer_id
                and record.get("source_instance_id") == producer_id
                and backend_report.get("source_instance_verified") is True
                and backend_report.get("regenerated_source_exact") is True
                and backend_report.get("status") == status
                and authenticated_exact
            )
            if not isolated_contract:
                raise VerificationError(
                    f"{cell['id']} seed {seed}: {backend} failed independent source authentication"
                )
            if record.get("returncode") != 0:
                raise VerificationError(f"{cell['id']} seed {seed}: {backend} result has nonzero exit status")
            if backend == "native-xor":
                if backend_report.get("conflict_budget") != protocol["solver_policy"]["native_conflict_budget"]:
                    raise VerificationError(f"{cell['id']} seed {seed}: native conflict budget differs")
                if status == "model_cap_inconclusive" and not (
                    backend_report.get("max_models") == 64
                    and backend_report.get("models_examined") == 64
                    and backend_report.get("nonlifting_models_blocked") == 64
                    and backend_report.get("source_model_valid") is True
                    and backend_report.get("source_witness_valid") is False
                ):
                    raise VerificationError(
                        f"{cell['id']} seed {seed}: native model-cap receipt is contradictory"
                    )
            elif status == "unsat" and backend_report.get("exhaustive") is not True:
                raise VerificationError(f"{cell['id']} seed {seed}: MITM UNSAT is not exhaustive")
            elif status == "not_run_resource_cap" and backend_report.get("exhaustive") is not False:
                raise VerificationError(f"{cell['id']} seed {seed}: MITM cap status is contradictory")
        if backend == "magma-f4" and status in {"unsat", "sat_basis_certificate_unverified_model"}:
            terminal = record.get("magma_terminal")
            expected_terminal = "unsat" if status == "unsat" else "sat"
            if (
                not isinstance(terminal, dict)
                or terminal.get("schema") != "koblitz_magma_f4_terminal.v1"
                or terminal.get("algorithm") != "direct-f4-sparse"
                or terminal.get("terminal_status") != expected_terminal
                or terminal.get("single_thread_requested") is not True
                or terminal.get("gpu_disabled") is not True
                or record.get("returncode") != 0
                or record.get("timed_out") is not False
                or (status == "unsat" and terminal.get("basis_size") != 1)
                or (status != "unsat" and not isinstance(terminal.get("basis_size"), int))
            ):
                raise VerificationError(f"{cell['id']} seed {seed}: invalid Magma F4 terminal certificate")
        point_validation_metrics = None
        if backend in {"wdsat", "cryptominisat"} and status in {
            "sat",
            "sat_nonlifting_model_inconclusive",
        }:
            point_validation_metrics = validate_external_point_witness(
                record,
                backend,
                manifest_path,
                recorded_manifest_path,
                panel,
                recorded_panel,
                producer_id,
                run["implementation"]["backend"]["path"],
                status == "sat",
                f"{cell['id']} seed {seed} {backend}",
            )
        if status == "unavailable_operational":
            metrics = {
                "wall_seconds": None,
                "user_seconds": None,
                "system_seconds": None,
                "total_core_seconds": None,
                "single_core_seconds": None,
                "peak_rss_bytes": None,
                "meter": None,
            }
            resource_complete = False
        else:
            metrics = validated_metrics(record, f"{cell['id']} seed {seed} {backend}")
            resource_complete = True
        conflicts = record.get("conflicts")
        if conflicts is not None and (
            not isinstance(conflicts, int) or isinstance(conflicts, bool) or conflicts < 0
        ):
            raise VerificationError(f"{cell['id']} seed {seed}: invalid {backend} conflicts")
        if (
            backend in {"native-xor", "wdsat", "cryptominisat"}
            and status
            in {
                "sat",
                "sat_nonlifting_model_inconclusive",
                "unsat",
                "unknown_inconclusive",
            }
            and conflicts is None
        ):
            raise VerificationError(f"{cell['id']} seed {seed}: {backend} omitted exposed conflicts")
        accounting_scope = "isolated_process" if resource_complete else "unavailable_operational"
        export_key = EXTERNAL_EXPORT.get(backend)
        command = record.get("command") or []
        input_path_verified: bool | None = None
        input_sha256: str | None = None
        if export_key is not None:
            input_sha256 = exports[export_key]["sha256"]
            if status == "unavailable_operational":
                input_path_verified = None
            else:
                expected_path = str(recorded_export_paths[export_key])
                input_path_verified = expected_path in command
                if not input_path_verified:
                    raise VerificationError(
                        f"{cell['id']} seed {seed}: {backend} command is not bound to {expected_path}"
                    )
        elif backend in {"native-xor", "direct-mitm"}:
            expected_path = str(recorded_manifest_path)
            input_path_verified = expected_path in command
            if not input_path_verified:
                raise VerificationError(
                    f"{cell['id']} seed {seed}: {backend} command is not bound to {expected_path}"
                )
        bindings.append(
            {
                "backend": backend,
                "status": status,
                "source_instance_sha256": source_id,
                "producer_source_instance_blake3": producer_id,
                "input_export_sha256": input_sha256,
                "input_path_verified": input_path_verified,
                "source_model_valid": model_valid,
                "source_witness_valid": witness_valid,
                "conflicts": record.get("conflicts"),
                "metrics": metrics,
                "point_witness_validation_metrics": point_validation_metrics,
                "accounting_scope": accounting_scope,
                "resource_complete": resource_complete,
            }
        )
    if {binding["source_instance_sha256"] for binding in bindings} != {source_id}:
        raise VerificationError(f"{cell['id']} seed {seed}: backend source identifiers differ")
    source_solver_bindings = [
        binding for binding in bindings if binding["backend"] != "direct-mitm"
    ]
    source_sat = {
        "sat",
        "sat_nonlifting_model_inconclusive",
        "model_cap_inconclusive",
        "sat_basis_certificate_unverified_model",
    }
    has_sat = any(binding["status"] in source_sat for binding in source_solver_bindings)
    has_unsat = any(binding["status"] == "unsat" for binding in source_solver_bindings)
    if has_sat and has_unsat:
        raise VerificationError(
            f"{cell['id']} seed {seed}: source-equivalent solvers report contradictory SAT and UNSAT"
        )

    outer_metrics = validated_metrics(outer, f"{cell['id']} seed {seed} whole-cell")
    wdsat_build = row.get("wdsat_build")
    wdsat_fixed_cost = {"present": False, "complete": False}
    if isinstance(wdsat_build, dict):
        source_copy = wdsat_build.get("source_copy")
        clean = wdsat_build.get("clean")
        if not isinstance(source_copy, dict) or not isinstance(clean, dict):
            raise VerificationError(f"{cell['id']} seed {seed}: malformed WDSat build receipt")
        copy_metrics = validated_metrics(source_copy, f"{cell['id']} seed {seed} WDSat source copy")
        clean_metrics = validated_metrics(clean, f"{cell['id']} seed {seed} WDSat clean")
        build_metrics = validated_metrics(wdsat_build, f"{cell['id']} seed {seed} WDSat build")
        config_wall = wdsat_build.get("configuration_wall_seconds")
        if not finite_nonnegative(config_wall):
            raise VerificationError(f"{cell['id']} seed {seed}: invalid WDSat configuration time")
        wdsat_fixed_cost = {
            "present": True,
            "complete": (
                wdsat_build.get("status") == "completed"
                and source_copy.get("returncode") == 0
                and source_copy.get("timed_out") is False
                and clean.get("returncode") == 0
                and clean.get("timed_out") is False
            ),
            "status": wdsat_build.get("status"),
            "source_copy": copy_metrics,
            "configuration_wall_seconds": config_wall,
            "clean": clean_metrics,
            "build": build_metrics,
            "binary_sha256": wdsat_build.get("binary_sha256"),
        }
    receipt = {
        "schema": "koblitz_pdp_source_binding.v1",
        "task": {"seed": seed, "cell": cell},
        "semantic_descriptor": semantic,
        "semantic_sha256": semantic_hash,
        "exports": exports,
        "source_instance_sha256": source_id,
        "producer_source_instance_blake3": producer_id,
        "execution_identity_sha256": execution_identity_sha256,
        "evidence_class": run["evidence_class"],
        "backend_bindings": bindings,
        "source_producer_process": producer_metrics,
        "source_producer_timing_ns": producer_timing_ns,
        "whole_cell_process": outer_metrics,
        "wdsat_fixed_cost": wdsat_fixed_cost,
        "resource_complete": all(binding["resource_complete"] for binding in bindings),
    }
    receipt_path = task_dir / "source-binding.json"
    if receipt_path.exists():
        previous = read_json(receipt_path)
        if previous != receipt:
            raise VerificationError(f"source artifacts differ from custody receipt {receipt_path}")
    elif write_receipt:
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    else:
        raise VerificationError(
            f"missing custody receipt {receipt_path}; create it immediately after execution with --write-receipts"
        )
    return receipt


def metric_distribution(values: list[float | int]) -> dict | None:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    if not finite:
        return None
    ordered = sorted(finite)
    return {
        "count": len(ordered),
        "min": ordered[0],
        "median": statistics.median(ordered),
        "max": ordered[-1],
    }


def validate_fixed_ggmp_discovery(
    protocol: dict, recorded_repo: Path, live_repo: Path
) -> dict:
    frozen = protocol.get("fixed_cost_receipts", {}).get("ggmp_public_discovery", {})
    path = HERE / str(frozen.get("path", ""))
    if not path.is_file() or sha256_file(path) != frozen.get("sha256"):
        raise VerificationError("frozen GGMP public-discovery receipt is missing or changed")
    report = read_json(path)
    if (
        report.get("schema") != "koblitz_ggmp_public_factor_discovery.v1"
        or report.get("search_space_complete") is not True
        or report.get("uses_target_scalar") is not False
        or report.get("uses_discrete_log_labels") is not False
        or report.get("enumerates_target_subgroup") is not False
        or report.get("selected", {}).get("curve_a") != frozen.get("selected_curve_a")
        or report.get("selected", {}).get("factor_index") != frozen.get("selected_factor_index")
    ):
        raise VerificationError("frozen GGMP discovery receipt violates its public-selection contract")
    rows = report.get("rows")
    if not isinstance(rows, list) or not rows:
        raise VerificationError("frozen GGMP discovery has no charged candidates")
    for index, row in enumerate(rows):
        for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes"):
            if not finite_nonnegative(row.get(field)):
                raise VerificationError(f"GGMP discovery row {index} has invalid {field}")
    core = report.get("total_discovery_core_seconds")
    peak = report.get("peak_discovery_rss_bytes")
    if not finite_nonnegative(core) or not finite_nonnegative(peak):
        raise VerificationError("GGMP discovery aggregate resources are invalid")
    if not math.isclose(float(core), sum(float(row["total_core_seconds"]) for row in rows), rel_tol=1e-9, abs_tol=1e-9):
        raise VerificationError("GGMP discovery aggregate core-seconds do not sum")
    if int(peak) != max(int(row["peak_rss_bytes"]) for row in rows):
        raise VerificationError("GGMP discovery aggregate peak RSS is wrong")
    return {
        "path": str(recorded_repo / path.relative_to(live_repo)),
        "sha256": frozen["sha256"],
        "candidate_processes": len(rows),
        "wall_seconds_sequential_sum": sum(float(row["wall_seconds"]) for row in rows),
        "total_core_seconds": core,
        "peak_rss_bytes": peak,
        "selected": report["selected"],
        "charged_once": frozen.get("charge_once_per_panel") is True,
    }


def summarize(protocol: dict, panel: Path, write_receipts: bool, allow_incomplete: bool) -> dict:
    run, execution_identity_sha256 = validate_panel_run(protocol, panel)
    recorded_panel = recorded_panel_root(protocol, panel)
    custody_summary = validate_relocated_custody(panel, recorded_panel, run)
    recorded_repo = Path(run["repo"])
    live_repo = repository_root(HERE)
    ggmp_discovery = validate_fixed_ggmp_discovery(protocol, recorded_repo, live_repo)
    receipts = []
    failures = []
    for seed, cell in tasks(protocol):
        try:
            receipts.append(
                verify_cell(
                    protocol,
                    panel,
                    seed,
                    cell,
                    write_receipts,
                    recorded_panel=recorded_panel,
                )
            )
        except VerificationError as error:
            failures.append({"seed": seed, "cell_id": cell["id"], "error": str(error)})
    if failures and not allow_incomplete:
        first = failures[0]
        raise VerificationError(
            f"panel incomplete or invalid ({len(failures)} task failures); first: "
            f"{first['cell_id']} seed {first['seed']}: {first['error']}"
        )

    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for receipt in receipts:
        cell_id = receipt["task"]["cell"]["id"]
        for binding in receipt["backend_bindings"]:
            grouped[(cell_id, binding["backend"])].append(binding)
    distributions = []
    for (cell_id, backend), rows in sorted(grouped.items()):
        statuses: dict[str, int] = defaultdict(int)
        for row in rows:
            statuses[row["status"]] += 1
        distributions.append(
            {
                "cell_id": cell_id,
                "backend": backend,
                "replications": len(rows),
                "statuses": dict(sorted(statuses.items())),
                "conflicts": metric_distribution(
                    [row["conflicts"] for row in rows if row["conflicts"] is not None]
                ),
                "wall_seconds": metric_distribution(
                    [row["metrics"]["wall_seconds"] for row in rows if row["metrics"]["wall_seconds"] is not None]
                ),
                "total_core_seconds": metric_distribution(
                    [row["metrics"]["total_core_seconds"] for row in rows if row["metrics"]["total_core_seconds"] is not None]
                ),
                "single_core_seconds": metric_distribution(
                    [row["metrics"]["single_core_seconds"] for row in rows if row["metrics"]["single_core_seconds"] is not None]
                ),
                "peak_rss_bytes": metric_distribution(
                    [row["metrics"]["peak_rss_bytes"] for row in rows if row["metrics"]["peak_rss_bytes"] is not None]
                ),
                "point_witness_validation_core_seconds": metric_distribution(
                    [
                        row["point_witness_validation_metrics"]["total_core_seconds"]
                        for row in rows
                        if row["point_witness_validation_metrics"] is not None
                    ]
                ),
                "point_witness_validation_wall_seconds": metric_distribution(
                    [
                        row["point_witness_validation_metrics"]["wall_seconds"]
                        for row in rows
                        if row["point_witness_validation_metrics"] is not None
                    ]
                ),
                "resource_complete_replications": sum(row["resource_complete"] for row in rows),
            }
        )
    whole_cell_distributions = []
    wdsat_build_distributions = []
    by_cell: dict[str, list[dict]] = defaultdict(list)
    for receipt in receipts:
        by_cell[receipt["task"]["cell"]["id"]].append(receipt)
    for cell_id, rows in sorted(by_cell.items()):
        whole_cell_distributions.append(
            {
                "cell_id": cell_id,
                "replications": len(rows),
                "wall_seconds": metric_distribution(
                    [row["whole_cell_process"]["wall_seconds"] for row in rows]
                ),
                "single_core_seconds": metric_distribution(
                    [row["whole_cell_process"]["single_core_seconds"] for row in rows]
                ),
                "total_core_seconds": metric_distribution(
                    [row["whole_cell_process"]["total_core_seconds"] for row in rows]
                ),
                "peak_rss_bytes": metric_distribution(
                    [row["whole_cell_process"]["peak_rss_bytes"] for row in rows]
                ),
                "source_producer_core_seconds": metric_distribution(
                    [row["source_producer_process"]["total_core_seconds"] for row in rows]
                ),
                "predicate_construction_ns": metric_distribution(
                    [row["source_producer_timing_ns"]["factor_base_predicate_construction"] for row in rows]
                ),
                "source_system_construction_ns": metric_distribution(
                    [row["source_producer_timing_ns"]["source_system_construction"] for row in rows]
                ),
            }
        )
        builds = [row["wdsat_fixed_cost"] for row in rows if row["wdsat_fixed_cost"]["present"]]
        wdsat_build_distributions.append(
            {
                "cell_id": cell_id,
                "receipts": len(builds),
                "complete_receipts": sum(build["complete"] for build in builds),
                "source_copy_core_seconds": metric_distribution(
                    [build["source_copy"]["total_core_seconds"] for build in builds]
                ),
                "configuration_wall_seconds": metric_distribution(
                    [build["configuration_wall_seconds"] for build in builds]
                ),
                "clean_core_seconds": metric_distribution(
                    [build["clean"]["total_core_seconds"] for build in builds]
                ),
                "build_core_seconds": metric_distribution(
                    [build["build"]["total_core_seconds"] for build in builds]
                ),
                "build_peak_rss_bytes": metric_distribution(
                    [build["build"]["peak_rss_bytes"] for build in builds]
                ),
            }
        )
    expected = protocol["acceptance"]["expected_tasks"]
    artifact_complete = len(receipts) == expected and not failures
    resource_complete = artifact_complete and all(receipt["resource_complete"] for receipt in receipts)
    wdsat_build_complete = artifact_complete and all(
        receipt["wdsat_fixed_cost"]["complete"] for receipt in receipts
    )
    explicit_terminal = {"sat", "unsat"}
    solver_matrix_terminal_complete = artifact_complete and all(
        binding["status"] in explicit_terminal
        for receipt in receipts
        for binding in receipt["backend_bindings"]
    )
    magma_executed = any(
        binding["backend"] == "magma-f4" and binding["status"] != "unavailable_operational"
        for receipt in receipts
        for binding in receipt["backend_bindings"]
    )
    outer_core = sum(float(receipt["whole_cell_process"]["total_core_seconds"]) for receipt in receipts)
    outer_wall = sum(float(receipt["whole_cell_process"]["wall_seconds"]) for receipt in receipts)
    outer_peak = max(
        [int(receipt["whole_cell_process"]["peak_rss_bytes"]) for receipt in receipts] or [0]
    )
    evidence_class = run["evidence_class"]
    scientific_admissible = evidence_class == "scientific_candidate" and artifact_complete
    full_solver_gate = (
        scientific_admissible
        and solver_matrix_terminal_complete
        and resource_complete
        and wdsat_build_complete
        and magma_executed
    )
    summary = {
        "schema": "koblitz_target_matched_pdp_panel_summary.v1",
        "protocol_sha256": canonical_sha256(protocol),
        "panel": str(recorded_panel),
        "execution_identity_sha256": execution_identity_sha256,
        "evidence_class": evidence_class,
        "scientific_evidence_admissible": scientific_admissible,
        "expected_tasks": expected,
        "verified_tasks": len(receipts),
        "failures": failures,
        "panel_artifact_complete": artifact_complete,
        "solver_matrix_terminal_complete": solver_matrix_terminal_complete,
        "per_backend_resources_complete": resource_complete,
        "wdsat_build_receipts_complete": wdsat_build_complete,
        "magma_executed": magma_executed,
        "full_solver_matrix_gate_passed": full_solver_gate,
        "charged_totals": {
            "task_processes_included": len(receipts),
            "task_outer_total_core_seconds": outer_core,
            "task_outer_sequential_wall_seconds": outer_wall,
            "task_outer_peak_rss_bytes": outer_peak,
            "ggmp_discovery": ggmp_discovery,
            "cold_total_core_seconds_including_ggmp_discovery": outer_core
            + float(ggmp_discovery["total_core_seconds"]),
            "cold_sequential_wall_seconds_including_ggmp_discovery": outer_wall
            + float(ggmp_discovery["wall_seconds_sequential_sum"]),
            "cold_peak_rss_bytes": max(outer_peak, int(ggmp_discovery["peak_rss_bytes"])),
            "complete_twenty_task_total": artifact_complete,
        },
        "source_instances": [
            {
                "seed": receipt["task"]["seed"],
                "cell_id": receipt["task"]["cell"]["id"],
                "semantic_sha256": receipt["semantic_sha256"],
                "source_instance_sha256": receipt["source_instance_sha256"],
                "resource_complete": receipt["resource_complete"],
            }
            for receipt in receipts
        ],
        "distributions": distributions,
        "whole_cell_distributions": whole_cell_distributions,
        "wdsat_build_distributions": wdsat_build_distributions,
        "claim": (
            "Complete scientific-candidate artifact set for the matched planted-PDP panel; no end-to-end index-calculus or SOTA inference"
            if scientific_admissible
            else "Operational smoke or incomplete artifact set; no scaling, end-to-end index-calculus, or SOTA inference"
        ),
    }
    if custody_summary is not None:
        recomputed = (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode()
        if recomputed != custody_summary:
            raise VerificationError("relocated panel does not reproduce its custody-bound summary")
    return summary


def self_test() -> None:
    protocol = read_json(DEFAULT_PROTOCOL)
    validate_protocol(protocol)
    assert len(tasks(protocol)) == 20
    first = tasks(protocol)[0]
    assert first[0] == 2026091301
    assert first[1]["id"] == "n31-l5-m3-standard-a1-f0"
    assert canonical_sha256({"b": 2, "a": 1}) == canonical_sha256({"a": 1, "b": 2})
    assert metric_distribution([3, 1, 2]) == {"count": 3, "min": 1.0, "median": 2.0, "max": 3.0}
    seed, cell = tasks(protocol)[-1]
    with tempfile.TemporaryDirectory(prefix="stage13-verifier-self-test-") as directory:
        temp_root = Path(directory)
        partial = temp_root / "stage-13-panel-20260909"
        invocation_path = partial / task_relpath(seed, cell) / "invocation.json"
        invocation_path.parent.mkdir(parents=True)
        recorded = Path("/recorded/panel")
        invocation_path.write_text(
            json.dumps(
                {
                    "runner_command": [
                        "python3",
                        "runner.py",
                        "--output",
                        str(recorded / task_relpath(seed, cell) / "matrix"),
                    ]
                }
            )
        )
        assert recorded_panel_root(protocol, partial) == recorded
        assert recorded_panel_root(protocol, partial, (seed, cell)) == recorded
        summary_bytes = b'{"frozen":"summary"}\n'
        (partial / "panel-summary.json").write_bytes(summary_bytes)
        run = {"source_revision": {"commit": "frozen-source-revision"}}
        custody_path = temp_root / "stage-13-custody.json"
        custody = {
            "schema": "koblitz_stage13_panel_custody.v1",
            "panel_archive": partial.name,
            "panel_original_root": str(recorded),
            "panel_source_revision": run["source_revision"]["commit"],
            "panel_summary_sha256": hashlib.sha256(summary_bytes).hexdigest(),
        }
        custody_path.write_text(json.dumps(custody))
        assert validate_relocated_custody(partial, recorded, run) == summary_bytes
        for field, wrong in (
            ("panel_archive", "wrong-archive"),
            ("panel_original_root", "/wrong/root"),
            ("panel_source_revision", "wrong-revision"),
            ("panel_summary_sha256", "0" * 64),
        ):
            changed = dict(custody)
            changed[field] = wrong
            custody_path.write_text(json.dumps(changed))
            try:
                validate_relocated_custody(partial, recorded, run)
            except VerificationError:
                pass
            else:
                raise AssertionError(f"bad custody field {field} was accepted")
        custody_path.unlink()
        try:
            validate_relocated_custody(partial, recorded, run)
        except VerificationError:
            pass
        else:
            raise AssertionError("relocated archive without custody was accepted")
        outside = temp_root / "outside.json"
        outside.write_text("{}")
        escape = partial / "escape.json"
        escape.symlink_to(outside)
        try:
            archive_relative_path(escape, partial, "self-test escape")
        except VerificationError:
            pass
        else:
            raise AssertionError("symlink escape was accepted")
    print(json.dumps({"self_test": "pass", "expected_tasks": len(tasks(protocol))}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--output", type=Path, help="write the summary JSON here")
    parser.add_argument("--write-receipts", action="store_true")
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.panel is None:
        parser.error("--panel is required unless --self-test is used")
    protocol = read_json(args.protocol)
    validate_protocol(protocol)
    summary = summarize(protocol, args.panel, args.write_receipts, args.allow_incomplete)
    text = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
