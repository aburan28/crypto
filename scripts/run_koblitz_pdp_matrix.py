#!/usr/bin/env python3
"""Run matched binary-Koblitz PDP instances through available solver backends.

Every subprocess is wrapped by ``process_meter.py``.  Missing tools and
watchdog expirations are recorded as operational outcomes; they are never
converted to UNSAT or scientific evidence.  The producer only exports the
source instance; native SAT and direct MITM run in separate metered processes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time


def run_timed(command: list[str], timeout: float, cwd: Path) -> dict:
    meter = Path(__file__).with_name("process_meter.py")
    with tempfile.TemporaryDirectory() as directory:
        scratch = Path(directory)
        stdout_path = scratch / "stdout"
        stderr_path = scratch / "stderr"
        metrics_path = scratch / "metrics.json"
        subprocess.run(
            [
                sys.executable,
                str(meter),
                "--cwd",
                str(cwd),
                "--timeout",
                str(timeout),
                "--stdout",
                str(stdout_path),
                "--stderr",
                str(stderr_path),
                "--metrics",
                str(metrics_path),
                "--",
                *command,
            ],
            check=True,
        )
        record = json.loads(metrics_path.read_text())
        record["stdout"] = stdout_path.read_text()
        record["stderr"] = stderr_path.read_text()
        return record


def parse_wdsat_model(stdout: str, n_vars: int) -> list[bool] | None:
    for line in stdout.splitlines():
        value = line.strip()
        if len(value) >= n_vars and set(value) <= {"0", "1"}:
            return [bit == "1" for bit in value[:n_vars]]
    return None


def validate_wdsat_anf(path: Path, model: list[bool]) -> bool:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    for line in lines[1:]:
        tokens = line.split()[1:-1]
        monomials: list[list[int]] = []
        has_t = False
        i = 0
        while i < len(tokens):
            token = tokens[i]
            if token == "T":
                has_t = True
                i += 1
            elif token.startswith("."):
                degree = int(token[1:])
                monomials.append([int(value) - 1 for value in tokens[i + 1 : i + 1 + degree]])
                i += degree + 1
            else:
                monomials.append([int(token) - 1])
                i += 1
        # Export convention: T is present exactly when the original
        # polynomial has no constant term.
        parity = not has_t
        for monomial in monomials:
            parity ^= all(model[index] for index in monomial)
        if parity:
            return False
    return True


def parse_cms_model(stdout: str, max_var: int) -> list[bool] | None:
    values: dict[int, bool] = {}
    for line in stdout.splitlines():
        if not line.startswith("v "):
            continue
        for token in line.split()[1:]:
            literal = int(token)
            if literal:
                values[abs(literal)] = literal > 0
    if not all(index in values for index in range(1, max_var + 1)):
        return None
    return [values[index] for index in range(1, max_var + 1)]


def validate_xor_dimacs(path: Path, model: list[bool]) -> bool:
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line[0] in {"c", "p"}:
            continue
        xor = line.startswith("x")
        tokens = line.split()[1:-1] if xor else line.split()[:-1]
        literals = [int(value) for value in tokens]
        values = [model[abs(literal) - 1] if literal > 0 else not model[abs(literal) - 1] for literal in literals]
        if xor:
            if sum(values) % 2 != 1:
                return False
        elif not any(values):
            return False
    return True


def source_artifact_snapshot(instance: Path, manifest: dict) -> dict:
    """Hash the three immutable source exports at one custody boundary."""
    exports = manifest.get("exports")
    expected = {"wdsat_anf", "cryptominisat_xor_dimacs", "magma_boolean_f4"}
    if not isinstance(exports, dict) or set(exports) != expected:
        raise ValueError("manifest must contain exactly the three source exports")
    receipt = {}
    for name in sorted(expected):
        descriptor = exports[name]
        relative = Path(descriptor.get("path", ""))
        if relative.is_absolute() or len(relative.parts) != 1 or relative.parts[0] in {"", ".", ".."}:
            raise ValueError(f"unsafe source export path for {name}")
        path = instance / relative
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"source export for {name} is missing or is a symlink")
        data = path.read_bytes()
        if descriptor.get("bytes") != len(data):
            raise ValueError(f"source export byte count changed for {name}")
        receipt[name] = {
            "path": str(relative),
            "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "manifest_blake3": descriptor.get("blake3"),
        }
    return receipt


def parse_magma_terminal(output: str) -> dict | None:
    """Parse the complete marker block emitted by ``instance.magma``."""
    marker_prefix = "KOBLITZ_MAGMA_"
    markers: dict[str, str] = {}
    for line in output.splitlines():
        line = line.strip()
        if not line.startswith(marker_prefix) or "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in markers:
            return None
        markers[key] = value.strip()
    required = {
        "KOBLITZ_MAGMA_SCHEMA",
        "KOBLITZ_MAGMA_ALGORITHM",
        "KOBLITZ_MAGMA_STATUS",
        "KOBLITZ_MAGMA_F4_DEGREES",
        "KOBLITZ_MAGMA_BASIS_SIZE",
        "KOBLITZ_MAGMA_CPU_SECONDS",
        "KOBLITZ_MAGMA_WALL_SECONDS",
    }
    if set(markers) != required:
        return None
    if markers["KOBLITZ_MAGMA_SCHEMA"] != "koblitz_magma_f4_terminal.v1":
        return None
    if markers["KOBLITZ_MAGMA_ALGORITHM"] != "direct-f4-sparse":
        return None
    status = markers["KOBLITZ_MAGMA_STATUS"]
    if status not in {"SAT", "UNSAT"}:
        return None
    degrees_text = markers["KOBLITZ_MAGMA_F4_DEGREES"]
    if not re.fullmatch(r"\[\s*(?:\d+(?:\s*,\s*\d+)*)?\s*\]", degrees_text):
        return None
    degrees = [int(value) for value in re.findall(r"\d+", degrees_text)]
    try:
        basis_size = int(markers["KOBLITZ_MAGMA_BASIS_SIZE"])
        cpu_seconds = float(markers["KOBLITZ_MAGMA_CPU_SECONDS"])
        wall_seconds = float(markers["KOBLITZ_MAGMA_WALL_SECONDS"])
    except ValueError:
        return None
    if (
        basis_size < 1
        or not math.isfinite(cpu_seconds)
        or not math.isfinite(wall_seconds)
        or cpu_seconds < 0
        or wall_seconds < 0
    ):
        return None
    if status == "UNSAT" and basis_size != 1:
        return None
    return {
        "schema": markers["KOBLITZ_MAGMA_SCHEMA"],
        "algorithm": markers["KOBLITZ_MAGMA_ALGORITHM"],
        "terminal_status": status.lower(),
        "f4_step_degrees": degrees,
        "basis_size": basis_size,
        "cpu_seconds": cpu_seconds,
        "wall_seconds": wall_seconds,
        "single_thread_requested": True,
        "gpu_disabled": True,
    }


def solver_status(run: dict, solver: str, instance: Path, manifest: dict) -> dict:
    terminal = None
    if run["timed_out"]:
        status = "timeout_inconclusive"
        conflicts = None
        model_valid = None
    elif solver == "wdsat":
        lines = [line.strip() for line in run["stdout"].splitlines() if line.strip()]
        model = parse_wdsat_model(run["stdout"], manifest["source_variables"])
        if model is not None:
            model_valid = validate_wdsat_anf(instance / "instance.anf", model)
            status = (
                "sat_source_model_unverified_point_witness"
                if model_valid
                else "sat_invalid_model"
            )
        elif any("UNSAT" in line for line in lines):
            status = "unsat"
            model_valid = None
        else:
            status = "solver_error"
            model_valid = None
        conflicts = next((int(line) for line in reversed(lines) if line.isdigit()), None)
    elif solver == "cryptominisat":
        if "s SATISFIABLE" in run["stdout"]:
            status = "sat"
            max_var = manifest["exports"]["cryptominisat_xor_dimacs"]["variables"]
            model = parse_cms_model(run["stdout"], max_var)
            model_valid = None if model is None else validate_xor_dimacs(instance / "instance.xor.cnf", model)
            if model_valid is not True:
                status = "sat_invalid_model"
            elif run["returncode"] != 10:
                status = "sat_invalid_terminal_status"
            else:
                status = "sat_source_model_unverified_point_witness"
        elif "s UNSATISFIABLE" in run["stdout"]:
            status = "unsat" if run["returncode"] == 20 else "unsat_invalid_terminal_status"
            model_valid = None
        elif run["returncode"] == 0:
            status = "unknown_inconclusive"
            model_valid = None
        else:
            status = "solver_error"
            model_valid = None
        matches = re.findall(r"(?:conflicts|Conflicts)\s*[:=]?\s*(\d+)", run["stdout"] + run["stderr"])
        conflicts = int(matches[-1]) if matches else None
    elif solver == "magma-f4":
        text = run["stdout"] + run["stderr"]
        conflicts = None
        model_valid = None
        terminal = parse_magma_terminal(run["stdout"])
        if run["returncode"] != 0 or "Runtime error" in text or "User error" in text:
            status = "solver_error"
        elif terminal is None:
            status = "terminal_certificate_missing"
        elif terminal["terminal_status"] == "unsat":
            status = "unsat"
        else:
            # A non-unit Boolean Gröbner basis establishes a proper ideal,
            # but this script deliberately does not run a second solver to
            # extract a model.  Keep the distinction visible and fail closed.
            status = "sat_basis_certificate_unverified_model"
    else:
        raise ValueError(f"unknown solver {solver}")
    return {
        "solver": solver,
        "status": status,
        "conflicts": conflicts,
        "source_model_valid": model_valid,
        "source_witness_valid": None,
        "returncode": run["returncode"],
        "timed_out": run["timed_out"],
        "metrics": run["metrics"],
        "command": run["command"],
        "magma_terminal": terminal,
    }


def parsed_source_assignment(run: dict, solver: str, manifest: dict) -> list[bool] | None:
    source_variables = manifest["source_variables"]
    if solver == "wdsat":
        model = parse_wdsat_model(run["stdout"], source_variables)
    elif solver == "cryptominisat":
        maximum = manifest["exports"]["cryptominisat_xor_dimacs"]["variables"]
        model = parse_cms_model(run["stdout"], maximum)
    else:
        raise ValueError(f"solver {solver} has no parsed assignment")
    if model is None or len(model) < source_variables:
        return None
    return model[:source_variables]


def isolated_backend_status(run: dict, backend: str, manifest: dict) -> dict:
    """Validate one isolated native-backend terminal record fail closed."""
    report = None
    try:
        report = json.loads(run["stdout"])
    except (json.JSONDecodeError, TypeError):
        pass
    if run["timed_out"]:
        status = "timeout_inconclusive"
    elif report is None:
        status = "backend_error" if run["returncode"] != 0 else "backend_contract_error"
    else:
        expected_id = manifest.get("source_instance", {}).get("id_blake3")
        artifacts = report.get("source_artifacts")
        artifacts_valid = (
            isinstance(artifacts, dict)
            and set(artifacts)
            == {"wdsat_anf", "cryptominisat_xor_dimacs", "magma_boolean_f4"}
            and all(
                isinstance(receipt, dict) and receipt.get("valid") is True
                for receipt in artifacts.values()
            )
        )
        contract_valid = (
            report.get("schema") == "koblitz_pdp_isolated_backend.v1"
            and report.get("backend") == backend
            and expected_id is not None
            and report.get("source_instance_id") == expected_id
            and report.get("source_instance_verified") is True
            and report.get("regenerated_source_exact") is True
            and artifacts_valid
        )
        result_status = report.get("status")
        if not contract_valid:
            status = "backend_contract_error"
        elif backend == "native-sat":
            if (
                run["returncode"] == 0
                and result_status == "sat"
                and report.get("source_model_valid") is True
                and report.get("source_witness_valid") is True
            ):
                status = "sat"
            elif run["returncode"] == 0 and result_status in {
                "unsat",
                "unknown_inconclusive",
                "model_cap_inconclusive",
            }:
                status = result_status
            elif run["returncode"] == 2 and result_status == "sat_invalid_model":
                status = "sat_invalid_model"
            else:
                status = "backend_contract_error"
        elif (
            run["returncode"] == 0
            and result_status == "sat"
            and report.get("source_witness_valid") is True
        ):
            status = "sat"
        elif (
            run["returncode"] == 0
            and result_status == "unsat"
            and report.get("exhaustive") is True
        ):
            status = "unsat"
        elif (
            run["returncode"] == 0
            and result_status == "not_run_resource_cap"
            and report.get("exhaustive") is False
        ):
            status = "not_run_resource_cap"
        elif run["returncode"] == 2 and result_status == "sat_invalid_witness":
            status = "sat_invalid_witness"
        else:
            status = "backend_contract_error"
    stats = report.get("stats", {}) if isinstance(report, dict) else {}
    return {
        "solver": backend,
        "status": status,
        "conflicts": stats.get("conflicts"),
        "source_model_valid": (
            report.get("source_model_valid") if isinstance(report, dict) else None
        ),
        "source_witness_valid": (
            report.get("source_witness_valid") if isinstance(report, dict) else None
        ),
        "source_instance_id": (
            report.get("source_instance_id") if isinstance(report, dict) else None
        ),
        "returncode": run["returncode"],
        "timed_out": run["timed_out"],
        "metrics": run["metrics"],
        "command": run["command"],
        "backend_report": report,
    }


def assignment_validation_status(run: dict, manifest: dict, assignment_path: Path) -> dict:
    """Validate a metered rational point-witness checker receipt."""
    report = None
    try:
        report = json.loads(run["stdout"])
    except (json.JSONDecodeError, TypeError):
        pass
    assignment_bytes = assignment_path.read_bytes()
    try:
        expected_assignment = json.loads(assignment_bytes)
    except json.JSONDecodeError:
        expected_assignment = None
    expected_id = manifest.get("source_instance", {}).get("id_blake3")
    contract_valid = (
        isinstance(report, dict)
        and report.get("schema") == "koblitz_pdp_assignment_validation.v1"
        and report.get("source_instance_id") == expected_id
        and report.get("source_instance_verified") is True
        and report.get("regenerated_source_exact") is True
        and report.get("assignment_values") == manifest.get("source_variables")
        and report.get("source_assignment") == expected_assignment
        and isinstance(report.get("assignment_blake3"), str)
        and re.fullmatch(r"[0-9a-f]{64}", report["assignment_blake3"])
    )
    if run["timed_out"]:
        status = "timeout_inconclusive"
    elif not contract_valid:
        status = "validation_contract_error"
    elif run["returncode"] == 0:
        status = (
            "valid_point_witness"
            if report.get("status") == "valid_point_witness"
            and report.get("source_model_valid") is True
            and report.get("source_witness_valid") is True
            else "validation_contract_error"
        )
    elif run["returncode"] == 2 and report.get("status") in {
        "invalid_source_model",
        "nonlifting_source_model",
    }:
        status = report["status"]
    else:
        status = "validation_backend_error"
    return {
        "status": status,
        "returncode": run["returncode"],
        "timed_out": run["timed_out"],
        "metrics": run["metrics"],
        "command": run["command"],
        "report": report,
    }
def run_assignment_validation(
    backend: Path,
    manifest_path: Path,
    manifest: dict,
    assignment: list[bool],
    label: str,
    timeout: float,
    instance: Path,
) -> dict:
    assignment_path = instance / f"{label}.source-model.json"
    assignment_bytes = (json.dumps(assignment, separators=(",", ":")) + "\n").encode()
    assignment_path.write_bytes(assignment_bytes)
    run = run_timed(
        [
            str(backend.resolve()),
            "validate-model",
            str(manifest_path),
            str(assignment_path),
        ],
        timeout,
        instance,
    )
    (instance / f"{label}.point-validation.stdout").write_text(run["stdout"])
    (instance / f"{label}.point-validation.stderr").write_text(run["stderr"])
    result = assignment_validation_status(run, manifest, assignment_path)
    result["assignment_path"] = assignment_path.name
    result["assignment_sha256"] = hashlib.sha256(assignment_bytes).hexdigest()
    if isinstance(result.get("report"), dict):
        # The backend's BLAKE3 binds the same file; its source identity and the
        # runner's independent SHA-256 custody are both retained.
        result["assignment_blake3"] = result["report"].get("assignment_blake3")
        if result["report"].get("assignment_values") != len(assignment):
            result["status"] = "validation_contract_error"
    return result


def version(binary: str | None) -> dict:
    if binary is None:
        return {"available": False, "path": None, "version": None}
    path = shutil.which(binary) if "/" not in binary else binary
    if not path or not Path(path).exists():
        return {"available": False, "path": path, "version": None}
    for option in ([path, "--version"], [path, "-h"]):
        try:
            result = subprocess.run(option, capture_output=True, text=True, timeout=10)
            text = (result.stdout + result.stderr).strip()
            if text:
                resolved = Path(path).resolve()
                return {
                    "available": True,
                    "path": str(resolved),
                    "version": text[:1000],
                    "sha256": hashlib.sha256(resolved.read_bytes()).hexdigest(),
                }
        except (OSError, subprocess.SubprocessError):
            pass
    resolved = Path(path).resolve()
    return {
        "available": True,
        "path": str(resolved),
        "version": "unreported",
        "sha256": hashlib.sha256(resolved.read_bytes()).hexdigest(),
    }


def git_commit(path: Path) -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=path,
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None


def wdsat_xor_atom_count(path: Path) -> int:
    """Count the global XOR atoms consumed by WDSat's ANF parser."""
    count = 0
    for line in path.read_text().splitlines():
        tokens = line.split()
        if not tokens or tokens[0] != "x":
            continue
        index = 1
        while index < len(tokens) and tokens[index] != "0":
            token = tokens[index]
            if token.startswith("."):
                try:
                    degree = int(token[1:])
                except ValueError as error:
                    raise ValueError(f"invalid ANF degree token {token!r}") from error
                if degree < 2 or index + degree >= len(tokens):
                    raise ValueError(f"invalid ANF monomial at token {index}")
                index += degree + 1
            else:
                index += 1
            count += 1
        if index >= len(tokens) or tokens[index] != "0":
            raise ValueError("ANF equation lacks its zero terminator")
    return count


def build_wdsat(source: Path, instance: Path, manifest: dict, timeout: float) -> tuple[dict, str | None]:
    """Build an instance-sized WDSat binary and return its charged receipt."""
    build_root = instance / "wdsat-build"
    build_root.mkdir()
    copy_receipt = run_timed(
        ["cp", "-R", str((source / "src").resolve()), str(build_root / "src")],
        timeout,
        instance,
    )
    cms = manifest["exports"]["cryptominisat_xor_dimacs"]
    n_source = int(manifest["source_variables"])
    max_degree = int(manifest["source_max_degree"])
    max_id = int(cms["variables"])
    max_eq = int(cms["cnf_clauses"])
    max_xeq = int(cms["xor_rows"])
    max_terms = int(manifest["source_max_monomials_per_equation"])
    source_xor_atoms = wdsat_xor_atom_count(instance / "instance.anf")
    config = "\n".join(
        [
            "#define __XG_ENHANCED__",
            f"#define __MAX_ANF_ID__ {n_source + 1}",
            f"#define __MAX_DEGREE__ {max_degree + 1}",
            f"#define __MAX_ID__ {max_id}",
            f"#define __MAX_BUFFER_SIZE__ {max(5000, max_terms * 8, max_id * 16, source_xor_atoms + 1)}",
            f"#define __MAX_EQ__ {max(64, max_eq + 16)}",
            f"#define __MAX_EQ_SIZE__ {max_degree + 2}",
            f"#define __MAX_XEQ__ {max(8, max_xeq + 2)}",
            f"#define __MAX_XEQ_SIZE__ {max(max_id + 1, max_terms + 2)}",
            "",
        ]
    )
    config_start = time.perf_counter()
    (build_root / "src" / "config.h").write_text(config)
    config_wall = time.perf_counter() - config_start
    clean = run_timed(["make", "-C", "src", "clean"], timeout, build_root)
    run = run_timed(["make", "-C", "src"], timeout, build_root)
    (build_root / "clean.stdout").write_text(clean["stdout"])
    (build_root / "clean.stderr").write_text(clean["stderr"])
    (build_root / "build.stdout").write_text(run["stdout"])
    (build_root / "build.stderr").write_text(run["stderr"])
    binary = build_root / "wdsat_solver"
    receipt = {
        "status": "completed" if run["returncode"] == 0 and binary.exists() else (
            "timeout_inconclusive" if run["timed_out"] else "failed_operational"
        ),
        "source_copy": {
            "returncode": copy_receipt["returncode"],
            "timed_out": copy_receipt["timed_out"],
            "metrics": copy_receipt["metrics"],
        },
        "configuration_wall_seconds": config_wall,
        "clean": {
            "returncode": clean["returncode"],
            "timed_out": clean["timed_out"],
            "metrics": clean["metrics"],
            "command": clean["command"],
        },
        "metrics": run["metrics"],
        "command": run["command"],
        "config": config,
        "source_xor_atoms": source_xor_atoms,
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest() if binary.exists() else None,
    }
    return receipt, str(binary.resolve()) if binary.exists() else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new result directory")
    parser.add_argument("--exporter", type=Path, required=True, help="built koblitz_pdp_export binary")
    parser.add_argument(
        "--backend", type=Path, required=True, help="built koblitz_pdp_backend binary"
    )
    parser.add_argument("--wdsat", help="WDSat binary")
    parser.add_argument("--wdsat-source", type=Path, help="WDSat checkout to right-size and build per cell")
    parser.add_argument("--cryptominisat", help="CryptoMiniSat binary")
    parser.add_argument("--magma", help="Magma binary")
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--conflicts", type=int, default=100_000)
    parser.add_argument(
        "--configs",
        default="15:5:standard,31:5:standard,31:5:ggmp,41:5:standard,59:9:standard,67:9:standard",
        help="comma-separated n:ell:basis[:curve_a:factor_index] cells",
    )
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output must be a new path")
    args.output.mkdir(parents=True)
    args.output = args.output.resolve()
    if not args.exporter.exists():
        parser.error("exporter does not exist; build the release example first")
    if not args.backend.exists():
        parser.error("backend does not exist; build the release example first")

    tools = {
        "wdsat": (
            {
                "available": args.wdsat_source is not None and (args.wdsat_source / "src" / "makefile").exists(),
                "path": str(args.wdsat_source.resolve()) if args.wdsat_source else None,
                "version": "per-instance source build",
                "source_commit": git_commit(args.wdsat_source) if args.wdsat_source else None,
            }
            if args.wdsat_source
            else version(args.wdsat)
        ),
        "cryptominisat": version(args.cryptominisat),
        "magma": version(args.magma),
        "isolated_backend": version(str(args.backend.resolve())),
    }
    report = {
        "schema": "koblitz_pdp_matched_matrix.v1",
        "scope": "PDP construction and solver matrix; not an end-to-end ECDLP result",
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "logical_cpus": os.cpu_count(),
        },
        "policy": {
            "unknown_is_unsat": False,
            "solver_subprocess_resources_charged": True,
            "python_orchestration_included_in_child_sum": False,
            "whole_panel_requires_outer_process_meter": True,
            "single_thread_requested": True,
            "missing_tool_is_negative_evidence": False,
            "producer_runs_solver_backends": False,
            "native_sat_and_mitm_isolated": True,
        },
        "tools": tools,
        "instances": [],
    }

    for cell_index, cell in enumerate(args.configs.split(",")):
        parts = cell.split(":")
        if len(parts) not in (3, 5):
            parser.error(f"bad cell {cell!r}; expected n:ell:basis[:curve_a:factor_index]")
        n_text, ell_text, basis = parts[:3]
        curve_a, factor_index = (int(parts[3]), int(parts[4])) if len(parts) == 5 else (1, 0)
        n, ell = int(n_text), int(ell_text)
        suffix = f"-a{curve_a}-f{factor_index}" if len(parts) == 5 else ""
        instance = args.output / f"n{n}-l{ell}-m3-{basis}{suffix}"
        command = [
            str(args.exporter.resolve()),
            str(n),
            str(ell),
            basis,
            str(args.seed + cell_index),
            str(args.conflicts),
            str(instance),
        ]
        if len(parts) == 5:
            command.extend([str(curve_a), str(factor_index)])
        command.append("--export-only")
        generated = run_timed(command, args.timeout, args.exporter.parent)
        (args.output / f"{instance.name}.export.stdout").write_text(
            generated["stdout"].rstrip() + "\n"
        )
        (args.output / f"{instance.name}.export.stderr").write_text(generated["stderr"])
        entry = {
            "cell": {
                "n": n,
                "ell": ell,
                "m": 3,
                "basis": basis,
                "curve_a": curve_a,
                "factor_index": factor_index,
            },
            "generator": {
                "status": "timeout_inconclusive" if generated["timed_out"] else (
                    "completed" if generated["returncode"] == 0 else "failed_operational"
                ),
                "returncode": generated["returncode"],
                "metrics": generated["metrics"],
                "command": command,
            },
            "solvers": [],
        }
        manifest_path = instance / "manifest.json"
        if not manifest_path.exists():
            entry["artifact_status"] = {
                "status": "missing_manifest_operational",
                "expected_path": str(manifest_path),
                "asserts_nothing_about": "the PDP, factor base, or solver performance",
            }
            report["instances"].append(entry)
            (args.output / "progress.json").write_text(json.dumps(report, indent=2) + "\n")
            continue
        manifest = json.loads(manifest_path.read_text())
        entry["manifest"] = manifest
        if (
            manifest.get("native_sat", {}).get("status") != "not_run_in_export_process"
            or manifest.get("direct_meet_in_the_middle", {}).get("status")
            != "not_run_in_export_process"
        ):
            entry["artifact_status"] = {
                "status": "producer_isolation_contract_error",
                "asserts_nothing_about": "isolated native-SAT or MITM resources",
            }
            report["instances"].append(entry)
            (args.output / "progress.json").write_text(json.dumps(report, indent=2) + "\n")
            continue
        try:
            entry["source_artifacts_before"] = source_artifact_snapshot(instance, manifest)
        except (OSError, ValueError) as error:
            entry["artifact_status"] = {
                "status": "source_artifact_snapshot_failed",
                "error": str(error),
                "asserts_nothing_about": "backend correctness or performance",
            }
            report["instances"].append(entry)
            (args.output / "progress.json").write_text(json.dumps(report, indent=2) + "\n")
            continue

        native_run = run_timed(
            [
                str(args.backend.resolve()),
                "native-sat",
                str(manifest_path),
                str(args.conflicts),
            ],
            args.timeout,
            instance,
        )
        (instance / "native-sat.stdout").write_text(native_run["stdout"])
        (instance / "native-sat.stderr").write_text(native_run["stderr"])
        entry["solvers"].append(
            isolated_backend_status(native_run, "native-sat", manifest)
        )

        mitm_run = run_timed(
            [str(args.backend.resolve()), "direct-mitm", str(manifest_path)],
            args.timeout,
            instance,
        )
        (instance / "direct-mitm.stdout").write_text(mitm_run["stdout"])
        (instance / "direct-mitm.stderr").write_text(mitm_run["stderr"])
        entry["solvers"].append(
            isolated_backend_status(mitm_run, "direct-mitm", manifest)
        )

        wdsat_binary = tools["wdsat"]["path"]
        if args.wdsat_source and tools["wdsat"]["available"]:
            entry["wdsat_build"], wdsat_binary = build_wdsat(
                args.wdsat_source, instance, manifest, args.timeout
            )
        if tools["wdsat"]["available"] and wdsat_binary:
            branch_variables = ",".join(str(index) for index in range(1, 3 * ell + 1))
            run = run_timed(
                [wdsat_binary, "-i", str(instance / "instance.anf"), "-g", branch_variables],
                args.timeout,
                instance,
            )
            (instance / "wdsat.stdout").write_text(run["stdout"])
            (instance / "wdsat.stderr").write_text(run["stderr"])
            row = solver_status(run, "wdsat", instance, manifest)
            if row["status"] == "sat_source_model_unverified_point_witness":
                assignment = parsed_source_assignment(run, "wdsat", manifest)
                if assignment is None:
                    row["status"] = "sat_invalid_model"
                else:
                    validation = run_assignment_validation(
                        args.backend,
                        manifest_path,
                        manifest,
                        assignment,
                        "wdsat",
                        args.timeout,
                        instance,
                    )
                    row["point_witness_validation"] = validation
                    row["source_witness_valid"] = validation["status"] == "valid_point_witness"
                    row["status"] = (
                        "sat"
                        if validation["status"] == "valid_point_witness"
                        else "sat_nonlifting_model_inconclusive"
                        if validation["status"] == "nonlifting_source_model"
                        else "sat_point_validation_error"
                    )
            entry["solvers"].append(row)
            if args.wdsat_source:
                build_root = instance / "wdsat-build"
                for name in ["clean.stdout", "clean.stderr", "build.stdout", "build.stderr"]:
                    source = build_root / name
                    if source.exists():
                        (instance / f"wdsat-{name}").write_bytes(source.read_bytes())
                (instance / "wdsat-config.h").write_text(entry["wdsat_build"]["config"])
                shutil.rmtree(build_root)
        else:
            entry["solvers"].append({"solver": "wdsat", "status": "unavailable_operational"})

        if tools["cryptominisat"]["available"]:
            run = run_timed(
                [tools["cryptominisat"]["path"], "--verb", "1", "--threads", "1", str(instance / "instance.xor.cnf")],
                args.timeout,
                instance,
            )
            (instance / "cryptominisat.stdout").write_text(run["stdout"])
            (instance / "cryptominisat.stderr").write_text(run["stderr"])
            row = solver_status(run, "cryptominisat", instance, manifest)
            if row["status"] == "sat_source_model_unverified_point_witness":
                assignment = parsed_source_assignment(run, "cryptominisat", manifest)
                if assignment is None:
                    row["status"] = "sat_invalid_model"
                else:
                    validation = run_assignment_validation(
                        args.backend,
                        manifest_path,
                        manifest,
                        assignment,
                        "cryptominisat",
                        args.timeout,
                        instance,
                    )
                    row["point_witness_validation"] = validation
                    row["source_witness_valid"] = validation["status"] == "valid_point_witness"
                    row["status"] = (
                        "sat"
                        if validation["status"] == "valid_point_witness"
                        else "sat_nonlifting_model_inconclusive"
                        if validation["status"] == "nonlifting_source_model"
                        else "sat_point_validation_error"
                    )
            entry["solvers"].append(row)
        else:
            entry["solvers"].append({"solver": "cryptominisat", "status": "unavailable_operational"})

        if tools["magma"]["available"]:
            run = run_timed(
                [
                    tools["magma"]["path"],
                    "-t",
                    "1",
                    "-b",
                    str(instance / "instance.magma"),
                ],
                args.timeout,
                instance,
            )
            (instance / "magma.stdout").write_text(run["stdout"])
            (instance / "magma.stderr").write_text(run["stderr"])
            entry["solvers"].append(solver_status(run, "magma-f4", instance, manifest))
        else:
            entry["solvers"].append({"solver": "magma-f4", "status": "unavailable_operational"})
        try:
            entry["source_artifacts_after"] = source_artifact_snapshot(instance, manifest)
            entry["source_artifacts_unchanged"] = (
                entry["source_artifacts_before"] == entry["source_artifacts_after"]
            )
        except (OSError, ValueError) as error:
            entry["source_artifacts_after"] = {
                "status": "source_artifact_snapshot_failed",
                "error": str(error),
            }
            entry["source_artifacts_unchanged"] = False
        if not entry["source_artifacts_unchanged"]:
            entry["artifact_status"] = {
                "status": "source_artifacts_changed_during_backend_matrix",
                "asserts_nothing_about": "cross-backend comparability",
            }
        report["instances"].append(entry)
        (args.output / "progress.json").write_text(json.dumps(report, indent=2) + "\n")

    (args.output / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
