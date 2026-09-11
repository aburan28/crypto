#!/usr/bin/env python3
"""Pure licensed-Magma execution helpers for the Phase-B replay packet."""

from __future__ import annotations

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
from typing import Any


MODEL_SCHEMA = "koblitz_magma_f4_model.v1"


class ExternalMagmaError(RuntimeError):
    """A Magma process or its terminal evidence violates the local contract."""


F4_SUFFIX = """\
cpu_start := Cputime();
wall_start := Realtime();
G, D := GroebnerBasis(I : Al := "Direct", Faugere := true, Dense := false, Nthreads := 1);
cpu_seconds := Cputime(cpu_start);
wall_seconds := Realtime(wall_start);
printf "KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1\\n";
printf "KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse\\n";
if #G eq 1 and G[1] eq R!1 then
  printf "KOBLITZ_MAGMA_STATUS=UNSAT\\n";
else
  printf "KOBLITZ_MAGMA_STATUS=SAT\\n";
end if;
printf "KOBLITZ_MAGMA_F4_DEGREES=%o\\n", D;
printf "KOBLITZ_MAGMA_BASIS_SIZE=%o\\n", #G;
printf "KOBLITZ_MAGMA_CPU_SECONDS=%o\\n", cpu_seconds;
printf "KOBLITZ_MAGMA_WALL_SECONDS=%o\\n", wall_seconds;
quit;
"""
F4_BODY_WITHOUT_QUIT = F4_SUFFIX.removesuffix("quit;\n")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ExternalMagmaError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ExternalMagmaError(f"expected a JSON object in {path}")
    return value


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def finite_nonnegative(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and value >= 0
    )


def validate_meter_record(record: dict[str, Any], expected_command: list[str], timeout: float) -> None:
    if record.get("command") != expected_command:
        raise ExternalMagmaError("process-meter command changed")
    if record.get("watchdog_seconds") != timeout:
        raise ExternalMagmaError("process-meter watchdog changed")
    if not isinstance(record.get("returncode"), int):
        raise ExternalMagmaError("process-meter returncode is missing")
    if type(record.get("timed_out")) is not bool or type(record.get("orphan_group_terminated")) is not bool:
        raise ExternalMagmaError("process-meter terminal booleans are invalid")
    metrics = record.get("metrics")
    if not isinstance(metrics, dict):
        raise ExternalMagmaError("process-meter metrics are missing")
    for name in ("wall_seconds", "user_seconds", "system_seconds", "total_core_seconds", "single_core_seconds"):
        if not finite_nonnegative(metrics.get(name)):
            raise ExternalMagmaError(f"invalid process metric {name}")
    peak = metrics.get("peak_rss_bytes")
    if not isinstance(peak, int) or isinstance(peak, bool) or peak < 0:
        raise ExternalMagmaError("invalid process metric peak_rss_bytes")
    if metrics.get("meter") != "fresh-process getrusage(RUSAGE_CHILDREN)":
        raise ExternalMagmaError("unexpected process meter")
    total = float(metrics["total_core_seconds"])
    if not math.isclose(total, float(metrics["user_seconds"]) + float(metrics["system_seconds"]), rel_tol=1e-9, abs_tol=1e-9):
        raise ExternalMagmaError("total core-seconds do not equal user plus system")
    if not math.isclose(total, float(metrics["single_core_seconds"]), rel_tol=1e-9, abs_tol=1e-9):
        raise ExternalMagmaError("single-core alias and total core-seconds disagree")


def run_metered(
    meter: Path,
    command: list[str],
    cwd: Path,
    timeout: float,
    attempt: Path,
    environment: dict[str, str] | None = None,
) -> dict[str, Any]:
    if attempt.exists() or attempt.is_symlink():
        raise ExternalMagmaError(f"refusing to overwrite attempt directory {attempt}")
    attempt.mkdir(parents=True)
    stdout_path = attempt / "stdout"
    stderr_path = attempt / "stderr"
    metrics_path = attempt / "metrics.json"
    receipt_path = attempt / "receipt.json"
    meter_command = [
        str(Path(sys.executable).resolve()), str(meter.resolve()),
        "--cwd", str(cwd.resolve()), "--timeout", str(timeout),
        "--stdout", str(stdout_path.resolve()), "--stderr", str(stderr_path.resolve()),
        "--metrics", str(metrics_path.resolve()), "--", *command,
    ]
    launcher = subprocess.run(meter_command, cwd=cwd, env=environment, check=False)
    if launcher.returncode != 0 or not metrics_path.is_file():
        raise ExternalMagmaError(
            f"process meter failed for {attempt}; launcher returncode {launcher.returncode}"
        )
    record = read_json(metrics_path)
    validate_meter_record(record, command, timeout)
    receipt = {
        "meter_launcher_command": meter_command,
        "meter_launcher_returncode": launcher.returncode,
        "process": record,
        "stdout": {
            "path": str(stdout_path.resolve()), "bytes": stdout_path.stat().st_size,
            "sha256": sha256_file(stdout_path),
        },
        "stderr": {
            "path": str(stderr_path.resolve()), "bytes": stderr_path.stat().st_size,
            "sha256": sha256_file(stderr_path),
        },
    }
    write_json(receipt_path, receipt)
    return receipt


def receipt_stream(receipt: dict[str, Any], stream: str) -> str:
    record = receipt.get(stream)
    if not isinstance(record, dict):
        raise ExternalMagmaError(f"metered {stream} record is missing")
    path = Path(str(record.get("path", "")))
    data = path.read_bytes()
    if len(data) != record.get("bytes") or sha256_bytes(data) != record.get("sha256"):
        raise ExternalMagmaError(f"metered {stream} changed: {path}")
    return data.decode(errors="replace")


def receipt_stdout(receipt: dict[str, Any]) -> str:
    return receipt_stream(receipt, "stdout")


def receipt_stderr(receipt: dict[str, Any]) -> str:
    return receipt_stream(receipt, "stderr")


def binary_identity(value: str | Path, version_options: list[list[str]]) -> dict[str, Any]:
    text = str(value)
    located = shutil.which(text) if "/" not in text else text
    if located is None:
        raise ExternalMagmaError(f"cannot resolve executable {text}")
    path = Path(located).resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ExternalMagmaError(f"executable is missing or not executable: {path}")
    version = None
    version_command = None
    version_returncode = None
    failed = []
    for suffix in version_options:
        command = [str(path), *suffix]
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=15, check=False)
        except (OSError, subprocess.SubprocessError) as error:
            failed.append({"command": command, "error": str(error)})
            continue
        output = (result.stdout + result.stderr).strip()
        if output and result.returncode == 0:
            version = output[:4000]
            version_command = command
            version_returncode = 0
            break
        failed.append({"command": command, "returncode": result.returncode, "output": output[:1000]})
    return {
        "path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path),
        "version": version or "unreported", "version_command": version_command,
        "version_returncode": version_returncode, "failed_version_probes": failed,
    }


def recognized_magma_version(output: str) -> bool:
    return bool(re.search(r"\bMagma V2\.\d+-\d+\b", output)) or any(
        re.fullmatch(r"2\.\d+-\d+", line.strip()) is not None for line in output.splitlines()
    )


def cpu_model() -> str:
    if platform.system() == "Darwin":
        result = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True, check=False)
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        match = re.search(r"^model name\s*:\s*(.+)$", cpuinfo.read_text(), re.MULTILINE)
        if match:
            return match.group(1).strip()
    return platform.processor() or "unreported"


def host_identity() -> dict[str, Any]:
    return {
        "node": platform.node(), "platform": platform.platform(),
        "system": platform.system(), "release": platform.release(),
        "machine": platform.machine(), "cpu_model": cpu_model(),
        "logical_cpus": os.cpu_count(), "python": platform.python_version(),
        "python_executable": str(Path(sys.executable).resolve()),
    }


def split_magma_source(path: Path) -> str:
    text = path.read_text()
    if not text.endswith(F4_SUFFIX):
        raise ExternalMagmaError(f"{path}: Magma input is not the frozen direct-F4 template")
    prefix = text[: -len(F4_SUFFIX)]
    required = ("SetNthreads(1);\n", "SetGPU(false);\n", "SetSeed(1);\n", "R := BooleanPolynomialRing(")
    if not all(item in prefix for item in required) or not prefix.endswith("I := ideal<R | F>;\n"):
        raise ExternalMagmaError(f"{path}: unsafe or incomplete Magma source prefix")
    return prefix


def render_witness_script(magma_path: Path, exclusions: list[list[bool]], source_variables: int) -> str:
    if not isinstance(source_variables, int) or source_variables <= 0:
        raise ExternalMagmaError("source_variables must be positive")
    for index, assignment in enumerate(exclusions):
        if len(assignment) != source_variables or any(type(bit) is not bool for bit in assignment):
            raise ExternalMagmaError(f"invalid excluded assignment {index}")
    prefix = split_magma_source(magma_path)
    rows = ["  [ " + ", ".join(f"K!{int(bit)}" for bit in assignment) + " ]" for assignment in exclusions]
    excluded = "[\n" + ",\n".join(rows) + "\n]" if rows else "[]"
    sat_call = "sat_model, S := SAT(G : Verbose := false);" if not exclusions else "sat_model, S := SAT(G : Exclude := Excluded, Verbose := false);"
    model = f"""\
K := GF(2);
Excluded := {excluded};
model_cpu_start := Cputime();
model_wall_start := Realtime();
{sat_call}
model_cpu_seconds := Cputime(model_cpu_start);
model_wall_seconds := Realtime(model_wall_start);
printf "KOBLITZ_MODEL_SCHEMA={MODEL_SCHEMA}\\n";
if sat_model then
  printf "KOBLITZ_MODEL_STATUS=SAT\\n";
  printf "KOBLITZ_MODEL_BITS=";
  for value in S do
    printf "%o", Integers()!value;
  end for;
  printf "\\n";
else
  printf "KOBLITZ_MODEL_STATUS=UNSAT\\n";
  printf "KOBLITZ_MODEL_BITS=-\\n";
end if;
printf "KOBLITZ_MODEL_VARIABLES={source_variables}\\n";
printf "KOBLITZ_MODEL_EXCLUDED_COUNT={len(exclusions)}\\n";
printf "KOBLITZ_MODEL_CPU_SECONDS=%o\\n", model_cpu_seconds;
printf "KOBLITZ_MODEL_WALL_SECONDS=%o\\n", model_wall_seconds;
quit;
"""
    return prefix + F4_BODY_WITHOUT_QUIT + model


def parse_model_terminal(output: str, source_variables: int, excluded_count: int) -> dict[str, Any] | None:
    markers: dict[str, str] = {}
    for line in output.splitlines():
        if not line.strip().startswith("KOBLITZ_MODEL_") or "=" not in line:
            continue
        key, value = line.strip().split("=", 1)
        if key in markers:
            return None
        markers[key] = value.strip()
    required = {
        "KOBLITZ_MODEL_SCHEMA", "KOBLITZ_MODEL_STATUS", "KOBLITZ_MODEL_BITS",
        "KOBLITZ_MODEL_VARIABLES", "KOBLITZ_MODEL_EXCLUDED_COUNT",
        "KOBLITZ_MODEL_CPU_SECONDS", "KOBLITZ_MODEL_WALL_SECONDS",
    }
    if set(markers) != required or markers["KOBLITZ_MODEL_SCHEMA"] != MODEL_SCHEMA:
        return None
    status = markers["KOBLITZ_MODEL_STATUS"]
    if status not in {"SAT", "UNSAT"}:
        return None
    try:
        variables = int(markers["KOBLITZ_MODEL_VARIABLES"])
        excluded = int(markers["KOBLITZ_MODEL_EXCLUDED_COUNT"])
        cpu = float(markers["KOBLITZ_MODEL_CPU_SECONDS"])
        wall = float(markers["KOBLITZ_MODEL_WALL_SECONDS"])
    except ValueError:
        return None
    if variables != source_variables or excluded != excluded_count or not finite_nonnegative(cpu) or not finite_nonnegative(wall):
        return None
    bits = markers["KOBLITZ_MODEL_BITS"]
    if status == "SAT":
        if len(bits) != source_variables or set(bits) - {"0", "1"}:
            return None
        assignment = [bit == "1" for bit in bits]
    else:
        if bits != "-":
            return None
        assignment = None
    return {
        "schema": MODEL_SCHEMA, "status": status.lower(), "source_variables": variables,
        "excluded_count": excluded, "assignment": assignment,
        "cpu_seconds": cpu, "wall_seconds": wall,
    }


def parse_magma_terminal(output: str) -> dict[str, Any] | None:
    markers: dict[str, str] = {}
    for line in output.splitlines():
        line = line.strip()
        if not line.startswith("KOBLITZ_MAGMA_") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in markers:
            return None
        markers[key] = value.strip()
    required = {
        "KOBLITZ_MAGMA_SCHEMA", "KOBLITZ_MAGMA_ALGORITHM", "KOBLITZ_MAGMA_STATUS",
        "KOBLITZ_MAGMA_F4_DEGREES", "KOBLITZ_MAGMA_BASIS_SIZE",
        "KOBLITZ_MAGMA_CPU_SECONDS", "KOBLITZ_MAGMA_WALL_SECONDS",
    }
    if set(markers) != required or markers["KOBLITZ_MAGMA_SCHEMA"] != "koblitz_magma_f4_terminal.v1" or markers["KOBLITZ_MAGMA_ALGORITHM"] != "direct-f4-sparse":
        return None
    status = markers["KOBLITZ_MAGMA_STATUS"]
    degrees_text = markers["KOBLITZ_MAGMA_F4_DEGREES"]
    if status not in {"SAT", "UNSAT"} or re.fullmatch(r"\[\s*(?:\d+(?:\s*,\s*\d+)*)?\s*\]", degrees_text) is None:
        return None
    try:
        basis_size = int(markers["KOBLITZ_MAGMA_BASIS_SIZE"])
        cpu = float(markers["KOBLITZ_MAGMA_CPU_SECONDS"])
        wall = float(markers["KOBLITZ_MAGMA_WALL_SECONDS"])
    except ValueError:
        return None
    if basis_size < 1 or not finite_nonnegative(cpu) or not finite_nonnegative(wall) or (status == "UNSAT" and basis_size != 1):
        return None
    return {
        "schema": markers["KOBLITZ_MAGMA_SCHEMA"], "algorithm": markers["KOBLITZ_MAGMA_ALGORITHM"],
        "terminal_status": status.lower(), "f4_step_degrees": [int(value) for value in re.findall(r"\d+", degrees_text)],
        "basis_size": basis_size, "cpu_seconds": cpu, "wall_seconds": wall,
        "single_thread_requested": True, "gpu_disabled": True,
    }


def assignment_validation_status(run: dict[str, Any], manifest: dict[str, Any], assignment_path: Path) -> dict[str, Any]:
    try:
        report = json.loads(run["stdout"])
    except (json.JSONDecodeError, TypeError):
        report = None
    try:
        assignment = json.loads(assignment_path.read_bytes())
    except json.JSONDecodeError:
        assignment = None
    expected_id = manifest.get("source_instance", {}).get("id_blake3")
    valid = (
        isinstance(report, dict)
        and report.get("schema") == "koblitz_pdp_assignment_validation.v1"
        and report.get("source_instance_id") == expected_id
        and report.get("source_instance_verified") is True
        and report.get("regenerated_source_exact") is True
        and report.get("assignment_values") == manifest.get("source_variables")
        and report.get("source_assignment") == assignment
        and isinstance(report.get("assignment_blake3"), str)
        and re.fullmatch(r"[0-9a-f]{64}", report["assignment_blake3"]) is not None
    )
    if run["timed_out"]:
        status = "timeout_inconclusive"
    elif not valid:
        status = "validation_contract_error"
    elif run["returncode"] == 0 and report.get("status") == "valid_point_witness" and report.get("source_model_valid") is True and report.get("source_witness_valid") is True:
        status = "valid_point_witness"
    elif run["returncode"] == 2 and report.get("status") in {"invalid_source_model", "nonlifting_source_model"}:
        status = report["status"]
    else:
        status = "validation_backend_error"
    return {"status": status, "report": report}


class _MatrixContract:
    parse_magma_terminal = staticmethod(parse_magma_terminal)
    assignment_validation_status = staticmethod(assignment_validation_status)


MATRIX = _MatrixContract()


def classify_f4(receipt: dict[str, Any]) -> tuple[str, dict[str, Any] | None]:
    process = receipt["process"]
    stdout = receipt_stdout(receipt)
    stderr = receipt_stderr(receipt)
    if process["timed_out"]:
        return "timeout_inconclusive", None
    terminal = parse_magma_terminal(stdout)
    if process["returncode"] != 0 or "Runtime error" in stdout + stderr or "User error" in stdout + stderr:
        return "solver_error", terminal
    if terminal is None:
        return "terminal_certificate_missing", None
    if terminal["terminal_status"] == "unsat":
        return "unsat", terminal
    return "sat_basis_certificate_unverified_model", terminal


def same_f4_terminal(primary: dict[str, Any], repeated: dict[str, Any]) -> bool:
    return all(primary.get(key) == repeated.get(key) for key in (
        "schema", "algorithm", "terminal_status", "f4_step_degrees",
        "basis_size", "single_thread_requested", "gpu_disabled",
    ))


def validation_status(
    receipt: dict[str, Any], manifest: dict[str, Any], assignment_path: Path
) -> tuple[str, dict[str, Any] | None]:
    process = receipt["process"]
    interpreted = assignment_validation_status({
        "stdout": receipt_stdout(receipt), "stderr": receipt_stderr(receipt),
        "returncode": process["returncode"], "timed_out": process["timed_out"],
        "metrics": process["metrics"], "command": process["command"],
    }, manifest, assignment_path)
    return interpreted["status"], interpreted.get("report")
