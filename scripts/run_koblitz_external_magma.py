#!/usr/bin/env python3
"""Run the frozen Stage 13 Koblitz PDP inputs on a licensed Magma host.

This runner consumes the retained inputs; it never regenerates or modifies
them.  Every experimental Magma solver command and witness-validation command
is executed through the repository process meter.  A proper Boolean ideal is
accepted as a PDP SAT result only after a separately metered source-model and
rational-point check.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import sys
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
STAGE = (
    REPO
    / "research"
    / "sat_factor_base_review_20260908"
    / "continuation-05-sota-gates"
)
DEFAULT_PROTOCOL = STAGE / "stage-13-pdp-panel-protocol.json"
DEFAULT_PANEL = STAGE / "stage-13-panel-20260909"
DEFAULT_METER = HERE / "process_meter.py"
DEFAULT_BACKEND = REPO / "target" / "release" / "examples" / "koblitz_pdp_backend"
RUN_SCHEMA = "koblitz_external_magma_run.v1"
TASK_SCHEMA = "koblitz_external_magma_task.v1"
MANIFEST_SCHEMA = "koblitz_external_magma_dry_run_manifest.v1"
MODEL_SCHEMA = "koblitz_magma_f4_model.v1"
MAX_MODEL_CAP = 64
STAGE13_ARCHIVE_COMMIT = "5dbe54dd468e914b5872d10d06e7c44d01b0a236"
PINNED_STAGE13_VERIFIER_SHA256 = "fed21d90bf21bb5eadd31b234980cc1b9b910c06480bfa65eb2a9244959b1749"
PINNED_MATRIX_CONTRACT_SHA256 = "794adfac8151df2887f600944ed899f826bd25adbfec12a10c61200c7e278bb6"
PINNED_METER_SHA256 = "8705343c10c129941ff3a1068be8211c0a495d2ea9a5200bc26b64ad26e3065f"
RELEVANT_RUNTIME_PATHS = (
    "scripts/run_koblitz_external_magma.py",
    "scripts/run_koblitz_pdp_matrix.py",
    "scripts/process_meter.py",
    "research/sat_factor_base_review_20260908/continuation-05-sota-gates/verify_stage13_pdp_panel.py",
    "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-13-pdp-panel-protocol.json",
    "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-13-custody.json",
    "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-13-panel-20260909",
)
ARCHIVE_PINNED_PATHS = RELEVANT_RUNTIME_PATHS[1:]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bootstrap_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


if bootstrap_sha256(HERE / "run_koblitz_pdp_matrix.py") != PINNED_MATRIX_CONTRACT_SHA256:
    raise RuntimeError("refusing to import an unpinned Koblitz matrix contract")
if (
    bootstrap_sha256(STAGE / "verify_stage13_pdp_panel.py")
    != PINNED_STAGE13_VERIFIER_SHA256
):
    raise RuntimeError("refusing to import an unpinned Stage 13 verifier")


MATRIX = load_module("koblitz_pdp_matrix_for_external_magma", HERE / "run_koblitz_pdp_matrix.py")
STAGE13 = load_module("koblitz_stage13_for_external_magma", STAGE / "verify_stage13_pdp_panel.py")


class ExternalMagmaError(RuntimeError):
    """The external run or one of its custody records is invalid."""


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


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def canonical_sha256(value: Any) -> str:
    return sha256_bytes(canonical_bytes(value))


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ExternalMagmaError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ExternalMagmaError(f"expected a JSON object in {path}")
    return value


def atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def safe_relative(path: Path, root: Path, label: str) -> str:
    if path.is_symlink() or not path.is_file():
        raise ExternalMagmaError(f"{label} must be a regular non-symlink file: {path}")
    try:
        relative = path.resolve(strict=True).relative_to(root.resolve(strict=True))
    except (OSError, ValueError) as error:
        raise ExternalMagmaError(f"{label} escapes the frozen panel: {path}") from error
    if ".." in relative.parts:
        raise ExternalMagmaError(f"unsafe relative path for {label}: {relative}")
    return str(relative)


def file_record(path: Path, root: Path, label: str) -> dict:
    return {
        "path": safe_relative(path, root, label),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def validate_stage13_archive(protocol_path: Path, panel: Path) -> tuple[dict, dict, dict]:
    protocol = STAGE13.read_json(protocol_path)
    STAGE13.validate_protocol(protocol)
    summary = STAGE13.summarize(
        protocol,
        panel,
        write_receipts=False,
        allow_incomplete=False,
    )
    if summary.get("expected_tasks") != 20 or summary.get("verified_tasks") != 20:
        raise ExternalMagmaError("Stage 13 archive is not the exact verified 20-task panel")
    if summary.get("panel_artifact_complete") is not True:
        raise ExternalMagmaError("Stage 13 archive is incomplete")
    custody_path = panel.parent / "stage-13-custody.json"
    custody = read_json(custody_path)
    if (
        custody.get("schema") != "koblitz_stage13_panel_custody.v1"
        or custody.get("panel_archive") != panel.name
        or custody.get("panel_files") != 599
        or custody.get("panel_summary_sha256") != sha256_file(panel / "panel-summary.json")
    ):
        raise ExternalMagmaError("Stage 13 custody record is inconsistent")
    try:
        panel_relative = panel.resolve(strict=True).relative_to(REPO.resolve(strict=True))
    except (OSError, ValueError) as error:
        raise ExternalMagmaError("Stage 13 panel is outside the repository") from error
    ancestry = subprocess.run(
        ["git", "merge-base", "--is-ancestor", STAGE13_ARCHIVE_COMMIT, "HEAD"],
        cwd=REPO,
        check=False,
    )
    if ancestry.returncode != 0:
        raise ExternalMagmaError("the pinned Stage 13 archive commit is not an ancestor of HEAD")
    archived_summary = subprocess.run(
        ["git", "show", f"{STAGE13_ARCHIVE_COMMIT}:{panel_relative}/panel-summary.json"],
        cwd=REPO,
        capture_output=True,
        check=False,
    )
    if (
        archived_summary.returncode != 0
        or sha256_bytes(archived_summary.stdout) != custody["panel_summary_sha256"]
    ):
        raise ExternalMagmaError("the pinned archive commit does not contain the custody-bound panel")
    return protocol, summary, custody


def task_id(seed: int, cell: dict) -> str:
    return f"seed-{seed}/{cell['id']}"


def build_dry_run_manifest(protocol_path: Path, panel: Path, meter: Path = DEFAULT_METER) -> dict:
    protocol_path = protocol_path.resolve()
    panel = panel.resolve()
    meter = meter.resolve()
    if sha256_file(meter) != PINNED_METER_SHA256:
        raise ExternalMagmaError("process meter does not match the pinned Stage 13 meter")
    if not archive_contracts_match():
        raise ExternalMagmaError(
            "Stage 13 verifier, matrix contract, protocol, custody, meter, or panel differs "
            "from the pinned archive commit"
        )
    protocol, summary, custody = validate_stage13_archive(protocol_path, panel)
    tasks = []
    for seed, cell in STAGE13.tasks(protocol):
        task_root = panel / STAGE13.task_relpath(seed, cell)
        instance = task_root / "matrix" / cell["id"]
        source_binding_path = task_root / "source-binding.json"
        source_binding = read_json(source_binding_path)
        manifest_path = instance / "manifest.json"
        manifest = read_json(manifest_path)
        if manifest.get("seed") != seed:
            raise ExternalMagmaError(f"{task_id(seed, cell)}: manifest seed changed")
        expected_cell = {
            name: cell[name] for name in ("n", "ell", "m", "basis", "curve_a", "factor_index")
        }
        observed_cell = {
            "n": manifest.get("n"),
            "ell": manifest.get("ell"),
            "m": manifest.get("m"),
            "basis": cell["basis"],
            "curve_a": manifest.get("curve_a"),
            "factor_index": manifest.get("factor_base_predicate", {}).get("factor_index", 0),
        }
        if observed_cell != expected_cell:
            raise ExternalMagmaError(f"{task_id(seed, cell)}: manifest cell changed")
        files = {
            "manifest": file_record(manifest_path, panel, "manifest"),
            "wdsat_anf": file_record(instance / "instance.anf", panel, "ANF export"),
            "cryptominisat_xor_dimacs": file_record(
                instance / "instance.xor.cnf", panel, "CNF-XOR export"
            ),
            "magma_boolean_f4": file_record(instance / "instance.magma", panel, "Magma export"),
            "source_binding": file_record(source_binding_path, panel, "source binding"),
        }
        for export_name in ("wdsat_anf", "cryptominisat_xor_dimacs", "magma_boolean_f4"):
            binding = source_binding.get("exports", {}).get(export_name)
            record = files[export_name]
            if not isinstance(binding, dict) or (
                binding.get("bytes") != record["bytes"]
                or binding.get("sha256") != record["sha256"]
            ):
                raise ExternalMagmaError(
                    f"{task_id(seed, cell)}: {export_name} disagrees with source binding"
                )
        source_id = source_binding.get("source_instance_sha256")
        producer_id = source_binding.get("producer_source_instance_blake3")
        if not isinstance(source_id, str) or not re.fullmatch(r"[0-9a-f]{64}", source_id):
            raise ExternalMagmaError(f"{task_id(seed, cell)}: invalid source SHA-256")
        if not isinstance(producer_id, str) or not re.fullmatch(r"[0-9a-f]{64}", producer_id):
            raise ExternalMagmaError(f"{task_id(seed, cell)}: invalid producer BLAKE3")
        if manifest.get("source_instance", {}).get("id_blake3") != producer_id:
            raise ExternalMagmaError(f"{task_id(seed, cell)}: producer identity changed")
        split_magma_source(instance / "instance.magma")
        tasks.append(
            {
                "id": task_id(seed, cell),
                "seed": seed,
                "cell": cell,
                "source_instance_sha256": source_id,
                "semantic_sha256": source_binding.get("semantic_sha256"),
                "producer_source_instance_blake3": producer_id,
                "source_variables": manifest.get("source_variables"),
                "files": files,
            }
        )
    expected_ids = [task_id(seed, cell) for seed, cell in STAGE13.tasks(protocol)]
    if len(tasks) != 20 or [task["id"] for task in tasks] != expected_ids:
        raise ExternalMagmaError("dry-run task inventory is not the frozen seed-major panel")
    body = {
        "schema": MANIFEST_SCHEMA,
        "archive_commit": STAGE13_ARCHIVE_COMMIT,
        "runner_repository_revision": git_revision(REPO),
        "panel_source_revision": custody["panel_source_revision"],
        "archive_root": str(panel),
        "protocol": {
            "path": str(protocol_path),
            "file_sha256": sha256_file(protocol_path),
            "canonical_sha256": canonical_sha256(protocol),
        },
        "custody": {
            "path": str(panel.parent / "stage-13-custody.json"),
            "sha256": sha256_file(panel.parent / "stage-13-custody.json"),
            "panel_summary_sha256": custody["panel_summary_sha256"],
            "corrected_summary_sha256": custody.get("corrected_summary_sha256"),
        },
        "meter": {
            "path": str(meter),
            "sha256": sha256_file(meter),
        },
        "contract_pins": {
            "stage13_verifier": {
                "path": str((STAGE / "verify_stage13_pdp_panel.py").resolve()),
                "sha256": sha256_file(STAGE / "verify_stage13_pdp_panel.py"),
                "expected_sha256": PINNED_STAGE13_VERIFIER_SHA256,
            },
            "matrix_contract": {
                "path": str((HERE / "run_koblitz_pdp_matrix.py").resolve()),
                "sha256": sha256_file(HERE / "run_koblitz_pdp_matrix.py"),
                "expected_sha256": PINNED_MATRIX_CONTRACT_SHA256,
            },
            "meter_expected_sha256": PINNED_METER_SHA256,
            "archive_tree_matches_pinned_commit": True,
        },
        "relevant_checkout": relevant_checkout_state(),
        "task_order": "seed-major",
        "expected_tasks": 20,
        "tasks": tasks,
        "claim_boundary": protocol["claim_boundary"],
        "archive_summary_claim": summary["claim"],
    }
    result = dict(body)
    result["manifest_sha256"] = canonical_sha256(body)
    return result


def resolve_manifest_path(panel: Path, record: dict, label: str) -> Path:
    relative = Path(str(record.get("path", "")))
    if relative.is_absolute() or ".." in relative.parts or not relative.parts:
        raise ExternalMagmaError(f"unsafe {label} path in manifest")
    path = panel / relative
    safe_relative(path, panel, label)
    return path


def verify_task_inputs(panel: Path, task: dict) -> dict:
    current = {}
    for name, frozen in task["files"].items():
        path = resolve_manifest_path(panel, frozen, f"{task['id']} {name}")
        observed = {
            "path": frozen["path"],
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        if observed != frozen:
            raise ExternalMagmaError(f"{task['id']}: frozen input changed: {name}")
        current[name] = observed
    binding = read_json(resolve_manifest_path(panel, task["files"]["source_binding"], "source binding"))
    if (
        binding.get("source_instance_sha256") != task["source_instance_sha256"]
        or binding.get("producer_source_instance_blake3")
        != task["producer_source_instance_blake3"]
    ):
        raise ExternalMagmaError(f"{task['id']}: source identity changed")
    return current


def split_magma_source(path: Path) -> str:
    text = path.read_text()
    if not text.endswith(F4_SUFFIX):
        raise ExternalMagmaError(f"{path}: Magma input is not the frozen direct-F4 template")
    prefix = text[: -len(F4_SUFFIX)]
    required_prefix = (
        "SetNthreads(1);\n",
        "SetGPU(false);\n",
        "SetSeed(1);\n",
        "R := BooleanPolynomialRing(",
    )
    if not all(item in prefix for item in required_prefix) or not prefix.endswith("I := ideal<R | F>;\n"):
        raise ExternalMagmaError(f"{path}: unsafe or incomplete Magma source prefix")
    return prefix


def render_witness_script(magma_path: Path, exclusions: list[list[bool]], source_variables: int) -> str:
    if not isinstance(source_variables, int) or source_variables <= 0:
        raise ExternalMagmaError("source_variables must be a positive integer")
    for index, assignment in enumerate(exclusions):
        if len(assignment) != source_variables or any(type(bit) is not bool for bit in assignment):
            raise ExternalMagmaError(f"invalid excluded assignment {index}")
    prefix = split_magma_source(magma_path)
    excluded_lines = []
    for assignment in exclusions:
        values = ", ".join(f"K!{int(bit)}" for bit in assignment)
        excluded_lines.append(f"  [ {values} ]")
    excluded_literal = "[\n" + ",\n".join(excluded_lines) + "\n]" if excluded_lines else "[]"
    sat_call = (
        "sat_model, S := SAT(G : Verbose := false);"
        if not exclusions
        else "sat_model, S := SAT(G : Exclude := Excluded, Verbose := false);"
    )
    model_body = f"""\
K := GF(2);
Excluded := {excluded_literal};
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
    return prefix + F4_BODY_WITHOUT_QUIT + model_body


def parse_model_terminal(output: str, source_variables: int, excluded_count: int) -> dict | None:
    markers: dict[str, str] = {}
    for raw in output.splitlines():
        line = raw.strip()
        if not line.startswith("KOBLITZ_MODEL_") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in markers:
            return None
        markers[key] = value.strip()
    required = {
        "KOBLITZ_MODEL_SCHEMA",
        "KOBLITZ_MODEL_STATUS",
        "KOBLITZ_MODEL_BITS",
        "KOBLITZ_MODEL_VARIABLES",
        "KOBLITZ_MODEL_EXCLUDED_COUNT",
        "KOBLITZ_MODEL_CPU_SECONDS",
        "KOBLITZ_MODEL_WALL_SECONDS",
    }
    if set(markers) != required or markers["KOBLITZ_MODEL_SCHEMA"] != MODEL_SCHEMA:
        return None
    status = markers["KOBLITZ_MODEL_STATUS"]
    if status not in {"SAT", "UNSAT"}:
        return None
    try:
        variables = int(markers["KOBLITZ_MODEL_VARIABLES"])
        excluded = int(markers["KOBLITZ_MODEL_EXCLUDED_COUNT"])
        cpu_seconds = float(markers["KOBLITZ_MODEL_CPU_SECONDS"])
        wall_seconds = float(markers["KOBLITZ_MODEL_WALL_SECONDS"])
    except ValueError:
        return None
    if (
        variables != source_variables
        or excluded != excluded_count
        or not math.isfinite(cpu_seconds)
        or not math.isfinite(wall_seconds)
        or cpu_seconds < 0
        or wall_seconds < 0
    ):
        return None
    bits_text = markers["KOBLITZ_MODEL_BITS"]
    if status == "SAT":
        if len(bits_text) != source_variables or set(bits_text) - {"0", "1"}:
            return None
        assignment = [bit == "1" for bit in bits_text]
    else:
        if bits_text != "-":
            return None
        assignment = None
    return {
        "schema": MODEL_SCHEMA,
        "status": status.lower(),
        "source_variables": variables,
        "excluded_count": excluded,
        "assignment": assignment,
        "cpu_seconds": cpu_seconds,
        "wall_seconds": wall_seconds,
    }


def finite_nonnegative(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and value >= 0
    )


def validate_meter_record(record: dict, expected_command: list[str], timeout: float) -> None:
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
    if not math.isclose(
        total,
        float(metrics["user_seconds"]) + float(metrics["system_seconds"]),
        rel_tol=1e-9,
        abs_tol=1e-9,
    ):
        raise ExternalMagmaError("total core-seconds do not equal user plus system")
    if not math.isclose(total, float(metrics["single_core_seconds"]), rel_tol=1e-9, abs_tol=1e-9):
        raise ExternalMagmaError("single-core and total core-seconds disagree")


def process_receipt_resource_complete(receipt: dict) -> bool:
    process = receipt.get("process", {})
    return (
        process.get("timed_out") is False
        and process.get("orphan_group_terminated") is False
    )


def run_metered(
    meter: Path,
    command: list[str],
    cwd: Path,
    timeout: float,
    attempt: Path,
    environment: dict[str, str] | None = None,
) -> dict:
    if attempt.exists():
        raise ExternalMagmaError(f"refusing to overwrite attempt directory {attempt}")
    attempt.mkdir(parents=True)
    stdout_path = attempt / "stdout"
    stderr_path = attempt / "stderr"
    metrics_path = attempt / "metrics.json"
    receipt_path = attempt / "receipt.json"
    meter_command = [
        str(Path(sys.executable).resolve()),
        str(meter.resolve()),
        "--cwd",
        str(cwd.resolve()),
        "--timeout",
        str(timeout),
        "--stdout",
        str(stdout_path.resolve()),
        "--stderr",
        str(stderr_path.resolve()),
        "--metrics",
        str(metrics_path.resolve()),
        "--",
        *command,
    ]
    launcher = subprocess.run(meter_command, cwd=REPO, env=environment, check=False)
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
            "path": str(stdout_path.resolve()),
            "bytes": stdout_path.stat().st_size,
            "sha256": sha256_file(stdout_path),
        },
        "stderr": {
            "path": str(stderr_path.resolve()),
            "bytes": stderr_path.stat().st_size,
            "sha256": sha256_file(stderr_path),
        },
    }
    atomic_json(receipt_path, receipt)
    return receipt


def attempt_path(root: Path, stem: str) -> Path:
    index = 1
    while (root / f"{stem}-{index:03d}").exists():
        index += 1
    return root / f"{stem}-{index:03d}"


def launcher_option(command: list, option: str) -> str | None:
    positions = [index for index, value in enumerate(command) if value == option]
    if len(positions) != 1 or positions[0] + 1 >= len(command):
        return None
    return str(command[positions[0] + 1])


def collect_interrupted_processes(
    task_root: Path,
    meter: Path,
    magma: dict,
    backend: dict,
    magma_input: Path,
    manifest_path: Path,
) -> tuple[list[dict], list[str]]:
    """Recover and charge complete attempt receipts left before task finalization."""
    recovered = []
    receipt_paths = sorted(task_root.rglob("receipt.json")) if task_root.exists() else []
    for receipt_path in receipt_paths:
        receipt = read_json(receipt_path)
        process = receipt.get("process")
        launcher = receipt.get("meter_launcher_command")
        if not isinstance(process, dict) or not isinstance(launcher, list):
            raise ExternalMagmaError(f"malformed interrupted process receipt {receipt_path}")
        command = process.get("command")
        timeout = process.get("watchdog_seconds")
        cwd_text = launcher_option(launcher, "--cwd")
        if not isinstance(command, list) or not finite_nonnegative(timeout) or cwd_text is None:
            raise ExternalMagmaError(f"incomplete interrupted process receipt {receipt_path}")
        cwd = Path(cwd_text)
        validate_saved_receipt(receipt, command, float(timeout), meter, cwd, task_root)
        if command and command[0] == magma["path"] and len(command) == 5 and command[1:4] == ["-t", "1", "-b"]:
            script = Path(command[4])
            if script.resolve() == magma_input.resolve():
                kind = "primary_f4"
            else:
                try:
                    script.resolve(strict=True).relative_to(task_root.resolve(strict=True))
                except (OSError, ValueError) as error:
                    raise ExternalMagmaError(
                        f"interrupted witness script escapes task directory: {script}"
                    ) from error
                if script.name != "witness.magma":
                    raise ExternalMagmaError(f"unexpected interrupted Magma script {script}")
                kind = "f4_plus_sat_g"
        elif command and command[0] == backend["path"] and len(command) == 4 and command[1] == "validate-model":
            if Path(command[2]).resolve() != manifest_path.resolve():
                raise ExternalMagmaError("interrupted validator used a different manifest")
            try:
                Path(command[3]).resolve(strict=True).relative_to(task_root.resolve(strict=True))
            except (OSError, ValueError) as error:
                raise ExternalMagmaError(
                    "interrupted validator assignment escapes task directory"
                ) from error
            kind = "validator"
        else:
            raise ExternalMagmaError(
                f"unrecognized interrupted command in {receipt_path}: {command}"
            )
        recovered_row = {
            "kind": kind,
            "receipt_path": str(receipt_path.resolve()),
            "receipt_sha256": sha256_file(receipt_path),
            "receipt": receipt,
        }
        if kind in {"primary_f4", "f4_plus_sat_g"}:
            prior_status, prior_terminal = classify_f4(receipt)
            recovered_row["f4_status"] = prior_status
            recovered_row["f4_terminal"] = prior_terminal
        recovered.append(recovered_row)
    attempt_directories = sorted(
        path
        for pattern in ("f4-attempt-*", "model-attempt-*")
        for path in task_root.glob(pattern)
        if path.is_dir()
    )
    resolved_task_root = task_root.resolve(strict=True)
    covered_roots = {
        resolved_task_root
        / Path(row["receipt_path"]).resolve(strict=True).relative_to(resolved_task_root).parts[0]
        for row in recovered
    }
    incomplete = [
        str(path.resolve()) for path in attempt_directories if path.resolve() not in covered_roots
    ]
    for attempt in attempt_directories:
        if not attempt.name.startswith("model-attempt-"):
            continue
        witness_receipt_path = attempt / "magma-process" / "receipt.json"
        if not witness_receipt_path.is_file():
            continue
        witness_output = receipt_stdout(read_json(witness_receipt_path))
        if "KOBLITZ_MODEL_STATUS=SAT" not in witness_output:
            continue
        # A crash may occur before or after the validator.  This runner charges
        # the recovered processes but does not silently substitute a rerun for
        # the abandoned model-level verdict.
        resolved = str(attempt.resolve())
        if resolved not in incomplete:
            incomplete.append(resolved)
    return recovered, incomplete


def receipt_stdout(receipt: dict) -> str:
    path = Path(receipt["stdout"]["path"])
    data = path.read_bytes()
    if len(data) != receipt["stdout"]["bytes"] or sha256_bytes(data) != receipt["stdout"]["sha256"]:
        raise ExternalMagmaError(f"metered stdout changed: {path}")
    return data.decode(errors="replace")


def receipt_stderr(receipt: dict) -> str:
    path = Path(receipt["stderr"]["path"])
    data = path.read_bytes()
    if len(data) != receipt["stderr"]["bytes"] or sha256_bytes(data) != receipt["stderr"]["sha256"]:
        raise ExternalMagmaError(f"metered stderr changed: {path}")
    return data.decode(errors="replace")


def binary_identity(value: str | Path, version_options: list[list[str]]) -> dict:
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
    failed_version_probes = []
    for suffix in version_options:
        command = [str(path), *suffix]
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=15, check=False)
        except (OSError, subprocess.SubprocessError) as error:
            failed_version_probes.append({"command": command, "error": str(error)})
            continue
        combined = (result.stdout + result.stderr).strip()
        if combined and result.returncode == 0:
            version = combined[:4000]
            version_command = command
            version_returncode = result.returncode
            break
        failed_version_probes.append(
            {
                "command": command,
                "returncode": result.returncode,
                "output": combined[:1000],
            }
        )
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "version": version or "unreported",
        "version_command": version_command,
        "version_returncode": version_returncode,
        "failed_version_probes": failed_version_probes,
    }


def recognized_magma_version(output: str) -> bool:
    return bool(re.search(r"\bMagma V2\.\d+-\d+\b", output)) or any(
        re.fullmatch(r"2\.\d+-\d+", line.strip()) is not None
        for line in output.splitlines()
    )


def cpu_model() -> str:
    if platform.system() == "Darwin":
        result = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        match = re.search(r"^model name\s*:\s*(.+)$", cpuinfo.read_text(), re.MULTILINE)
        if match:
            return match.group(1).strip()
    return platform.processor() or "unreported"


def host_identity() -> dict:
    return {
        "node": platform.node(),
        "platform": platform.platform(),
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "cpu_model": cpu_model(),
        "logical_cpus": os.cpu_count(),
        "python": platform.python_version(),
        "python_executable": str(Path(sys.executable).resolve()),
    }


def git_revision(path: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=path,
            capture_output=True,
            text=True,
            check=True,
            timeout=15,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise ExternalMagmaError(f"cannot resolve repository revision at {path}: {error}") from error
    revision = result.stdout.strip()
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ExternalMagmaError(f"invalid repository revision at {path}: {revision!r}")
    return revision


def relevant_checkout_state() -> dict:
    result = subprocess.run(
        ["git", "status", "--porcelain=v1", "--", *RELEVANT_RUNTIME_PATHS],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
    )
    porcelain = result.stdout.splitlines()
    return {
        "commit": git_revision(REPO),
        "dirty": bool(porcelain),
        "porcelain": porcelain,
        "paths": list(RELEVANT_RUNTIME_PATHS),
    }


def archive_contracts_match() -> bool:
    result = subprocess.run(
        ["git", "diff", "--quiet", STAGE13_ARCHIVE_COMMIT, "--", *ARCHIVE_PINNED_PATHS],
        cwd=REPO,
        check=False,
    )
    if result.returncode not in {0, 1}:
        raise ExternalMagmaError("cannot compare current archive contracts with pinned commit")
    return result.returncode == 0


def classify_f4(receipt: dict) -> tuple[str, dict | None]:
    process = receipt["process"]
    stdout = receipt_stdout(receipt)
    stderr = receipt_stderr(receipt)
    if process["timed_out"]:
        return "timeout_inconclusive", None
    terminal = MATRIX.parse_magma_terminal(stdout)
    if process["returncode"] != 0 or "Runtime error" in stdout + stderr or "User error" in stdout + stderr:
        return "solver_error", terminal
    if terminal is None:
        return "terminal_certificate_missing", None
    if terminal["terminal_status"] == "unsat":
        return "unsat", terminal
    return "sat_basis_certificate_unverified_model", terminal


def same_f4_terminal(primary: dict, repeated: dict) -> bool:
    keys = (
        "schema",
        "algorithm",
        "terminal_status",
        "f4_step_degrees",
        "basis_size",
        "single_thread_requested",
        "gpu_disabled",
    )
    return all(primary.get(key) == repeated.get(key) for key in keys)


def validation_status(receipt: dict, manifest: dict, assignment_path: Path) -> tuple[str, dict | None]:
    process = receipt["process"]
    run = {
        "stdout": receipt_stdout(receipt),
        "stderr": receipt_stderr(receipt),
        "returncode": process["returncode"],
        "timed_out": process["timed_out"],
        "metrics": process["metrics"],
        "command": process["command"],
    }
    parsed = MATRIX.assignment_validation_status(run, manifest, assignment_path)
    return parsed["status"], parsed.get("report")


def execute_task(
    packet: dict,
    task: dict,
    output: Path,
    meter: Path,
    magma: dict,
    backend: dict,
    environment: dict[str, str],
    policy: dict,
) -> dict:
    panel = Path(packet["archive_root"])
    before = verify_task_inputs(panel, task)
    task_root = output / "tasks" / task["id"]
    task_root.mkdir(parents=True, exist_ok=True)
    manifest_path = resolve_manifest_path(panel, task["files"]["manifest"], "manifest")
    magma_path = resolve_manifest_path(panel, task["files"]["magma_boolean_f4"], "Magma input")
    manifest = read_json(manifest_path)
    interrupted_processes, incomplete_prior_attempts = collect_interrupted_processes(
        task_root,
        meter,
        magma,
        backend,
        magma_path,
        manifest_path,
    )
    task_progress_path = task_root / "task-progress.json"

    f4_command = [magma["path"], "-t", "1", "-b", str(magma_path)]
    f4_receipt = run_metered(
        meter,
        f4_command,
        task_root,
        policy["f4_watchdog_seconds"],
        attempt_path(task_root, "f4-attempt"),
        environment,
    )
    f4_status, f4_terminal = classify_f4(f4_receipt)
    result = {
        "schema": TASK_SCHEMA,
        "id": task["id"],
        "seed": task["seed"],
        "cell": task["cell"],
        "source_instance_sha256": task["source_instance_sha256"],
        "producer_source_instance_blake3": task["producer_source_instance_blake3"],
        "inputs_before": before,
        "f4": {
            "status": f4_status,
            "terminal": f4_terminal,
            "receipt": f4_receipt,
            "resource_complete": (
                f4_status in {"unsat", "sat_basis_certificate_unverified_model"}
                and process_receipt_resource_complete(f4_receipt)
            ),
        },
        "model_attempts": [],
        "interrupted_processes": interrupted_processes,
        "incomplete_prior_attempts": incomplete_prior_attempts,
        "contradictions": (
            ["an interrupted prior attempt has no complete process receipt"]
            if incomplete_prior_attempts
            else []
        ),
        "started_at": now(),
    }
    if any(
        row["receipt"]["process"].get("orphan_group_terminated") is True
        for row in interrupted_processes
    ):
        result["contradictions"].append(
            "an interrupted process left a descendant group requiring forced termination"
        )
    if any(row["kind"] == "validator" for row in interrupted_processes):
        result["contradictions"].append(
            "an interrupted validator receipt is charged but not admitted as replacement evidence"
        )
    for prior in interrupted_processes:
        if prior["kind"] not in {"primary_f4", "f4_plus_sat_g"}:
            continue
        prior_status = prior.get("f4_status")
        prior_terminal = prior.get("f4_terminal")
        if prior_status == "unsat":
            result["contradictions"].append(
                "an interrupted clean Magma receipt reported UNSAT for a frozen planted instance"
            )
        elif prior_status == "sat_basis_certificate_unverified_model" and (
            f4_terminal is None
            or prior_terminal is None
            or not same_f4_terminal(f4_terminal, prior_terminal)
        ):
            result["contradictions"].append(
                "an interrupted clean Magma terminal disagrees with the rerun primary F4 terminal"
            )
    atomic_json(task_progress_path, result)

    if f4_receipt["process"].get("orphan_group_terminated") is True:
        result["final_status"] = "f4_orphan_group_inconclusive"
    elif f4_status != "sat_basis_certificate_unverified_model":
        result["final_status"] = f4_status
        if f4_status == "unsat":
            result["contradictions"].append(
                "frozen task is a planted SAT instance but Magma reported UNSAT"
            )
            result["final_status"] = "contradictory_unsat"
    else:
        exclusions: list[list[bool]] = []
        for ordinal in range(1, policy["model_cap"] + 1):
            attempt = attempt_path(task_root, "model-attempt")
            attempt.mkdir(parents=True)
            witness_path = attempt / "witness.magma"
            witness_text = render_witness_script(
                magma_path,
                exclusions,
                int(task["source_variables"]),
            )
            witness_path.write_text(witness_text)
            witness_hash = sha256_bytes(witness_text.encode())
            witness_command = [magma["path"], "-t", "1", "-b", str(witness_path.resolve())]
            witness_receipt = run_metered(
                meter,
                witness_command,
                attempt,
                policy["model_watchdog_seconds"],
                attempt / "magma-process",
                environment,
            )
            repeated_status, repeated_terminal = classify_f4(witness_receipt)
            model_terminal = parse_model_terminal(
                receipt_stdout(witness_receipt),
                int(task["source_variables"]),
                len(exclusions),
            )
            model_row = {
                "ordinal": ordinal,
                "excluded_count": len(exclusions),
                "witness_script": {
                    "path": str(witness_path.resolve()),
                    "bytes": len(witness_text.encode()),
                    "sha256": witness_hash,
                    "source_magma_sha256": task["files"]["magma_boolean_f4"]["sha256"],
                    "recomputes_direct_sparse_f4": True,
                    "witness_process_includes_f4_recomputation": True,
                },
                "f4_status": repeated_status,
                "f4_terminal": repeated_terminal,
                "model_terminal": model_terminal,
                "witness_receipt": witness_receipt,
            }
            result["model_attempts"].append(model_row)
            atomic_json(task_progress_path, result)
            if repeated_terminal is None or not same_f4_terminal(f4_terminal, repeated_terminal):
                result["contradictions"].append(
                    f"model attempt {ordinal} did not reproduce the measured F4 terminal"
                )
                result["final_status"] = "witness_f4_contradiction"
                break
            if witness_receipt["process"]["timed_out"]:
                result["final_status"] = "witness_timeout_inconclusive"
                break
            if witness_receipt["process"].get("orphan_group_terminated") is True:
                result["final_status"] = "witness_orphan_group_inconclusive"
                break
            if repeated_status != "sat_basis_certificate_unverified_model":
                result["final_status"] = "witness_solver_error_inconclusive"
                break
            if model_terminal is None:
                result["final_status"] = "witness_terminal_missing_inconclusive"
                break
            if model_terminal["status"] == "unsat":
                if not exclusions:
                    result["contradictions"].append(
                        "proper F4 basis but SAT(G) returned UNSAT before any exclusion"
                    )
                    result["final_status"] = "witness_search_contradiction"
                else:
                    result["contradictions"].append(
                        "SAT(G) exhausted after only nonlifting models for a frozen planted PDP instance"
                    )
                    result["final_status"] = "nonlifting_models_exhausted_contradiction"
                break

            assignment = model_terminal["assignment"]
            if assignment in exclusions:
                result["contradictions"].append(
                    f"SAT(G) returned an assignment already excluded at model attempt {ordinal}"
                )
                result["final_status"] = "duplicate_excluded_model"
                break
            assignment_path = attempt / "assignment.json"
            assignment_bytes = (json.dumps(assignment, separators=(",", ":")) + "\n").encode()
            assignment_path.write_bytes(assignment_bytes)
            validation_command = [
                backend["path"],
                "validate-model",
                str(manifest_path),
                str(assignment_path.resolve()),
            ]
            validation_receipt = run_metered(
                meter,
                validation_command,
                attempt,
                policy["validation_watchdog_seconds"],
                attempt / "validation-process",
            )
            checked_status, report = validation_status(validation_receipt, manifest, assignment_path)
            if validation_receipt["process"].get("orphan_group_terminated") is True:
                checked_status = "validation_orphan_group_inconclusive"
            model_row["assignment"] = {
                "path": str(assignment_path.resolve()),
                "bytes": len(assignment_bytes),
                "sha256": sha256_bytes(assignment_bytes),
            }
            model_row["validation"] = {
                "status": checked_status,
                "report": report,
                "receipt": validation_receipt,
            }
            atomic_json(task_progress_path, result)
            if checked_status == "valid_point_witness":
                result["final_status"] = "sat"
                result["valid_point_witness_ordinal"] = ordinal
                break
            if checked_status == "nonlifting_source_model":
                exclusions.append(assignment)
                continue
            if checked_status == "timeout_inconclusive":
                result["final_status"] = "validation_timeout_inconclusive"
                break
            if checked_status == "validation_orphan_group_inconclusive":
                result["final_status"] = "validation_orphan_group_inconclusive"
                break
            result["contradictions"].append(
                f"SAT(G) assignment failed the source-bound validator: {checked_status}"
            )
            result["final_status"] = "invalid_model_evidence"
            break
        else:
            result["final_status"] = "model_cap_inconclusive"

    after = verify_task_inputs(panel, task)
    result["inputs_after"] = after
    result["inputs_unchanged"] = before == after
    if not result["inputs_unchanged"]:
        result["contradictions"].append("frozen source inputs changed during execution")
        result["final_status"] = "source_custody_failure"
    result["resource_complete"] = (
        result["f4"]["resource_complete"]
        and all(
            process_receipt_resource_complete(row["witness_receipt"])
            and (
                "validation" not in row
                or process_receipt_resource_complete(row["validation"]["receipt"])
            )
            for row in result["model_attempts"]
        )
    )
    result["finished_at"] = now()
    atomic_json(task_progress_path, result)
    atomic_json(task_root / "task-result.json", result)
    return result


def validate_saved_receipt(
    receipt: dict,
    expected_command: list[str],
    timeout: float,
    meter: Path,
    cwd: Path,
    allowed_root: Path,
) -> tuple[str, str]:
    process = receipt.get("process")
    if not isinstance(process, dict):
        raise ExternalMagmaError("saved process receipt is missing its process record")
    validate_meter_record(process, expected_command, timeout)
    stdout_record = receipt.get("stdout")
    stderr_record = receipt.get("stderr")
    if not isinstance(stdout_record, dict) or not isinstance(stderr_record, dict):
        raise ExternalMagmaError("saved process receipt lacks raw-output custody")
    stdout_path = Path(str(stdout_record.get("path", "")))
    stderr_path = Path(str(stderr_record.get("path", "")))
    if stdout_path.name != "stdout" or stderr_path.name != "stderr":
        raise ExternalMagmaError("saved process receipt uses unexpected raw-output names")
    if stdout_path.parent != stderr_path.parent:
        raise ExternalMagmaError("saved stdout and stderr do not share an attempt directory")
    attempt = stdout_path.parent
    try:
        attempt.resolve(strict=True).relative_to(allowed_root.resolve(strict=True))
    except (OSError, ValueError) as error:
        raise ExternalMagmaError("saved process attempt escapes its task directory") from error
    for artifact in (stdout_path, stderr_path, attempt / "metrics.json", attempt / "receipt.json"):
        if artifact.is_symlink() or not artifact.is_file():
            raise ExternalMagmaError(f"saved process artifact is missing or a symlink: {artifact}")
    metrics = read_json(attempt / "metrics.json")
    if metrics != process:
        raise ExternalMagmaError("embedded process metrics differ from metrics.json")
    on_disk_receipt = read_json(attempt / "receipt.json")
    if on_disk_receipt != receipt:
        raise ExternalMagmaError("embedded process receipt differs from receipt.json")
    expected_launcher = [
        str(Path(sys.executable).resolve()),
        str(meter.resolve()),
        "--cwd",
        str(cwd.resolve()),
        "--timeout",
        str(timeout),
        "--stdout",
        str(stdout_path.resolve()),
        "--stderr",
        str(stderr_path.resolve()),
        "--metrics",
        str((attempt / "metrics.json").resolve()),
        "--",
        *expected_command,
    ]
    if receipt.get("meter_launcher_command") != expected_launcher:
        raise ExternalMagmaError("saved meter launcher command changed")
    if receipt.get("meter_launcher_returncode") != 0:
        raise ExternalMagmaError("saved meter launcher did not complete")
    return receipt_stdout(receipt), receipt_stderr(receipt)


def validate_completed_task(
    path: Path,
    expected: dict,
    packet: dict,
    tools: dict,
    policy: dict,
) -> dict:
    result = read_json(path)
    task_root = path.parent
    panel = Path(packet["archive_root"])
    current_inputs = verify_task_inputs(panel, expected)
    if (
        result.get("schema") != TASK_SCHEMA
        or result.get("id") != expected["id"]
        or result.get("seed") != expected["seed"]
        or result.get("cell") != expected["cell"]
        or result.get("source_instance_sha256") != expected["source_instance_sha256"]
        or result.get("producer_source_instance_blake3")
        != expected["producer_source_instance_blake3"]
        or result.get("inputs_before") != current_inputs
        or result.get("inputs_after") != current_inputs
        or result.get("inputs_unchanged") is not True
    ):
        raise ExternalMagmaError(f"invalid completed task identity or custody in {path}")
    meter = Path(tools["meter"]["path"])
    manifest_path = resolve_manifest_path(panel, expected["files"]["manifest"], "manifest")
    magma_path = resolve_manifest_path(panel, expected["files"]["magma_boolean_f4"], "Magma input")
    manifest = read_json(manifest_path)
    interrupted = result.get("interrupted_processes", [])
    if not isinstance(interrupted, list):
        raise ExternalMagmaError(f"{path}: interrupted-process inventory is invalid")
    for prior in interrupted:
        if not isinstance(prior, dict) or prior.get("kind") not in {
            "primary_f4",
            "f4_plus_sat_g",
            "validator",
        }:
            raise ExternalMagmaError(f"{path}: interrupted-process classification is invalid")
        receipt_path = Path(str(prior.get("receipt_path", "")))
        try:
            receipt_path.resolve(strict=True).relative_to(task_root.resolve(strict=True))
        except (OSError, ValueError) as error:
            raise ExternalMagmaError(f"{path}: interrupted receipt escapes task directory") from error
        if (
            receipt_path.is_symlink()
            or not receipt_path.is_file()
            or prior.get("receipt_sha256") != sha256_file(receipt_path)
            or prior.get("receipt") != read_json(receipt_path)
        ):
            raise ExternalMagmaError(f"{path}: interrupted receipt custody changed")
        prior_receipt = prior["receipt"]
        prior_process = prior_receipt.get("process", {})
        prior_command = prior_process.get("command")
        prior_timeout = prior_process.get("watchdog_seconds")
        prior_launcher = prior_receipt.get("meter_launcher_command", [])
        prior_cwd = launcher_option(prior_launcher, "--cwd")
        if (
            not isinstance(prior_command, list)
            or not finite_nonnegative(prior_timeout)
            or prior_cwd is None
        ):
            raise ExternalMagmaError(f"{path}: interrupted receipt is incomplete")
        validate_saved_receipt(
            prior_receipt,
            prior_command,
            float(prior_timeout),
            meter,
            Path(prior_cwd),
            task_root,
        )
        if prior["kind"] in {"primary_f4", "f4_plus_sat_g"}:
            if not (
                prior_command
                and prior_command[0] == tools["magma"]["path"]
                and len(prior_command) == 5
                and prior_command[1:4] == ["-t", "1", "-b"]
            ):
                raise ExternalMagmaError(f"{path}: interrupted Magma command changed")
            prior_status, prior_terminal = classify_f4(prior_receipt)
            if (
                prior.get("f4_status") != prior_status
                or prior.get("f4_terminal") != prior_terminal
            ):
                raise ExternalMagmaError(f"{path}: interrupted Magma interpretation changed")
        elif not (
            prior_command
            and prior_command[0] == tools["backend"]["path"]
            and len(prior_command) == 4
            and prior_command[1] == "validate-model"
        ):
            raise ExternalMagmaError(f"{path}: interrupted validator command changed")
    incomplete_prior = result.get("incomplete_prior_attempts", [])
    if not isinstance(incomplete_prior, list):
        raise ExternalMagmaError(f"{path}: incomplete-prior-attempt inventory is invalid")
    for item in incomplete_prior:
        prior_path = Path(str(item))
        try:
            prior_path.resolve(strict=True).relative_to(task_root.resolve(strict=True))
        except (OSError, ValueError) as error:
            raise ExternalMagmaError(f"{path}: incomplete prior attempt path changed") from error

    f4 = result.get("f4")
    if not isinstance(f4, dict) or not isinstance(f4.get("receipt"), dict):
        raise ExternalMagmaError(f"{path}: missing primary F4 receipt")
    f4_command = [tools["magma"]["path"], "-t", "1", "-b", str(magma_path)]
    validate_saved_receipt(
        f4["receipt"],
        f4_command,
        policy["f4_watchdog_seconds"],
        meter,
        task_root,
        task_root,
    )
    observed_f4_status, observed_f4_terminal = classify_f4(f4["receipt"])
    if f4.get("status") != observed_f4_status or f4.get("terminal") != observed_f4_terminal:
        raise ExternalMagmaError(f"{path}: stored F4 interpretation changed")
    if f4.get("resource_complete") is not (
        observed_f4_status in {"unsat", "sat_basis_certificate_unverified_model"}
        and process_receipt_resource_complete(f4["receipt"])
    ):
        raise ExternalMagmaError(f"{path}: stored F4 resource classification changed")
    expected_prior_contradictions = []
    if incomplete_prior:
        expected_prior_contradictions.append(
            "an interrupted prior attempt has no complete process receipt"
        )
    if any(
        prior["receipt"]["process"].get("orphan_group_terminated") is True
        for prior in interrupted
    ):
        expected_prior_contradictions.append(
            "an interrupted process left a descendant group requiring forced termination"
        )
    if any(prior["kind"] == "validator" for prior in interrupted):
        expected_prior_contradictions.append(
            "an interrupted validator receipt is charged but not admitted as replacement evidence"
        )
    for prior in interrupted:
        if prior["kind"] not in {"primary_f4", "f4_plus_sat_g"}:
            continue
        if prior.get("f4_status") == "unsat":
            expected_prior_contradictions.append(
                "an interrupted clean Magma receipt reported UNSAT for a frozen planted instance"
            )
        elif prior.get("f4_status") == "sat_basis_certificate_unverified_model" and (
            observed_f4_terminal is None
            or prior.get("f4_terminal") is None
            or not same_f4_terminal(observed_f4_terminal, prior["f4_terminal"])
        ):
            expected_prior_contradictions.append(
                "an interrupted clean Magma terminal disagrees with the rerun primary F4 terminal"
            )
    for contradiction in expected_prior_contradictions:
        if contradiction not in result.get("contradictions", []):
            raise ExternalMagmaError(f"{path}: interrupted Magma contradiction was hidden")

    exclusions: list[list[bool]] = []
    valid_witness_ordinals: list[int] = []
    model_attempts = result.get("model_attempts")
    if not isinstance(model_attempts, list) or len(model_attempts) > policy["model_cap"]:
        raise ExternalMagmaError(f"{path}: invalid model-attempt inventory")
    for expected_ordinal, row in enumerate(model_attempts, 1):
        if not isinstance(row, dict) or row.get("ordinal") != expected_ordinal:
            raise ExternalMagmaError(f"{path}: nonconsecutive model attempt")
        if row.get("excluded_count") != len(exclusions):
            raise ExternalMagmaError(f"{path}: model exclusion count changed")
        script_record = row.get("witness_script")
        if not isinstance(script_record, dict):
            raise ExternalMagmaError(f"{path}: missing witness-script record")
        script_path = Path(str(script_record.get("path", "")))
        try:
            script_path.resolve(strict=True).relative_to(task_root.resolve(strict=True))
        except (OSError, ValueError) as error:
            raise ExternalMagmaError(f"{path}: witness script escapes task directory") from error
        if script_path.is_symlink() or not script_path.is_file():
            raise ExternalMagmaError(f"{path}: witness script is missing or a symlink")
        expected_script = render_witness_script(
            magma_path,
            exclusions,
            int(expected["source_variables"]),
        ).encode()
        actual_script = script_path.read_bytes()
        if (
            actual_script != expected_script
            or script_record.get("bytes") != len(actual_script)
            or script_record.get("sha256") != sha256_bytes(actual_script)
            or script_record.get("source_magma_sha256")
            != expected["files"]["magma_boolean_f4"]["sha256"]
            or script_record.get("recomputes_direct_sparse_f4") is not True
            or script_record.get("witness_process_includes_f4_recomputation") is not True
        ):
            raise ExternalMagmaError(f"{path}: witness script or its derivation record changed")
        witness_receipt = row.get("witness_receipt")
        if not isinstance(witness_receipt, dict):
            raise ExternalMagmaError(f"{path}: missing witness process receipt")
        witness_command = [tools["magma"]["path"], "-t", "1", "-b", str(script_path.resolve())]
        validate_saved_receipt(
            witness_receipt,
            witness_command,
            policy["model_watchdog_seconds"],
            meter,
            script_path.parent,
            task_root,
        )
        repeated_status, repeated_terminal = classify_f4(witness_receipt)
        parsed_model = parse_model_terminal(
            receipt_stdout(witness_receipt),
            int(expected["source_variables"]),
            len(exclusions),
        )
        if (
            row.get("f4_status") != repeated_status
            or row.get("f4_terminal") != repeated_terminal
            or row.get("model_terminal") != parsed_model
            or observed_f4_terminal is None
            or repeated_terminal is None
            or not same_f4_terminal(observed_f4_terminal, repeated_terminal)
        ):
            raise ExternalMagmaError(f"{path}: repeated F4/model interpretation changed")
        if parsed_model is None or parsed_model["status"] != "sat":
            continue
        assignment = parsed_model["assignment"]
        if assignment in exclusions:
            raise ExternalMagmaError(f"{path}: model search repeated an excluded assignment")
        assignment_record = row.get("assignment")
        validation = row.get("validation")
        if not isinstance(assignment_record, dict) or not isinstance(validation, dict):
            raise ExternalMagmaError(f"{path}: SAT model lacks assignment validation")
        assignment_path = Path(str(assignment_record.get("path", "")))
        if assignment_path != script_path.parent / "assignment.json":
            raise ExternalMagmaError(f"{path}: assignment path changed")
        if assignment_path.is_symlink() or not assignment_path.is_file():
            raise ExternalMagmaError(f"{path}: assignment file is missing or a symlink")
        assignment_bytes = assignment_path.read_bytes()
        try:
            decoded_assignment = json.loads(assignment_bytes)
        except json.JSONDecodeError as error:
            raise ExternalMagmaError(f"{path}: assignment JSON is invalid") from error
        if (
            decoded_assignment != assignment
            or assignment_record.get("bytes") != len(assignment_bytes)
            or assignment_record.get("sha256") != sha256_bytes(assignment_bytes)
        ):
            raise ExternalMagmaError(f"{path}: assignment custody changed")
        validation_receipt = validation.get("receipt")
        if not isinstance(validation_receipt, dict):
            raise ExternalMagmaError(f"{path}: validation process receipt is missing")
        validation_command = [
            tools["backend"]["path"],
            "validate-model",
            str(manifest_path),
            str(assignment_path.resolve()),
        ]
        validate_saved_receipt(
            validation_receipt,
            validation_command,
            policy["validation_watchdog_seconds"],
            meter,
            script_path.parent,
            task_root,
        )
        checked_status, checked_report = validation_status(
            validation_receipt, manifest, assignment_path
        )
        if validation_receipt["process"].get("orphan_group_terminated") is True:
            checked_status = "validation_orphan_group_inconclusive"
        if validation.get("status") != checked_status or validation.get("report") != checked_report:
            raise ExternalMagmaError(f"{path}: validator interpretation changed")
        if checked_status == "valid_point_witness":
            valid_witness_ordinals.append(expected_ordinal)
        elif checked_status == "nonlifting_source_model":
            exclusions.append(assignment)

    final_status = result.get("final_status")
    if final_status == "sat":
        if (
            observed_f4_status != "sat_basis_certificate_unverified_model"
            or len(valid_witness_ordinals) != 1
            or result.get("valid_point_witness_ordinal") != valid_witness_ordinals[0]
            or result.get("contradictions") != []
        ):
            raise ExternalMagmaError(f"{path}: SAT leaf lacks its unique valid point witness")
    if final_status == "model_cap_inconclusive" and len(exclusions) != policy["model_cap"]:
        raise ExternalMagmaError(f"{path}: model-cap leaf did not charge every exclusion")
    observed_resource_complete = (
        f4.get("resource_complete") is True
        and all(
            process_receipt_resource_complete(row["witness_receipt"])
            and (
                "validation" not in row
                or process_receipt_resource_complete(row["validation"]["receipt"])
            )
            for row in model_attempts
        )
    )
    if result.get("resource_complete") is not observed_resource_complete:
        raise ExternalMagmaError(f"{path}: task resource-completeness flag changed")
    return result


def process_metrics(receipt: dict) -> dict:
    metrics = receipt["process"]["metrics"]
    return {
        "wall_seconds": float(metrics["wall_seconds"]),
        "user_seconds": float(metrics["user_seconds"]),
        "system_seconds": float(metrics["system_seconds"]),
        "total_core_seconds": float(metrics["total_core_seconds"]),
        "single_core_seconds": float(metrics["single_core_seconds"]),
        "peak_rss_bytes": int(metrics["peak_rss_bytes"]),
    }


def charged_metrics(rows: list[dict]) -> dict:
    return {
        "processes": len(rows),
        "total_core_seconds": sum(row["total_core_seconds"] for row in rows),
        "user_seconds": sum(row["user_seconds"] for row in rows),
        "system_seconds": sum(row["system_seconds"] for row in rows),
        "sequential_wall_seconds": sum(row["wall_seconds"] for row in rows),
        "peak_rss_bytes": max((row["peak_rss_bytes"] for row in rows), default=0),
        "peak_rss_scope": "maximum platform ru_maxrss among separately metered child processes",
    }


def distribution(values: list[float | int]) -> dict | None:
    if not values:
        return None
    numeric = [float(value) for value in values]
    return {
        "count": len(numeric),
        "min": min(numeric),
        "median": statistics.median(numeric),
        "max": max(numeric),
        "sum": sum(numeric),
    }


def phase_distribution(rows: list[dict]) -> dict:
    return {
        "processes": len(rows),
        "total_core_seconds": distribution([row["total_core_seconds"] for row in rows]),
        "wall_seconds": distribution([row["wall_seconds"] for row in rows]),
        "peak_rss_bytes": distribution([row["peak_rss_bytes"] for row in rows]),
    }


def summarize_run(
    packet: dict,
    selection: list[dict],
    results: list[dict],
    policy: dict,
) -> dict:
    expected_ids = {task["id"] for task in packet["tasks"]}
    selected_ids = {task["id"] for task in selection}
    result_ids = {result["id"] for result in results}
    full_selection = selected_ids == expected_ids and len(selection) == 20
    all_attempted = result_ids == selected_ids and len(results) == len(selection)
    all_f4_terminal = all(
        result.get("f4", {}).get("status")
        in {"sat_basis_certificate_unverified_model", "unsat"}
        for result in results
    )
    all_resources = all(result.get("resource_complete") is True for result in results)
    all_witnesses = all(result.get("final_status") == "sat" for result in results)
    no_contradictions = all(not result.get("contradictions") for result in results)
    synthetic_test_mode = policy.get("synthetic_test_mode") is True
    archive_contracts_pinned = policy.get("archive_contracts_pinned") is True
    relevant_checkout_clean = policy.get("relevant_checkout_clean") is True
    status_counts: dict[str, int] = {}
    primary_f4_metrics: list[dict] = []
    witness_metrics: list[dict] = []
    validator_metrics: list[dict] = []
    interrupted_metrics: list[dict] = []
    by_cell: dict[str, dict[str, list[dict]]] = {}
    for result in results:
        status = str(result.get("final_status", "missing"))
        status_counts[status] = status_counts.get(status, 0) + 1
        cell_id = result["cell"]["id"]
        cell_phases = by_cell.setdefault(
            cell_id,
            {
                "primary_f4": [],
                "f4_plus_sat_g": [],
                "validator": [],
                "interrupted_prior": [],
            },
        )
        for prior in result.get("interrupted_processes", []):
            recovered = process_metrics(prior["receipt"])
            interrupted_metrics.append(recovered)
            cell_phases["interrupted_prior"].append(recovered)
        primary = process_metrics(result["f4"]["receipt"])
        primary_f4_metrics.append(primary)
        cell_phases["primary_f4"].append(primary)
        for model_attempt in result["model_attempts"]:
            witness = process_metrics(model_attempt["witness_receipt"])
            witness_metrics.append(witness)
            cell_phases["f4_plus_sat_g"].append(witness)
            validation = model_attempt.get("validation")
            if isinstance(validation, dict):
                validator = process_metrics(validation["receipt"])
                validator_metrics.append(validator)
                cell_phases["validator"].append(validator)
    every_process = (
        primary_f4_metrics + witness_metrics + validator_metrics + interrupted_metrics
    )
    per_cell_distributions = []
    for cell_id in sorted(by_cell):
        phases = by_cell[cell_id]
        cell_results = [result for result in results if result["cell"]["id"] == cell_id]
        cell_statuses: dict[str, int] = {}
        for result in cell_results:
            status = str(result["final_status"])
            cell_statuses[status] = cell_statuses.get(status, 0) + 1
        per_cell_distributions.append(
            {
                "cell_id": cell_id,
                "tasks": len(cell_results),
                "statuses": dict(sorted(cell_statuses.items())),
                "primary_f4": phase_distribution(phases["primary_f4"]),
                "f4_plus_sat_g": phase_distribution(phases["f4_plus_sat_g"]),
                "validator": phase_distribution(phases["validator"]),
                "interrupted_prior": phase_distribution(phases["interrupted_prior"]),
                "all_processes": phase_distribution(
                    phases["primary_f4"]
                    + phases["f4_plus_sat_g"]
                    + phases["validator"]
                    + phases["interrupted_prior"]
                ),
            }
        )
    return {
        "schema": "koblitz_external_magma_summary.v1",
        "dry_run_manifest_sha256": packet["manifest_sha256"],
        "expected_tasks": 20,
        "selected_tasks": len(selection),
        "completed_tasks": len(results),
        "full_selection": full_selection,
        "all_inputs_attempted": all_attempted,
        "all_f4_terminals_resource_complete": all_f4_terminal and all_resources,
        "all_sat_point_witnesses_validated": all_witnesses,
        "no_contradictions": no_contradictions,
        "synthetic_test_mode": synthetic_test_mode,
        "archive_contracts_pinned": archive_contracts_pinned,
        "relevant_checkout_clean": relevant_checkout_clean,
        "tool_identity_boundary": (
            "Executable hashes and successful version probes bind the observed tools; "
            "they do not authenticate Magma license entitlement or reviewer independence"
        ),
        "status_counts": dict(sorted(status_counts.items())),
        "accounting_scope": (
            "Sum of separately metered primary F4, F4-plus-SAT(G), validator, and "
            "recovered interrupted processes; Python orchestration and version probes excluded"
        ),
        "charged_totals": {
            "primary_f4": charged_metrics(primary_f4_metrics),
            "f4_plus_sat_g": charged_metrics(witness_metrics),
            "validator": charged_metrics(validator_metrics),
            "interrupted_prior": charged_metrics(interrupted_metrics),
            "all_metered_processes": charged_metrics(every_process),
        },
        "per_cell_distributions": per_cell_distributions,
        "full_gate_passed": (
            full_selection
            and all_attempted
            and all_f4_terminal
            and all_resources
            and all_witnesses
            and no_contradictions
            and not synthetic_test_mode
            and archive_contracts_pinned
            and relevant_checkout_clean
        ),
        "claim_boundary": (
            "Tool-identity-bound Magma execution of the frozen planted-PDP panel; "
            "license entitlement and reviewer independence require external attestation, and "
            "this is not an end-to-end index-calculus, novelty, or SOTA result"
        ),
    }


def execute(args: argparse.Namespace) -> dict:
    packet = build_dry_run_manifest(args.protocol, args.panel, args.meter)
    tasks_by_id = {task["id"]: task for task in packet["tasks"]}
    if args.only_task:
        unknown = sorted(set(args.only_task) - set(tasks_by_id))
        if unknown:
            raise ExternalMagmaError(f"unknown --only-task values: {unknown}")
        selection = [task for task in packet["tasks"] if task["id"] in set(args.only_task)]
    else:
        selection = packet["tasks"]
    if not 1 <= args.model_cap <= MAX_MODEL_CAP:
        raise ExternalMagmaError(f"--model-cap must be between 1 and {MAX_MODEL_CAP}")
    output = args.output.resolve()
    panel_root = Path(packet["archive_root"]).resolve()
    try:
        output.relative_to(panel_root)
    except ValueError:
        pass
    else:
        raise ExternalMagmaError("output must not equal or reside inside the frozen Stage 13 panel")
    progress_path = output / "run.json"
    if output.exists() and not args.resume:
        raise ExternalMagmaError("output exists; choose a new path or use --resume")

    magma = binary_identity(args.magma, [["--version"]])
    backend = binary_identity(args.backend, [["--version"]])
    minisat = binary_identity(args.minisat, [["--version"], ["-h"]])
    for name, identity in (("magma", magma), ("backend", backend), ("minisat", minisat)):
        if identity.get("version_returncode") != 0 or identity.get("version") == "unreported":
            raise ExternalMagmaError(f"{name} did not return a successful version probe")
    magma_version_recognized = recognized_magma_version(magma["version"])
    magma["recognized_version_banner"] = magma_version_recognized
    if not args.synthetic_test_mode and not magma_version_recognized:
        raise ExternalMagmaError(
            "Magma --version must contain a 'Magma V2.<digits>-<digits>' banner or "
            "a standalone '2.<digits>-<digits>' line"
        )
    if Path(minisat["path"]).name != "minisat":
        raise ExternalMagmaError("--minisat must resolve to an executable named 'minisat'")
    meter = args.meter.resolve()
    tools = {
        "magma": magma,
        "backend": backend,
        "minisat": minisat,
        "meter": {"path": str(meter), "sha256": sha256_file(meter)},
        "runner": {"path": str(Path(__file__).resolve()), "sha256": sha256_file(Path(__file__))},
        "matrix_contract": {
            "path": str((HERE / "run_koblitz_pdp_matrix.py").resolve()),
            "sha256": sha256_file(HERE / "run_koblitz_pdp_matrix.py"),
        },
        "stage13_verifier": {
            "path": str((STAGE / "verify_stage13_pdp_panel.py").resolve()),
            "sha256": sha256_file(STAGE / "verify_stage13_pdp_panel.py"),
            "expected_sha256": PINNED_STAGE13_VERIFIER_SHA256,
        },
    }
    policy = {
        "parallel_workers": 1,
        "magma_threads": 1,
        "gpu_disabled_by_input": True,
        "f4_watchdog_seconds": args.f4_timeout,
        "model_watchdog_seconds": args.model_timeout,
        "validation_watchdog_seconds": args.validation_timeout,
        "model_cap": args.model_cap,
        "unknown_and_timeout_are_inconclusive": True,
        "witness_process_includes_f4_recomputation": True,
        "synthetic_test_mode": args.synthetic_test_mode,
        "tool_identity_bound": True,
        "license_entitlement_authenticated": False,
        "archive_contracts_pinned": packet["contract_pins"][
            "archive_tree_matches_pinned_commit"
        ],
        "relevant_checkout_clean": packet["relevant_checkout"]["dirty"] is False,
        "relevant_checkout_state": packet["relevant_checkout"],
    }
    frozen = {
        "dry_run_manifest_sha256": packet["manifest_sha256"],
        "selection": [task["id"] for task in selection],
        "tools": tools,
        "host": host_identity(),
        "policy": policy,
    }
    if not output.exists():
        output.mkdir(parents=True)
        progress = {
            "schema": RUN_SCHEMA,
            **frozen,
            "started_at": now(),
            "status": "running",
            "tasks": {},
        }
        atomic_json(progress_path, progress)
    else:
        progress = read_json(progress_path)
        changed = [key for key, value in frozen.items() if progress.get(key) != value]
        if progress.get("schema") != RUN_SCHEMA:
            changed.append("schema")
        if changed:
            raise ExternalMagmaError(f"resume changes frozen run identity: {', '.join(changed)}")
        progress["status"] = "running"
        progress.setdefault("resumed_at", []).append(now())
        atomic_json(progress_path, progress)

    environment = os.environ.copy()
    environment["PATH"] = str(Path(minisat["path"]).parent) + os.pathsep + environment.get("PATH", "")
    results = []
    for task in selection:
        result_path = output / "tasks" / task["id"] / "task-result.json"
        if result_path.is_file():
            indexed = progress.get("tasks", {}).get(task["id"])
            if indexed is not None and (
                not isinstance(indexed, dict)
                or indexed.get("receipt") != str(result_path.resolve())
                or indexed.get("receipt_sha256") != sha256_file(result_path)
            ):
                raise ExternalMagmaError(
                    f"resume task index disagrees with its immutable leaf: {task['id']}"
                )
            result = validate_completed_task(result_path, task, packet, tools, policy)
        else:
            result = execute_task(
                packet,
                task,
                output,
                meter,
                magma,
                backend,
                environment,
                policy,
            )
        results.append(result)
        progress["tasks"][task["id"]] = {
            "final_status": result["final_status"],
            "resource_complete": result["resource_complete"],
            "receipt": str(result_path.resolve()),
            "receipt_sha256": sha256_file(result_path),
        }
        atomic_json(progress_path, progress)

    summary = summarize_run(packet, selection, results, policy)
    summary["finished_at"] = now()
    atomic_json(output / "summary.json", summary)
    progress["status"] = "complete" if summary["all_inputs_attempted"] else "partial"
    progress["finished_at"] = summary["finished_at"]
    progress["summary"] = summary
    progress["summary_sha256"] = sha256_file(output / "summary.json")
    atomic_json(progress_path, progress)
    return summary


def self_test() -> dict:
    packet = build_dry_run_manifest(DEFAULT_PROTOCOL, DEFAULT_PANEL, DEFAULT_METER)
    if packet["expected_tasks"] != 20 or len(packet["tasks"]) != 20:
        raise AssertionError("wrong dry-run inventory")
    first = packet["tasks"][0]
    if first["id"] != "seed-2026091301/n31-l5-m3-standard-a1-f0":
        raise AssertionError("wrong first frozen task")
    panel = Path(packet["archive_root"])
    magma_path = resolve_manifest_path(panel, first["files"]["magma_boolean_f4"], "Magma input")
    script = render_witness_script(magma_path, [], int(first["source_variables"]))
    if "GroebnerBasis(I : Al := \"Direct\"" not in script or "SAT(G : Verbose := false)" not in script:
        raise AssertionError("witness script lacks direct F4 or SAT(G)")
    sample = "\n".join(
        [
            f"KOBLITZ_MODEL_SCHEMA={MODEL_SCHEMA}",
            "KOBLITZ_MODEL_STATUS=SAT",
            "KOBLITZ_MODEL_BITS=" + "0" * int(first["source_variables"]),
            f"KOBLITZ_MODEL_VARIABLES={first['source_variables']}",
            "KOBLITZ_MODEL_EXCLUDED_COUNT=0",
            "KOBLITZ_MODEL_CPU_SECONDS=0.1",
            "KOBLITZ_MODEL_WALL_SECONDS=0.2",
        ]
    )
    if parse_model_terminal(sample, int(first["source_variables"]), 0) is None:
        raise AssertionError("strict model parser rejected a valid record")
    if parse_model_terminal(sample + "\nKOBLITZ_MODEL_STATUS=UNSAT", int(first["source_variables"]), 0) is not None:
        raise AssertionError("strict model parser accepted duplicate markers")
    return {
        "self_test": "pass",
        "expected_tasks": len(packet["tasks"]),
        "manifest_sha256": packet["manifest_sha256"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--meter", type=Path, default=DEFAULT_METER)
    parser.add_argument("--magma")
    parser.add_argument("--backend", type=Path, default=DEFAULT_BACKEND)
    parser.add_argument("--minisat")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--only-task", action="append", default=[])
    parser.add_argument("--f4-timeout", type=float, default=120.0)
    parser.add_argument("--model-timeout", type=float, default=120.0)
    parser.add_argument("--validation-timeout", type=float, default=120.0)
    parser.add_argument("--model-cap", type=int, default=MAX_MODEL_CAP)
    parser.add_argument(
        "--synthetic-test-mode",
        action="store_true",
        help="permit fake tool banners for tests; this permanently disables full_gate_passed",
    )
    parser.add_argument("--dry-run-manifest", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        if args.self_test:
            print(json.dumps(self_test(), indent=2, sort_keys=True))
            return
        if args.dry_run_manifest:
            print(
                json.dumps(
                    build_dry_run_manifest(args.protocol, args.panel, args.meter),
                    indent=2,
                    sort_keys=True,
                )
            )
            return
        if args.magma is None or args.minisat is None or args.output is None:
            parser.error("normal execution requires --magma, --minisat, and --output")
        for timeout_name in ("f4_timeout", "model_timeout", "validation_timeout"):
            value = getattr(args, timeout_name)
            if not math.isfinite(value) or value <= 0:
                parser.error(f"--{timeout_name.replace('_', '-')} must be finite and positive")
        summary = execute(args)
        print(json.dumps(summary, indent=2, sort_keys=True))
    except (ExternalMagmaError, STAGE13.VerificationError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
