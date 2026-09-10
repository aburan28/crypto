#!/usr/bin/env python3
"""Run the frozen Stage 18 degree-23 IC versus automorphism-rho panel."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any

import verify_stage18_degree23_panel as verify


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_PROTOCOL = HERE / "stage-18-degree23-panel-protocol.json"
DEFAULT_AMENDMENT = HERE / "stage-18-amendment-01-lock-correction.json"
DEFAULT_OUTPUT = HERE / "stage-18-degree23-panel-lock-corrected-20260909"
FAILED_V1_OUTPUT = HERE / "stage-18-degree23-panel-20260909"
DEFAULT_METER = REPO / "scripts" / "process_meter.py"
DEFAULT_LOCK = HERE / "stage-18-corrected-Cargo.lock"
DEFAULT_IC_BINARY = REPO / "target" / "release" / "examples" / "koblitz_algebraic_e2e"
DEFAULT_DISCOVERY_BINARY = REPO / "target" / "release" / "examples" / "koblitz_public_factor_base_discovery"
LOCK_ARCHIVE_RELATIVE = Path("dependency-lock/Cargo.lock")
BUILD_WATCHDOG_SECONDS = 1800.0


class RunnerError(RuntimeError):
    pass


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w") as output:
        json.dump(value, output, indent=2, sort_keys=True)
        output.write("\n")
        output.flush()
        os.fsync(output.fileno())
    temporary.replace(path)


def command_output(command: list[str]) -> dict:
    completed = subprocess.run(command, cwd=REPO, text=True, capture_output=True, check=False)
    return {
        "command": command,
        "returncode": completed.returncode,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
    }


def file_identity(path: Path) -> dict:
    resolved = path.resolve(strict=True)
    if path.is_symlink() or not resolved.is_file():
        raise RunnerError(f"required tool is not a regular file: {path}")
    return {"path": str(resolved), "bytes": resolved.stat().st_size, "sha256": verify.sha256_file(resolved)}


def protocol_source(path: Path) -> Path:
    resolved = path.resolve(strict=True)
    if path.is_symlink() or not resolved.is_file():
        raise RunnerError("Stage 18 protocol is not a regular file")
    if resolved != DEFAULT_PROTOCOL.resolve(strict=True):
        raise RunnerError("scientific execution requires the committed Stage 18 protocol path")
    return resolved


def amendment_source(path: Path, protocol: dict) -> tuple[Path, dict]:
    resolved = path.resolve(strict=True)
    if path.is_symlink() or not resolved.is_file():
        raise RunnerError("Stage 18 lock amendment is not a regular file")
    if resolved != DEFAULT_AMENDMENT.resolve(strict=True):
        raise RunnerError("corrected execution requires the committed Stage 18 lock amendment")
    amendment = verify.read_json(resolved)
    verify.validate_amendment(amendment, protocol)
    if verify.sha256_file(resolved) != verify.EXPECTED_AMENDMENT_FILE_SHA256:
        raise RunnerError("Stage 18 lock amendment bytes changed")
    return resolved, amendment


def host_identity() -> dict:
    memory_bytes = None
    try:
        memory_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (AttributeError, OSError, ValueError):
        pass
    cpu_model = platform.processor() or None
    if platform.system() == "Darwin":
        probe = command_output(["sysctl", "-n", "machdep.cpu.brand_string"])
        if probe["returncode"] == 0 and probe["stdout"]:
            cpu_model = probe["stdout"]
        memory_probe = command_output(["sysctl", "-n", "hw.memsize"])
        if memory_probe["returncode"] == 0 and memory_probe["stdout"].isdigit():
            memory_bytes = int(memory_probe["stdout"])
    return {
        "platform": platform.platform(),
        "uname": list(platform.uname()),
        "cpu_model": cpu_model,
        "logical_cpus": os.cpu_count(),
        "physical_memory_bytes": memory_bytes,
        "python": command_output([str(Path(sys.executable).resolve()), "--version"]),
        "rustc": command_output(["rustc", "-Vv"]),
        "cargo": command_output(["cargo", "-V"]),
    }


def git_state(output: Path | None = None) -> dict:
    commit = command_output(["git", "rev-parse", "HEAD"])
    if commit["returncode"] != 0:
        raise RunnerError("checkout is not a Git repository")
    branch = command_output(["git", "branch", "--show-current"])
    ancestry = subprocess.run(
        ["git", "merge-base", "--is-ancestor", verify.ALGORITHM_BASE_COMMIT, "HEAD"],
        cwd=REPO, check=False,
    ).returncode == 0
    status_command = ["git", "status", "--porcelain=v1", "--untracked-files=all", "--", "."]
    if output is not None:
        try:
            relative = output.resolve().relative_to(REPO.resolve())
        except ValueError:
            pass
        else:
            status_command.extend([f":(exclude){relative}", f":(exclude){relative}/**"])
    status = command_output(status_command)
    if status["returncode"] != 0:
        raise RunnerError("cannot inspect checkout status")
    lines = status["stdout"].splitlines() if status["stdout"] else []
    tracked_delta = verify.git_tracked_delta(REPO, commit["stdout"])
    return {
        "commit": commit["stdout"], "branch": branch["stdout"] or None,
        "algorithm_base_commit": verify.ALGORITHM_BASE_COMMIT,
        "algorithm_base_is_ancestor": ancestry, "tracked_delta": tracked_delta,
        "dirty": bool(lines), "porcelain": lines,
    }


def source_identities() -> dict:
    identities = {}
    for relative, expected in verify.SOURCE_SHA256.items():
        path = REPO / relative
        identity = file_identity(path)
        if identity["sha256"] != expected:
            raise RunnerError(f"source hash differs from frozen Stage 18 input: {relative}")
        identity["expected_sha256"] = expected
        identities[relative] = identity
    lock = file_identity(DEFAULT_LOCK)
    if lock["sha256"] != verify.LOCK_SHA256:
        raise RunnerError("tracked dependency lock snapshot changed")
    lock["expected_sha256"] = verify.LOCK_SHA256
    identities["dependency_lock_snapshot"] = lock
    return identities


def assert_source_unchanged(frozen: dict, output: Path, allow_dirty: bool) -> None:
    current = source_identities()
    if current != frozen["sources"]:
        raise RunnerError("frozen source identities changed during execution")
    state = git_state(output)
    if not state["algorithm_base_is_ancestor"]:
        raise RunnerError("frozen algorithm commit is no longer an ancestor")
    if state["commit"] != frozen["source_revision"]["commit"]:
        raise RunnerError("checkout commit changed during execution")
    if state["tracked_delta"] != frozen["source_revision"].get("tracked_delta"):
        raise RunnerError("tracked source delta changed during execution")
    if state["dirty"] and not allow_dirty:
        raise RunnerError("scientific checkout became dirty during execution")
    if not allow_dirty:
        verify.require_expected_control_delta(state["tracked_delta"])


def safe_output(path: Path) -> Path:
    resolved = path.resolve()
    if resolved in {Path("/").resolve(), REPO.resolve(), (REPO / ".git").resolve()}:
        raise RunnerError(f"unsafe output directory: {resolved}")
    if path.exists() and path.is_symlink():
        raise RunnerError("output directory cannot be a symlink")
    return resolved


def corrected_output(path: Path) -> Path:
    resolved = safe_output(path)
    if resolved == FAILED_V1_OUTPUT.resolve():
        raise RunnerError("the retained v1 failure cannot be resumed or overwritten")
    if resolved != DEFAULT_OUTPUT.resolve():
        raise RunnerError("corrected Stage 18 requires the amendment-bound fresh output path")
    return resolved


def install_lockfile(panel: Path) -> dict:
    source = DEFAULT_LOCK.resolve(strict=True)
    destination = REPO / "Cargo.lock"
    archive = panel / LOCK_ARCHIVE_RELATIVE
    archive.parent.mkdir(parents=True, exist_ok=True)
    if archive.exists():
        if archive.is_symlink() or not archive.is_file():
            raise RunnerError("archived Cargo.lock is not a regular file")
        if verify.sha256_file(archive) != verify.LOCK_SHA256:
            raise RunnerError("archived Cargo.lock differs from the frozen snapshot")
    else:
        shutil.copyfile(source, archive)
    if verify.sha256_file(archive) != verify.LOCK_SHA256:
        raise RunnerError("archived Cargo.lock failed its custody check")
    if destination.exists():
        if destination.is_symlink() or not destination.is_file():
            raise RunnerError("workspace Cargo.lock is not a regular file")
        if verify.sha256_file(destination) != verify.LOCK_SHA256:
            raise RunnerError("workspace Cargo.lock differs from the frozen snapshot")
        installed = False
    else:
        shutil.copyfile(source, destination)
        installed = True
    if verify.sha256_file(destination) != verify.LOCK_SHA256:
        raise RunnerError("installed Cargo.lock failed its custody check")
    return {
        "source": file_identity(source), "destination": file_identity(destination),
        "archive": file_identity(archive),
        "installed_by_runner": installed,
    }


def meter_command(meter: Path, cwd: Path, timeout: float, stdout: Path, stderr: Path,
                  metrics: Path, command: list[str]) -> list[str]:
    return [
        str(Path(sys.executable).resolve()), str(meter.resolve()),
        "--cwd", str(cwd.resolve()), "--timeout", str(float(timeout)),
        "--stdout", str(stdout.resolve()), "--stderr", str(stderr.resolve()),
        "--metrics", str(metrics.resolve()), "--", *command,
    ]


def run_meter(meter: Path, cwd: Path, timeout: float, directory: Path,
              command: list[str], stdout_name: str = "result.json") -> int:
    invocation = {
        "schema": "koblitz_degree23_process_invocation.v1", "cwd": str(cwd.resolve()),
        "command": command, "watchdog_seconds": timeout, "started_at": now(),
    }
    atomic_json(directory / "process-invocation.json", invocation)
    wrapped = meter_command(
        meter, cwd, timeout, directory / stdout_name, directory / "stderr.txt",
        directory / "metrics.json", command,
    )
    started = time.perf_counter()
    completed = subprocess.run(wrapped, cwd=REPO, check=False)
    invocation.update({
        "meter_command": wrapped, "meter_launcher_returncode": completed.returncode,
        "driver_observed_wall_seconds": time.perf_counter() - started, "finished_at": now(),
    })
    atomic_json(directory / "process-invocation.json", invocation)
    return completed.returncode


def build_command() -> list[str]:
    cargo = shutil.which("cargo")
    if cargo is None:
        raise RunnerError("cargo is unavailable")
    return [
        str(Path(cargo).resolve()), "build", "--release", "--locked",
        "--example", "koblitz_algebraic_e2e",
        "--example", "koblitz_public_factor_base_discovery",
    ]


def validate_build(panel: Path, expected_command: list[str]) -> dict:
    directory = panel / "build"
    expected = {"process-invocation.json", "stdout.txt", "stderr.txt", "metrics.json", "receipt.json"}
    if not directory.is_dir() or directory.is_symlink():
        raise RunnerError("build receipt is missing")
    observed = {path.name for path in directory.iterdir()}
    if observed != expected or any(path.is_symlink() or not path.is_file() for path in directory.iterdir()):
        raise RunnerError("build artifact inventory changed")
    meter = verify.read_json(directory / "metrics.json")
    _, success = verify.validate_meter(meter, expected_command, BUILD_WATCHDOG_SECONDS)
    if not success:
        raise RunnerError("locked release build did not complete")
    receipt = verify.read_json(directory / "receipt.json")
    if receipt.get("command") != expected_command or receipt.get("status") != "verified":
        raise RunnerError("build receipt changed")
    for key, path in (("ic_binary", DEFAULT_IC_BINARY), ("discovery_binary", DEFAULT_DISCOVERY_BINARY)):
        if receipt.get("binaries", {}).get(key) != file_identity(path):
            raise RunnerError(f"built executable identity changed: {key}")
    return receipt


def ensure_build(panel: Path, meter: Path) -> dict:
    command = build_command()
    final = panel / "build"
    staging = panel / ".staging-build"
    if final.exists():
        if staging.exists():
            raise RunnerError("both completed and interrupted build leaves exist")
        return validate_build(panel, command)
    if staging.exists():
        raise RunnerError("interrupted build retained; use a new panel path")
    staging.mkdir()
    run_meter(meter, REPO, BUILD_WATCHDOG_SECONDS, staging, command, "stdout.txt")
    meter_record = verify.read_json(staging / "metrics.json")
    _, success = verify.validate_meter(meter_record, command, BUILD_WATCHDOG_SECONDS)
    receipt = {
        "schema": "koblitz_degree23_build_receipt.v1", "status": "verified" if success else "failed",
        "command": command, "dependency_lock_sha256": verify.sha256_file(REPO / "Cargo.lock"),
        "binaries": {},
    }
    if success:
        receipt["binaries"] = {
            "ic_binary": file_identity(DEFAULT_IC_BINARY),
            "discovery_binary": file_identity(DEFAULT_DISCOVERY_BINARY),
        }
    atomic_json(staging / "receipt.json", receipt)
    staging.replace(final)
    if not success:
        raise RunnerError("locked release build failed; retained build receipt")
    return validate_build(panel, command)


def initial_run(protocol: dict, protocol_path: Path, amendment: dict,
                amendment_path: Path, failed_v1_custody: dict, panel: Path,
                meter: Path, allow_dirty: bool) -> dict:
    state = git_state(panel)
    if not state["algorithm_base_is_ancestor"]:
        raise RunnerError("PR105/head is not an ancestor of the execution checkout")
    if state["dirty"] and not allow_dirty:
        raise RunnerError("scientific execution requires a clean checkout")
    if not allow_dirty:
        verify.require_expected_control_delta(state["tracked_delta"])
    sources = source_identities()
    frozen_protocol = panel / "protocol.json"
    frozen_protocol.write_bytes(protocol_path.read_bytes())
    frozen_amendment = panel / "amendment.json"
    frozen_amendment.write_bytes(amendment_path.read_bytes())
    return {
        "schema": verify.RUN_SCHEMA,
        "protocol_sha256": verify.canonical_sha256(protocol),
        "protocol_file_sha256": verify.sha256_file(frozen_protocol),
        "amendment_sha256": verify.canonical_sha256(amendment),
        "amendment_file_sha256": verify.sha256_file(frozen_amendment),
        "amendment_source": file_identity(amendment_path),
        "failed_v1_custody": failed_v1_custody,
        "execution_roots": {
            "repository": str(REPO.resolve()), "panel": str(panel.resolve()),
        },
        "source_revision": state,
        "sources": sources,
        "host": host_identity(),
        "evidence_class": "operational_smoke" if allow_dirty else "scientific_candidate",
        "implementation": {
            "runner": file_identity(Path(__file__)),
            "verifier": file_identity(HERE / "verify_stage18_degree23_panel.py"),
            "meter": file_identity(meter),
        },
        "tools": {
            "ic_binary": {"path": str(DEFAULT_IC_BINARY.resolve()), "sha256": None},
            "discovery_binary": {"path": str(DEFAULT_DISCOVERY_BINARY.resolve()), "sha256": None},
        },
        "started_at": now(), "status": "initializing", "tasks": {}, "outer_attempts": [],
    }


def immutable_run_identity(run: dict, protocol: dict, meter: Path, panel: Path,
                           allow_dirty: bool) -> None:
    if run.get("schema") != verify.RUN_SCHEMA:
        raise RunnerError("resume run schema changed")
    if run.get("protocol_sha256") != verify.canonical_sha256(protocol):
        raise RunnerError("resume protocol changed")
    if run.get("execution_roots") != {
        "repository": str(REPO.resolve()), "panel": str(panel.resolve()),
    }:
        raise RunnerError("resume moved the execution repository or panel")
    implementation = {
        "runner": file_identity(Path(__file__)),
        "verifier": file_identity(HERE / "verify_stage18_degree23_panel.py"),
        "meter": file_identity(meter),
    }
    if run.get("implementation") != implementation:
        raise RunnerError("runner, verifier, or meter changed before resume")
    if run.get("evidence_class") != ("operational_smoke" if allow_dirty else "scientific_candidate"):
        raise RunnerError("resume changed evidence class")
    assert_source_unchanged(run, panel, allow_dirty)


def invalid_receipt(task: dict, leaf: Path, error: Exception) -> dict:
    receipt = {
        "schema": verify.RECEIPT_SCHEMA, "task_id": task["id"], "attempt": 1,
        "command": task["command"], "status": "invalid_artifact",
        "verified_terminal": False, "error": str(error),
        "result_sha256": verify.sha256_file(leaf / "result.json") if (leaf / "result.json").is_file() else None,
        "stderr_sha256": verify.sha256_file(leaf / "stderr.txt") if (leaf / "stderr.txt").is_file() else None,
        "metrics_sha256": verify.sha256_file(leaf / "metrics.json") if (leaf / "metrics.json").is_file() else None,
    }
    if (leaf / "metrics.json").is_file():
        try:
            meter = verify.read_json(leaf / "metrics.json")
            metrics, _ = verify.validate_meter(meter, task["command"], 600)
        except verify.VerificationError:
            pass
        else:
            receipt["metrics"] = metrics
    return receipt


def execute_task(protocol: dict, panel: Path, meter: Path, task: dict) -> dict:
    final = panel / "tasks" / task["id"]
    staging = panel / ".staging" / task["id"]
    if final.exists():
        if staging.exists():
            raise RunnerError(f"completed and interrupted leaves coexist: {task['id']}")
        return verify.verify_task_leaf(protocol, panel, task)
    if staging.exists():
        raise RunnerError(f"interrupted attempt retained for {task['id']}; use a new panel path")
    staging.mkdir(parents=True)
    invocation = {
        "schema": "koblitz_degree23_task_invocation.v1", "task": task,
        "cwd": str(REPO.resolve()), "command": task["command"],
        "watchdog_seconds": 600.0, "attempt": 1, "started_at": now(),
    }
    atomic_json(staging / "invocation.json", invocation)
    wrapped = meter_command(
        meter, REPO, 600.0, staging / "result.json", staging / "stderr.txt",
        staging / "metrics.json", task["command"],
    )
    started = time.perf_counter()
    completed = subprocess.run(wrapped, cwd=REPO, check=False)
    invocation.update({
        "meter_command": wrapped, "meter_launcher_returncode": completed.returncode,
        "driver_observed_wall_seconds": time.perf_counter() - started,
        "finished_at": now(),
    })
    atomic_json(staging / "invocation.json", invocation)
    final.parent.mkdir(parents=True, exist_ok=True)
    staging.replace(final)
    try:
        return verify.verify_task_leaf(protocol, panel, task, write_receipt=True)
    except verify.VerificationError as error:
        receipt = invalid_receipt(task, final, error)
        atomic_json(final / "receipt.json", receipt)
        return receipt


def write_preliminary_summary(protocol: dict, panel: Path) -> dict:
    summary = verify.summarize(
        protocol, panel, allow_incomplete=True, validate_controls=False,
    )
    atomic_json(panel / "summary.json", summary)
    return summary


def post_execution_custody(run: dict) -> dict:
    binaries = {}
    for key, identity in run["tools"].items():
        binaries[key] = file_identity(Path(identity["path"]))
        if binaries[key] != identity:
            raise RunnerError(f"binary changed during execution: {key}")
    return {
        "sources": source_identities(),
        "binaries": binaries,
        "root_lock_sha256": verify.sha256_file(REPO / "Cargo.lock"),
    }


def inner_execute(args: argparse.Namespace) -> int:
    if args.resume:
        raise RunnerError("corrected Stage 18 cannot resume or retry any prior panel")
    protocol_path = protocol_source(args.protocol)
    protocol = verify.read_json(protocol_path)
    verify.validate_protocol(protocol)
    amendment_path, amendment = amendment_source(args.amendment, protocol)
    failed_v1_custody = verify.validate_failed_v1_custody(amendment, protocol)
    panel = corrected_output(args.output)
    meter = args.meter.resolve(strict=True)
    if verify.sha256_file(meter) != verify.SOURCE_SHA256["scripts/process_meter.py"]:
        raise RunnerError("process meter differs from the frozen implementation")
    run_path = panel / "run.json"
    if run_path.exists():
        raise RunnerError("corrected Stage 18 output already contains a run and cannot be retried")
    panel.mkdir(parents=True, exist_ok=True)
    run = initial_run(
        protocol, protocol_path, amendment, amendment_path, failed_v1_custody,
        panel, meter, args.allow_dirty,
    )
    atomic_json(run_path, run)
    run["status"] = "running"
    atomic_json(run_path, run)

    run["dependency_lock"] = install_lockfile(panel)
    atomic_json(run_path, run)
    try:
        build = ensure_build(panel, meter)
    except RunnerError:
        build_receipt = panel / "build" / "receipt.json"
        if build_receipt.is_file():
            run["build_receipt_sha256"] = verify.sha256_file(build_receipt)
            run["status"] = "inconclusive_build_failure"
            atomic_json(run_path, run)
        raise
    run["build_receipt_sha256"] = verify.sha256_file(panel / "build" / "receipt.json")
    run["tools"] = build["binaries"]
    atomic_json(run_path, run)
    tasks = verify.task_plan(
        protocol, run["tools"]["discovery_binary"]["path"], run["tools"]["ic_binary"]["path"]
    )
    for position, task in enumerate(tasks):
        if position == 2:
            discovery_ok = all(
                run.get("tasks", {}).get(key, {}).get("status") == "verified"
                for key in ("discovery-a0", "discovery-a1")
            )
            if not discovery_ok:
                run["status"] = "inconclusive_public_discovery"
                run["post_execution"] = post_execution_custody(run)
                atomic_json(run_path, run)
                write_preliminary_summary(protocol, panel)
                return 1
        assert_source_unchanged(run, panel, args.allow_dirty)
        receipt = execute_task(protocol, panel, meter, task)
        run["tasks"][task["id"]] = {
            "status": receipt["status"],
            "receipt": str((panel / "tasks" / task["id"] / "receipt.json").resolve()),
            "receipt_sha256": verify.sha256_file(panel / "tasks" / task["id"] / "receipt.json"),
        }
        atomic_json(run_path, run)
    assert_source_unchanged(run, panel, args.allow_dirty)
    run["post_execution"] = post_execution_custody(run)
    atomic_json(run_path, run)
    summary = write_preliminary_summary(protocol, panel)
    run["status"] = summary["status"]
    run["finished_tasks_at"] = now()
    run["summary_sha256"] = verify.sha256_file(panel / "summary.json")
    atomic_json(run_path, run)
    return 0 if summary.get("task_panel_complete") and run["evidence_class"] == "scientific_candidate" else 1


def next_outer_attempt(panel: Path) -> tuple[int, Path, Path]:
    root = panel / "outer-attempts"
    root.mkdir(parents=True, exist_ok=True)
    if any(path.name.startswith(".staging-") for path in root.iterdir()):
        raise RunnerError("interrupted outer driver attempt is retained; inspect it before continuing")
    indices = [int(path.name) for path in root.iterdir() if path.is_dir() and path.name.isdigit()]
    index = max(indices, default=0) + 1
    return index, root / f".staging-{index:04d}", root / f"{index:04d}"


def artifact_manifest(panel: Path) -> dict:
    excluded = {"run.json", "artifact-manifest.json", "verification.json"}
    files = []
    for path in sorted(panel.rglob("*")):
        if path.is_dir():
            if path.is_symlink():
                raise RunnerError(f"symlink directory in panel: {path}")
            continue
        if path.is_symlink() or not path.is_file():
            raise RunnerError(f"non-regular panel artifact: {path}")
        relative = str(path.relative_to(panel))
        if relative in excluded or relative.endswith(".tmp"):
            continue
        files.append({"path": relative, "bytes": path.stat().st_size, "sha256": verify.sha256_file(path)})
    return {
        "schema": "koblitz_degree23_artifact_manifest.v1",
        "excluded_mutable_controls": sorted(excluded),
        "files": files, "file_count": len(files),
    }


def finalize_outer(protocol: dict, panel: Path, outer_index: int) -> dict:
    try:
        summary = verify.summarize(
            protocol, panel, allow_incomplete=False, validate_controls=False
        )
    except verify.VerificationError:
        summary = verify.summarize(
            protocol, panel, allow_incomplete=True, validate_controls=False,
        )
    atomic_json(panel / "summary.json", summary)
    manifest = artifact_manifest(panel)
    atomic_json(panel / "artifact-manifest.json", manifest)
    verification = {
        "schema": "koblitz_degree23_panel_verification.v1",
        "status": summary["status"], "outer_attempt_finalized": outer_index,
        "protocol_sha256": verify.canonical_sha256(protocol),
        "amendment_sha256": summary["amendment_sha256"],
        "summary_sha256": verify.sha256_file(panel / "summary.json"),
        "artifact_manifest_sha256": verify.sha256_file(panel / "artifact-manifest.json"),
        "exact_frozen_task_count": 12,
        "claim_boundary": summary["claim_boundary"],
    }
    atomic_json(panel / "verification.json", verification)
    run_path = panel / "run.json"
    if run_path.is_file():
        run = verify.read_json(run_path)
        run["status"] = summary["status"]
        run["finished_at"] = now()
        run["summary_sha256"] = verification["summary_sha256"]
        run["artifact_manifest_sha256"] = verification["artifact_manifest_sha256"]
        run["verification_sha256"] = verify.sha256_file(panel / "verification.json")
        atomic_json(run_path, run)
    return summary


def outer_execute(args: argparse.Namespace) -> int:
    if args.resume:
        raise RunnerError("corrected Stage 18 cannot resume or retry any prior panel")
    protocol_path = protocol_source(args.protocol)
    protocol = verify.read_json(protocol_path)
    verify.validate_protocol(protocol)
    amendment_path, amendment = amendment_source(args.amendment, protocol)
    verify.validate_failed_v1_custody(amendment, protocol)
    panel = corrected_output(args.output)
    existed = panel.exists()
    if existed:
        raise RunnerError("corrected Stage 18 output already exists; retries require a new amendment")
    panel.mkdir(parents=True, exist_ok=True)
    index, staging, final = next_outer_attempt(panel)
    staging.mkdir()
    child = [
        str(Path(sys.executable).resolve()), str(Path(__file__).resolve()), "--inner",
        "--protocol", str(protocol_path), "--amendment", str(amendment_path),
        "--output", str(panel),
        "--meter", str(args.meter.resolve()),
    ]
    if args.allow_dirty:
        child.append("--allow-dirty")
    invocation = {
        "schema": "koblitz_degree23_outer_invocation.v1", "attempt": index,
        "command": child, "cwd": str(REPO.resolve()),
        "watchdog_seconds": 7200.0, "started_at": now(),
    }
    atomic_json(staging / "invocation.json", invocation)
    wrapped = meter_command(
        args.meter, REPO, 7200.0, staging / "stdout.txt", staging / "stderr.txt",
        staging / "metrics.json", child,
    )
    started = time.perf_counter()
    completed = subprocess.run(wrapped, cwd=REPO, check=False)
    invocation.update({
        "meter_command": wrapped, "meter_launcher_returncode": completed.returncode,
        "driver_observed_wall_seconds": time.perf_counter() - started, "finished_at": now(),
    })
    atomic_json(staging / "invocation.json", invocation)
    metrics = verify.read_json(staging / "metrics.json")
    receipt = {
        "schema": "koblitz_degree23_outer_receipt.v1", "attempt": index,
        "returncode": metrics["returncode"], "timed_out": metrics["timed_out"],
        "orphan_group_terminated": metrics["orphan_group_terminated"],
        "metrics": metrics["metrics"], "command": child,
    }
    atomic_json(staging / "receipt.json", receipt)
    staging.replace(final)
    run_path = panel / "run.json"
    if run_path.is_file():
        run = verify.read_json(run_path)
        run.setdefault("outer_attempts", []).append({
            "attempt": index,
            "receipt": str((final / "receipt.json").resolve()),
            "receipt_sha256": verify.sha256_file(final / "receipt.json"),
        })
        atomic_json(run_path, run)
    summary = finalize_outer(protocol, panel, index) if (panel / "run.json").is_file() else None
    print(json.dumps(summary or receipt, indent=2, sort_keys=True))
    return 0 if summary and summary["status"] == "complete_verified_panel" and completed.returncode == 0 else 1


def plan(args: argparse.Namespace) -> dict:
    protocol_path = protocol_source(args.protocol)
    protocol = verify.read_json(protocol_path)
    verify.validate_protocol(protocol)
    amendment_path, amendment = amendment_source(args.amendment, protocol)
    failed_v1_custody = verify.validate_failed_v1_custody(amendment, protocol)
    output = corrected_output(args.output)
    if output.exists():
        raise RunnerError("corrected Stage 18 output already exists and cannot be retried")
    sources = source_identities()
    tasks = verify.task_plan(protocol, str(DEFAULT_DISCOVERY_BINARY.resolve()), str(DEFAULT_IC_BINARY.resolve()))
    state = git_state(args.output)
    if not args.allow_dirty:
        if state["dirty"]:
            raise RunnerError("scientific plan requires a clean checkout")
        verify.require_expected_control_delta(state["tracked_delta"])
    return {
        "schema": "koblitz_degree23_replication_panel_plan.v1",
        "protocol_sha256": verify.canonical_sha256(protocol),
        "amendment_sha256": verify.canonical_sha256(amendment),
        "amendment_file_sha256": verify.sha256_file(amendment_path),
        "failed_v1_custody": failed_v1_custody,
        "corrected_output": str(output),
        "source_revision": state, "tracked_delta": state["tracked_delta"], "sources": sources,
        "dependency_lock_destination": str((REPO / "Cargo.lock").resolve()),
        "build_command": build_command(), "build_watchdog_seconds": BUILD_WATCHDOG_SECONDS,
        "tasks": tasks, "inner_watchdog_seconds": 600,
        "outer_watchdog_seconds": 7200, "scientific_task_attempts": 1,
    }


def self_test() -> dict:
    protocol = verify.read_json(DEFAULT_PROTOCOL)
    verify.validate_protocol(protocol)
    amendment_path, amendment = amendment_source(DEFAULT_AMENDMENT, protocol)
    failed_v1_custody = verify.validate_failed_v1_custody(amendment, protocol)
    tasks = verify.task_plan(protocol, str(DEFAULT_DISCOVERY_BINARY.resolve()), str(DEFAULT_IC_BINARY.resolve()))
    if [task["id"] for task in tasks[:4]] != ["discovery-a0", "discovery-a1", "row-01-ic", "row-01-rho-auto"]:
        raise AssertionError("wrong frozen task ordering")
    if tasks[-1]["command"][-3:] != ["100000", "1000", "divisor:0,2"]:
        raise AssertionError("wrong frozen command tail")
    source_identities()
    checks = 6
    try:
        safe_output(REPO)
    except RunnerError:
        checks += 1
    else:
        raise AssertionError("repository root was accepted as output")
    command = ["/bin/true", "frozen"]
    with tempfile.TemporaryDirectory(prefix="stage18-runner-self-test-") as temporary:
        root = Path(temporary)
        path = root / "atomic.json"
        atomic_json(path, {"ok": True})
        if verify.read_json(path) != {"ok": True} or path.with_name("atomic.json.tmp").exists():
            raise AssertionError("atomic JSON helper failed")
        checks += 1
        wrapped = meter_command(DEFAULT_METER, REPO, 600.0, root / "out", root / "err", root / "metrics", command)
        if wrapped[-3:] != ["--", "/bin/true", "frozen"] or wrapped[2:6] != ["--cwd", str(REPO.resolve()), "--timeout", "600.0"]:
            raise AssertionError("meter command construction changed")
        checks += 1
        panel = root / "panel"
        panel.mkdir()
        index, staging, final = next_outer_attempt(panel)
        if (index, final.name) != (1, "0001"):
            raise AssertionError("outer attempt numbering changed")
        staging.mkdir()
        try:
            next_outer_attempt(panel)
        except RunnerError:
            checks += 1
        else:
            raise AssertionError("interrupted outer attempt was ignored")
        initial_panel = root / "initial-panel"
        initial_panel.mkdir()
        initial = initial_run(
            protocol, DEFAULT_PROTOCOL, amendment, amendment_path, failed_v1_custody,
            initial_panel, DEFAULT_METER, True,
        )
        if (
            initial.get("amendment_sha256") != verify.EXPECTED_AMENDMENT_SHA256
            or initial.get("failed_v1_custody") != failed_v1_custody
            or initial.get("evidence_class") != "operational_smoke"
        ):
            raise AssertionError("corrected initial-run binding changed")
        checks += 1
        try:
            corrected_output(FAILED_V1_OUTPUT)
        except RunnerError:
            checks += 1
        else:
            raise AssertionError("retained v1 output was accepted for corrected execution")
    state = git_state()
    if not state["algorithm_base_is_ancestor"]:
        raise AssertionError("frozen base ancestry check failed")
    delta = state.get("tracked_delta")
    if not isinstance(delta, dict) or delta.get("schema") != "koblitz_stage18_tracked_delta.v1":
        raise AssertionError("tracked source delta was not recorded")
    checks += 1
    synthetic = {
        "schema": "koblitz_stage18_tracked_delta.v1",
        "base_commit": verify.ALGORITHM_BASE_COMMIT,
        "head_commit": "2" * 40,
        "entries": [dict(entry, object_id="1" * 40) for entry in verify.EXPECTED_CONTROL_DELTA],
        "matches_expected_control_only_delta": True,
    }
    verify.require_expected_control_delta(synthetic)
    checks += 1
    return {
        "self_test": "pass", "checks": checks, "frozen_tasks": len(tasks),
        "protocol_sha256": verify.canonical_sha256(protocol),
        "amendment_sha256": verify.canonical_sha256(amendment),
        "failed_v1_scientific_tasks": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--amendment", type=Path, default=DEFAULT_AMENDMENT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--meter", type=Path, default=DEFAULT_METER)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true", help="permanently label the run operational_smoke")
    parser.add_argument("--plan", action="store_true", help="print the exact frozen plan without building or running")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--inner", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        if args.self_test:
            print(json.dumps(self_test(), indent=2, sort_keys=True))
            return
        if args.plan:
            print(json.dumps(plan(args), indent=2, sort_keys=True))
            return
        status = inner_execute(args) if args.inner else outer_execute(args)
        raise SystemExit(status)
    except (RunnerError, verify.VerificationError, OSError, subprocess.SubprocessError) as error:
        parser.exit(1, f"stage18 runner failed: {error}\n")


if __name__ == "__main__":
    main()
