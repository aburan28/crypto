#!/usr/bin/env python3
"""Create a charged clean Rust build capsule for the Phase-B native-F4 arm.

The historical Phase-B builder intentionally preserves its 2026-09-10 source
closure.  Current main also embeds ``docs/ic/calibration.json`` at compile
time, so this additive builder records that file without changing the legacy
receipt validator or invalidating archived Stage-20 build capsules.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import platform
import shutil
import stat
import subprocess
from typing import Any

import build_koblitz_phase_b_tools as legacy
import run_koblitz_blind_pdp_phase_b as phase_b


SCHEMA = "koblitz_phase_b_native_f4_tool_build_receipt.v1"
RUST_LOCK_PATH = legacy.RUST_LOCK_PATH
RUST_PATHS = [*legacy.RUST_PATHS, "docs/ic/calibration.json"]
ROLES = [
    "archive-source",
    "extract-source",
    "install-lock",
    "version-cargo",
    "version-rustc",
    "vendor",
    "build",
    "copy-exporter",
    "copy-backend",
]
BOUNDARY = {
    "cost_scope": (
        "Current source archival, extraction, frozen-lock installation, tool-version probes, "
        "vendoring, compilation and binary copies; dependency acquisition is separately charged"
    ),
    "cpu_scope": (
        "RUSAGE_CHILDREN includes every waited build child; one compilation job is requested"
    ),
    "memory_scope": (
        "largest child high-water RSS, not simultaneous aggregate process-tree memory"
    ),
    "single_core_elapsed_seconds": None,
    "full_cost_gate_passed": False,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise phase_b.PhaseBError(message)


def git(source: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=source, text=True, capture_output=True, check=True
    ).stdout.strip()


def clean_state(source: Path, expected: str | None = None) -> dict[str, Any]:
    commit = git(source, "rev-parse", "HEAD")
    phase_b.require_hex40(commit, "native-F4 build source revision")
    porcelain = git(source, "status", "--porcelain=v1", "--untracked-files=all").splitlines()
    require(not porcelain, f"native-F4 build source must be clean: {source}")
    if expected is not None:
        require(commit == expected, "native-F4 build source revision changed")
    links = []
    for line in git(source, "ls-files", "--stage").splitlines():
        metadata, path = line.split("\t", 1)
        mode, _, _ = metadata.split()
        if mode == "160000":
            links.append(path)
    require(not links, "native-F4 build source contains gitlinks")
    return {"commit": commit, "dirty": False, "porcelain": []}


def source_objects(source: Path, revision: str = "HEAD") -> dict[str, str]:
    if revision != "HEAD":
        phase_b.require_hex40(revision, "native-F4 Rust source revision")
    return {path: git(source, "rev-parse", f"{revision}:{path}") for path in RUST_PATHS}


def tool_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    metadata = resolved.stat()
    require(
        stat.S_ISREG(metadata.st_mode) and os.access(path, os.X_OK),
        f"build tool is not an executable regular file: {path}",
    )
    data = resolved.read_bytes()
    return {
        "path": str(path.absolute()),
        "resolved_path": str(resolved),
        "bytes": len(data),
        "sha256": phase_b.sha256_bytes(data),
    }


def resource_summary(processes: list[dict[str, Any]]) -> dict[str, Any]:
    metrics = [row["metrics"] for row in processes]
    return {
        "total_core_seconds": round(
            math.fsum(row["total_core_seconds"] for row in metrics), 12
        ),
        "summed_process_wall_seconds": round(
            math.fsum(row["wall_seconds"] for row in metrics), 12
        ),
        "largest_child_peak_rss_bytes": max(
            (row["peak_rss_bytes"] for row in metrics), default=0
        ),
        "single_core_elapsed_seconds": None,
        "aggregate_process_tree_peak_rss_bytes": None,
    }


def expected_commands(receipt: dict[str, Any]) -> dict[str, list[str]]:
    commands = receipt["toolchain"]
    output = Path(receipt["source_archive"]["path"]).parent
    source_root = output / "source"
    build_root = output / "build"
    jobs = str(receipt["requested_parallel_jobs"])
    return {
        "archive-source": [
            commands["git"]["path"],
            "archive",
            "--format=tar",
            f"--output={output / 'source.tar'}",
            receipt["source_commit"],
            *RUST_PATHS,
        ],
        "extract-source": [
            commands["tar"]["path"],
            "-xf",
            str(output / "source.tar"),
            "-C",
            str(source_root),
        ],
        "install-lock": [
            commands["cp"]["path"],
            str(source_root / RUST_LOCK_PATH),
            str(source_root / "Cargo.lock"),
        ],
        "version-cargo": [commands["cargo"]["path"], "--version", "--verbose"],
        "version-rustc": [commands["rustc"]["path"], "--version", "--verbose"],
        "vendor": [
            commands["cargo"]["path"],
            "vendor",
            "--locked",
            "--offline",
            "--versioned-dirs",
            str(output / "vendor"),
        ],
        "build": [
            commands["cargo"]["path"],
            "build",
            "--release",
            "--locked",
            "--offline",
            "--jobs",
            jobs,
            "--target-dir",
            str(build_root),
            "--example",
            "koblitz_pdp_export",
            "--example",
            "koblitz_pdp_backend",
        ],
        "copy-exporter": [
            commands["cp"]["path"],
            str(build_root / "release/examples/koblitz_pdp_export"),
            str(output / "bin/koblitz_pdp_export"),
        ],
        "copy-backend": [
            commands["cp"]["path"],
            str(build_root / "release/examples/koblitz_pdp_backend"),
            str(output / "bin/koblitz_pdp_backend"),
        ],
    }


def build(source: Path, output: Path, jobs: int, timeout: int) -> dict[str, Any]:
    source = source.resolve()
    output = output.resolve()
    require(type(jobs) is int and 1 <= jobs <= 64, "invalid native-F4 build job count")
    require(timeout > 0, "invalid native-F4 build timeout")
    require(not output.exists() and not output.is_symlink(), "build output must be new")
    state = clean_state(source)
    require(source == phase_b.REPO.resolve(), "build source must be this checkout")
    objects = source_objects(source, state["commit"])
    helper_hash = phase_b.sha256_file(Path(phase_b.__file__), "Phase-B helper")
    names = ("git", "tar", "cp", "cargo", "rustc")
    commands = {}
    for name in names:
        found = shutil.which(name)
        require(found is not None, f"required build tool missing: {name}")
        commands[name] = tool_identity(Path(found))

    output.mkdir()
    evidence = output / "evidence"
    evidence.mkdir()
    source_root = output / "source"
    source_root.mkdir()
    build_root = output / "build"
    environment = phase_b.safe_child_environment()
    environment["PATH"] = os.pathsep.join(
        dict.fromkeys(
            [str(Path(row["path"]).parent) for row in commands.values()]
            + environment["PATH"].split(os.pathsep)
        )
    )
    environment.update(
        {
            "HOME": str(Path.home()),
            "CARGO_HOME": str(
                Path(os.environ.get("CARGO_HOME", Path.home() / ".cargo")).resolve()
            ),
            "RUSTC": commands["rustc"]["path"],
            "CARGO_INCREMENTAL": "0",
        }
    )
    processes: list[dict[str, Any]] = []

    def measured(role: str, command: list[str], cwd: Path) -> dict[str, Any]:
        record = phase_b.run_metered(
            role=role,
            command=command,
            cwd=cwd,
            task_root=evidence,
            input_paths=[],
            timeout=timeout,
            meter=phase_b.DEFAULT_METER,
            environment=environment,
            markers=set(),
        )
        require(
            record["returncode"] == 0 and not record["timed_out"],
            f"native-F4 {role} failed; retain {evidence} and use a new output",
        )
        record["environment"] = dict(environment)
        processes.append(phase_b.compact_process(record))
        return record

    source_tar = output / "source.tar"
    measured(
        "archive-source",
        [
            commands["git"]["path"],
            "archive",
            "--format=tar",
            f"--output={source_tar}",
            state["commit"],
            *RUST_PATHS,
        ],
        source,
    )
    source_archive = phase_b.executable_identity(
        source_tar, "native-F4 source archive", executable=False
    )
    measured(
        "extract-source",
        [commands["tar"]["path"], "-xf", str(source_tar), "-C", str(source_root)],
        output,
    )
    measured(
        "install-lock",
        [
            commands["cp"]["path"],
            str(source_root / RUST_LOCK_PATH),
            str(source_root / "Cargo.lock"),
        ],
        output,
    )
    shutil.copyfile(source_root / "Cargo.lock", evidence / "Cargo.lock")
    measured("version-cargo", [commands["cargo"]["path"], "--version", "--verbose"], source_root)
    measured("version-rustc", [commands["rustc"]["path"], "--version", "--verbose"], source_root)
    vendor = measured(
        "vendor",
        [
            commands["cargo"]["path"],
            "vendor",
            "--locked",
            "--offline",
            "--versioned-dirs",
            str(output / "vendor"),
        ],
        source_root,
    )
    (source_root / ".cargo").mkdir()
    (source_root / ".cargo/config.toml").write_text(vendor["stdout_text"])
    (evidence / "cargo-config.toml").write_text(vendor["stdout_text"])
    cargo_home = output / "cargo-home"
    cargo_home.mkdir()
    environment["CARGO_HOME"] = str(cargo_home)
    measured(
        "build",
        [
            commands["cargo"]["path"],
            "build",
            "--release",
            "--locked",
            "--offline",
            "--jobs",
            str(jobs),
            "--target-dir",
            str(build_root),
            "--example",
            "koblitz_pdp_export",
            "--example",
            "koblitz_pdp_backend",
        ],
        source_root,
    )
    require(
        (source_root / "Cargo.lock").read_bytes() == (evidence / "Cargo.lock").read_bytes(),
        "frozen Cargo lock changed during native-F4 compilation",
    )
    packaged = output / "bin"
    packaged.mkdir()
    binaries = {}
    for name in ("koblitz_pdp_export", "koblitz_pdp_backend"):
        role = "exporter" if name.endswith("export") else "backend"
        source_binary = build_root / "release/examples" / name
        destination = packaged / name
        measured(
            f"copy-{role}",
            [commands["cp"]["path"], str(source_binary), str(destination)],
            output,
        )
        os.chmod(destination, 0o755)
        binaries[role] = phase_b.executable_identity(destination, role)

    require(clean_state(source, state["commit"]) == state, "source changed during native-F4 build")
    for name, identity in commands.items():
        require(tool_identity(Path(identity["path"])) == identity, f"build tool changed: {name}")
    require(
        phase_b.sha256_file(Path(phase_b.__file__), "Phase-B helper") == helper_hash,
        "Phase-B helper changed during native-F4 build",
    )
    receipt = {
        "schema": SCHEMA,
        "status": "completed",
        "source_commit": state["commit"],
        "source_clean": True,
        "source_state": state,
        "rust_source_paths": RUST_PATHS,
        "rust_source_objects": objects,
        "source_archive": source_archive,
        "requested_parallel_jobs": jobs,
        "toolchain": commands,
        "platform": {
            "system": platform.system(),
            "machine": platform.machine(),
            "release": platform.release(),
        },
        "environment": environment,
        "build_processes": processes,
        "binaries": binaries,
        "resources": resource_summary(processes),
        "accounting_boundary": BOUNDARY,
        "evidence_inventory": phase_b.all_regular_inventory(evidence),
        "build_driver_sha256": phase_b.sha256_file(Path(__file__), "native-F4 build driver"),
        "process_meter_sha256": phase_b.sha256_file(phase_b.DEFAULT_METER, "process meter"),
        "phase_b_helper_sha256": helper_hash,
    }
    receipt["receipt_payload_sha256"] = phase_b.canonical_sha256(receipt)
    phase_b.write_json_new(output / "receipt.json", receipt)
    validate_receipt(output / "receipt.json", binaries, objects, state)
    return receipt


def validate_receipt(
    path: Path,
    identities: dict[str, Any],
    objects: dict[str, str] | None = None,
    expected_state: dict[str, Any] | None = None,
) -> dict[str, Any]:
    value, raw = phase_b.read_json(path, "native-F4 build receipt")
    payload = dict(value)
    self_hash = payload.pop("receipt_payload_sha256", None)
    require(phase_b.canonical_sha256(payload) == self_hash, "build receipt self-hash differs")
    require(
        value.get("schema") == SCHEMA
        and value.get("status") == "completed"
        and value.get("source_clean") is True,
        "native-F4 build receipt is not complete and clean",
    )
    state = value.get("source_state")
    require(
        isinstance(state, dict)
        and set(state) == {"commit", "dirty", "porcelain"}
        and state.get("dirty") is False
        and state.get("porcelain") == [],
        "native-F4 build receipt has invalid source state",
    )
    phase_b.require_hex40(state.get("commit"), "native-F4 build receipt commit")
    require(value.get("source_commit") == state["commit"], "build commit fields disagree")
    if expected_state is not None:
        require(state == expected_state, "build receipt uses a different implementation state")
    require(value.get("rust_source_paths") == RUST_PATHS, "build source paths changed")
    require(set(value.get("rust_source_objects", {})) == set(RUST_PATHS), "build objects changed")
    if objects is not None:
        require(value["rust_source_objects"] == objects, "build source object ids changed")
    require(set(value.get("binaries", {})) == {"exporter", "backend"}, "build binaries changed")
    for name, identity in identities.items():
        require(
            all(value["binaries"][name].get(key) == identity.get(key) for key in ("sha256", "bytes")),
            f"native-F4 build does not bind {name}",
        )
    require(value.get("accounting_boundary") == BOUNDARY, "build accounting boundary changed")
    trusted = {
        "build_driver_sha256": phase_b.sha256_file(Path(__file__), "native-F4 build driver"),
        "process_meter_sha256": phase_b.sha256_file(phase_b.DEFAULT_METER, "process meter"),
        "phase_b_helper_sha256": phase_b.sha256_file(Path(phase_b.__file__), "Phase-B helper"),
    }
    for key, expected in trusted.items():
        require(value.get(key) == expected, f"native-F4 build {key} changed")
    evidence = path.parent / "evidence"
    require(
        phase_b.all_regular_inventory(evidence) == value.get("evidence_inventory"),
        "native-F4 build evidence inventory changed",
    )
    processes = value.get("build_processes")
    require(
        isinstance(processes, list) and [row.get("role") for row in processes] == ROLES,
        "native-F4 build process roles changed",
    )
    expected = expected_commands(value)
    expected_files = {"Cargo.lock", "cargo-config.toml"}
    for process in processes:
        role = process["role"]
        require(process.get("command") == expected[role], f"native-F4 {role} command changed")
        for suffix in ("intent.json", "metrics.json", "stdout", "stderr"):
            expected_files.add(f"{role}.{suffix}")
        metrics, _ = phase_b.read_json(evidence / f"{role}.metrics.json", f"{role} metrics")
        intent, _ = phase_b.read_json(evidence / f"{role}.intent.json", f"{role} intent")
        require(
            all(
                metrics.get(key) == process.get(key)
                for key in ("command", "returncode", "timed_out", "orphan_group_terminated", "metrics")
            )
            and intent.get("command") == process.get("command")
            and intent.get("environment") == process.get("environment"),
            f"native-F4 {role} raw receipt changed",
        )
        require(
            process.get("returncode") == 0
            and process.get("timed_out") is False
            and process.get("orphan_group_terminated") is False,
            f"native-F4 {role} did not terminate cleanly",
        )
        resources = process["metrics"]
        for key in (
            "wall_seconds",
            "user_seconds",
            "system_seconds",
            "total_core_seconds",
            "single_core_seconds",
            "peak_rss_bytes",
        ):
            number = resources.get(key)
            require(
                type(number) in (int, float) and math.isfinite(number) and number >= 0,
                f"native-F4 {role} has invalid {key}",
            )
        require(
            math.isclose(
                resources["total_core_seconds"],
                resources["user_seconds"] + resources["system_seconds"],
                abs_tol=1e-9,
                rel_tol=0,
            )
            and resources["single_core_seconds"] == resources["total_core_seconds"]
            and resources.get("meter") == "fresh-process getrusage(RUSAGE_CHILDREN)",
            f"native-F4 {role} resource accounting changed",
        )
        for suffix in ("stdout", "stderr"):
            require(
                phase_b.sha256_file(evidence / f"{role}.{suffix}", f"{role} {suffix}")
                == process.get(f"{suffix}_sha256"),
                f"native-F4 {role} {suffix} changed",
            )
    require(
        {row["path"] for row in value["evidence_inventory"]} == expected_files,
        "native-F4 build evidence file set changed",
    )
    require(
        (evidence / "cargo-config.toml").read_bytes()
        == (evidence / "vendor.stdout").read_bytes(),
        "native-F4 vendor configuration changed",
    )
    require(value.get("resources") == resource_summary(processes), "build resources changed")
    return {"receipt_sha256": phase_b.sha256_bytes(raw), "receipt": value}


def copy_capsule(path: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(path, destination / "receipt.json")
    shutil.copytree(path.parent / "evidence", destination / "evidence")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()
    try:
        result = build(
            args.source,
            args.output,
            jobs=args.jobs,
            timeout=args.timeout,
        )
        print(
            phase_b.pretty_bytes(
                {
                    "receipt": str(args.output / "receipt.json"),
                    "resources": result["resources"],
                    "binaries": result["binaries"],
                }
            ).decode(),
            end="",
        )
    except (OSError, subprocess.CalledProcessError, phase_b.PhaseBError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
