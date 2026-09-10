#!/usr/bin/env python3
"""Build and meter the one frozen, overprovisioned WDSat Phase-B binary."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess

import run_koblitz_blind_pdp_phase_b as phase_b


BUILD_RECEIPT_SCHEMA = "koblitz_pdp_phase_b_wdsat_build_receipt.v1"


def git_state(source: Path) -> tuple[str, list[str]]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=source,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"],
        cwd=source,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.splitlines()
    return commit, status


def parse_config(config: str) -> dict[str, int]:
    macros: dict[str, int] = {}
    for line in config.splitlines():
        match = re.fullmatch(r"\s*#define\s+(\S+)\s+(\d+)\s*", line)
        if match is not None:
            macros[match.group(1)] = int(match.group(2))
    names = {
        "max_anf_id": "__MAX_ANF_ID__",
        "max_degree": "__MAX_DEGREE__",
        "max_id": "__MAX_ID__",
        "max_buffer_size": "__MAX_BUFFER_SIZE__",
        "max_eq": "__MAX_EQ__",
        "max_eq_size": "__MAX_EQ_SIZE__",
        "max_xeq": "__MAX_XEQ__",
        "max_xeq_size": "__MAX_XEQ_SIZE__",
    }
    try:
        return {field: macros[macro] for field, macro in names.items()}
    except KeyError as error:
        raise phase_b.PhaseBError(f"frozen WDSat config lacks {error.args[0]}") from error


def successful_process(record: dict, role: str) -> dict:
    if record["timed_out"] or record["returncode"] != 0:
        raise phase_b.PhaseBError(
            f"WDSat {role} failed: returncode={record['returncode']} timed_out={record['timed_out']}"
        )
    metrics = record["metrics"]
    return {
        "role": role,
        "command": record["command"],
        "returncode": record["returncode"],
        "timed_out": record["timed_out"],
        "orphan_group_terminated": record["orphan_group_terminated"],
        "stdout_sha256": record["stdout_sha256"],
        "stderr_sha256": record["stderr_sha256"],
        "metrics": {
            name: metrics[name]
            for name in (
                "wall_seconds",
                "user_seconds",
                "system_seconds",
                "total_core_seconds",
                "single_core_seconds",
                "peak_rss_bytes",
                "meter",
            )
        },
    }


def build(
    protocol_path: Path,
    source: Path,
    config_path: Path,
    output: Path,
    meter: Path,
) -> dict:
    protocol, _ = phase_b.read_json(protocol_path, "Phase-B protocol")
    phase_b.validate_protocol(protocol)
    policy = protocol["wdsat_build"]
    if output.exists() or output.is_symlink():
        raise phase_b.PhaseBError(f"WDSat build output must be new: {output}")
    source_metadata = source.lstat()
    if stat.S_ISLNK(source_metadata.st_mode) or not stat.S_ISDIR(source_metadata.st_mode):
        raise phase_b.PhaseBError("WDSat source must be a real directory")
    if not (source / "src/makefile").is_file():
        raise phase_b.PhaseBError("WDSat source lacks src/makefile")
    commit, status = git_state(source)
    if commit != policy["source_commit"] or status:
        raise phase_b.PhaseBError(
            f"WDSat source must be clean commit {policy['source_commit']}; got {commit} with {status}"
        )
    config_bytes = phase_b.regular_file_bytes(config_path, "frozen WDSat config")
    if phase_b.sha256_bytes(config_bytes) != policy["frozen_config_sha256"]:
        raise phase_b.PhaseBError("WDSat config bytes differ from the frozen protocol")
    try:
        config_text = config_bytes.decode()
    except UnicodeDecodeError as error:
        raise phase_b.PhaseBError("WDSat config is not UTF-8 text") from error
    if parse_config(config_text) != policy["limits"]:
        raise phase_b.PhaseBError("WDSat config macros differ from the frozen protocol limits")
    output.mkdir(parents=False)
    build_root = output / "build"
    build_root.mkdir()
    environment = phase_b.safe_child_environment()
    markers = phase_b.forbidden_material(protocol)
    timeout = protocol["execution"]["per_process_watchdog_seconds"]
    cp = Path(shutil.which("cp") or "/bin/cp").resolve()
    make = Path(shutil.which("make") or "/usr/bin/make").resolve()
    processes = []
    copied = phase_b.run_metered(
        role="source-copy",
        command=[str(cp), "-R", str((source / "src").resolve()), str((build_root / "src").resolve())],
        cwd=output,
        task_root=output,
        input_paths=[],
        timeout=timeout,
        meter=meter,
        environment=environment,
        markers=markers,
    )
    processes.append(successful_process(copied, "source-copy"))
    installed = phase_b.run_metered(
        role="config-install",
        command=[str(cp), str(config_path.resolve()), str((build_root / "src/config.h").resolve())],
        cwd=output,
        task_root=output,
        input_paths=[],
        timeout=timeout,
        meter=meter,
        environment=environment,
        markers=markers,
    )
    processes.append(successful_process(installed, "config-install"))
    if phase_b.regular_file_bytes(build_root / "src/config.h", "installed WDSat config") != config_bytes:
        raise phase_b.PhaseBError("installed WDSat config differs from its frozen source")
    source_inventory = phase_b.all_regular_inventory(build_root / "src")
    phase_b.write_json_new(
        output / "configured-source-inventory.json",
        {
            "schema": "koblitz_pdp_phase_b_wdsat_source_inventory.v1",
            "source_commit": commit,
            "config_sha256": policy["frozen_config_sha256"],
            "inventory": source_inventory,
            "inventory_sha256": phase_b.canonical_sha256(source_inventory),
        },
    )
    cleaned = phase_b.run_metered(
        role="clean",
        command=[str(make), "-C", "src", "clean"],
        cwd=build_root,
        task_root=output,
        input_paths=[build_root / "src/config.h"],
        timeout=timeout,
        meter=meter,
        environment=environment,
        markers=markers,
    )
    processes.append(successful_process(cleaned, "clean"))
    compiled = phase_b.run_metered(
        role="build",
        command=[str(make), "-C", "src"],
        cwd=build_root,
        task_root=output,
        input_paths=[build_root / "src/config.h"],
        timeout=timeout,
        meter=meter,
        environment=environment,
        markers=markers,
    )
    processes.append(successful_process(compiled, "build"))
    built_binary = build_root / "wdsat_solver"
    phase_b.regular_file_bytes(built_binary, "built WDSat binary")
    final_binary = output / "wdsat_solver"
    copied_binary = phase_b.run_metered(
        role="binary-copy",
        command=[str(cp), str(built_binary.resolve()), str(final_binary.resolve())],
        cwd=output,
        task_root=output,
        input_paths=[built_binary],
        timeout=timeout,
        meter=meter,
        environment=environment,
        markers=markers,
    )
    processes.append(successful_process(copied_binary, "binary-copy"))
    os.chmod(final_binary, 0o755)
    binary_hash = phase_b.sha256_file(final_binary, "final WDSat binary")
    receipt = {
        "schema": BUILD_RECEIPT_SCHEMA,
        "status": "completed",
        "source_commit": commit,
        "source_dirty": False,
        "config": config_text,
        "config_sha256": policy["frozen_config_sha256"],
        "binary_sha256": binary_hash,
        "limits": policy["limits"],
        "build_processes": processes,
        "claim_boundary": "Clean source-copy, config-install, clean, build, and binary-copy costs only; exact blind-target export sizing is charged and checked separately before solving",
    }
    phase_b.validate_wdsat_build_receipt(
        receipt,
        phase_b.pretty_bytes(receipt),
        phase_b.executable_identity(final_binary, "final WDSat binary"),
        protocol,
    )
    phase_b.write_json_new(output / "receipt.json", receipt)
    inventory = phase_b.all_regular_inventory(output, {"build-seal.json"})
    build_seal = {
        "schema": "koblitz_pdp_phase_b_wdsat_build_seal.v1",
        "status": "build_frozen",
        "source_commit": commit,
        "config_sha256": policy["frozen_config_sha256"],
        "binary_path": "wdsat_solver",
        "binary_sha256": binary_hash,
        "receipt_path": "receipt.json",
        "receipt_sha256": phase_b.sha256_file(output / "receipt.json", "WDSat build receipt"),
        "inventory": inventory,
        "inventory_sha256": phase_b.canonical_sha256(inventory),
    }
    phase_b.write_json_new(output / "build-seal.json", build_seal)
    return build_seal


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=phase_b.DEFAULT_PROTOCOL)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--meter", type=Path, default=phase_b.DEFAULT_METER)
    args = parser.parse_args()
    try:
        protocol, _ = phase_b.read_json(args.protocol.resolve(), "Phase-B protocol")
        phase_b.validate_protocol(protocol)
        config = (
            args.config.resolve()
            if args.config is not None
            else (phase_b.REPO / protocol["wdsat_build"]["frozen_config_path"]).resolve()
        )
        result = build(
            args.protocol.resolve(),
            args.source.resolve(),
            config,
            args.output.resolve(),
            args.meter.resolve(),
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (OSError, phase_b.PhaseBError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
