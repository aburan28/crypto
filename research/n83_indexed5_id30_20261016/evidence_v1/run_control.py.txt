#!/usr/bin/env python3
"""Frozen same-point old/ID30 N83 K1024 correctness control."""

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DYNAMIC = ROOT / "research/n83_root11_dynamic_20261016/run.py"
SOURCE = ROOT / "examples/koblitz_orbit_dlp_fast.rs"
SOURCE_SHA256 = "f14d9747cd24627d39d7cae4c9bef69ed2819544125ef15f31306ce5d75f3fc7"
K = 1024
CAP_SECONDS = 900
ARMS = ("legacy_id29", "opt_in_id30")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def utc():
    return datetime.now(timezone.utc).isoformat()


def dynamic_runner():
    spec = importlib.util.spec_from_file_location("n83_id30_dynamic", DYNAMIC)
    dynamic = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dynamic)
    dynamic.SOURCE_SHA256 = SOURCE_SHA256
    dynamic.parent.EXPECTED[dynamic.parent.SOURCE] = SOURCE_SHA256
    return dynamic


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    output = args.output.resolve()
    dynamic = dynamic_runner()
    hashes = dynamic.frozen_inputs()
    if digest(SOURCE) != SOURCE_SHA256 or output.exists():
        raise ValueError("source hash differs or output already exists")
    output.mkdir(parents=True)
    copied_binary = output / "koblitz_orbit_dlp_fast"
    shutil.copy2(binary, copied_binary)
    write(output / "source.json", {
        "schema": "n83.id30-control-source/v1", "hashes": hashes,
        "source_sha256": SOURCE_SHA256,
        "binary_sha256": digest(copied_binary),
        "runner_sha256": digest(Path(__file__)),
        "public_input_sha256": digest(dynamic.cold.TARGETS),
        "arm_order": ARMS, "k": K, "cap_seconds": CAP_SECONDS,
    })
    for name in ARMS:
        folder = output / "runs" / name
        folder.mkdir(parents=True)
        environment = os.environ.copy()
        environment.update(dynamic.environment(folder, "root11_dynamic"))
        environment["KIC_INDEXED5_ID30"] = "1" if name == "opt_in_id30" else "0"
        command = [str(copied_binary), f"construct:83:1:{K}",
                   str(dynamic.cold.TARGETS), "29",
                   str(folder / "targets_result.jsonl")]
        started = utc()
        t0 = time.monotonic()
        with (folder / "stdout.jsonl").open("wb") as stdout, (
                folder / "stderr.log").open("wb") as stderr:
            try:
                result = subprocess.run(command, cwd=ROOT, env=environment,
                                        stdout=stdout, stderr=stderr,
                                        timeout=CAP_SECONDS, check=False)
                code = result.returncode
                status = "PASS_PRODUCER" if code == 0 else "FAILED"
            except subprocess.TimeoutExpired:
                code = None
                status = "TIMEOUT"
        write(folder / "usage.json", {
            "status": status, "exit_code": code,
            "started_utc": started, "finished_utc": utc(),
            "wall_seconds": time.monotonic() - t0,
            "command": command,
            "algorithm_flags": {"KIC_INDEXED5_ID30": environment["KIC_INDEXED5_ID30"]},
            "hardware_tier": "local_correctness_control",
        })
        print(f"{name}: {status}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
