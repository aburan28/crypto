#!/usr/bin/env python3
"""Bounded macOS feasibility diagnostic for the current-source N83 K4096 base."""

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import psutil
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DYNAMIC = ROOT / "research/n83_root11_dynamic_20261016/run.py"
SOURCE = ROOT / "examples/koblitz_orbit_dlp_fast.rs"
SOURCE_SHA256 = "3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03"
K = 4096
CAP_SECONDS = 900
RSS_CAP_BYTES = 21 * 1024**3
MIN_AVAILABLE_BEFORE_BYTES = 21 * 1024**3
MIN_AVAILABLE_DURING_BYTES = 5 * 1024**3

spec = importlib.util.spec_from_file_location("n83_large_dynamic", DYNAMIC)
dynamic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dynamic)
dynamic.SOURCE_SHA256 = SOURCE_SHA256
dynamic.parent.EXPECTED[dynamic.parent.SOURCE] = SOURCE_SHA256


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def now():
    return datetime.now(timezone.utc).isoformat()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    output = args.output.resolve()
    hashes = dynamic.frozen_inputs()
    if digest(SOURCE) != SOURCE_SHA256:
        raise ValueError("current-source producer hash differs")
    if output.exists():
        raise ValueError("output already exists")
    output.mkdir(parents=True)
    folder = output / "runs" / "k4096root11dynamic"
    folder.mkdir(parents=True)
    memory = psutil.virtual_memory()
    write(output / "host.json", {
        "schema": "n83.large-local-host/v1", "tier": "local_diagnostic",
        "utc": now(), "platform": dict(zip(
            ("sysname", "nodename", "release", "version", "machine"), os.uname())),
        "physical_cores": psutil.cpu_count(logical=False),
        "logical_cpus": psutil.cpu_count(logical=True),
        "memory_total_bytes": memory.total,
        "memory_available_before_bytes": memory.available,
        "memory_available_floor_bytes": MIN_AVAILABLE_DURING_BYTES,
        "rss_cap_bytes": RSS_CAP_BYTES, "cap_seconds": CAP_SECONDS,
    })
    write(output / "source.json", {
        "schema": "n83.large-local-source/v1", "hashes": hashes,
        "source_sha256": SOURCE_SHA256, "binary_sha256": digest(binary),
        "runner_sha256": digest(Path(__file__)), "k": K,
        "layout": "root11_dynamic", "target_sha256": digest(dynamic.cold.TARGETS),
    })
    if memory.available < MIN_AVAILABLE_BEFORE_BYTES:
        write(output / "preflight.json", {"status": "REJECTED_MEMORY_HEADROOM"})
        return 2
    write(output / "preflight.json", {"status": "PASS_MEMORY_HEADROOM"})
    command = [str(binary), f"construct:83:1:{K}",
               str(dynamic.cold.TARGETS), "29",
               str(folder / "targets_result.jsonl")]
    environment = os.environ.copy()
    environment.update(dynamic.environment(folder, "root11_dynamic"))
    started = now()
    t0 = time.monotonic()
    peak = 0
    reason = None
    with (folder / "stdout.jsonl").open("wb") as stdout, (
            folder / "stderr.log").open("wb") as stderr:
        child = psutil.Popen(command, cwd=ROOT, env=environment,
                             stdout=stdout, stderr=stderr)
        while child.poll() is None:
            try:
                peak = max(peak, child.memory_info().rss)
            except psutil.NoSuchProcess:
                break
            available = psutil.virtual_memory().available
            elapsed = time.monotonic() - t0
            if peak > RSS_CAP_BYTES:
                reason = "STOPPED_RSS_GUARD"
            elif available < MIN_AVAILABLE_DURING_BYTES:
                reason = "STOPPED_HOST_HEADROOM_GUARD"
            elif elapsed > CAP_SECONDS:
                reason = "TIMEOUT"
            if reason:
                child.send_signal(signal.SIGTERM)
                try:
                    child.wait(timeout=5)
                except psutil.TimeoutExpired:
                    child.kill()
                break
            time.sleep(0.25)
        code = child.wait()
    usage = {
        "status": reason or ("PASS_PRODUCER" if code == 0 else "FAILED"),
        "exit_code": code, "started_utc": started, "finished_utc": now(),
        "wall_seconds": time.monotonic() - t0, "peak_rss_bytes": peak,
        "memory_available_after_bytes": psutil.virtual_memory().available,
        "major_faults": None, "cgroup_before": {}, "cgroup_after": {},
        "command": command, "resource_tier": "local_diagnostic",
    }
    write(folder / "usage.json", usage)
    print(json.dumps(usage, sort_keys=True))
    return 0 if usage["status"] == "PASS_PRODUCER" else 1


if __name__ == "__main__":
    raise SystemExit(main())
