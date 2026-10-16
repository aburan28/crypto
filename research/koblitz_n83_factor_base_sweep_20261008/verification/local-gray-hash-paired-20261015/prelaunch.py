#!/usr/bin/env python3
"""Retain a fail-closed admission receipt before the local N83 paired gate."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
PAIR = HERE.parent / "k600-gray-hash-paired-protocol-20261014"
sys.path.insert(0, str(PAIR))
import preflight  # noqa: E402
import supervisor  # noqa: E402


IMAGE = "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31"
MAX_LOAD_PER_LOGICAL_CPU = 1.0
MIN_IDLE_PERCENT = 30.0


def command(argv: list[str], timeout: int = 30) -> dict:
    started = time.monotonic()
    try:
        result = subprocess.run(argv, capture_output=True, text=True,
                                timeout=timeout, check=False)
        return {"argv": argv, "exit_code": result.returncode,
                "stdout": result.stdout, "stderr": result.stderr,
                "wall_seconds_L0": time.monotonic() - started}
    except (OSError, subprocess.SubprocessError) as error:
        return {"argv": argv, "exit_code": None, "stdout": "",
                "stderr": f"{type(error).__name__}: {error}",
                "wall_seconds_L0": time.monotonic() - started}


def classify(top: dict, load: dict, cores: dict, containers: dict,
             image: dict, cgroup: dict, system_free: int) -> tuple[str, dict]:
    observations = (top, load, cores, containers, image, cgroup)
    if any(row["exit_code"] != 0 for row in observations):
        return "PRODUCER_FAILURE_PRELAUNCH_PROBE", {}
    cpu = re.search(r"CPU usage: .*?([0-9]+(?:\.[0-9]+)?)% idle", top["stdout"])
    load_numbers = re.findall(r"[0-9]+(?:\.[0-9]+)?", load["stdout"])
    try:
        idle = float(cpu.group(1)) if cpu else None
        load_1m = float(load_numbers[0])
        logical_cpus = int(cores["stdout"].strip())
    except (ValueError, IndexError, AttributeError):
        return "PRODUCER_FAILURE_PRELAUNCH_PARSE", {}
    if (idle is None or logical_cpus <= 0
            or image["stdout"].strip() != f"{IMAGE} arm64 linux"
            or cgroup["stdout"].splitlines() != ["4294967296", "0"]):
        return "PRODUCER_FAILURE_PRELAUNCH_IDENTITY", {}
    other_containers = [row for row in containers["stdout"].splitlines()
                        if row and row != "conductor-db"]
    measured = {"idle_percent": idle, "load_1m": load_1m,
                "logical_cpus": logical_cpus,
                "other_containers": other_containers,
                "system_free_bytes": system_free}
    if (idle < MIN_IDLE_PERCENT or load_1m > MAX_LOAD_PER_LOGICAL_CPU * logical_cpus
            or other_containers or system_free < 10 * (1 << 30)):
        return "HOST_BUSY_NO_SOLVER_LAUNCH", measured
    return "PASS_LOCAL_LINUX_ADMISSION", measured


def run(auth: Path, archives: Path, output: Path) -> dict:
    if output.exists():
        raise ValueError("prelaunch receipt path must be fresh")
    plan = preflight.make_plan(archives)
    plan_sha256 = hashlib.sha256(preflight.PLAN.read_bytes()).hexdigest()
    if preflight.PLAN.read_bytes() != (json.dumps(plan, indent=2, sort_keys=True) + "\n").encode():
        raise ValueError("frozen paired plan replay changed")
    record = json.loads(auth.read_text())
    run_output = Path(record["output_dir"])
    error = supervisor.authorization_error(record, run_output, plan_sha256)
    if error:
        raise ValueError(error)
    if run_output.exists():
        raise ValueError("authorized run output already exists")
    probes = {
        "top": command(["top", "-l", "1", "-n", "0"], 30),
        "load": command(["sysctl", "-n", "vm.loadavg"]),
        "cores": command(["sysctl", "-n", "hw.logicalcpu"]),
        "containers": command(["docker", "ps", "--format", "{{.Names}}"]),
        "image": command(["docker", "image", "inspect", IMAGE,
                          "--format", "{{.Id}} {{.Architecture}} {{.Os}}"]),
        "cgroup": command(["docker", "run", "--rm", "--pull", "never",
                           "--platform", "linux/arm64", "--network", "none",
                           "--cpus", "1", "--memory", "4096m",
                           "--memory-swap", "4096m", IMAGE, "sh", "-c",
                           "cat /sys/fs/cgroup/memory.max /sys/fs/cgroup/memory.swap.max"]),
    }
    free_bytes = shutil.disk_usage("/private/tmp").free
    status, measured = classify(*(probes[key] for key in
                                 ("top", "load", "cores", "containers", "image", "cgroup")),
                                free_bytes)
    result = {"schema": "n83.local-paired-prelaunch/v1", "status": status,
              "observed_at_utc": datetime.now(timezone.utc).isoformat(),
              "plan_sha256": plan_sha256,
              "authorization_sha256": hashlib.sha256(auth.read_bytes()).hexdigest(),
              "archive_sha256": {p.name: preflight.sha256(p) for p in sorted(archives.glob("*.jsonl.gz"))},
              "source_commit": command(["git", "rev-parse", "HEAD"])["stdout"].strip(),
              "probes": probes, "measured": measured,
              "n83_solver_launched": False, "conservative_charge_seconds": 0}
    with output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--archives", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.authorization, args.archives, args.output)
    print(json.dumps({"status": result["status"], "measured": result["measured"]}, sort_keys=True))
    return 0 if result["status"] == "PASS_LOCAL_LINUX_ADMISSION" else 1


if __name__ == "__main__":
    raise SystemExit(main())
