#!/usr/bin/env python3
"""Run one fresh-process control and retain its resource/provenance receipt."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--env", action="append", default=[])
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("missing command after --")

    overrides: dict[str, str] = {}
    for item in args.env:
        key, separator, value = item.partition("=")
        if not separator or not key:
            parser.error(f"invalid --env value: {item!r}")
        overrides[key] = value

    args.out.mkdir(parents=True, exist_ok=False)
    stdout_path = args.out / "stdout.txt"
    stderr_path = args.out / "stderr.txt"
    executable = Path(command[0]).resolve()
    if not executable.is_file():
        parser.error(f"executable does not exist: {executable}")

    environment = os.environ.copy()
    environment.update(overrides)
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started_unix_ns = time.time_ns()
    started = time.perf_counter_ns()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        completed = subprocess.run(
            command,
            check=False,
            env=environment,
            stdout=stdout,
            stderr=stderr,
        )
    wall_ms = (time.perf_counter_ns() - started) / 1_000_000
    ended_unix_ns = time.time_ns()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)

    receipt = {
        "schema_version": "1.0",
        "task_id": args.task_id,
        "command": command,
        "environment_overrides": dict(sorted(overrides.items())),
        "executable": str(executable),
        "executable_sha256": sha256(executable),
        "started_unix_ns": started_unix_ns,
        "ended_unix_ns": ended_unix_ns,
        "wall_ms": wall_ms,
        "child_user_ms": (after.ru_utime - before.ru_utime) * 1000,
        "child_system_ms": (after.ru_stime - before.ru_stime) * 1000,
        "peak_rss_bytes": after.ru_maxrss,
        "exit_code": completed.returncode,
        "stdout_bytes": stdout_path.stat().st_size,
        "stdout_sha256": sha256(stdout_path),
        "stderr_bytes": stderr_path.stat().st_size,
        "stderr_sha256": sha256(stderr_path),
    }
    (args.out / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
