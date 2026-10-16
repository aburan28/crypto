#!/usr/bin/env python3
"""Retain an unambiguous exit receipt for the required release library suite."""

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
COMMAND = ["cargo", "test", "--release", "--lib", "--offline", "--", "--test-threads=1"]


def utc():
    return datetime.now(timezone.utc).isoformat()


def main():
    started = utc()
    with (HERE / "release-lib.log").open("x") as output:
        result = subprocess.run(COMMAND, cwd=ROOT, stdout=output,
                                stderr=subprocess.STDOUT, check=False)
    receipt = {
        "schema": "n83.release-library-validation/v1",
        "command": COMMAND,
        "cwd": str(ROOT),
        "runner_pid": os.getpid(),
        "started_utc": started,
        "finished_utc": utc(),
        "exit_code": result.returncode,
        "log_bytes": (HERE / "release-lib.log").stat().st_size,
    }
    with (HERE / "release-lib-exit.json").open("x") as output:
        json.dump(receipt, output, sort_keys=True, indent=2)
        output.write("\n")
    print(json.dumps(receipt, sort_keys=True), flush=True)
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
