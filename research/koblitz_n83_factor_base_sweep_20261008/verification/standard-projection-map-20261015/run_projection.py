#!/usr/bin/env python3
"""Run one exact N83 projection-map construction under a hard outer wall."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import resource
import shutil
import subprocess
import time


HERE = Path(__file__).resolve().parent
BINARY_SHA256 = "d25688592fda9adfb71eb1da545daaecc160fb0cdf9c14fbcd7d80cab31934e0"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dimension", type=int)
    parser.add_argument("scratch", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--wall", type=int, default=300)
    args = parser.parse_args()
    if not 2 <= args.dimension <= 18 or not 1 <= args.wall <= 600:
        raise ValueError("dimension or wall cap out of bounds")
    binary = args.binary.resolve(strict=True)
    if hashlib.sha256(binary.read_bytes()).hexdigest() != BINARY_SHA256:
        raise ValueError("binary identity changed")
    scratch = args.scratch.resolve(strict=True)
    output = scratch / f"d{args.dimension}.jsonl"
    receipt_path = HERE / f"d{args.dimension}-launch.json"
    if output.exists() or receipt_path.exists():
        raise ValueError("dimension output or receipt already exists")
    commit = subprocess.run(["git", "rev-parse", "HEAD"], check=True,
                            capture_output=True, text=True).stdout.strip()
    command = [str(binary), "construct", str(output), str(args.dimension)]
    free_before = shutil.disk_usage(scratch).free
    started = time.monotonic()
    try:
        result = subprocess.run(command, capture_output=True, timeout=args.wall)
        timeout = False
        code, stdout, stderr = result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired as error:
        timeout = True
        code, stdout, stderr = None, error.stdout or b"", error.stderr or b""
    wall = time.monotonic() - started
    use = resource.getrusage(resource.RUSAGE_CHILDREN)
    (scratch / f"d{args.dimension}.stdout").write_bytes(stdout)
    (scratch / f"d{args.dimension}.stderr").write_bytes(stderr)
    receipt = {
        "schema": "n83.standard-projection-map-launch/v1",
        "source_commit": commit,
        "binary_sha256": BINARY_SHA256,
        "argv": command,
        "wall_cap_seconds": args.wall,
        "wall_seconds_L0": wall,
        "child_user_seconds_L0": use.ru_utime,
        "child_system_seconds_L0": use.ru_stime,
        "child_maxrss_bytes_L0": use.ru_maxrss,
        "scratch_free_before_bytes": free_before,
        "scratch_free_after_bytes": shutil.disk_usage(scratch).free,
        "timeout": timeout,
        "exit_code": code,
        "stdout_sha256": hashlib.sha256(stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(stderr).hexdigest(),
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest() if output.exists() else None,
        "result": json.loads(stdout) if code == 0 else None,
        "stderr_text": stderr.decode(errors="replace") if code != 0 else "",
        "selected_best_total_runtime": None,
    }
    with receipt_path.open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
