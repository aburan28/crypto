#!/usr/bin/env python3
"""Construct each smaller N83 standard-subspace base once with exact native replay."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
BINARY_SHA256 = "56003eed67a803eda956fe6003c5f64ccdfda96345eac2b55d13bc0b58497cb8"
MIN_FREE_BYTES = 8 * 1024**3


def write_new_json(path: Path, value: dict) -> None:
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: construct_frontier.py BINARY NEW_SCRATCH_DIRECTORY")
    binary = Path(sys.argv[1]).resolve(strict=True)
    scratch = Path(sys.argv[2]).resolve()
    if hashlib.sha256(binary.read_bytes()).hexdigest() != BINARY_SHA256:
        raise ValueError("native exporter binary identity changed")
    source_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    ).stdout.strip()
    if len(source_commit) != 40:
        raise ValueError("source commit is not pinned")
    if shutil.disk_usage(scratch.parent).free < MIN_FREE_BYTES:
        raise ValueError("less than 8 GiB of scratch is free before construction")
    scratch.mkdir()
    progress = HERE / "construction.jsonl"
    rows = []
    env = dict(os.environ, RAYON_NUM_THREADS="1")
    with progress.open("x") as log:
        for dimension in range(1, 18):
            if shutil.disk_usage(scratch).free < MIN_FREE_BYTES:
                row = {"dimension": dimension, "status": "DEFER_SCRATCH_BELOW_8_GIB"}
                rows.append(row)
                log.write(json.dumps(row, sort_keys=True) + "\n")
                log.flush()
                os.fsync(log.fileno())
                break
            directory = scratch / f"d{dimension}"
            command = [str(binary), "construct", str(directory), str(dimension)]
            started = time.monotonic()
            try:
                result = subprocess.run(command, env=env, capture_output=True, timeout=180)
                status = "PASS_NATIVE_CONSTRUCT_AND_EXACT_REPLAY" if result.returncode == 0 else "CONSTRUCT_FAILURE"
                stdout, stderr, exit_code = result.stdout, result.stderr, result.returncode
            except subprocess.TimeoutExpired as error:
                status = "CONSTRUCT_TIMEOUT_180S"
                stdout, stderr, exit_code = error.stdout or b"", error.stderr or b"", None
            wall = time.monotonic() - started
            (scratch / f"d{dimension}.stdout").write_bytes(stdout)
            (scratch / f"d{dimension}.stderr").write_bytes(stderr)
            row = {
                "dimension": dimension,
                "status": status,
                "source_commit": source_commit,
                "binary_sha256": BINARY_SHA256,
                "argv": command,
                "rayon_num_threads": 1,
                "wall_seconds_L0": wall,
                "exit_code": exit_code,
                "stdout_sha256": hashlib.sha256(stdout).hexdigest(),
                "stderr_sha256": hashlib.sha256(stderr).hexdigest(),
                "stderr_text": stderr.decode("utf-8", errors="replace") if status != "PASS_NATIVE_CONSTRUCT_AND_EXACT_REPLAY" else "",
            }
            if status == "PASS_NATIVE_CONSTRUCT_AND_EXACT_REPLAY":
                manifest = json.loads((directory / "manifest.json").read_text())
                if json.loads(stdout) != manifest or manifest["source_dimension"] != dimension:
                    raise ValueError(f"dimension {dimension} stdout/manifest mismatch")
                row["manifest"] = manifest
            rows.append(row)
            log.write(json.dumps(row, sort_keys=True) + "\n")
            log.flush()
            os.fsync(log.fileno())
    write_new_json(HERE / "construction.json", {
        "schema": "n83.standard-subspace-frontier-construction/v1",
        "source_commit": source_commit,
        "binary_sha256": BINARY_SHA256,
        "requested_dimensions": list(range(1, 18)),
        "rows": rows,
        "selected_best_total_runtime": None,
    })
    print(json.dumps({"rows": len(rows), "passed": sum(row["status"].startswith("PASS") for row in rows)}, sort_keys=True))


if __name__ == "__main__":
    main()
