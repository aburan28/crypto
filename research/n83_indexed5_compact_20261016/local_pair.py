#!/usr/bin/env python3
"""Four-arm macOS cold diagnostic for opt-in compact exact-root slots."""

import argparse
import hashlib
import json
import os
import re
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/koblitz_orbit_dlp_fast.rs"
SOURCE_SHA256 = "e1b143ec8ba80425a1102b65c3e7c0b67a4b4e9a7e0bc116bcaf13e76c5651f9"
TARGETS = ROOT / "research/n83_cold_factor_base_20261016/targets.jsonl"
ANSWERS = ROOT / "research/n83_cold_factor_base_20261016/answers.jsonl"
K = 1152
ARMS = ("off_a", "compact_a", "compact_b", "off_b")
BASE_FLAGS = {
    "KIC_PACKED_ROOT_TABLE": "1", "KIC_PACKED_ROOT_TABLE_DENSE": "0",
    "KIC_PACKED_ROOT_TABLE_LOAD75": "0", "KIC_STATE28": "0",
    "KIC_STATE26": "0", "KIC_ROOT_ONLY22": "0", "KIC_ROOT_ONLY21": "0",
    "KIC_ROOT_ONLY11": "1", "KIC_INDEXED_ROOT5": "1",
    "KIC_INDEXED_ROOT4": "0", "KIC_ROOT_BLOOM": "0",
    "KIC_SYMMETRIC_S3": "1", "KIC_INDEXED5_LOAD75": "1",
    "KIC_INDEXED5_FIXED30": "0", "KIC_INDEXED4_FIXED31": "0",
    "KIC_SPARSE_CANONICAL": "1", "KIC_RANK_THREADS": "4",
    "KIC_TARGET_THREADS": "1",
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, item):
    path.write_text(json.dumps(item, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    output = args.output.resolve()
    if output.exists():
        raise ValueError("output directory already exists")
    if (sha256(SOURCE) != SOURCE_SHA256
            or sha256(TARGETS) !=
            "2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d"):
        raise ValueError("source or exposed public target differs")
    output.mkdir(parents=True)
    save(output / "source.json", {
        "schema": "n83.indexed5-compact-local-source/v1",
        "source_sha256": SOURCE_SHA256,
        "binary_sha256": sha256(binary),
        "target_sha256": sha256(TARGETS),
        "answer_sidecar_sha256": sha256(ANSWERS),
        "runner_sha256": sha256(Path(__file__)),
        "k": K, "rank_seed": 29, "arm_order": ARMS,
        "hardware_tier": "local_diagnostic",
    })
    for arm in ARMS:
        folder = output / "runs" / arm
        folder.mkdir(parents=True)
        flags = {**BASE_FLAGS, "KIC_DUMP_BASE": str(folder / "base_dump.jsonl")}
        compact = arm.startswith("compact")
        flags["KIC_INDEXED5_COMPACT27"] = "1" if compact else "0"
        env = os.environ.copy()
        env.update(flags)
        command = [str(binary), f"construct:83:1:{K}", str(TARGETS),
                   "29", str(folder / "targets_result.jsonl")]
        save(folder / "command.json", {
            "argv": command, "env_flags": {key: value for key, value in flags.items()
                                          if key.startswith("KIC_")},
            "compact27": compact,
        })
        start = time.monotonic()
        try:
            completed = subprocess.run(["/usr/bin/time", "-l", *command],
                                       env=env, cwd=ROOT, capture_output=True,
                                       text=True, timeout=900, check=False)
            (folder / "stdout.jsonl").write_text(completed.stdout)
            (folder / "stderr.txt").write_text(completed.stderr)
            match = re.search(r"(\d+)\s+maximum resident set size",
                              completed.stderr)
            telemetry_denied = (completed.returncode != 0
                                and "sysctl kern.clockrate: Operation not permitted"
                                in completed.stderr
                                and (folder / "targets_result.jsonl").exists())
            usage = {"status": ("PASS_PRODUCER" if completed.returncode == 0
                                else "INFRA_TELEMETRY_FAILURE" if telemetry_denied
                                else "PRODUCER_FAILURE"),
                     "exit_code": completed.returncode,
                     "outer_wall_seconds": time.monotonic() - start,
                     "peak_rss_bytes": int(match.group(1)) if match else None,
                     "tier": "local_diagnostic"}
        except subprocess.TimeoutExpired as error:
            (folder / "stdout.jsonl").write_text(
                error.stdout.decode() if isinstance(error.stdout, bytes)
                else error.stdout or "")
            (folder / "stderr.txt").write_text(
                error.stderr.decode() if isinstance(error.stderr, bytes)
                else error.stderr or "")
            usage = {"status": "TIMEOUT", "cap_seconds": 900,
                     "outer_wall_seconds": time.monotonic() - start,
                     "tier": "local_diagnostic"}
        save(folder / "usage.json", usage)
        print(f"{arm}: {usage['status']} after {usage['outer_wall_seconds']:.3f}s",
              flush=True)


if __name__ == "__main__":
    main()
