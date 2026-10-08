#!/usr/bin/env python3
"""Matched msolve thread scaling on frozen radical-image inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import time


INPUTS = {
    "prime-b11-t0": "prime-holdout-b11-t0-radical_symmetric.msolve",
    "prime-b11-t2": "prime-holdout-b11-t2-radical_symmetric.msolve",
    "binary-k3-t0": "binary-k3-t0-radical_symmetric.msolve",
}
THREADS = [1, 2, 4, 8, 16]


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def leading_hash(stdout):
    start, end = stdout.rfind("["), stdout.rfind("]")
    if start < 0 or end < start:
        raise ValueError("missing leading ideal")
    monomials = sorted(item.strip() for item in stdout[start + 1:end].split(",")
                       if item.strip())
    return hashlib.sha256("\n".join(monomials).encode()).hexdigest(), len(monomials)


def run(args):
    frozen = Path(args.run).resolve()
    output = Path(args.output).resolve()
    if output.exists():
        raise SystemExit(f"refusing to overwrite {output}")
    raw = output / "raw"
    raw.mkdir(parents=True)
    rows = []
    for input_name, filename in INPUTS.items():
        input_path = frozen / "inputs" / filename
        for repetition in range(args.repetitions):
            order = THREADS[repetition:] + THREADS[:repetition]
            for threads in order:
                command = ["msolve", "-f", str(input_path), "-g", "1", "-v", "0",
                           "-l", "2", "-t", str(threads), "--random-seed", "0"]
                started = time.perf_counter()
                completed = subprocess.run(command, capture_output=True, text=True,
                                           timeout=args.timeout)
                wall = time.perf_counter() - started
                completed.check_returncode()
                digest, terms = leading_hash(completed.stdout)
                stem = f"{input_name}-t{threads}-r{repetition}"
                stdout_path, stderr_path = raw / f"{stem}.stdout", raw / f"{stem}.stderr"
                stdout_path.write_text(completed.stdout)
                stderr_path.write_text(completed.stderr)
                row = {
                    "input": input_name,
                    "input_path": str(input_path.relative_to(frozen)),
                    "input_sha256": sha256_file(input_path),
                    "threads": threads,
                    "repetition": repetition,
                    "wall_seconds": wall,
                    "leading_ideal_sha256": digest,
                    "leading_terms": terms,
                    "stdout": str(stdout_path.relative_to(output)),
                    "stderr": str(stderr_path.relative_to(output)),
                    "stdout_sha256": sha256_file(stdout_path),
                    "stderr_sha256": sha256_file(stderr_path),
                }
                rows.append(row)
                print(input_name, threads, repetition, f"{wall:.4f}", flush=True)

    summary = []
    for input_name in INPUTS:
        hashes = {row["leading_ideal_sha256"] for row in rows if row["input"] == input_name}
        assert len(hashes) == 1
        baseline = statistics.median(row["wall_seconds"] for row in rows
                                     if row["input"] == input_name and row["threads"] == 1)
        for threads in THREADS:
            median = statistics.median(row["wall_seconds"] for row in rows
                                       if row["input"] == input_name and row["threads"] == threads)
            summary.append({
                "input": input_name,
                "threads": threads,
                "median_seconds": median,
                "speedup_over_one_thread": baseline / median,
                "leading_ideal_match": True,
            })
    (output / "processes.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))
    result = {
        "status": "complete",
        "classification": "wall-time engineering diagnostic",
        "frozen_solver_run": os.path.relpath(frozen, output),
        "repetitions": args.repetitions,
        "threads": THREADS,
        "msolve_options": ["-g", "1", "-v", "0", "-l", "2", "--random-seed", "0"],
        "summary": summary,
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    manifest = {
        "status": "complete",
        "source_sha256": sha256_file(Path(__file__)),
        "processes": len(rows),
        "all_leading_ideals_match": True,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def verify(args):
    output = Path(args.output).resolve()
    manifest = json.loads((output / "manifest.json").read_text())
    rows = [json.loads(line) for line in (output / "processes.jsonl").read_text().splitlines()]
    summary = json.loads((output / "summary.json").read_text())
    assert manifest["source_sha256"] == sha256_file(Path(__file__))
    assert manifest["processes"] == len(rows) == len(INPUTS) * len(THREADS) * 3
    for row in rows:
        frozen = (output / summary["frozen_solver_run"]).resolve()
        assert sha256_file(frozen / row["input_path"]) == row["input_sha256"]
        assert sha256_file(output / row["stdout"]) == row["stdout_sha256"]
        assert sha256_file(output / row["stderr"]) == row["stderr_sha256"]
    for input_name in INPUTS:
        assert len({row["leading_ideal_sha256"] for row in rows
                    if row["input"] == input_name}) == 1
    print(json.dumps({"status": "pass", "processes": len(rows)}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    runner = subparsers.add_parser("run")
    runner.add_argument("--run", required=True)
    runner.add_argument("--output", required=True)
    runner.add_argument("--repetitions", type=int, default=3)
    runner.add_argument("--timeout", type=float, default=120)
    checker = subparsers.add_parser("verify")
    checker.add_argument("--output", required=True)
    args = parser.parse_args()
    (run if args.command == "run" else verify)(args)


if __name__ == "__main__":
    main()
