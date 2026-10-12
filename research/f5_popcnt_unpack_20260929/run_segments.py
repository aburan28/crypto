#!/usr/bin/env python3
"""Replay one-thread or two-thread F5 seeds in short isolated reservations.

The first qualified attempt for each seed is selected without reading its
timing ratio. All failures and partial outputs remain in the artifact tree.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


SEEDS = ("frozen", "holdout_a", "holdout_b", "holdout_c")
MAX_ATTEMPTS = 3


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def read_last_json_line(path):
    try:
        return json.loads(path.read_text().splitlines()[-1])
    except (OSError, IndexError, json.JSONDecodeError):
        return None


def read_json(path):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def qualifies(code, isolation, paired):
    if code != 0 or not isolation or not paired:
        return False
    left = isolation.get("left_on_reserved", {}).get("user_threads")
    runs = paired.get("runs")
    return bool(
        isolation.get("exit_status") == 0
        and isinstance(left, list)
        and not left
        and isolation.get("contended_samples") == 0
        and paired.get("status") == "complete"
        and isinstance(runs, list)
        and len(runs) == 22
        and all(run.get("status") == "ok" for run in runs)
    )


def run_attempt(args, root, seed, attempt):
    directory = root / seed / f"attempt{attempt}"
    directory.mkdir(parents=True, exist_ok=True)
    isolation_path = directory / "isolation.jsonl"
    output_path = directory / "paired.json"
    lock = str(Path(os.environ["RUNNER_TEMP"]) / "f5-popcnt-unpack-bench.lock")
    home = os.environ["HOME"]
    command = [
        "sudo", "-n", "env",
        f"PATH={os.environ['PATH']}",
        f"HOME={home}",
        f"RUSTUP_HOME={os.environ.get('RUSTUP_HOME', home + '/.rustup')}",
        f"CARGO_HOME={os.environ.get('CARGO_HOME', home + '/.cargo')}",
        f"GITHUB_SHA={os.environ.get('GITHUB_SHA', '')}",
        f"GITHUB_HEAD_REF={os.environ.get('GITHUB_HEAD_REF', '')}",
        "python3", "tools/isolated_bench.py", "reserve",
        "--lock", lock,
        "--cpus", args.reserve,
        "--out", str(isolation_path),
        "--settle", "10", "--period", "2",
        "--label", f"f5-popcnt-unpack-t{args.threads}-{seed}-attempt{attempt}",
        "--", "taskset", "-c", args.pin,
        "python3", "research/f5_popcnt_unpack_20260929/confirm_segment.py",
        str(args.binary), str(output_path),
        "--pairs", "5", "--threads", str(args.threads), "--workload", seed,
    ]
    try:
        process = subprocess.run(command, capture_output=True, text=True, timeout=180, check=False)
        code = process.returncode
        stdout, stderr = process.stdout, process.stderr
    except subprocess.TimeoutExpired as exc:
        code = "timeout"
        stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else exc.stdout or ""
        stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else exc.stderr or ""
    (directory / "stdout.txt").write_text(stdout)
    (directory / "stderr.txt").write_text(stderr)
    isolation = read_last_json_line(isolation_path)
    paired = read_json(output_path)
    qualified = qualifies(code, isolation, paired)
    status = {
        "attempt": attempt,
        "seed": seed,
        "threads": args.threads,
        "command_exit": code,
        "isolation_exit": isolation.get("exit_status") if isolation else None,
        "left_user_threads": isolation.get("left_on_reserved", {}).get("user_threads") if isolation else None,
        "contended_samples": isolation.get("contended_samples") if isolation else None,
        "paired_status": paired.get("status") if paired else None,
        "calls": len(paired.get("runs", [])) if paired else 0,
        "qualified": qualified,
    }
    save(directory / "STATUS.json", status)
    return status


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--reserve", required=True)
    parser.add_argument("--pin", required=True)
    args = parser.parse_args()
    args.binary = args.binary.resolve()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    report = {"status": "running", "threads": args.threads, "seeds": {}}
    save(root / "INDEX.json", report)
    for seed in SEEDS:
        attempts = []
        selected = None
        for attempt in range(1, MAX_ATTEMPTS + 1):
            status = run_attempt(args, root, seed, attempt)
            attempts.append(status)
            report["seeds"][seed] = {"attempts": attempts, "selected_attempt": selected}
            save(root / "INDEX.json", report)
            if status["qualified"]:
                selected = attempt
                break
        report["seeds"][seed]["selected_attempt"] = selected
        save(root / "INDEX.json", report)
    report["status"] = (
        "qualified" if all(report["seeds"][seed]["selected_attempt"] for seed in SEEDS)
        else "incomplete_qualification"
    )
    save(root / "INDEX.json", report)
    return 0 if report["status"] == "qualified" else 1


if __name__ == "__main__":
    sys.exit(main())
