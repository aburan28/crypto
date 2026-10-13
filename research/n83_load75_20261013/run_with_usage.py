#!/usr/bin/env python3
"""Bound one direct child and record its kernel-reported peak resident set."""

import argparse
import json
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--limit-seconds", type=int, default=900)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or args.limit_seconds <= 0:
        parser.error("a command and positive time cap are required")

    started = time.monotonic()
    child = subprocess.Popen(command)
    timed_out = False
    try:
        try:
            exit_status = child.wait(timeout=args.limit_seconds)
        except subprocess.TimeoutExpired:
            timed_out = True
            child.kill()
            child.wait()
            exit_status = 124
    except BaseException:
        child.kill()
        child.wait()
        raise
    wall_seconds = time.monotonic() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    system = platform.system()
    peak_unit = "bytes" if system == "Darwin" else "KiB"
    receipt = {
        "schema": "n83.child-usage/v1",
        "command": command,
        "exit_status": exit_status,
        "timed_out": timed_out,
        "limit_seconds": args.limit_seconds,
        "wall_seconds": wall_seconds,
        "peak_rss": usage.ru_maxrss,
        "peak_rss_unit": peak_unit,
        "user_cpu_seconds": usage.ru_utime,
        "system_cpu_seconds": usage.ru_stime,
        "involuntary_context_switches": usage.ru_nivcsw,
        "major_faults": usage.ru_majflt,
    }
    Path(args.receipt).write_text(json.dumps(receipt, sort_keys=True) + "\n")
    return exit_status


if __name__ == "__main__":
    sys.exit(main())
