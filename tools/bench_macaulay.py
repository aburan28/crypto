#!/usr/bin/env python3
"""Repeat the existing Rust Macaulay probe, retaining machine-readable timings.

Example:
  python3 tools/bench_macaulay.py --repeats 3 --n 31 --dim 16 --targets 4
"""
import argparse
import json
import os
import platform
import subprocess
import time
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--n", type=int, default=31)
    p.add_argument("--dim", type=int, default=16)
    p.add_argument("--targets", type=int, default=4)
    p.add_argument("--timeout", type=int, default=1800)
    p.add_argument("--output", default="macaulay-benchmark.jsonl")
    a = p.parse_args()
    if a.repeats < 1 or a.timeout < 1:
        p.error("repeats and timeout must be positive")
    root = Path(__file__).resolve().parents[1]
    cmd = ["cargo", "run", "--release", "--example", "macaulay_syzygy_probe", "--",
           "--n", str(a.n), "--dim", str(a.dim), "--targets", str(a.targets)]
    with open(a.output, "a", encoding="utf-8") as output:
        for iteration in range(a.repeats):
            started = time.perf_counter()
            try:
                proc = subprocess.run(cmd, cwd=root, text=True, capture_output=True,
                                      timeout=a.timeout, check=False)
                record = {"status": "ok" if proc.returncode == 0 else "error",
                          "exit_code": proc.returncode, "stdout": proc.stdout,
                          "stderr": proc.stderr}
            except subprocess.TimeoutExpired as exc:
                record = {"status": "timeout", "exit_code": None,
                          "stdout": str(exc.stdout or ""), "stderr": str(exc.stderr or "")}
            record.update({"iteration": iteration, "command": cmd, "elapsed_wall_s": time.perf_counter()-started,
                           "host": platform.node(), "platform": platform.platform(),
                           "cpu_count": os.cpu_count()})
            output.write(json.dumps(record)+"\n")
            output.flush()
            if record["status"] != "ok":
                raise SystemExit(f"probe {record['status']} at iteration {iteration}")
if __name__ == "__main__":
    main()
