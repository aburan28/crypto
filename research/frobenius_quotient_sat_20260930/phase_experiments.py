#!/usr/bin/env python3
"""Phase-aware SAT experiment matrix, reusing the independently verified runner.

Dry run by default. Pass --execute to run bounded jobs.
Example: python3 phase_experiments.py --execute --fields 13 --budget 30 --targets 2
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
FORMS = ("subgroup_xor", "line_subgroup", "specialized_xor")
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--execute", action="store_true")
    p.add_argument("--fields", default="13,19")
    p.add_argument("--budget", type=float, default=30.0)
    p.add_argument("--targets", type=int, default=2)
    p.add_argument("--seed", type=int, default=20261010)
    p.add_argument("--procs", type=int, default=1)
    p.add_argument("--outdir", type=Path, default=HERE / "phase_results")
    a = p.parse_args()
    if a.budget <= 0 or a.targets < 1 or a.procs < 1:
        p.error("budget, targets and procs must be positive")
    fields = [int(s) for s in a.fields.split(",")]
    if any(n not in (13,19,29) for n in fields):
        p.error("supported fields: 13,19,29 (confirm Curve supports field before running)")
    # Preserve the same target generation seed for each pair of formulations.
    jobs = []
    for n in fields:
        for phases in (False, True):
            for form in FORMS:
                cmd = [sys.executable, str(HERE / "run.py"), "--n", str(n),
                       "--s", str({13:4,19:5,29:6}[n]), "--seed", str(a.seed),
                       "--targets", str(a.targets), "--budget", str(a.budget),
                       "--procs", str(a.procs), "--forms", form,
                       "--out", str(a.outdir / f"n{n}_{'phase' if phases else 'direct'}_{form}.jsonl")]
                if phases:
                    cmd.append("--phases")
                jobs.append((n, phases, form, cmd))
    for n, phases, form, cmd in jobs:
        print(" ".join(cmd), flush=True)
        if not a.execute:
            continue
        a.outdir.mkdir(parents=True, exist_ok=True)
        t0 = time.perf_counter()
        # run.py enforces per-solver budget. This outer limit includes model
        # construction and all target/formulation attempts.
        timeout = max(60, int(2 * a.targets * (a.budget + 30) / max(1, a.procs) + 120))
        try:
            r = subprocess.run(cmd, cwd=HERE, text=True, capture_output=True, timeout=timeout)
            status = "ok" if r.returncode == 0 else "error"
            stdout, stderr, code = r.stdout, r.stderr, r.returncode
        except subprocess.TimeoutExpired as exc:
            status, code = "timeout", None
            stdout, stderr = str(exc.stdout or ""), str(exc.stderr or "")
        receipt = {"n":n,"phases":phases,"form":form,"status":status,
                   "returncode":code,"elapsed_wall_s":time.perf_counter()-t0,
                   "command":cmd,"platform":platform.platform(),
                   "machine":platform.machine(),
                   "stdout_sha256":hashlib.sha256(stdout.encode()).hexdigest(),
                   "stderr":stderr[-2000:]}
        with (a.outdir / "receipts.jsonl").open("a") as f:
            f.write(json.dumps(receipt)+"\n")
        if status != "ok":
            print(f"WARNING {status}: {n} {phases} {form}", file=sys.stderr)
if __name__ == "__main__":
    main()
