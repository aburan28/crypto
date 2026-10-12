#!/usr/bin/env python3
"""Run the registered surplus control, one process per job, restartably.

    python3 run_grid.py BINARY WORKDIR [--parallel 2]

BINARY is the frozen `dreg_ladder` (PREREGISTRATION.md, "Harness").  Each
job writes WORKDIR/<job>.jsonl and .log.  A finished job is skipped, so after
a container restart the same command resumes; only jobs in flight re-run.
Order is the registration's: (7, 4) and (8, 4) whole, then (11, 5)'s four
unsatisfiable draws one process each.
"""
import argparse
import json
import os
import subprocess
import time
from pathlib import Path

SEED = "20260926"
JOBS = ([(c, "all") for c in ("7:4:7", "8:4:7")]
        + [("11:5:6", f"u{k}") for k in range(4)])


def finished(stem, job):
    try:
        if job == "all":  # the whole-cell mode prints its summary table last
            return "| n | ℓ | S |" in Path(f"{stem}.log").read_text()
        lines = Path(f"{stem}.jsonl").read_text().splitlines()
        return any("outcome" in json.loads(l) for l in lines if l.strip())
    except (OSError, ValueError):
        return False


def command(binary, cell, job):
    base = ["timeout", "96h", binary, "--cells", cell, "--ffd-max", "5", "--seed", SEED]
    if job == "all":
        return base + ["--unsat", "4", "--controls", "0"]
    return base + ["--unsat-index", job[1:]]


def stamp():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("binary")
    ap.add_argument("workdir")
    ap.add_argument("--parallel", type=int, default=2)
    args = ap.parse_args()
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, F4_F2_MAX_ROWS="50000000", F4_F2_MAX_COLS="50000000")
    stem = lambda cell, job: work / f"{cell.replace(':', '-')}.{job}"
    pending = [(c, j) for c, j in JOBS if not finished(stem(c, j), j)]
    running = []
    with open(work / "queue.log", "a") as qlog:
        qlog.write(f"{stamp()} start, pending {pending}\n")
        qlog.flush()
        while pending or running:
            while pending and len(running) < args.parallel:
                cell, job = pending.pop(0)
                s = stem(cell, job)
                out, err = open(f"{s}.jsonl", "w"), open(f"{s}.log", "w")
                p = subprocess.Popen(command(args.binary, cell, job), stdout=out, stderr=err, env=env)
                running.append((p, cell, job, time.time()))
                qlog.write(f"{stamp()} launch {cell} {job} pid {p.pid}\n")
                qlog.flush()
            time.sleep(10)
            for item in list(running):
                p, cell, job, t0 = item
                if p.poll() is not None:
                    running.remove(item)
                    qlog.write(f"{stamp()} done {cell} {job} exit {p.returncode} wall {time.time() - t0:.0f}s\n")
                    qlog.flush()
        qlog.write(f"{stamp()} queue empty\n")


if __name__ == "__main__":
    main()
