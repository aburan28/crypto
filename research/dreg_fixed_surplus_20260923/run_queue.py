#!/usr/bin/env python3
"""Run the remaining ladder work one process per draw, restartably.

    python3 run_queue.py BINARY WORKDIR [--parallel 2]

Each job writes WORKDIR/<job>.jsonl and .log.  A job whose .jsonl already
holds a finished line (one with an "outcome") is skipped, so after a
container restart the same command resumes where the finished work ends;
only jobs in flight are re-run.  Order is the addendum's: the primary cell's
draws first, then (15, 5)'s, then the controls.
"""
import argparse
import json
import os
import subprocess
import time
from pathlib import Path

SEED = "20260928"
JOBS = ([("13:5:6", f"u{k}") for k in range(4)]
        + [("15:5:6", f"u{k}") for k in range(4)]
        + [("13:5:6", "control"), ("15:5:6", "control"), ("13:4:6", "control")])


def finished(path):
    try:
        return any("outcome" in json.loads(l) for l in path.read_text().splitlines() if l.strip())
    except (OSError, ValueError):
        return False


def command(binary, cell, job):
    base = ["timeout", "96h", binary, "--cells", cell, "--ffd-max", "5", "--seed", SEED]
    if job == "control":
        return base + ["--control-only", "--controls", "1"]
    return base + ["--unsat-index", job[1:]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("binary")
    ap.add_argument("workdir")
    ap.add_argument("--parallel", type=int, default=2)
    args = ap.parse_args()
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, F4_F2_MAX_ROWS="50000000", F4_F2_MAX_COLS="50000000")
    pending = [(c, j) for c, j in JOBS if not finished(work / f"{c.replace(':', '-')}.{j}.jsonl")]
    running = []
    with open(work / "queue.log", "a") as qlog:
        qlog.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} start, pending {pending}\n")
        qlog.flush()
        while pending or running:
            while pending and len(running) < args.parallel:
                cell, job = pending.pop(0)
                stem = work / f"{cell.replace(':', '-')}.{job}"
                out, err = open(f"{stem}.jsonl", "w"), open(f"{stem}.log", "w")
                p = subprocess.Popen(command(args.binary, cell, job), stdout=out, stderr=err, env=env)
                running.append((p, cell, job, time.time()))
                qlog.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} launch {cell} {job} pid {p.pid}\n")
                qlog.flush()
            time.sleep(30)
            for item in list(running):
                p, cell, job, t0 = item
                if p.poll() is not None:
                    running.remove(item)
                    qlog.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} done {cell} {job} "
                               f"exit {p.returncode} wall {time.time() - t0:.0f}s\n")
                    qlog.flush()
        qlog.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} queue empty\n")


if __name__ == "__main__":
    main()
