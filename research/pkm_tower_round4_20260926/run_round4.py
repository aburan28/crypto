#!/usr/bin/env python3
"""Round 4 of the PKM tower measurement (note section 13.4): the signature
engine (`--engine sig`, its default options) on systems that F4 finished in
round 2, with F4 re-run beside it under the same conditions.

    cargo build --release --example pkm_tower_pilot
    python3 research/pkm_tower_round4_20260926/run_round4.py

Each size is one process per engine, F4 first and then the signature engine,
each measuring the size's two targets. A process runs under the rounds' 14 GB
address-space cap, and its peak resident memory (`ru_maxrss` through `wait4`)
and exit status go to `memory.jsonl`. The rows go to
`runs/<engine cell>-N<N>.jsonl` and the example's stderr, with the trace of
every step, to the matching `.log`. `progress.txt` records when each process
started and ended. Once a signature system times out, the larger sizes of its
`m` are not run (section 13.4).
"""

import json
import os
import resource
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "runs")
BIN = os.environ.get(
    "BIN", os.path.join(HERE, "..", "..", "target", "release", "examples", "pkm_tower_pilot"))
P1 = "2013265921"
COMMON = ["--p", P1, "--kinds", "kummer", "--controls", "tower", "--planted", "0",
          "--random", "2", "--ladder-t", "none", "--budget", "1800",
          "--max-nnz", "2000000000", "--trace"]
# (signature cell, F4 cell (T for the tower engine), m, the flag bounding t at this m, sizes N)
CELLS = [
    ("G2", "T2", 2, "--max-t", [12, 14, 16, 18, 20]),
    ("G3", "T3", 3, "--max-t-m3", [9, 12, 15]),
    ("G4", "T4", 4, "--max-t-m4", [8, 12]),
]


def limit():
    cap = 14_000_000 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))


def progress(line):
    with open(os.path.join(HERE, "progress.txt"), "a") as f:
        f.write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " " + line + "\n")


def run(name, engine, m, t_flag, n):
    t = n // m
    args = [BIN, "--engine", engine, *COMMON, "--m", str(m), "--t-min", str(t), t_flag, str(t),
            "--out", os.path.join(OUT, name + ".jsonl")]
    progress(f"start {name}: {' '.join(args[1:])}")
    started = time.time()
    with open(os.path.join(OUT, name + ".log"), "w") as log:
        p = subprocess.Popen(args, stderr=log, stdout=subprocess.DEVNULL, preexec_fn=limit)
        _, status, usage = os.wait4(p.pid, 0)
    code = os.waitstatus_to_exitcode(status)
    progress(f"end {name} exit {code}")
    path = os.path.join(OUT, name + ".jsonl")
    rows = [json.loads(l) for l in open(path) if l.strip()] if os.path.exists(path) else []
    rows = [r for r in rows if "N" in r]
    rec = {
        "run": name, "engine": engine, "m": m, "N": n, "exit": code,
        "peak_rss_mb": round(usage.ru_maxrss / 1024, 1),
        "wall_s": round(time.time() - started, 2),
        "rows": len(rows),
        "timed_out": [r["timed_out"] for r in rows],
    }
    with open(os.path.join(HERE, "memory.jsonl"), "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)
    return code, rows


def main():
    os.makedirs(OUT, exist_ok=True)
    for sig_cell, f4_cell, m, t_flag, sizes in CELLS:
        for n in sizes:
            run(f"{f4_cell}-kummer-m{m}-p1-N{n}", "tower", m, t_flag, n)
            code, rows = run(f"{sig_cell}-kummer-m{m}-p1-N{n}", "sig", m, t_flag, n)
            if code != 0 or len(rows) != 2 or any(r["timed_out"] for r in rows):
                progress(f"{sig_cell} stops at N = {n}")
                break
    progress("round 4 done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
