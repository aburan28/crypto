#!/usr/bin/env python3
"""The gate of note section 14.2, with the context runs of section 14.3.

    cargo build --release --example pkm_tower_pilot
    python3 research/pkm_tower_round5_20260926/run_gate.py

Three signature-engine variants on round 4's systems, one process per variant
and size, each measuring the size's two targets:
- GP: steps by polynomial degree, position over term (the gate);
- GD: steps by polynomial degree, signature degree first (the gate);
- S4: round 4's engine, steps by signature degree, position over term
  (context: its D_lm).

The gate's systems (m = 3, N = 9 and 12) run first, then the context sizes.
Every process runs under the rounds' 14 GB address-space cap with a 600 s
budget per system. Its peak resident memory (`ru_maxrss` through `wait4`) and
exit status go to `memory.jsonl`, its rows to `runs/<cell>-...jsonl` and its
stderr, with the trace of every step, to the matching `.log`. `progress.txt`
records when each process started and ended. F4's rows are round 4's T cells.
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
COMMON = ["--engine", "sig", "--sig-rewrite", "ratio", "--p", P1, "--kinds", "kummer",
          "--controls", "tower", "--planted", "0", "--random", "2", "--ladder-t", "none",
          "--budget", "600", "--max-nnz", "2000000000", "--trace"]
VARIANTS = [
    ("GP", ["--sig-steps", "degree", "--sig-order", "pot"]),
    ("GD", ["--sig-steps", "degree", "--sig-order", "degree"]),
    ("S4", ["--sig-steps", "sugar", "--sig-order", "pot"]),
]
T_FLAG = {2: "--max-t", 3: "--max-t-m3", 4: "--max-t-m4"}
GATE = [(3, 9), (3, 12)]
CONTEXT = [(2, 12), (2, 14), (2, 16), (4, 8)]


def limit():
    cap = 14_000_000 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))


def progress(line):
    with open(os.path.join(HERE, "progress.txt"), "a") as f:
        f.write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " " + line + "\n")


def run(cell, flags, m, n):
    name = f"{cell}-kummer-m{m}-p1-N{n}"
    t = n // m
    args = [BIN, *COMMON, *flags, "--m", str(m), "--t-min", str(t), T_FLAG[m], str(t),
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
        "run": name, "cell": cell, "m": m, "N": n, "exit": code,
        "peak_rss_mb": round(usage.ru_maxrss / 1024, 1),
        "wall_s": round(time.time() - started, 2),
        "rows": len(rows),
        "timed_out": [r["timed_out"] for r in rows],
    }
    with open(os.path.join(HERE, "memory.jsonl"), "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    for m, n in GATE + CONTEXT:
        for cell, flags in VARIANTS:
            run(cell, flags, m, n)
    progress("gate done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
