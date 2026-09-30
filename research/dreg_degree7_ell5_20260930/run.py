#!/usr/bin/env python3
"""Runner for the degree-7 study at ℓ = 5: one draw a process, one at a time.

    python3 research/dreg_degree7_ell5_20260930/run.py BINARY [JOB ...]

A JOB is `N-L-7.uK`, e.g. `10-5-7.u0`; with none, the registered order runs.
A job whose `runs/cell-<JOB>.jsonl` already holds a line is skipped, so a
restart loses only the draw in flight.  Each draw is `--unsat-index K
--d-min 7` at its committed cell's seed (PREREGISTRATION.md), so it is the
same subspace and target as the committed degree-6 row.

Writes `runs/cell-<JOB>.{jsonl,log,hwm}` and appends to `runs/queue.log`.
`.hwm` is the process's peak resident set (VmHWM, kB), sampled every 5 s.
`KIC_SPARSE_DENSE_BUDGET_MB` is the cell's registered budget (BUDGET_MB)
unless set in the environment, and is recorded per job; it moves where the
dense finish starts, never an outcome.
"""
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
# cell -> (--cells value, seed of the committed degree-6 run)
CELLS = {
    "10-5-7": ("10:5:7", "20260925"),  # the (n, ℓ) grid
    "13-5-7": ("13:5:7", "20260928"),  # the fixed-surplus ladder
    "11-5-7": ("11:5:7", "20260926"),  # the surplus control (secondary)
}
# Switch right after band 7 where it fits (about 10.2 GiB at (10, 5));
# elsewhere that block is 17-42 GiB, so band 6 is eliminated sparsely.
BUDGET_MB = {"10-5-7": "11000", "13-5-7": "6000", "11-5-7": "6000"}
ORDER = [f"{c}.u{k}" for c in ("10-5-7", "13-5-7", "11-5-7") for k in range(4)]
LIMIT_SECS = 96 * 3600


def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
    print(line, flush=True)
    with open(RUNS / "queue.log", "a") as f:
        f.write(line + "\n")


def hwm_kb(pid):
    try:
        for line in open(f"/proc/{pid}/status"):
            if line.startswith("VmHWM:"):
                return int(line.split()[1])
    except OSError:
        pass
    return None


def run(binary, job):
    cell, k = job.split(".u")
    cells, seed = CELLS[cell]
    out = RUNS / f"cell-{job}.jsonl"
    if out.exists() and out.read_text().strip():
        log(f"{job}: done already, skipped")
        return
    budget = os.environ.get("KIC_SPARSE_DENSE_BUDGET_MB", BUDGET_MB[cell])
    env = dict(
        os.environ,
        F4_F2_MAX_ROWS="50000000",
        F4_F2_MAX_COLS="50000000",
        KIC_SPARSE_DENSE_FINISH="1",
        KIC_SPARSE_DENSE_BUDGET_MB=budget,
        RAYON_NUM_THREADS=os.environ.get("RAYON_NUM_THREADS", "4"),
    )
    cmd = [binary, "--cells", cells, "--ffd-max", "5", "--seed", seed,
           "--unsat-index", k, "--d-min", "7"]
    log(f"{job}: start budget_mb={budget} threads={env['RAYON_NUM_THREADS']} cmd={' '.join(cmd)}")
    started = time.time()
    peak = 0
    with open(out, "w") as fo, open(RUNS / f"cell-{job}.log", "w") as fe:
        p = subprocess.Popen(cmd, stdout=fo, stderr=fe, env=env)
        while p.poll() is None:
            peak = max(peak, hwm_kb(p.pid) or 0)
            (RUNS / f"cell-{job}.hwm").write_text(f"{peak}\n")
            if time.time() - started > LIMIT_SECS:
                p.kill()
                log(f"{job}: watchdog at {LIMIT_SECS} s, killed (a resource limit, not evidence)")
            time.sleep(5)
    log(f"{job}: exit {p.returncode} after {time.time() - started:.0f} s, peak_rss_kB {peak}")


def main():
    binary = sys.argv[1]
    for job in sys.argv[2:] or ORDER:
        run(binary, job)


if __name__ == "__main__":
    main()
