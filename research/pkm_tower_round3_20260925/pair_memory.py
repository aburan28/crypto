#!/usr/bin/env python3
"""Paired baseline/candidate runs of note section 12.2: the round-2 build and
this round's build measure the same systems one after the other, and each run
records its peak resident memory (the kernel's `ru_maxrss` for that process
alone, through `wait4`) and its wall time. Wall time is a practicality note
(`AGENTS.md` section 2), and one run per build is not a timing claim.

    python3 pair_memory.py BASELINE_BIN CANDIDATE_BIN

writes `pairs/<system>-<build>.jsonl`, `.log`, and one summary line per run to
`pairs/memory.jsonl`. The deterministic fields of each pair must agree, which
`compare_builds.py`'s rule checks again here.
"""

import json
import os
import resource
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "pairs")
P0, P1 = "786433", "2013265921"
COMMON = ["--engine", "tower", "--planted", "0", "--random", "1", "--ladder-t", "none",
          "--budget", "7200", "--max-nnz", "2000000000"]
# Old sizes only: each system here was finished by round 2.
SYSTEMS = [
    ("K0-kummer-m2-p0-N20-t0", ["--p", P0, "--kinds", "kummer", "--m", "2", "--controls", "tower",
                                "--t-min", "10", "--max-t", "10"]),
    ("D1-kummer-m4-p1-N12-t0", ["--p", P1, "--kinds", "kummer", "--m", "4", "--controls", "tower",
                                "--t-min", "3", "--max-t-m4", "3"]),
    ("M3-kummer-m3-p1-N15-t0", ["--p", P1, "--kinds", "kummer", "--m", "3", "--controls", "tower",
                                "--t-min", "5", "--max-t-m3", "5"]),
    ("K1-kummer-m2-p1-N22-t0", ["--p", P1, "--kinds", "kummer", "--m", "2", "--controls", "tower",
                                "--t-min", "11", "--max-t", "11"]),
]
TIMING = {"ms", "wall_s"}


def limit():
    # The rounds' 14 GB address-space cap.
    cap = 14_000_000 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))


def run(binary, build, name, args):
    stem = os.path.join(OUT, f"{name}-{build}")
    started = time.time()
    with open(stem + ".log", "w") as log:
        p = subprocess.Popen([binary, *COMMON, *args, "--out", stem + ".jsonl"],
                             stderr=log, stdout=subprocess.DEVNULL, preexec_fn=limit)
        _, status, usage = os.wait4(p.pid, 0)
    wall = time.time() - started
    rows = [json.loads(l) for l in open(stem + ".jsonl") if l.strip()]
    rows = [r for r in rows if "N" in r]
    rec = {
        "system": name, "build": build, "binary": binary,
        "exit": os.waitstatus_to_exitcode(status),
        "peak_rss_mb": round(usage.ru_maxrss / 1024, 1),
        "wall_s": round(wall, 2),
        "f4_ms": rows[0]["ms"] if rows else None,
        "solving_degree_max": rows[0]["solving_degree_max"] if rows else None,
    }
    with open(os.path.join(OUT, "memory.jsonl"), "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)
    return rows


def main():
    baseline, candidate = sys.argv[1], sys.argv[2]
    os.makedirs(OUT, exist_ok=True)
    bad = 0
    for name, args in SYSTEMS:
        a = run(baseline, "baseline", name, args)
        b = run(candidate, "candidate", name, args)
        if len(a) != 1 or len(b) != 1:
            print(f"{name}: expected one row per build, got {len(a)} and {len(b)}")
            bad += 1
            continue
        diff = sorted(f for f in (set(a[0]) | set(b[0])) - TIMING if a[0].get(f) != b[0].get(f))
        if diff:
            print(f"{name}: the builds differ on {diff}")
            bad += 1
    print("pairs agree" if bad == 0 else f"{bad} pairs differ")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
