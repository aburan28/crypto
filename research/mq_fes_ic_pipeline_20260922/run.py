#!/usr/bin/env python3
"""Matched full-DLP runs of `ic run` across named arms.

Every arm sees identical curves, factor bases (legacy index 0), target
seeds and relation batch size 1.  Repetitions are interleaved across arms
(arm 1, arm 2, ..., arm 1, ...) so host drift hits every arm alike.  Each
record keeps the full `ic` JSON report; nothing is dropped.

    python3 run.py --arm baseline=BIN:mq-fes:REV --arm candidate=BIN:mq-fes:REV \
        --arm enumerate=BIN:enumerate:REV --output DIR
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time

INSTANCES = [(13, 0), (17, 1), (23, 0), (23, 1)]
TRAIN_SEEDS = list(range(1, 9))
HOLDOUT_SEEDS = list(range(1001, 1009))
REPS = 3
TIMEOUT_S = 600


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_one(binary, solver, degree, curve_a, seed):
    cmd = [binary, "run", "--json", "--degree", str(degree), "--curve-a", str(curve_a),
           "--solver", solver, "--summands", "2", "--random-target", "--seed", str(seed),
           "--batch", "1", "--max-trials", "200000"]
    start = time.perf_counter()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return {"command": cmd, "status": "timed_out", "process_seconds": time.perf_counter() - start}
    elapsed = time.perf_counter() - start
    try:
        report = json.loads(proc.stdout)
    except json.JSONDecodeError:
        report = {"status": "failed", "stderr": proc.stderr[-4000:]}
    return {"command": cmd, "exit": proc.returncode, "process_seconds": elapsed, "report": report}


def parse_arm(spec):
    name, rest = spec.split("=", 1)
    binary, solver, rev = rest.split(":", 2)
    return name, binary, solver, rev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", action="append", required=True, help="NAME=BINARY:SOLVER:REVISION")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    arms = [parse_arm(a) for a in args.arm]
    if os.path.exists(args.output):
        sys.exit(f"refusing to overwrite {args.output}")
    os.makedirs(args.output)
    manifest = {
        "schema_version": 2,
        "experiment": "mq_fes_ic_pipeline_20260922",
        "arms": {name: {"binary_sha256": sha256(b), "solver": s, "source_revision": rev}
                 for name, b, s, rev in arms},
        "instances": [{"degree": n, "curve_a": a, "factor_index": 0, "summands": 2} for n, a in INSTANCES],
        "train_seeds": TRAIN_SEEDS,
        "holdout_seeds": HOLDOUT_SEEDS,
        "repetitions": REPS,
        "batch": 1,
        "order": "per (instance, seed, rep): arms in the listed order",
        "host": {"platform": platform.platform(), "cpus": os.cpu_count()},
    }
    with open(os.path.join(args.output, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    with open(os.path.join(args.output, "raw.jsonl"), "w") as raw:
        for n, a in INSTANCES:
            for split, seeds in (("train", TRAIN_SEEDS), ("holdout", HOLDOUT_SEEDS)):
                for seed in seeds:
                    for rep in range(REPS):
                        for name, binary, solver, _ in arms:
                            rec = run_one(binary, solver, n, a, seed)
                            rec.update({"arm": name, "degree": n, "curve_a": a, "seed": seed,
                                        "split": split, "rep": rep})
                            raw.write(json.dumps(rec) + "\n")
                            raw.flush()
                    print(f"n={n} a={a} {split} seed={seed} done", flush=True)


if __name__ == "__main__":
    main()
