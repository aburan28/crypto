#!/usr/bin/env python3
"""Matched baseline/candidate full-DLP runs of `ic run --solver mq-fes`.

Both binaries see identical curves, factor bases (legacy index 0), target
seeds and relation batch size 1.  Repetitions are interleaved
(baseline, candidate, baseline, ...) so host drift hits both arms alike.
Every record keeps the full `ic` JSON report; nothing is dropped.

    python3 run.py --baseline BIN --candidate BIN --output DIR
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


def run_one(binary, degree, curve_a, seed):
    cmd = [binary, "run", "--json", "--degree", str(degree), "--curve-a", str(curve_a),
           "--solver", "mq-fes", "--summands", "2", "--random-target", "--seed", str(seed),
           "--batch", "1", "--max-trials", "200000"]
    start = time.perf_counter()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return {"status": "timed_out", "process_seconds": time.perf_counter() - start, "command": cmd}
    elapsed = time.perf_counter() - start
    try:
        report = json.loads(proc.stdout)
    except json.JSONDecodeError:
        report = {"status": "failed", "stderr": proc.stderr[-4000:]}
    return {"command": cmd, "exit": proc.returncode, "process_seconds": elapsed, "report": report}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--baseline-rev", required=True)
    ap.add_argument("--candidate-rev", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    if os.path.exists(args.output):
        sys.exit(f"refusing to overwrite {args.output}")
    os.makedirs(args.output)
    manifest = {
        "schema_version": 1,
        "experiment": "mq_fes_ic_pipeline_20260922",
        "arms": {
            "baseline": {"binary_sha256": sha256(args.baseline), "source_revision": args.baseline_rev},
            "candidate": {"binary_sha256": sha256(args.candidate), "source_revision": args.candidate_rev},
        },
        "instances": [{"degree": n, "curve_a": a, "factor_index": 0, "summands": 2} for n, a in INSTANCES],
        "train_seeds": TRAIN_SEEDS,
        "holdout_seeds": HOLDOUT_SEEDS,
        "repetitions": REPS,
        "batch": 1,
        "order": "per (instance, seed, rep): baseline then candidate",
        "host": {"platform": platform.platform(), "cpus": os.cpu_count()},
    }
    with open(os.path.join(args.output, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    raw = open(os.path.join(args.output, "raw.jsonl"), "w")
    for n, a in INSTANCES:
        for split, seeds in (("train", TRAIN_SEEDS), ("holdout", HOLDOUT_SEEDS)):
            for seed in seeds:
                for rep in range(REPS):
                    for arm, binary in (("baseline", args.baseline), ("candidate", args.candidate)):
                        rec = run_one(binary, n, a, seed)
                        rec.update({"arm": arm, "degree": n, "curve_a": a, "seed": seed,
                                    "split": split, "rep": rep})
                        raw.write(json.dumps(rec) + "\n")
                        raw.flush()
                print(f"n={n} a={a} {split} seed={seed} done", flush=True)
    raw.close()


if __name__ == "__main__":
    main()
