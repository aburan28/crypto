#!/usr/bin/env python3
"""Run every (n, encoding) cell of the Hamming-ideal PDP experiment as its
own process, four at a time, largest cells first. Raw output goes to
results/<tag>/n<n>_<enc>.jsonl with a log beside it; nothing is overwritten.

    python3 run.py --tag main            # the campaign
    python3 run.py --tag tau --tau-sweep # the budget-sensitivity check
"""
import argparse, hashlib, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "hamming_pdp", "target", "release", "hamming_pdp")

# Amended before these cells ran (see RESULT.md, "Deviations"): w = 2 at
# every n, 8 targets at n >= 17, budget 2^28 row-word XORs, matrix cap
# 2^26 words. Small cells run first so that partial evidence can be frozen.
CELLS = [
    (7, 2, 24), (11, 2, 24), (13, 2, 24), (17, 2, 8), (19, 2, 8),
]
ENCODINGS = ["C", "FC", "QFC", "MONO", "SUB"]

def source_hashes():
    out = {}
    src = os.path.join(HERE, "hamming_pdp", "src")
    for f in sorted(os.listdir(src)):
        with open(os.path.join(src, f), "rb") as fh:
            out[f] = hashlib.sha256(fh.read()).hexdigest()
    with open(BIN, "rb") as fh:
        out["binary"] = hashlib.sha256(fh.read()).hexdigest()
    return out

def run_cell(args):
    outdir, n, w, targets, enc, budget, cap = args
    out = os.path.join(outdir, f"n{n}_{enc}.jsonl")
    log = os.path.join(outdir, f"n{n}_{enc}.log")
    if os.path.exists(out):
        return (n, enc, "exists")
    cmd = [BIN, "--n", str(n), "--w", str(w), "--targets", str(targets), "--encodings", enc,
           "--budget-log2", str(budget), "--matrix-cap-log2", str(cap), "--out", out + ".part"]
    # Resume: a cell interrupted by a lost container keeps its finished
    # targets and runs only the missing seeds, appending to the same file.
    done = set()
    if os.path.exists(out + ".part"):
        with open(out + ".part") as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if r.get("kind") == "cell":
                    done.add(r["seed"])
    missing = [s for s in range(targets) if s not in done]
    if done and missing:
        cmd += ["--seeds", ",".join(map(str, missing)), "--append"]
    t0 = time.time()
    rc = 0
    if missing:
        with open(log, "a") as lf:
            lf.write(" ".join(cmd) + "\n")
            lf.flush()
            rc = subprocess.call(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if rc == 0:
        os.rename(out + ".part", out)
    return (n, enc, f"rc={rc} {time.time()-t0:.0f}s")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--tau-sweep", action="store_true")
    a = ap.parse_args()
    outdir = os.path.join(HERE, "results", a.tag)
    os.makedirs(outdir, exist_ok=True)
    jobs = []
    if a.tau_sweep:
        # n = 11, 4 targets, FC and MONO, budgets 2^26, 2^28, 2^30 (cap 2^27).
        for enc in ["FC", "MONO", "QFC"]:
            for b in [26, 28, 30]:
                jobs.append((outdir, 11, 2, 4, f"{enc}", b, 27))
        # distinct output names per budget
        jobs = [(o, n, w, t, e, b, c) for (o, n, w, t, e, b, c) in jobs]
        def run_tau(j):
            o, n, w, t, e, b, c = j
            out = os.path.join(o, f"n{n}_{e}_tau{b}.jsonl")
            log = out.replace(".jsonl", ".log")
            if os.path.exists(out):
                return (n, e, b, "exists")
            cmd = [BIN, "--n", str(n), "--w", str(w), "--targets", str(t), "--encodings", e,
                   "--budget-log2", str(b), "--matrix-cap-log2", str(c), "--out", out + ".part"]
            with open(log, "w") as lf:
                lf.write(" ".join(cmd) + "\n"); lf.flush()
                rc = subprocess.call(cmd, stdout=lf, stderr=subprocess.STDOUT)
            os.rename(out + ".part", out)
            return (n, e, b, f"rc={rc}")
        with ThreadPoolExecutor(a.jobs) as ex:
            for r in ex.map(run_tau, jobs):
                print(r, flush=True)
    else:
        for (n, w, t) in CELLS:
            for enc in ENCODINGS:
                jobs.append((outdir, n, w, t, enc, 28, 26))
        with ThreadPoolExecutor(a.jobs) as ex:
            for r in ex.map(run_cell, jobs):
                print(r, flush=True)
    with open(os.path.join(outdir, "manifest.json"), "w") as f:
        json.dump({"tag": a.tag, "cells": CELLS, "encodings": ENCODINGS, "sources": source_hashes(),
                   "host": os.uname().machine, "cpus": os.cpu_count(),
                   "finished": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, f, indent=1)

if __name__ == "__main__":
    main()
