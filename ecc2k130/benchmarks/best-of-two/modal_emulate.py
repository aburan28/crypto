"""Run the walk-constant emulator (examples/ecc2k130_walk_constant.rs) on Modal CPU
containers, one container per cell, and keep each cell's JSON line unchanged.

    cargo build --release --example ecc2k130_walk_constant
    modal run modal_emulate.py --out emulation-<tag>.jsonl

Launcher only: the binary built in this tree is uploaded (Ubuntu 24.04 to
match the build host), its sha256 and the source commit go in the receipt.
"""
import hashlib
import json
import os
import pathlib
import subprocess
import time

import modal

HERE = pathlib.Path(__file__).resolve().parent
CPU = 32

image = modal.Image.from_registry("ubuntu:24.04", add_python="3.12")
if modal.is_local():
    ROOT = HERE.parents[2]
    BIN = ROOT / "target/release/examples/ecc2k130_walk_constant"
    image = image.add_local_file(BIN, "/usr/local/bin/walkc", copy=True)
app = modal.App("ecc2k130-best-of-two-emulation")

TRIALS = {23: 200000, 31: 200000, 37: 40000, 41: 20000, 53: 10000, 59: 100000}
WALKS = [("sigma", "hw"), ("sigma2", "hw"), ("sigma2", "key")]
if os.environ.get("WALKS") == "v2":
    WALKS = [("sigmapm", "hw"), ("sigma2t", "hw"), ("sigma2t", "key")]
DISTS = ["ecc2k130", "native"]


@app.function(image=image, cpu=CPU, memory=32768, timeout=4 * 60 * 60)
def cell(argv: list) -> dict:
    t0 = time.time()
    p = subprocess.run(["walkc", *argv, "--threads", str(CPU)], capture_output=True, text=True, check=True)
    return {"argv": argv, "json": p.stdout.strip(), "stderr": p.stderr[-4000:], "wall_s": round(time.time() - t0, 1)}


@app.local_entrypoint()
def main(out: str = "emulation.jsonl", seed: int = 20261014):
    jobs = [["--n", str(n), "--walk", w, "--score", sc, "--dist", d, "--walks", "16",
             "--trials", str(TRIALS[n]), "--seed", str(seed)]
            for n in TRIALS for (w, sc) in WALKS for d in DISTS]
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    sha = hashlib.sha256(BIN.read_bytes()).hexdigest()  # the binary the image was built from
    res = list(cell.map(jobs, order_outputs=True))
    with open(HERE / out, "w") as f:
        for r in res:
            f.write(r["json"] + "\n")
    with open(HERE / out.replace(".jsonl", ".log"), "w") as f:
        for r in res:
            f.write(" ".join(r["argv"]) + f"  ({r['wall_s']} s)\n" + r["stderr"] + "\n")
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip()
    receipt = {"started_utc": started, "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "commit": commit, "binary_sha256": sha,
               "image": "ubuntu:24.04", "cpu_per_container": CPU, "cells": len(jobs)}
    (HERE / out.replace(".jsonl", ".launch.json")).write_text(json.dumps(receipt, indent=1))
    print(f"{len(jobs)} cells -> {out}")
