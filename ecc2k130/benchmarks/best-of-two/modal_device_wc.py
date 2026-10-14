"""The device walks' own iteration constant (src/walkconstant.cpp, Ref::step on
the type-II normal basis) for SIGMA_SIGN = 0, 1, 2, one Modal container per cell.

    g++ <WALK_CONSTANT_FLAGS> -DECC_SIGMA_SIGN={0,1,2} src/walkconstant.cpp -o device-wc/wc-sign{0,1,2}
    modal run modal_device_wc.py --out device-wc/device-wc-v1.jsonl

Launcher only; each cell's JSON line is kept unchanged.  The binaries are
built in this tree (Ubuntu 24.04 to match), uploaded, and hashed in the receipt.
"""
import hashlib
import json
import pathlib
import subprocess
import time

import modal

HERE = pathlib.Path(__file__).resolve().parent
CPU = 64
SIGNS = (0, 1, 2)

image = modal.Image.from_registry("ubuntu:24.04", add_python="3.12")
if modal.is_local():
    for s in SIGNS:
        image = image.add_local_file(HERE / f"device-wc/wc-sign{s}", f"/usr/local/bin/wc-sign{s}", copy=True)
app = modal.App("ecc2k130-best-of-two-device-wc")

CELLS = [(23, 8, 200000), (41, 16, 4000)]
DISTS = ["native", "ecc2k130"]


@app.function(image=image, cpu=CPU, memory=32768, timeout=6 * 60 * 60)
def cell(sign: int, argv: list) -> dict:
    t0 = time.time()
    p = subprocess.run([f"wc-sign{sign}", *argv, "--threads", str(CPU)], capture_output=True, text=True, check=True)
    line = p.stdout.strip().splitlines()[-1]
    d = json.loads(line)
    d["sigma_sign"] = sign
    return {"json": json.dumps(d), "stderr": p.stderr[-2000:], "argv": argv, "sign": sign,
            "wall_s": round(time.time() - t0, 1)}


@app.local_entrypoint()
def main(out: str = "device-wc/device-wc-v1.jsonl", seed: int = 20261014):
    jobs = [(s, ["--n", str(n), "--walk", "sigma", "--dist", d, "--walks", str(w),
                 "--trials", str(t), "--seed", str(seed)])
            for (n, w, t) in CELLS for d in DISTS for s in SIGNS]
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    shas = {s: hashlib.sha256((HERE / f"device-wc/wc-sign{s}").read_bytes()).hexdigest() for s in SIGNS}
    res = list(cell.starmap(jobs, order_outputs=True))
    with open(HERE / out, "w") as f:
        for r in res:
            f.write(r["json"] + "\n")
    with open(HERE / out.replace(".jsonl", ".log"), "w") as f:
        for r in res:
            f.write(f"sign {r['sign']} " + " ".join(r["argv"]) + f"  ({r['wall_s']} s)\n" + r["stderr"] + "\n")
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=HERE).stdout.strip()
    receipt = {"started_utc": started, "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "commit": commit, "binary_sha256": shas, "image": "ubuntu:24.04", "cpu_per_container": CPU,
               "build": "g++ -O2 -std=c++17 -DECC_WALK_TABLE=1 -DECC_TABLE_PIVOT_BYTES=0 -pthread "
                        "-DECC_TABLE_BRANCHES=8 -DECC_WITNESS=0 -DECC_SIGMA_SIGN=<s> src/walkconstant.cpp",
               "cells": len(jobs)}
    (HERE / out.replace(".jsonl", ".launch.json")).write_text(json.dumps(receipt, indent=1))
    print(f"{len(jobs)} cells -> {out}")
