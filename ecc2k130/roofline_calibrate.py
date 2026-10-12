#!/usr/bin/env python3
"""Which pipe peaks are consistent with every measured build?  roofline.py's check.

roofline.py turns a build into dynamic lane-instructions per scalar update on
each pipe; a pipe's peak rate turns that into SM-clocks per update, and the
largest of those is a floor on the time an update can take.  A peak is
therefore falsifiable: any measured build that ran faster than the floor a
model puts under it proves that model's peak too low.  This script rebuilds
every arm of the frozen automatic sweeps (benchmarks/autosweep: one session,
one card, arms alternating with their base), prices each under several pipe
models, and prints, per model, each build's efficiency

    eta = floor(model) / measured time        (<= 1 for a valid model)

A model whose largest eta exceeds 1 by more than the rate noise is refuted.
Among the survivors, the one whose eta is flattest across builds of the same
geometry is the one that explains the knob-to-knob differences -- the sweep
arms were chosen to move one pipe's work at a time.

    ./roofline_calibrate.py                    # both GPUs, every sweep arm
    ./roofline_calibrate.py --gpu rtx-pro-6000 --json out.json

Needs what roofline.py needs (nvcc and nvdisasm, CUDA 13.3+); about 3 s a build.
"""
import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from roofline import MACHINES, analyse  # noqa: E402

B200_BASE = "PACKED_TOP_CLMAD=1 PACKED_ALU_SQUARE=0 PACKED_ONB_INV=1"

# Pipe models: name -> {pipe: lanes per SM-clock}.  "int" is alu + fma, for the
# hypothesis that every integer op shares one datapath (the hardware-limits
# probe measured a LOP3+IMAD mix at 69 lanes/SM-clock, one op alone at 62-64).
MODELS = {
    "rtx-pro-6000": [
        ("clmad 1.62 only (ONE-BLOCK-GEOMETRY)", {"clmad": 1.62}),
        ("clmad 2.00 only", {"clmad": 2.00}),
        ("alu 64, clmad 1.62", {"alu": 64.0, "clmad": 1.62, "issue": 128.0}),
        ("alu 64, clmad 2.00", {"alu": 64.0, "clmad": 2.00, "issue": 128.0}),
        ("int 69 (alu+fma), clmad 2.00", {"int": 69.0, "clmad": 2.00, "issue": 128.0}),
        ("alu 128 (Nsight's peak), clmad 2.00", {"alu": 128.0, "clmad": 2.00, "issue": 128.0}),
    ],
    "b200": [
        ("alu 63.7, clmad 29.1", {"alu": 63.7, "clmad": 29.1, "issue": 128.0}),
        ("int 69 (alu+fma), clmad 29.1", {"int": 69.0, "clmad": 29.1, "issue": 128.0}),
        ("alu 128 (Nsight's peak), clmad 29.1", {"alu": 128.0, "clmad": 29.1, "issue": 128.0}),
    ],
}


def receipts(gpus):
    rows = []
    if "rtx-pro-6000" in gpus:
        path = "benchmarks/autosweep/rtx-pro-6000/summary.json"
        d = json.load(open(os.path.join(HERE, path)))
        rows.append({"gpu": "rtx-pro-6000", "arm": "base (20 B/s knob set)", "knobs": "",
                     "bps": d["baseScreenedBps"], "source": path})
        for a in d["star"]:
            rows.append({"gpu": "rtx-pro-6000", "arm": a["arm"], "knobs": a["knobs"], "bps": a["screenedBps"],
                         "regs": a.get("registers"), "source": path})
    if "b200" in gpus:
        path = "benchmarks/autosweep/b200/summary.json"
        d = json.load(open(os.path.join(HERE, path)))
        rows.append({"gpu": "b200", "arm": "base (" + B200_BASE + ")", "knobs": B200_BASE,
                     "bps": d["baseScreenedBps"], "source": path})
        for a in d["star"]:
            rows.append({"gpu": "b200", "arm": a["arm"], "knobs": B200_BASE + " " + a["knobs"],
                         "bps": a["screenedBps"], "regs": a.get("registers"), "source": path})
        for g in d["greedy"]:
            if g.get("decision") == "keep":
                rows.append({"gpu": "b200", "arm": "greedy " + g["name"], "knobs": B200_BASE + " " + g["knobs"],
                             "bps": g["screenedBps"], "regs": g.get("registers"), "source": path})
    return rows


def geometry(defs):
    return "%sx%s b%s" % (defs.get("ECC_THREADS"), defs.get("ECC_MINBLOCKS"), defs.get("ECC_BATCH"))


def work(by_pipe, pipe):
    if pipe == "int":
        return float(by_pipe.get("alu", 0)) + float(by_pipe.get("fma", 0))
    return float(by_pipe.get(pipe, 0))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gpu", action="append", choices=sorted(MACHINES), help="default: both")
    ap.add_argument("--json", help="write every row and model here")
    a = ap.parse_args()
    gpus = a.gpu or ["rtx-pro-6000", "b200"]
    out = []
    for r in receipts(gpus):
        m = MACHINES[r["gpu"]]
        res = analyse("gpu-preset", r["knobs"], r["gpu"])
        bp = res["by_pipe"]
        t_meas = m["sms"] * m["clock_ghz"] * 1e9 / (r["bps"] * 1e9)
        row = dict(r)
        row.update({"geometry": geometry(res["defs"]), "regsBuilt": res["regs"], "spills": res["spills"],
                    "smClocksPerUpdate": t_meas, "clmad": work(bp, "clmad"), "expectedClmad": res["expectedClmad"],
                    "alu": work(bp, "alu"), "fma": work(bp, "fma"), "issue": work(bp, "issue"),
                    "lsu": work(bp, "lsu"), "eta": {}, "binding": {}})
        for name, model in MODELS[r["gpu"]]:
            floors = {p: work(bp, p) / rate for p, rate in model.items()}
            bind = max(floors, key=floors.get)
            row["eta"][name] = floors[bind] / t_meas
            row["binding"][name] = bind
        out.append(row)
        chk = "" if row["expectedClmad"] is None or abs(row["clmad"] - row["expectedClmad"]) < 0.01 else \
            "  CLMAD MISMATCH (expected %.3f)" % row["expectedClmad"]
        print("%-12s %-30s %-12s %7.3f B/s  %6.2f clk  alu %7.1f fma %6.1f clmad %6.2f issue %7.1f%s" % (
            r["gpu"], r["arm"][:30], row["geometry"], r["bps"], t_meas, row["alu"], row["fma"], row["clmad"],
            row["issue"], chk), flush=True)
    for gpu in gpus:
        rows = [o for o in out if o["gpu"] == gpu]
        if not rows:
            continue
        base_geo = rows[0]["geometry"]
        print("\n%s: efficiency eta = model floor / measured time, per model (<= 1 unless the model is refuted)"
              % MACHINES[gpu]["name"])
        names = [n for n, _ in MODELS[gpu]]
        print("  %-30s %-12s %s" % ("arm", "geometry", "  ".join("%13s" % ("M%d" % (i + 1)) for i in range(len(names)))))
        for o in rows:
            print("  %-30s %-12s %s" % (o["arm"][:30], o["geometry"],
                                        "  ".join("%7.3f %-5s" % (o["eta"][n], o["binding"][n][:5]) for n in names)))
        print("  models:")
        for i, n in enumerate(names):
            etas = [o["eta"][n] for o in rows]
            same = [o["eta"][n] for o in rows if o["geometry"] == base_geo]
            mean = sum(same) / len(same)
            sd = math.sqrt(sum((e - mean) ** 2 for e in same) / max(1, len(same) - 1))
            worst = max(rows, key=lambda o: o["eta"][n])
            verdict = "REFUTED by %s (eta %.3f)" % (worst["arm"], worst["eta"][n]) if worst["eta"][n] > 1.02 \
                else "consistent"
            print("    M%d %-40s max eta %.3f; %s geometry: mean %.3f sd %.3f over %d builds -- %s" % (
                i + 1, n, max(etas), base_geo, mean, sd, len(same), verdict))
    if a.json:
        with open(a.json, "w") as f:
            json.dump({"models": {g: [[n, m] for n, m in MODELS[g]] for g in gpus}, "rows": out}, f, indent=1)


if __name__ == "__main__":
    main()
