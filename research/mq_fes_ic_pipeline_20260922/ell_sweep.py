#!/usr/bin/env python3
"""Per-call mq-fes oracle cost against ℓ = dim V, both arms, and a fit.

Each (curve, arm, seed) is one cold `ic run`; the oracle counters are
divided by that run's oracle calls.  Fits are least squares of log2(cost)
on ℓ, reported as the base-2 exponent per unit of ℓ.

    python3 ell_sweep.py --arm baseline=BIN:REV --arm mqfes-linear=BIN:REV --output DIR
"""
import argparse
import json
import math
import os
import statistics
import subprocess
import sys

CURVES = [(7, 1), (9, 1), (11, 1), (17, 1), (23, 0), (23, 1), (13, 0)]
SEEDS = [1, 2, 3, 4, 5]


def fit(xs, ys):
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return slope, my - slope * mx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", action="append", required=True, help="NAME=BINARY:REVISION")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    if os.path.exists(args.output):
        sys.exit(f"refusing to overwrite {args.output}")
    os.makedirs(args.output)
    arms = [(a.split("=", 1)[0], *a.split("=", 1)[1].split(":", 1)) for a in args.arm]
    records = []
    with open(f"{args.output}/raw.jsonl", "w") as raw:
        for n, a in CURVES:
            for seed in SEEDS:
                for name, binary, rev in arms:
                    cmd = [binary, "run", "--json", "--degree", str(n), "--curve-a", str(a),
                           "--solver", "mq-fes", "--random-target", "--seed", str(seed), "--batch", "1"]
                    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
                    report = json.loads(proc.stdout)
                    rec = {"arm": name, "revision": rev, "degree": n, "curve_a": a, "seed": seed,
                           "command": cmd, "report": report}
                    raw.write(json.dumps(rec) + "\n")
                    records.append(rec)
    per = {}
    for r in records:
        rep = r["report"]
        ok = rep.get("status") == "complete" and rep.get("result", {}).get("verified") is True
        p = rep["counts"]["mq_fes_profile"]
        key = (r["arm"], rep["factor_base"]["dimension"], r["degree"], r["curve_a"])
        e = per.setdefault(key, {"verified": 0, "runs": 0, "word_ops": [], "walk_us": [], "oracle_us": []})
        e["runs"] += 1
        e["verified"] += ok
        if p["calls"]:
            e["word_ops"].append(p.get("word_ops", 2 * p["points"]) / p["calls"])
            e["walk_us"].append(p["walk_ns"] / p["calls"] / 1e3)
            e["oracle_us"].append((p["build_ns"] + p["walk_ns"] + p["lift_ns"]) / p["calls"] / 1e3)
    rows = []
    for (arm, ell, n, a), e in sorted(per.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])):
        rows.append({"arm": arm, "ell": ell, "degree": n, "curve_a": a, "runs": e["runs"],
                     "verified": e["verified"],
                     "word_ops_per_call": statistics.median(e["word_ops"]),
                     "walk_us_per_call": statistics.median(e["walk_us"]),
                     "oracle_us_per_call": statistics.median(e["oracle_us"])})
    fits = {}
    for arm in dict.fromkeys(r["arm"] for r in rows):
        pts = [r for r in rows if r["arm"] == arm and r["ell"] >= 3]
        xs = [r["ell"] for r in pts]
        fits[arm] = {
            "sizes": sorted(set(xs)),
            "word_ops_log2_slope": fit(xs, [math.log2(r["word_ops_per_call"]) for r in pts])[0],
            "walk_us_log2_slope": fit(xs, [math.log2(max(r["walk_us_per_call"], 1e-3)) for r in pts])[0],
        }
    summary = {"schema_version": 1, "rows": rows, "fits": fits,
               "note": "slope is d log2(cost) / d ell; a 2^{c·ell} cost has slope c"}
    with open(f"{args.output}/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    for r in rows:
        print(f"{r['arm']:>13} ell={r['ell']:>2} n={r['degree']:>2} a={r['curve_a']} "
              f"ver={r['verified']}/{r['runs']} wordops/call={r['word_ops_per_call']:.4g} "
              f"walk µs/call={r['walk_us_per_call']:.2f} oracle µs/call={r['oracle_us_per_call']:.1f}")
    print(json.dumps(fits, indent=2))


if __name__ == "__main__":
    main()
