#!/usr/bin/env python3
"""Paired instruction counts: the triple-table collector against `scaled`.

    python3 check.py SCALED_WORKER TRIPLE_WORKER --fixtures 32 --seed 20260924 --out confirm.json

Both workers come from `build_arms.sh`.  Every fixture is drawn once and run on
both arms with the same target and algorithm seeds; each run is counted by
callgrind (`Collected:`, the tournament's own unit) and its report goes through
`oracle.py`, the tournament's independent checker, with the arm's own summand
count.  A report that fails the checker, or recovers a different logarithm from
its pair, voids the run: the script exits non-zero and prints no ratio.

The summary statistic is fixed in PREREGISTRATION.md: per cell, the geometric
mean of the paired triple/scaled ratios with a 95% percentile bootstrap interval,
the median beside it.
"""
import argparse
import json
import math
import random
import re
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "research/ic_candidate_tournament_20260915"))
import oracle  # noqa: E402

# Round 0023's configuration: round 0020's promoted config, max_trials raised to 65,536.
PAIR = {"batch_trials": 1, "linear_algebra": "sparse", "max_trials": 65536,
        "solver": "pair_table", "summands": 3}
TRIPLE = dict(PAIR, solver="triple_table", summands=4)
CELLS = {"n23a1": (23, 1), "n37a0": (37, 0), "n43a1": (43, 1)}


def job(cell, tseed, aseed, config):
    n, a = CELLS[cell]
    return {"degree": n, "curve_a": a, "target_seeds": [tseed], "algorithm_seed": aseed,
            "config": config,
            "factor_base": {"kind": "subgroup_orbits", "seed": 43, "points": 222}}


def plain(worker, payload):
    out = subprocess.run([worker], input=json.dumps(payload), text=True,
                         capture_output=True, timeout=1800)
    return json.loads(out.stdout)


def counted(worker, payload):
    out = subprocess.run(["valgrind", "--tool=callgrind", "--callgrind-out-file=/dev/null", worker],
                         input=json.dumps(payload), text=True, capture_output=True, timeout=3600)
    found = re.search(r"Collected\s*:\s*(\d+)", out.stderr)
    if not found:
        raise RuntimeError("callgrind reported no instruction count")
    return int(found.group(1)), json.loads(out.stdout)


def one(scaled, triple, cell, tseed, aseed):
    fixture = plain(scaled, dict(job(cell, tseed, aseed, PAIR), mode="fixture"))["fixture"]
    ir_s, rep_s = counted(scaled, dict(job(cell, tseed, aseed, PAIR), mode="ic"))
    ir_t, rep_t = counted(triple, dict(job(cell, tseed, aseed, TRIPLE), mode="ic"))
    oracle.verify(rep_s, fixture, summands=3)
    oracle.verify(rep_t, fixture, summands=4)
    logs_s = [s["recovered"] for s in sorted(rep_s["solutions"], key=lambda s: s["index"])]
    logs_t = [s["recovered"] for s in sorted(rep_t["solutions"], key=lambda s: s["index"])]
    if logs_s != logs_t:
        raise RuntimeError(f"{cell}: arms recovered different logarithms")
    return {"cell": cell, "target_seed": tseed, "algorithm_seed": aseed,
            "ir_scaled": ir_s, "ir_triple": ir_t, "ratio": ir_t / ir_s,
            "trials_scaled": rep_s.get("trials"), "trials_triple": rep_t.get("trials"),
            "columns_scaled": rep_s.get("columns"), "columns_triple": rep_t.get("columns")}


def summarise(ratios, resamples=10_000, seed=0):
    logs = [math.log(x) for x in ratios]
    gmean = math.exp(statistics.fmean(logs))
    rng = random.Random(seed)
    boots = sorted(math.exp(statistics.fmean(rng.choices(logs, k=len(logs))))
                   for _ in range(resamples))
    return {"n": len(ratios), "geometric_mean": gmean,
            "ci95": [boots[int(0.025 * resamples)], boots[int(0.975 * resamples) - 1]],
            "median": statistics.median(ratios), "min": min(ratios), "max": max(ratios)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scaled")
    ap.add_argument("triple")
    ap.add_argument("--cells", default="n23a1,n37a0,n43a1")
    ap.add_argument("--fixtures", type=int, default=32)
    ap.add_argument("--seed", type=int, default=20260924)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    rows = []
    for cell in args.cells.split(","):
        for _ in range(args.fixtures):
            rows.append(one(args.scaled, args.triple, cell, rng.getrandbits(64), rng.getrandbits(64)))
            print(json.dumps({k: rows[-1][k] for k in ("cell", "ratio", "trials_scaled", "trials_triple")}),
                  flush=True)
    summary = {cell: summarise([r["ratio"] for r in rows if r["cell"] == cell])
               for cell in args.cells.split(",")}
    Path(args.out).write_text(json.dumps({"seed": args.seed, "fixtures_per_cell": args.fixtures,
                                          "summary": summary, "rows": rows}, indent=1) + "\n")
    for cell, s in summary.items():
        print(f"{cell}: triple/scaled geometric mean {s['geometric_mean']:.3f} "
              f"[{s['ci95'][0]:.3f}, {s['ci95'][1]:.3f}], median {s['median']:.3f}, n={s['n']}")


if __name__ == "__main__":
    main()
