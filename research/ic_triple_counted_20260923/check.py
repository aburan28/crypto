#!/usr/bin/env python3
"""Paired instruction counts: counted sizing against the triple arm and against rho.

    python3 check.py TRIPLE_WORKER COUNTED_WORKER --seed 20260926 --fixtures 32 --out confirm.json

Both workers come from `build_arms.sh`; rho is the triple worker's rho mode.
Every fixture is drawn once and run on all three arms with the same target and
algorithm seeds; each run is counted by callgrind (`Collected:`) and its report
goes through `oracle.py` -- summands 4 for the IC arms, rho mode for rho.

Per PREREGISTRATION.md:

- a finished report that fails the checker, or three arms that disagree on a
  logarithm, is a defect: the script stops and prints no ratio;
- a run that does not finish (a report that is not `complete`, or a timeout) is
  recorded as unfinished, and that fixture's cell then carries no comparison;
- at every cell but `n43a1`, the two IC arms must report the same base, trials,
  columns and relations (the identity check).

Fixtures run four at a time; callgrind's count does not depend on that.
"""
import argparse
import importlib.util
import json
import random
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "research/ic_candidate_tournament_20260915"))
import oracle  # noqa: E402

# The triple study's statistic, unchanged: loaded by path, since this file shares its name.
_spec = importlib.util.spec_from_file_location("triple_check", REPO / "research/ic_triple_table_20260923/check.py")
_triple_check = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_triple_check)
summarise = _triple_check.summarise

# Round 0023's configuration: round 0020's promoted config, max_trials raised to 65,536.
PAIR = {"batch_trials": 1, "linear_algebra": "sparse", "max_trials": 65536,
        "solver": "pair_table", "summands": 3}
ARMS = {"triple": dict(PAIR, solver="triple_table", summands=4),
        "counted": dict(PAIR, solver="triple_counted", summands=4)}
CELLS = {"n23a1": (23, 1), "n37a0": (37, 0), "n43a1": (43, 1), "n59a0": (59, 0), "n61a1": (61, 1)}
CHANGED = {"n43a1"}  # the one cell where the count chooses a different base
SAME = ("factor_base_orbits", "trials", "columns", "accepted_relations", "duplicate_relations")


def job(cell, tseed, aseed, config):
    n, a = CELLS[cell]
    return {"degree": n, "curve_a": a, "target_seeds": [tseed], "algorithm_seed": aseed,
            "config": config, "factor_base": {"kind": "subgroup_orbits", "seed": 43, "points": 222}}


def counted(worker, payload):
    try:
        out = subprocess.run(["valgrind", "--tool=callgrind", "--callgrind-out-file=/dev/null", worker],
                             input=json.dumps(payload), text=True, capture_output=True, timeout=3600)
    except subprocess.TimeoutExpired:
        return None, {"status": "timeout"}
    found = re.search(r"Collected\s*:\s*(\d+)", out.stderr)
    if not found:
        raise RuntimeError("callgrind reported no instruction count")
    return int(found.group(1)), json.loads(out.stdout)


def logs(report):
    return [s["recovered"] for s in sorted(report["solutions"], key=lambda s: s["index"])]


def one(args):
    triple, counted_worker, cell, tseed, aseed, stream = args
    fixture = json.loads(subprocess.run([triple], input=json.dumps(dict(job(cell, tseed, aseed, PAIR), mode="fixture")),
                                        text=True, capture_output=True, timeout=1800).stdout)["fixture"]
    runs = {"triple": counted(triple, dict(job(cell, tseed, aseed, ARMS["triple"]), mode="ic")),
            "counted": counted(counted_worker, dict(job(cell, tseed, aseed, ARMS["counted"]), mode="ic")),
            "rho": counted(triple, dict(job(cell, tseed, aseed, PAIR), mode="rho"))}
    row = {"cell": cell, "stream": stream, "target_seed": tseed, "algorithm_seed": aseed, "unfinished": []}
    for arm, (ir, rep) in runs.items():
        row[f"ir_{arm}"] = ir
        if rep.get("status") != "complete":
            row["unfinished"].append(arm)
            continue
        if arm == "rho":
            oracle.verify(rep, fixture, expected_mode="rho")
        else:
            oracle.verify(rep, fixture, summands=4)
            row[f"trials_{arm}"] = rep.get("trials")
            row[f"columns_{arm}"] = rep.get("columns")
            row[f"accepted_{arm}"] = rep.get("accepted_relations")
    finished = [arm for arm in runs if arm not in row["unfinished"]]
    if len({tuple(logs(runs[arm][1])) for arm in finished}) > 1:
        raise RuntimeError(f"{cell} {tseed}: arms recovered different logarithms")
    if cell not in CHANGED and not {"triple", "counted"} & set(row["unfinished"]):
        t, c = runs["triple"][1], runs["counted"][1]
        differ = [k for k in SAME if t.get(k) != c.get(k)]
        if differ:
            raise RuntimeError(f"{cell} {tseed}: identity broken on {differ}")
    if not row["unfinished"]:
        row["counted_over_triple"] = row["ir_counted"] / row["ir_triple"]
        row["counted_over_rho"] = row["ir_counted"] / row["ir_rho"]
        row["triple_over_rho"] = row["ir_triple"] / row["ir_rho"]
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("triple")
    ap.add_argument("counted")
    ap.add_argument("--cells", default=",".join(CELLS))
    ap.add_argument("--fixtures", type=int, default=32)
    ap.add_argument("--seed", type=int, default=20260926)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--extend", action="append", default=[], metavar="CELL:COUNT:SEED",
                    help="more fixtures at one cell from a separate stream, appended after the registered ones")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    cells = args.cells.split(",")
    work = [(args.triple, args.counted, cell, rng.getrandbits(64), rng.getrandbits(64), "registered")
            for cell in cells for _ in range(args.fixtures)]
    for spec in args.extend:
        cell, count, seed = spec.split(":")
        more = random.Random(int(seed))
        work += [(args.triple, args.counted, cell, more.getrandbits(64), more.getrandbits(64), f"extension:{seed}")
                 for _ in range(int(count))]
    rows = []
    with ProcessPoolExecutor(args.jobs) as pool:
        for row in pool.map(one, work):
            rows.append(row)
            print(json.dumps({k: row.get(k) for k in ("cell", "counted_over_triple", "counted_over_rho",
                                                      "trials_triple", "trials_counted", "unfinished")}),
                  flush=True)
    summary = {}
    extended = [spec.split(":")[0] for spec in args.extend]
    scopes = [(cell, cell) for cell in cells] + [(f"{cell}/registered-only", cell) for cell in extended]
    for name, cell in scopes:
        mine = [r for r in rows if r["cell"] == cell and (name == cell or r["stream"] == "registered")]
        summary_key = name
        if any(r["unfinished"] for r in mine):
            summary[summary_key] = {"comparison": None,
                                    "unfinished": [(r["target_seed"], r["unfinished"]) for r in mine if r["unfinished"]]}
            continue
        mean = lambda k: sum(r[k] for r in mine) / len(mine)
        summary[summary_key] = {key: summarise([r[key] for r in mine])
                                for key in ("counted_over_triple", "counted_over_rho", "triple_over_rho")}
        summary[summary_key]["fixtures"] = len(mine)
        summary[summary_key]["mean_ir"] = {arm: mean(f"ir_{arm}") for arm in ("triple", "counted", "rho")}
        summary[summary_key]["mean_trials"] = {arm: mean(f"trials_{arm}") for arm in ("triple", "counted")}
    Path(args.out).write_text(json.dumps({"seed": args.seed, "fixtures_per_cell": args.fixtures, "extend": args.extend,
                                          "summary": summary, "rows": rows}, indent=1) + "\n")
    for cell, s in summary.items():
        if s.get("comparison", 1) is None:
            print(f"{cell}: NO COMPARISON, unfinished runs {s['unfinished']}")
            continue
        parts = [f"{k} {v['geometric_mean']:.3f} [{v['ci95'][0]:.3f}, {v['ci95'][1]:.3f}]"
                 for k, v in s.items() if isinstance(v, dict) and "geometric_mean" in v]
        print(f"{cell}: " + "; ".join(parts) + f"; mean trials {s['mean_trials']}")


if __name__ == "__main__":
    main()
