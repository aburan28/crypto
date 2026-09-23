#!/usr/bin/env python3
"""The matched-canonical-form rho on the confirmation's own fixtures.

    python3 rho_control.py TRIPLE_WORKER CONTROL_WORKER --confirm confirm.json --out rho-control.json

CONTROL_WORKER is the triple tree with `rho-normal-basis.patch`; see
PREREGISTRATION-rho-control.md.  For every fixture in `confirm.json` the
control runs once in rho mode under callgrind and its report goes through
`oracle.py`'s rho mode; a correct scalar is the same logarithm the other arms
recovered, since the logarithm is unique.  The worker's current rho is re-run
natively (its walk is deterministic in the seeds) only to read its walk length
for the pre-registered check -- its instruction count is the committed one.
"""
import argparse
import importlib.util
import json
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "research/ic_candidate_tournament_20260915"))
import oracle  # noqa: E402

_spec = importlib.util.spec_from_file_location("triple_check", REPO / "research/ic_triple_table_20260923/check.py")
_triple_check = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_triple_check)
summarise = _triple_check.summarise

PAIR = {"batch_trials": 1, "linear_algebra": "sparse", "max_trials": 65536,
        "solver": "pair_table", "summands": 3}
CELLS = {"n23a1": (23, 1), "n37a0": (37, 0), "n43a1": (43, 1), "n59a0": (59, 0), "n61a1": (61, 1)}


def job(cell, tseed, aseed, mode):
    n, a = CELLS[cell]
    return {"degree": n, "curve_a": a, "target_seeds": [tseed], "algorithm_seed": aseed, "config": PAIR,
            "factor_base": {"kind": "subgroup_orbits", "seed": 43, "points": 222}, "mode": mode}


def run(worker, payload):
    return json.loads(subprocess.run([worker], input=json.dumps(payload), text=True,
                                     capture_output=True, timeout=1800).stdout)


def one(args):
    triple, control, row = args
    cell, t, a = row["cell"], row["target_seed"], row["algorithm_seed"]
    fixture = run(triple, job(cell, t, a, "fixture"))["fixture"]
    out = subprocess.run(["valgrind", "--tool=callgrind", "--callgrind-out-file=/dev/null", control],
                         input=json.dumps(job(cell, t, a, "rho")), text=True, capture_output=True, timeout=3600)
    ir = int(re.search(r"Collected\s*:\s*(\d+)", out.stderr).group(1))
    rep = json.loads(out.stdout)
    res = {"cell": cell, "stream": row["stream"], "target_seed": t, "algorithm_seed": a, "ir_control": ir,
           "status": rep.get("status")}
    if rep.get("status") != "complete":
        return res
    oracle.verify(rep, fixture, expected_mode="rho")
    old = run(triple, job(cell, t, a, "rho"))
    res.update(additions_control=rep["solutions"][0]["walk_group_additions"],
               additions_rho=old["solutions"][0]["walk_group_additions"],
               ir_rho=row["ir_rho"], ir_counted=row["ir_counted"],
               control_over_rho=ir / row["ir_rho"], counted_over_control=row["ir_counted"] / ir)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("triple")
    ap.add_argument("control")
    ap.add_argument("--confirm", default=str(HERE / "confirm.json"))
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rows = json.loads(Path(args.confirm).read_text())["rows"]
    results = []
    with ProcessPoolExecutor(args.jobs) as pool:
        for res in pool.map(one, [(args.triple, args.control, r) for r in rows]):
            results.append(res)
            print(json.dumps({k: res.get(k) for k in ("cell", "status", "counted_over_control",
                                                      "additions_control", "additions_rho")}), flush=True)
    summary = {}
    for cell in CELLS:
        mine = [r for r in results if r["cell"] == cell]
        if any(r["status"] != "complete" for r in mine):
            summary[cell] = {"comparison": None, "unfinished": [r["target_seed"] for r in mine if r["status"] != "complete"]}
            continue
        adds = sum(r["additions_control"] for r in mine) / sum(r["additions_rho"] for r in mine)
        summary[cell] = {"fixtures": len(mine),
                         "counted_over_control": summarise([r["counted_over_control"] for r in mine]),
                         "control_over_rho": summarise([r["control_over_rho"] for r in mine]),
                         "mean_additions_control_over_rho": adds,
                         "additions_ratio_by_fixture": summarise([r["additions_control"] / r["additions_rho"] for r in mine]),
                         "mean_ir_control": sum(r["ir_control"] for r in mine) / len(mine)}
    Path(args.out).write_text(json.dumps({"summary": summary, "rows": results}, indent=1) + "\n")
    for cell, s in summary.items():
        if s.get("comparison", 1) is None:
            print(f"{cell}: NO COMPARISON {s}")
            continue
        c, k = s["counted_over_control"], s["control_over_rho"]
        print(f"{cell}: counted/control {c['geometric_mean']:.3f} [{c['ci95'][0]:.3f}, {c['ci95'][1]:.3f}]; "
              f"control/rho {k['geometric_mean']:.3f} [{k['ci95'][0]:.3f}, {k['ci95'][1]:.3f}]; "
              f"mean additions control/rho {s['mean_additions_control_over_rho']:.3f}; n={s['fixtures']}")


if __name__ == "__main__":
    main()
