#!/usr/bin/env python3
"""Curve-vs-presentation ICC on measured full solves (results.jsonl).

Grid per (n, V family): curves x subspaces; families geom (g1, g2) and rand
(1, 2), plus all five V together.  Metrics: log total seconds (yield phase +
closure + linear algebra), log groebner seconds, closure trials to solve.
Model and bootstrap as presentation_icc.py.
"""
import importlib.util, json, math, random, sys
from collections import defaultdict
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "picc", Path(__file__).resolve().parents[2] / "research/notes/koblitz-isogeny/presentation_icc.py")
picc = importlib.util.module_from_spec(spec); spec.loader.exec_module(picc)

M = {
    "log_seconds": lambda d: math.log(d["seconds"]),
    "log_groebner_seconds": lambda d: math.log(d["groebner_seconds"]),
    "closure_trials": lambda d: d["closure_trials"],
}
rows = [json.loads(l) for l in open(Path(__file__).with_name("results.jsonl"))]
cells = defaultdict(dict)
for r in rows:
    n, a2, a6, l, v = r["key"].split()
    assert r["row"]["verified"] and r["row"]["log"] == r["row"]["planted"]
    cells[int(n)].setdefault(int(a6), {})[v] = r["row"]
out = {}
for n in sorted(cells):
    for fam, vs in [("geom", ["g1", "g2"]), ("rand", ["1", "2"]), ("all5", ["mono", "g1", "g2", "1", "2"])]:
        curves = [c for c in cells[n] if all(v in cells[n][c] for v in vs)]
        for name, f in M.items():
            grid = [[f(cells[n][c][v]) for v in vs] for c in curves]
            est = picc.icc(grid)[0]
            rng = random.Random(1)
            boots = sorted(picc.icc([grid[rng.randrange(len(grid))] for _ in grid])[0] for _ in range(2000))
            ci = (boots[50], boots[1949])
            out[f"{n}/{fam}/{name}"] = {"curves": len(curves), "icc": est, "ci95": ci}
            print(f"n={n} {fam:5s} {name:22s} curves={len(curves)} ICC={est:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]")
json.dump(out, open(Path(__file__).with_name("icc_full_solves.json"), "w"), indent=1)
