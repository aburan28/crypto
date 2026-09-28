#!/usr/bin/env python3
"""Walk lengths of the control against the worker's rho (PREREGISTRATION-rho-control.md, item 1).

    python3 research/ic_triple_counted_20260923/walk_lengths.py

Ratio of mean walk additions, control over rho, with a 95% bootstrap interval
over fixtures (10,000 resamples, seed 0), beside rho's expected steps
sqrt(pi r / 2) / sqrt(2n).  Added after the control ran: the registered check
is the point ratio, and this is the spread that ratio carries at 32 fixtures.
"""
import json
import math
import random
from pathlib import Path

CELLS = {"n23a1": (23, 4_196_903), "n37a0": (37, 230_603_167), "n43a1": (43, 4_644_189_029),
         "n59a0": (59, 10_063_074_221), "n61a1": (61, 11_514_943_771)}
rows = json.loads((Path(__file__).resolve().parent / "rho-control.json").read_text())["rows"]
rng = random.Random(0)
for cell, (n, r) in CELLS.items():
    mine = [x for x in rows if x["cell"] == cell]
    c = [x["additions_control"] for x in mine]
    o = [x["additions_rho"] for x in mine]
    boots = sorted(sum(c[i] for i in idx) / sum(o[i] for i in idx)
                   for idx in ([rng.randrange(len(mine)) for _ in mine] for _ in range(10_000)))
    print(f"{cell}: n={len(mine):3d}  mean additions control {sum(c)/len(c):7.0f}  rho {sum(o)/len(o):7.0f}  "
          f"expected {math.sqrt(math.pi * r / 2) / math.sqrt(2 * n):7.0f}  ratio {sum(c)/sum(o):.3f} "
          f"[{boots[250]:.3f}, {boots[9749]:.3f}]")
