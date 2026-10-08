#!/usr/bin/env python3
"""Curve-vs-presentation variance split for the presentation sweep.

Two-way crossed random-effects model without replication,
    Y[E,V] = mu + u_E + w_V + e,
estimated by ANOVA mean squares on the complete curves x subspaces grid.
ICC_E = s2_u / (s2_u + s2_w + s2_e), 95 % CI by bootstrapping curves.
Pre-registered thresholds (design note): H0 if the CI upper bound < 0.05 on
every metric; H1 if ICC_E >= 0.10 with CI lower bound > 0 on any cost metric.

Usage: presentation_icc.py pres_*.jsonl
"""
import json, math, random, sys
from collections import defaultdict

METRICS = {
    "yield_per_F2": lambda r: (r["relations"] / r["probes"]) / r["factor_base_points"] ** 2,
    "log_us_per_call": lambda r: math.log(r["groebner_ns"] / r["groebner_calls"] / 1e3),
    "reductions_per_call": lambda r: r["reductions"] / r["groebner_calls"],
    "log_projected_solve": lambda r: math.log(
        (r["unknowns"] + 4) / max(r["relations"] / r["probes"], 1e-9)
        * r["groebner_ns"] / r["groebner_calls"]),
}

def icc(grid):
    # grid: list of rows (curves) of equal-length lists (subspaces)
    a, b = len(grid), len(grid[0])
    gm = sum(map(sum, grid)) / (a * b)
    rm = [sum(r) / b for r in grid]
    cm = [sum(grid[i][j] for i in range(a)) / a for j in range(b)]
    ss_r = b * sum((x - gm) ** 2 for x in rm)
    ss_c = a * sum((x - gm) ** 2 for x in cm)
    ss_t = sum((grid[i][j] - gm) ** 2 for i in range(a) for j in range(b))
    ms_r, ms_c = ss_r / (a - 1), ss_c / (b - 1)
    ms_e = (ss_t - ss_r - ss_c) / ((a - 1) * (b - 1))
    su = max((ms_r - ms_e) / b, 0.0)
    sw = max((ms_c - ms_e) / a, 0.0)
    tot = su + sw + ms_e
    return (su / tot if tot > 0 else 0.0), su, sw, ms_e

def main(paths):
    rows = [json.loads(l) for p in paths for l in open(p)]
    groups = defaultdict(list)
    for r in rows:
        if "skipped" in r or r.get("inconsistent", 0) or not r.get("bsgs_ok", True):
            continue
        groups[(r["n"], r["a2"], r["l"])].append(r)
    out = {}
    for key, rs in sorted(groups.items()):
        seeds = sorted({r["basis_seed"] for r in rs})
        by = defaultdict(dict)
        for r in rs:
            by[r["a6"]][r["basis_seed"]] = r
        curves = [c for c in by if all(s in by[c] for s in seeds)]
        res = {"curves": len(curves), "subspaces": len(seeds)}
        for name, f in METRICS.items():
            grid = [[f(by[c][s]) for s in seeds] for c in curves]
            point, su, sw, se = icc(grid)
            rng = random.Random(7)
            boots = sorted(icc([grid[rng.randrange(len(grid))] for _ in grid])[0] for _ in range(2000))
            res[name] = {"icc_curve": point, "ci95": [boots[50], boots[1949]],
                         "var_curve": su, "var_subspace": sw, "var_residual": se,
                         "share_subspace": sw / (su + sw + se) if su + sw + se else 0.0}
        out["n=%d a2=%d l=%d" % key] = res
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main(sys.argv[1:])
