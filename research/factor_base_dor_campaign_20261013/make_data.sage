#!/usr/bin/env sage
# Consolidates recorded raw-cell aggregates into figures/data.json for the
# plot layout script. Values are copied from the raw JSON reports; the only
# arithmetic is averaging over seeds, as analysis.txt reports.
import json, glob

grid = {}
for p in sorted(glob.glob("raw/m2_n*_d*_s*.json")):
    d = json.load(open(p))
    for r in d["results"]:
        if "skipped" in r or "d_star_mean" not in r:
            continue
        n = d["field_degree"]; dim = r["dim"]; st = r["structure"]
        works = []
        for t in r["targets"]:
            if t["d_star"] is None: continue
            sh = t.get("macaulay_shape", {}).get(str(t["d_star"]))
            if sh: works.append(sh[0] * sh[1])
        grid.setdefault((n, st, dim), []).append({
            "d_mean": r["d_star_mean"], "d_min": r["d_star_min"], "d_max": r["d_star_max"],
            "frac": r["decomposable_fraction"],
            "work": sum(works) / len(works) if works else None,
        })

cells = []
for (n, st, dim), rs in sorted(grid.items()):
    def m(k): return sum(r[k] for r in rs) / len(rs)
    w = [r["work"] for r in rs if r["work"] is not None]
    cells.append({"n": n, "structure": st, "dim": dim, "seeds": len(rs),
                  "d_star_mean": m("d_mean"), "d_star_min": min(r["d_min"] for r in rs),
                  "d_star_max": max(r["d_max"] for r in rs),
                  "decomposable_fraction": m("frac"),
                  "macaulay_work_at_dstar": (sum(w) / len(w)) if w else None,
                  "work_per_relation_proxy": ((sum(w) / len(w)) / m("frac")) if (w and m("frac") > 0) else None})

deltas = {}
for p in sorted(glob.glob("raw/delta_n*_m2.json")):
    d = json.load(open(p))
    n = int(p.split("_")[1][1:])
    deltas["n%d" % n] = [{"delta": r["delta"], "dim": r.get("dim"),
                          "d_star_mean": r.get("d_star_mean"),
                          "decomposable_fraction": r.get("decomposable_fraction"),
                          "work_per_relation_proxy": r.get("work_per_relation_proxy")}
                         for r in d["rows"] if "skipped" not in r]

sym4 = {}
for p in sorted(glob.glob("raw/sym4_*.json")):
    d = json.load(open(p))
    sym4[p.split("/")[-1][:-5]] = [
        {"structure": r["structure"], "nvars": r.get("nvars"),
         "d_star_mean": r.get("d_star_mean"),
         "point_decomposable_fraction": r.get("point_decomposable_fraction"),
         "boolean_unsat_fraction": r.get("boolean_unsat_fraction")}
        for r in d["results"] if "skipped" not in r]

stack = {}
for p in sorted(glob.glob("raw/stack_n*.json")):
    stack[p.split("/")[-1][:-5]] = json.load(open(p))["rows"]

out = {"schema": "factor-base-dor-campaign-20261013-v1",
       "provenance": "aggregates of raw/m2_* (seeds 1-3), raw/delta_*, raw/sym4_*, raw/stack_*; engine sage/factor_base_dor.sage",
       "m2_grid": cells, "delta_sweeps_m2": deltas, "symmetric_s4": sym4, "stacking": stack}
with open("figures/data.json", "w") as f:
    json.dump(out, f, indent=1)
print("wrote figures/data.json:", len(cells), "m2 cells")
