"""Lay out a recorded volcano-rho run: no computation beyond arranging the
values in results.json / trials.csv (AGENTS.md "Implementation language").

usage: python plot.py RUN_DIR
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "volcano-diagram"))
from volcano_diagram import render  # noqa: E402

run = Path(sys.argv[1])
R = json.loads((run / "results.json").read_text())
rows = list(csv.DictReader(open(run / "trials.csv")))
depth = R["problem"]["volcano_depth"]
ell = R["params"]["ell"]
exp = R["rho"]["expected_iters_sqrt_pi_n_over_2"]
rid = run.name
LEVEL_C = ["#2a6fdb", "#e07b00", "#2e9e5b", "#b03a8c", "#7a7a7a"]


def save(fig, name):
    fig.savefig(run / f"{name}.svg", bbox_inches="tight")
    fig.savefig(run / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# 1. The volcano, in the repository's standard diagram (tools/volcano-diagram).
arms = {a["vertex"]: a for a in R["arms"] if a["vertex"] is not None}
spec = {
    "schema": "volcano-diagram/v1",
    "ell": ell,
    "depth": depth,
    "title": f"{ell}-isogeny volcano, depth {depth}, p = {R['problem']['p']}",
    "subtitle": rid,
    "value_label": "mean rho iterations (measured curves)",
    "vertices": [
        {
            "id": v["id"],
            "label": f"v{v['id']}",
            "level": v["level"],
            "nbrs": v["nbrs"],
            "value": arms[v["id"]]["summary"]["iters"]["mean"] if v["id"] in arms else None,
            "highlight": "root" if v["id"] == 0 else None,
        }
        for v in R["vertices"]
    ],
}
(run / "volcano.diagram.json").write_text(json.dumps(spec, indent=1))
render(spec, str(run / "volcano"))

# 2. Per-curve min / mean / max (and every trial) by level.
order = sorted(R["arms"], key=lambda a: (a["kind"] != "isogenous", a["level"] if a["level"] is not None else 9, a["arm"]))
by_arm = defaultdict(list)
for r in rows:
    by_arm[r["arm"]].append(int(r["iters"]))
fig, ax = plt.subplots(figsize=(max(9, len(order) * 0.28), 4.6))
for i, a in enumerate(order):
    s = a["summary"]["iters"]
    c = LEVEL_C[a["level"] % len(LEVEL_C)] if a["kind"] == "isogenous" else "#000000"
    ys = by_arm[a["arm"]]
    ax.scatter([i] * len(ys), ys, s=3, alpha=0.25, color=c)
    ax.vlines(i, s["min"], s["max"], color=c, lw=1)
    ax.errorbar(i, s["mean"], yerr=s["ci95_half"], fmt="o", color=c, ms=4, capsize=2)
ax.axhline(exp, color="red", ls="--", lw=1, label=f"sqrt(pi n/2) = {exp:,.0f}")
ax.set_xticks(range(len(order)), [a["arm"] for a in order], rotation=90, fontsize=6)
ax.set_ylabel("group operations to solve (iterations)")
for lv in range(depth + 1):
    ax.scatter([], [], color=LEVEL_C[lv], label=f"level {lv}" + (" (crater)" if lv == 0 else " (floor)" if lv == depth else ""))
ax.scatter([], [], color="k", label="isomorphic control (root rescaled)")
ax.legend(fontsize=7, ncol=3, loc="upper right")
ax.set_title(f"Same DLP on every curve: per-curve min-max bar, mean ± 95% CI, {R['params']['trials']} trials each\n{rid}", fontsize=9)
save(fig, "per_curve")

# 3. Pooled per level: min / mean / max of iterations and of ms.
fig, axs = plt.subplots(1, 2, figsize=(9, 3.8))
groups = [(f"L{l['level']}\n({l['curves']} curves)", l["summary"]) for l in R["levels"]]
groups.append(("all\nisogenous", R["all_isogenous"]))
if R["controls"]:
    groups.append(("isomorphic\ncontrol", R["controls"]))
for ax, key, lab in [(axs[0], "iters", "iterations"), (axs[1], "ms", "wall time per solve (ms)")]:
    for i, (g, s) in enumerate(groups):
        v = s[key]
        ax.vlines(i, v["min"], v["max"], color="#555555")
        ax.scatter([i, i], [v["min"], v["max"]], marker="_", s=120, color="#555555")
        ax.errorbar(i, v["mean"], yerr=v["ci95_half"], fmt="o", color="#2a6fdb", capsize=3)
        ax.annotate(f"{v['mean']:,.0f}" if key == "iters" else f"{v['mean']:.0f}", (i, v["mean"]), xytext=(5, 0), textcoords="offset points", fontsize=7)
    if key == "iters":
        ax.axhline(exp, color="red", ls="--", lw=1)
    ax.set_xticks(range(len(groups)), [g for g, _ in groups], fontsize=7)
    ax.set_ylabel(lab)
fig.suptitle(f"Pooled by volcano level: min, max, mean ± 95% CI   ({rid})", fontsize=9)
save(fig, "per_level")
print("wrote", run)
