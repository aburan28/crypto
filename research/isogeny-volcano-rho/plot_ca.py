"""Lay out a recorded `ca volcano --automorphisms` run (results.json,
trials.csv, volcano.diagram.json from aburan28/cryptanalysis): no
computation beyond arranging recorded values (AGENTS.md "Implementation
language").

usage: python plot_ca.py RUN_DIR
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
rid = run.name


def save(fig, name):
    fig.savefig(run / f"{name}.svg", bbox_inches="tight")
    fig.savefig(run / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# 1. The standard volcano diagram from ca's --diagram output.
spec = json.loads((run / "volcano.diagram.json").read_text())
# Vertex values from results.json's plain (isogenous) arms: a ca build
# before 31b72d2 wrote the folded arm's mean over the crater's plain one.
plain = {c["vertex"]: c["ops"]["mean"] for c in R["curves"] if c["kind"] == "isogenous"}
for v in spec["vertices"]:
    v["value"] = plain.get(v["id"])
spec["subtitle"] = f"{rid}  (colour: mean plain negation-map rho operations, ca)"
spec["vertices"][0]["highlight"] = "root, j = 0 (order-6 automorphism)"
render(spec, str(run / "volcano"))

# 2. Mean ± 95% CI by level, the controls, and the folded crater arm.
fig, ax = plt.subplots(figsize=(8, 4))
groups = [(f"L{l['level']}\n({l['curves']} curve{'s' if l['curves'] > 1 else ''})", l["ops"], "#2a6fdb")
          for l in R["levels"]]
groups.append(("isomorphic\ncontrols", R["controls"]["ops"], "#000000"))
groups.append(("j = 0 crater,\nautomorphism-\nfolded", R["automorphism"]["ops"], "#d1495b"))
for i, (g, s, c) in enumerate(groups):
    ax.vlines(i, s["min"], s["max"], color=c, alpha=0.35)
    ax.errorbar(i, s["mean"], yerr=s["ci95_half"], fmt="o", color=c, capsize=3)
    ax.annotate(f"{s['mean']:,.0f}", (i, s["mean"]), xytext=(6, 0), textcoords="offset points", fontsize=7)
ax.axhline(R["rho"]["expected_ops"], color="#2a6fdb", ls="--", lw=1, label=f"sqrt(pi n/4) = {R['rho']['expected_ops']:,.0f} (negation map)")
ax.axhline(R["rho"]["expected_ops"] / 3 ** 0.5, color="#d1495b", ls="--", lw=1, label=f"sqrt(pi n/12) = {R['rho']['expected_ops'] / 3 ** 0.5:,.0f} (order-6 folding)")
ax.set_xticks(range(len(groups)), [g for g, _, _ in groups], fontsize=7)
ax.set_ylabel("group operations to solve (ca ops)")
a = R["automorphism"]
ax.set_title(f"Same DLP on every curve: plain rho by level vs the folded walk on the j = 0 crater\n"
             f"speed-up {a['ops_speedup_vs_isogenous']:.2f} ± {a['ops_speedup_ci95_half']:.2f} (predicted {a['predicted']:.2f}); "
             f"{R['rho']['trials']} solves per curve; {rid}", fontsize=8)
ax.legend(fontsize=7, loc="upper right")
save(fig, "per_level")

# 3. Distribution on the crater curve: plain vs folded, same instance, same seeds.
by = defaultdict(list)
for r in csv.DictReader(open(run / "trials.csv")):
    by[r["arm"]].append(int(r["ops"]))
fig, ax = plt.subplots(figsize=(7, 3.6))
for arm, c, lab in [("v0", "#2a6fdb", "j = 0 crater, plain negation-map rho"),
                    ("v0+aut", "#d1495b", "j = 0 crater, automorphism-folded rho")]:
    xs = sorted(by[arm])
    ax.step(xs, [(i + 1) / len(xs) for i in range(len(xs))], where="post", color=c, label=f"{lab} (mean {sum(xs) / len(xs):,.0f})")
ax.set_xlabel("group operations to solve")
ax.set_ylabel("fraction of solves finished")
ax.set_title(f"Crater curve, same instance: cumulative distribution of {len(by['v0'])} solves per walk ({rid})", fontsize=8)
ax.legend(fontsize=7, loc="lower right")
save(fig, "crater_cdf")
print("wrote", run)
