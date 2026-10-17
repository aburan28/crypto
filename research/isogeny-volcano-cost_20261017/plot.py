"""Draw the recorded results of PROTOCOL v1 (analysis.json from volcano-stats,
classes/*/curves.jsonl): no computation beyond arranging recorded values
(AGENTS.md "Implementation language").

usage: python plot.py STUDY_DIR [ANALYSIS_JSON]   (default analysis-v1.1.json)
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

study = Path(sys.argv[1])
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "volcano-diagram"))
from volcano_diagram import render  # noqa: E402

A = json.loads((study / (sys.argv[2] if len(sys.argv) > 2 else "analysis-v1.1.json")).read_text())
out = study / "figures"
out.mkdir(exist_ok=True)
classes = {c["class"]: Path(c["curves"]).parent.name for c in A["classes"]}
metric = {(m["arm"], m["metric"]): m for m in A["metrics"]}


def save(fig, name):
    fig.savefig(out / f"{name}.svg", bbox_inches="tight")
    fig.savefig(out / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# 1. Volcano heat maps per class: rho cost and IC yield, relative to the class mean.
rho = metric.get(("rho-neg", "S"))
icy = metric.get(("ic-mitm", "relation_yield"))
for ci, name in classes.items():
    rows = [json.loads(l) for l in open(study / "classes" / name / "curves.jsonl")]
    m0 = [r for r in rows if r["model"] == 0]
    depth = max(r["level"] for r in m0)
    ell = 2 if name.startswith(("l2", "j")) else 3
    val = {}
    for key, m in (("rho", rho), ("ic", icy)):
        for c in (m or {}).get("curves", []):
            if c["class"] == ci:
                val.setdefault(c["vertex"], {})[key] = 100 * (c["rel_mean"] - 1)
    spec = {
        "schema": "volcano-diagram/v2",
        "ell": ell,
        "depth": depth,
        "title": f"{name}: {ell}-isogeny volcano over F_{m0[0]['p']} · every curve × 4 isomorphic models × 16 paired targets",
        "subtitle": "PROTOCOL v1 (research/isogeny-volcano-cost_20261017) · % from the class mean",
        "metrics": [
            {"key": "rho", "title": "rho.negation · mean S", "label": "% from class mean S"},
            {"key": "ic", "title": "ic.pipeline mitm · relation yield", "label": "% from class mean yield"},
        ],
        "vertices": [
            {"id": r["vertex"], "label": f"v{r['vertex']}", "level": r["level"], "nbrs": r["nbrs"],
             "values": val.get(r["vertex"], {}), "highlight": "root" if r["vertex"] == 0 else None}
            for r in m0
        ],
        "footnote": "Colours span each panel's observed range; compare them with the per-curve 95% CIs "
                    "(about ±5.6% rho, ±3.5% yield) and the horizontal test in analysis.json before reading a pattern.",
    }
    (out / f"volcano-{name}.diagram.json").write_text(json.dumps(spec, indent=1))
    render(spec, str(out / f"volcano-{name}"))

# 2. Forest plot: vertical ratio per tested metric, with the A/A arm, against the ±margin band.
rows = [("rho-neg", "S", "rho S"), ("rho-neg-aa", "S", "rho S (A/A seeds)"),
        ("ic-mitm", "relation_yield", "IC relation yield"), ("ic-mitm", "S", "IC S (not tested)"),
        ("bsgs-neg", "S", "BSGS S (exact)")]
rows = [r for r in rows if (r[0], r[1]) in metric and metric[(r[0], r[1])]["vertical"]]
fig, ax = plt.subplots(figsize=(7, 0.6 * len(rows) + 1.4))
M = A["margin"]
ax.axvspan(1 - M, 1 + M, color="#dfe8f5", label=f"equivalence margin ±{M:.0%}")
ax.axvline(1, color="#888", lw=0.8)
for i, (arm, met, lab) in enumerate(rows):
    v = metric[(arm, met)]["vertical"]
    y = len(rows) - 1 - i
    ax.plot(v["ci95"], [y, y], color="#999", lw=1)
    ax.plot(v["ci90"], [y, y], color="#2b3445", lw=3)
    ax.plot(v["floor_over_crater"], y, "o", color="#d1495b")
    ax.annotate(f"{v['floor_over_crater']:.3f}  p={v['perm_p']:.3g}", (v["ci95"][1], y), xytext=(6, -3),
                textcoords="offset points", fontsize=7)
ax.set_yticks(range(len(rows)), [r[2] for r in reversed(rows)], fontsize=8)
ax.set_xlabel("floor / crater (class-normalised; thick 90%, thin 95% CI)")
ax.set_title("Vertical effect, pooled over six classes", fontsize=9)
ax.legend(fontsize=7, loc="lower right")
save(fig, "vertical_forest")

# 3. Horizontal and model heterogeneity: observed τ against the permutation null's 95th percentile.
rows = [(a, m) for (a, m) in [("rho-neg", "S"), ("rho-neg-aa", "S"), ("ic-mitm", "relation_yield"), ("ic-mitm", "S")]
        if (a, m) in metric]
fig, axs = plt.subplots(1, 2, figsize=(9, 0.55 * len(rows) + 1.4), sharey=True)
for ax, key, title in [(axs[0], "horizontal", "between curves of one level (horizontal)"),
                       (axs[1], "model", "between isomorphic models of one curve")]:
    for i, (arm, met) in enumerate(rows):
        h = metric[(arm, met)][key]
        y = len(rows) - 1 - i
        ax.barh(y, 100 * h["tau_null_p95"], color="#e2e2e2", height=0.6)
        ax.plot(100 * h["tau_rel"], y, "o", color="#d1495b")
        ax.annotate(f"Q p={h['perm_p']:.3g}, max|z| p={h['max_z_perm_p']:.3g}", (100 * h["tau_null_p95"], y),
                    xytext=(4, -3), textcoords="offset points", fontsize=6.5)
    ax.set_title(title, fontsize=8)
    ax.set_xlabel("τ, % of mean (dot); grey: null 95th percentile")
axs[0].set_yticks(range(len(rows)), [f"{a} · {m}" for a, m in reversed(rows)], fontsize=7)
save(fig, "heterogeneity")
print("wrote", out)
