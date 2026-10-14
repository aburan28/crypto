"""Lay out recorded values as SVG graphs. Reads summary.jsonl (written by
`rho_walk_leads summary`); computes nothing beyond placing those numbers."""
import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

run = pathlib.Path(sys.argv[1])
rows = [json.loads(l) for l in (run / "summary.jsonl").read_text().splitlines() if l.strip()]
plt.rcParams.update({"font.size": 9, "svg.fonttype": "none"})

LABEL = {
    ("e1", "pm_coin"): "P±R, sign by hash of P (control)",
    ("e1", "pm_score"): "best of P+R, P−R",
    ("e1", "two_score"): "best of P+R_i, P+R_j",
    ("e1", "best4_score"): "best of P±R_i, P±R_j",
    ("rf", "k2"): "random function, best of 2",
    ("rf", "k4"): "random function, best of 4",
    ("e2", "lazy_dp"): "DP test on lazy x",
    ("e2", "lazy_index"): "table index on lazy x",
    ("e2", "pair_table"): "pair table R_i+R_j",
    ("e2", "shared_index"): "index shared across walks",
    ("e4", "neg_canon"): "negation map, fixed rule",
    ("e4", "neg_hashrule"): "negation map, rule per point",
    ("e4", "neg_epoch"): "negation map, rule per epoch",
}
for k, v in list(LABEL.items()):
    if k[0] == "e1":
        LABEL[("e1dp", "sel_" + k[1])] = v
        LABEL[("e1shift", k[1])] = v
REF = {"e1": "plain r-adding walk, doubling escape", "e1shift": "plain r-adding walk, shift escape",
       "e1dp": "plain walk, distinguished points, 64 walks", "rf": "random function", "e2": "canonical DP walk", "e4": "no negation map"}

panels = [("e1", 36), ("e1", 40), ("e1shift", 36), ("e1dp", 36), ("rf", 36), ("e2", 36), ("e4", 36)]
fig, axes = plt.subplots(len(panels), 1, figsize=(7.2, 12.5), sharex=True)
for ax, (exp, bits) in zip(axes, panels):
    sel = [r for r in rows if r.get("of") == exp and r["bits"] == bits and "speedup_ref_over_arm" in r]
    for y, r in enumerate(sel):
        lo, hi = r["ci95"]
        v = r["speedup_ref_over_arm"]
        ax.errorbar(v, y, xerr=[[v - lo], [hi - v]], fmt="o", color="#2f6db5", ms=4, capsize=2)
        if "predicted" in r:
            ax.plot(r["predicted"], y, marker="|", color="#c0392b", ms=12, mew=2)
    ax.set_yticks(range(len(sel)), [f'{LABEL.get((exp, r["arm"]), r["arm"])} (n={r["paired_n"]})' for r in sel])
    ax.axvline(1.0, color="#888", lw=0.8)
    ax.set_title(f"{exp}, N ≈ 2^{bits}: reference = {REF[exp]}", loc="left", fontsize=9)
    ax.set_xscale("log")
axes[-1].set_xlabel("steps speedup, reference mean / arm mean (paired, 95% bootstrap); red tick = derived prediction")
fig.tight_layout()
fig.savefig(run / "speedups.svg")

# cost side: measured ns per batched step and the derived op counts
bench = [r for r in rows if r.get("of") == "bench"]
fig, ax = plt.subplots(figsize=(7.2, 2.6))
arms = ["base", "pm_score", "best4"]
opmodel = {"base": 6.0, "pm_score": 8.0, "best4": 15.0}  # M-equivalents with S = M (README §Cost model)
xs = range(len(arms))
for k, bits in enumerate(sorted({r["bits"] for r in bench})):
    vals = [next((r["cost_ratio_vs_base"] for r in bench if r["bits"] == bits and r["arm"] == a), None) for a in arms]
    name = {64: "62-bit p (1 limb)", 256: "P-256 p (4 limbs)"}[bits]
    ax.bar([x + 0.25 * k - 0.25 for x in xs], vals, width=0.25, label=f"measured, {name}")
ax.bar([x + 0.25 for x in xs], [opmodel[a] / opmodel["base"] for a in arms], width=0.25, label="op-count model (S = M)", color="#999")
ax.set_xticks(list(xs), ["plain step", "best of P±R", "best of P±R_i, P±R_j"])
ax.set_ylabel("cost per step / plain step")
ax.legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig(run / "cost.svg")
print("wrote", run / "speedups.svg", run / "cost.svg")
