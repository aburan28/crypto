# Lays out the sealed session (sessions/k4-p29/records.jsonl, written by
# ecbench) as one figure: S per arm on a log axis with the IC arm split by
# phase, and §11.16's model prediction for the IC arm marked.  Recorded
# values only; the model point is the protocol's registered E3 figure.
import json, os
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
rows = [json.loads(l) for l in open(os.path.join(here, "..", "sessions", "k4-p29", "records.jsonl"))]
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "ink": "#0b0b0b", "ink2": "#52514e",
     "grid": "#e6e5e1", "surface": "#fcfcfb", "wash": "#9aa3ad", "plum": "#8e5bd1"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"],
                     "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.4, 3.6), gridspec_kw={"width_ratios": [1, 1.1]})
names = [r["arm"] for r in rows]
s = [r["cost"]["s"] for r in rows]
cols = [C["wash"] if "rho" in n else C["blue"] for n in names]
a1.bar(range(len(rows)), s, color=cols, width=0.6)
for i, v in enumerate(s):
    if v > 1e3:
        a1.text(i, v / 8, f"{v:.3g}", ha="center", fontsize=8, color="white")
    else:
        a1.text(i, v * 1.6, f"{v:.3g}", ha="center", fontsize=8, color=C["ink"])
rho = next(r for r in rows if r["arm"] == "rho-neg")
model = 1.3e9 * rho["cost"]["s"]
ic_i = names.index("ic-k4-d4")
a1.plot([ic_i - 0.35, ic_i + 0.35], [model, model], color=C["orange"], lw=1.6)
a1.text(ic_i + 0.38, model, "model (E3)", fontsize=7.5, color=C["orange"], va="center")
a1.set_yscale("log")
a1.set_ylim(0.3, 3e10)
a1.set_xticks(range(len(rows)), names, fontsize=8)
a1.set_ylabel("S = group-addition equivalents / √r")
a1.grid(axis="y", color=C["grid"], lw=0.6)
a1.set_title("E(F_29⁴), r = 44,263, one target", fontsize=9)
ic = rows[ic_i]
ph = [(p["name"], p["gae"]) for p in ic["phases"]]
tot = sum(g for _, g in ph)
a2.barh(range(len(ph)), [max(g, 1e-1) for _, g in ph], color=C["blue"])
for i, (n, g) in enumerate(ph):
    a2.text(max(g, 1e-1) * 1.5, i, f"{g:.3g} ({100 * g / tot:.4g} %)", va="center", fontsize=7.5, color=C["ink"])
a2.set_xscale("log")
a2.set_xlim(0.1, 1e14)
a2.set_yticks(range(len(ph)), [n for n, _ in ph], fontsize=8)
a2.invert_yaxis()
a2.set_xlabel("group-addition equivalents (log)")
a2.grid(axis="x", color=C["grid"], lw=0.6)
a2.set_title("the IC arm by phase", fontsize=9)
fig.tight_layout()
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(os.path.join(here, f"s_by_arm.{ext}"), **kw)
