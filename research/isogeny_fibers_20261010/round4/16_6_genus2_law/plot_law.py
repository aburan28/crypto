"""Plot recorded census summaries (no computation beyond reading them).

Reads raw/census_*.json (written by examples/r4_16_6_genus2_census.rs),
copies the plotted summary rows to plot_data.json, and draws
law_measured_vs_predicted.svg: measured/predicted ratio for the class
factor c (m = 2, 3) and the nominal-column cost (m = 2, 3) per arm and
curve. A value of 1 means the measurement equals Prop. 15.2's prediction.
"""
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in sorted(glob.glob(os.path.join(HERE, "raw", "census_*.json"))):
    d = json.load(open(f))
    for a in d["summary"]:
        rows.append(dict(file=os.path.basename(f), p=d["p"], b=d["b"], h=d["h"], r=d["r"],
                         size=d["size"], arm=a["arm"], n_reps=a["n_reps"],
                         predicted=a["predicted"],
                         **{k: dict(value=a[k]["value"], ci95=a[k]["ci95"], smallx=a[k]["smallx"]) for k in
                            ("c2", "c3", "c2_strict", "c3_strict", "cost2_nominal",
                             "cost3_nominal", "columns", "hit_rate3")}))
json.dump(rows, open(os.path.join(HERE, "plot_data.json"), "w"), indent=1)
scan = [r for r in rows if r["arm"] == "A"]
rows = [r for r in rows if "Aonly" not in r["file"]]

panels = [("c2_strict", 1, "class factor c, m = 2 (distinct summands)"),
          ("c3", 2, "class factor c, m = 3"),
          ("cost2_nominal", 3, "mitm cost, m = 2 (nominal plain columns)"),
          ("cost3_nominal", 4, "mitm cost, m = 3 (nominal plain columns)")]
fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharey=True)
colors = {1051: "#4C78A8", 2131: "#F58518", 4111: "#54A24B"}
markers = {50: "v", 80: "o", 100: "s", 160: "D", 200: "^"}
labels = [f"{r['arm']}\n{'C1' if 'c1' in r['file'] else 'C2'} |F|={r['size']}" for r in rows]
order = []
for lab in labels:
    key = lab.split("\n")[0] + "\n" + lab.split("\n")[1]
    if key not in order:
        order.append(key)
for ax, (key, idx, title) in zip(axes.flat, panels):
    seen = set()
    for r, lab in zip(rows, labels):
        x = order.index(lab)
        pred = r["predicted"][idx]
        m = r[key]["value"] / pred
        lo, hi = (v / pred for v in r[key]["ci95"])
        off = {1051: -0.22, 2131: 0.0, 4111: 0.22}[r["p"]]
        lbl = f"p = {r['p']}" if r["p"] not in seen else None
        seen.add(r["p"])
        ax.errorbar(x + off, m, yerr=[[m - lo], [hi - m]], fmt=markers[r["size"]],
                    color=colors[r["p"]], ms=5, capsize=2, label=lbl)
    ax.axhline(1.0, color="black", lw=0.8)
    ax.set_yscale("log")
    ax.set_ylim(0.3, 1.6)
    ax.set_yticks([0.3, 0.5, 0.7, 1.0, 1.4])
    ax.set_yticklabels(["0.3", "0.5", "0.7", "1", "1.4"])
    ax.set_xticks(range(len(order)))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xticklabels([o.replace("\n", " ") for o in order], fontsize=7, rotation=60, ha="right")
    ax.set_title(title, fontsize=9)
    ax.grid(axis="y", alpha=0.3)
axes[0][0].set_ylabel("measured / predicted (ratio to plain base)")
axes[1][0].set_ylabel("measured / predicted (ratio to plain base)")
axes[0][0].legend(fontsize=7, loc="lower left")
fig.suptitle("Genus-2 census (y² = x⁵ + b over F_p): measured ÷ Prop. 15.2 prediction; "
             "ratio of means over 40 random bases, 95% bootstrap CI", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "law_measured_vs_predicted.svg"))

# Figure 2: the A-arm class factor at m = 2 against |F|, with the derived
# finite-size line 1 + (h-1)(|A|-2)/(2|F|) recorded in finite_size_A.json.
fs = json.load(open(os.path.join(HERE, "finite_size_A.json")))
fig2, ax2 = plt.subplots(figsize=(6.4, 4.2))
for curve, col in (("C1", "#E45756"), ("C2", "#4C78A8")):
    pts = [r for r in scan if ("c1" in r["file"]) == (curve == "C1")]
    for r in pts:
        v = r["c2"]["value"]
        lo, hi = r["c2"]["ci95"]
        ax2.errorbar(r["size"] * {1051: 0.97, 2131: 1.0, 4111: 1.03}[r["p"]], v,
                     yerr=[[v - lo], [hi - v]], fmt=markers.get(r["size"], "o"),
                     color=col, capsize=2, ms=5)
    line = fs[curve]
    ax2.plot(line["size"], line["derived_c2_ratio"], color=col, lw=1, ls="--",
             label=f"{curve}: derived {line['formula']}")
ax2.axhline(1.0, color="black", lw=0.8)
ax2.set_xscale("log")
ax2.xaxis.set_minor_formatter(NullFormatter())
ax2.set_xticks([50, 80, 100, 160, 200])
ax2.set_xticklabels(["50", "80", "100", "160", "200"])
ax2.set_xlabel("|F| (base size)")
ax2.set_ylabel("c(A-closed) / c(plain) at m = 2")
ax2.set_title("A-closed base, m = 2: measured (points, 95% bootstrap CI, 40 bases)\n"
              "vs derived structured-pair correction (dashed)", fontsize=9)
ax2.legend(fontsize=7)
fig2.tight_layout()
fig2.savefig(os.path.join(HERE, "A_pair_correction.svg"))
print("rows", len(rows), "scan", len(scan))
