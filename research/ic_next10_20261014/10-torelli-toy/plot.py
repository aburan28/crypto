"""Figure: observed carrier mass vs Torelli-codimension prediction."""
import os, json, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
rows = R["rows"]

# one point per (target class, arm); signs E0/E1 shown with different markers
groups = [("A3", "positive control A3 (dim 2)", "#2a7ab9"),
          ("B7", "B7 factors of A7 (dim 3)", "#8a5cc2"),
          ("A5", "target A5 (dim 4)", "#d1495b"),
          ("W5", "W5 (dim 5)", "#e39b2d"),
          ("A7", "A7 (dim 6)", "#555555")]
FLOOR = -16

fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
ax = axes[0]
for pre, lab, col in groups:
    xs, ys, zx, zy = [], [], [], []
    for r in rows:
        if not r["target"].startswith(pre) or r["arm"] == "ci5_sub" and pre != "A3":
            continue
        if r["arm"] == "ci5" and pre == "A3":
            continue
        pred = math.log2(r["pred_mass_h1"])
        if r["mass"] > 0:
            xs.append(pred); ys.append(math.log2(r["mass"]))
        else:
            zx.append(pred); zy.append(FLOOR)
    ax.scatter(xs, ys, c=col, s=28, label=lab, edgecolor="white", linewidth=0.5)
    ax.scatter(zx, zy, c=col, s=28, marker="v", alpha=0.7)
lo, hi = -40, 12
ax.plot([lo, hi], [lo, hi], color="#999", lw=1, ls="--", label="observed = predicted (h=1)")
ax.fill_between([lo, hi], [lo, hi], [lo + 3, hi + 3], color="#ccc", alpha=0.3, label="h in [1, 8]")
ax.axhline(FLOOR + 0.6, color="#bbb", lw=0.6)
ax.text(lo + 1, FLOOR + 0.9, "zero carriers (triangles)", fontsize=8, color="#666")
ax.set_xlim(lo, hi); ax.set_ylim(FLOOR - 1, hi)
ax.set_xlabel("log2 predicted carrier mass  h * M_family * 2^-c(G,d)")
ax.set_ylabel("log2 observed carrier mass")
ax.set_title("Carrier mass per (family, target, sign)")
ax.legend(fontsize=7.5, loc="upper left")

# right: log2(observed / predicted) for every nonzero (family, target, sign)
ax = axes[1]
offs = {"A3": -0.2, "B7": -0.07, "W3": 0.07, "A5": 0.0, "W5": 0.0, "A7": 0.2, "W7": 0.3}
cols = {"A3": "#2a7ab9", "B7": "#8a5cc2", "W3": "#2a9d8f", "A7": "#555555", "W7": "#111111"}
labs = {"A3": "A3 (dim 2, n=3 control)", "W3": "W3 (dim 3, n=3)", "B7": "B7 factors (dim 3, n=7)",
        "A7": "A7 (dim 6, n=7)", "W7": "W7 (dim 7, n=7)"}
seen = set()
for r in rows:
    pre = r["target"][:2]
    if r["mass"] <= 0 or pre not in cols or r["arm"] == "ci5_sub" and pre != "A3":
        continue
    if r["arm"] == "ci5" and pre == "A3":
        continue
    y = math.log2(r["mass"] / r["pred_mass_h1"])
    ax.scatter(r["genus"] + offs[pre], y, c=cols[pre], s=30, edgecolor="white", linewidth=0.5,
               label=labs[pre] if pre not in seen else None)
    seen.add(pre)
ax.axhline(0, color="#999", lw=1, ls="--")
ax.axhspan(0, 3, color="#ccc", alpha=0.3)
ax.text(2.8, 0.5, "heuristic band, h in [1, 8]", fontsize=8, color="#666")
ax.text(5.1, -4.5, "A5, W5: zero carriers in every family\n(predicted mass 2^-1.3 .. 2^-15)", fontsize=8,
        color="#d1495b")
ax.set_xlim(2.6, 8.6)
ax.set_xlabel("genus G")
ax.set_ylabel("log2(observed / predicted carrier mass)")
ax.set_title("Where carriers exist, the heuristic under-predicts")
ax.legend(fontsize=7.5, loc="upper left")
fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(os.path.join(HERE, f"torelli_toy.{ext}"), dpi=150)
print("ok")
