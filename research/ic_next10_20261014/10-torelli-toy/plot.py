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

# right: A5 / W5 / A7 predicted mass per genus vs observation, cumulative over families
ax = axes[1]
for pre, lab, col in groups[2:]:
    by_g = {}
    for r in rows:
        if r["target"].startswith(pre) and r["arm"] != "ci5_sub":
            by_g.setdefault(r["genus"], [0.0, 0.0])
            by_g[r["genus"]][0] += r["pred_mass_h1"] / 2  # two signs -> per sign
            by_g[r["genus"]][1] += r["mass"]
    gs = sorted(by_g)
    ax.plot(gs, [math.log2(by_g[g][0]) for g in gs], "-o", color=col, label=f"{lab}: predicted (per sign, h=1)")
    ax.plot(gs, [math.log2(by_g[g][1]) if by_g[g][1] > 0 else FLOOR for g in gs], "v", color=col, ms=9,
            alpha=0.6)
pc = {}
for r in rows:
    if r["target"].startswith("A3") and r["arm"] != "ci5":
        pc.setdefault(r["genus"], [0.0, 0.0])
        pc[r["genus"]][0] += r["pred_mass_h1"] / 2
        pc[r["genus"]][1] += r["mass"] / 2
gs = sorted(pc)
ax.plot(gs, [math.log2(pc[g][0]) for g in gs], "-o", color="#2a7ab9", label="A3 control: predicted")
ax.plot(gs, [math.log2(pc[g][1]) for g in gs], "s", color="#2a7ab9", ms=8, label="A3 control: observed")
ax.axhline(0, color="#999", lw=0.8, ls=":")
ax.text(3.05, 0.3, "one curve class", fontsize=8, color="#666")
ax.set_ylim(FLOOR - 1, 8)
ax.set_xlabel("genus G (enumerated families only)")
ax.set_ylabel("log2 carrier mass (sum over enumerated families)")
ax.set_title("Toy ladder: triangles at bottom = zero observed")
ax.legend(fontsize=7.5, loc="lower left")
fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(os.path.join(HERE, f"torelli_toy.{ext}"), dpi=150)
print("ok")
