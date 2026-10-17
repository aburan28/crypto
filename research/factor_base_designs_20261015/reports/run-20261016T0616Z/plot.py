"""Plot recorded values for Amendment 4 (C1-C4) from runs/*.jsonl only: field_ops,
D and max_cols against b. The dashed lines are reference slopes (m - 1 bits per
bit), drawn through each panel's first point; no fit is computed here.
Rendering aid only (AGENTS.md "no Python")."""
import glob, json, os
from collections import defaultdict
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
runs = os.path.join(here, "..", "..", "runs")


def load(pattern):
    out = []
    for f in sorted(glob.glob(os.path.join(runs, pattern))):
        out += [json.loads(l) for l in open(f) if l.strip()]
    return out


c1, c2, c4 = load("C1-rescore.jsonl"), load("C2-*.jsonl"), load("C4-*.jsonl")
COL = {"A0": "#2a78d6", "N": "#eb6834", "A1z": "#1baf7a", "A1g": "#eda100", "A4": "#e87ba4",
       "A11": "#008300", "N11": "#7a5195", "H1": "#c0392b", "H2": "#555555"}
MK = {"A0": "o", "N": "s", "A1z": "^", "A1g": "v", "A4": "D", "A11": "P", "N11": "X", "H1": "<", "H2": ">"}


def src(r):
    return os.path.basename(r.get("replay_of", ""))[:2]


def cells(rows, key):
    """(b, mean of key) per cell where every system is ok (and trace-matched)."""
    by = defaultdict(list)
    for r in rows:
        by[r["t"]].append(r)
    pts = []
    for t, rs in sorted(by.items()):
        if all(r["status"] == "ok" and r.get("trace_matches", True) and r.get(key) is not None for r in rs):
            pts.append((rs[0]["b"], sum(r[key] for r in rs) / len(rs)))
    return pts


plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": "#888", "xtick.color": "#444", "ytick.color": "#444"})

# Figure 1: cost per arm, p1 series, m = 2, 3, 4.
fig, axes = plt.subplots(1, 3, figsize=(11, 4.0))
for ax, m in zip(axes, [2, 3, 4]):
    first = None
    for a in ["A0", "N", "A1z", "A1g", "A4", "A11", "N11", "H1", "H2"]:
        if m == 4:
            rows = [r for r in c1 if r["arm"] == a and r["m"] == 4] + [r for r in c2 if r["arm"] == a]
        else:
            rows = [r for r in c1 if r["arm"] == a and r["m"] == m and r.get("replay_of_prime") == "p1"]
            rows += [r for r in c4 if r["arm"] == a and r["m"] == m]
        pts = cells(rows, "field_ops")
        if not pts:
            continue
        ax.plot([b for b, _ in pts], [v for _, v in pts], color=COL[a], marker=MK[a], ms=5, lw=1.3,
                mec="white", mew=0.7, label=a)
        if first is None or pts[0][1] > first[1]:
            first = pts[0]
    if first:
        b0, v0 = first
        xs = [b0, b0 + 4]
        ax.plot(xs, [v0, v0 * 2 ** ((m - 1) * 4)], ls="--", color="#999", lw=1, label=f"slope m − 1 = {m - 1}")
    ax.set_yscale("log", base=2)
    ax.set_title(f"m = {m},  p1", fontsize=9)
    ax.set_xlabel("b = log2|V| per summand (bits)")
    ax.set_ylabel("field_ops (cell mean)")
    ax.grid(color="#e5e5e5", lw=0.6)
    ax.legend(frameon=False, fontsize=6.5, ncol=2, loc="upper left")
fig.suptitle("F4 row-reduction field multiplications vs factor-base bits, completed cells (run-20261016T0616Z)",
             fontsize=9.5)
fig.text(0.5, 0.01, "Rows: runs/C1-rescore.jsonl (replays whose trace matched the committed row), runs/C2-*, runs/C4-*. "
         "Dashed: reference slope m − 1 bits per bit through the highest first point (C1's break-even).",
         ha="center", fontsize=7, color="#555")
fig.tight_layout(rect=(0, 0.04, 1, 0.95))
fig.savefig(os.path.join(here, "ops_vs_b.svg"))
fig.savefig(os.path.join(here, "ops_vs_b.png"), dpi=150)

# Figure 2: the ablation, D / max_cols / field_ops at m = 2, 3.
fig, axes = plt.subplots(3, 2, figsize=(9.5, 9.5))
for j, m in enumerate([2, 3]):
    for i, (key, lab, log) in enumerate([("D", "D = solving_degree_max", False),
                                          ("max_cols", "max_cols", True),
                                          ("field_ops", "field_ops", True)]):
        ax = axes[i][j]
        for a in ["A0", "A11", "H1", "H2"]:
            rows = ([r for r in c4 if r["arm"] == a and r["m"] == m] if a in ("H1", "H2") else
                    [r for r in c1 if r["arm"] == a and r["m"] == m and os.path.basename(r["replay_of"]).startswith(f"{a}-m")])
            pts = cells(rows, key)
            if pts:
                ax.plot([b for b, _ in pts], [v for _, v in pts], color=COL[a], marker=MK[a], ms=5, lw=1.3,
                        mec="white", mew=0.7, label=a)
            for r in rows:
                if r["status"] != "ok" and r.get(key) is not None:
                    ax.plot([r["b"]], [r[key]], marker=MK[a], ms=7, mfc="none", mec=COL[a], ls="none")
        if log:
            ax.set_yscale("log", base=2)
        ax.set_title(f"m = {m}, p1: {key}", fontsize=9)
        ax.set_xlabel("b (bits)")
        ax.set_ylabel(lab)
        ax.grid(color="#e5e5e5", lw=0.6)
        ax.legend(frameon=False, fontsize=7)
fig.suptitle("Ablation of A11 (C4): H1 = A11's curve relation + Kummer towers; H2 = A0's S3 + A11's chain "
             "(run-20261016T0616Z)", fontsize=9)
fig.text(0.5, 0.005, "Filled: completed cells (mean of 2 systems). Open (H1/H2 only): budget-stopped systems at the value reached. "
         "A0/A11 values from the C1 replays of the first run.", ha="center", fontsize=7, color="#555")
fig.tight_layout(rect=(0, 0.02, 1, 0.96))
fig.savefig(os.path.join(here, "ablation.svg"))
fig.savefig(os.path.join(here, "ablation.png"), dpi=150)
print("wrote ops_vs_b and ablation (svg, png)")
