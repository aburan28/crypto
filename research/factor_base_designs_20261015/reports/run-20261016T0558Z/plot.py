"""Plot recorded D and max_cols against b at m = 4 (A0, A11): F4 rows at
t = 2, 3 and the G1 rows (Amendment 4). Reads committed rows only; a
rendering aid (AGENTS.md "no Python")."""
import json, os
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
runs = os.path.join(here, "..", "..", "runs")
COL = {"A0": "#2a78d6", "A11": "#008300"}
MK = {"A0": "o", "A11": "P"}
OFF = {"A0": -0.05, "A11": 0.05}


def load(name):
    path = os.path.join(runs, name)
    if not os.path.exists(path):
        return []
    return [r for r in (json.loads(l) for l in open(path) if l.strip()) if "summary" not in r]


rows = {}
for a in ["A0", "A11"]:
    f4 = [r for r in load(f"F4-{a}-m4.jsonl") if r["t"] <= 3]
    rows[a] = [dict(r, src="F4") for r in f4] + [dict(r, src="G1") for r in load(f"G1-{a}-m4.jsonl")]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#888"})
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
for a in ["A0", "A11"]:
    ok = sorted((r["b"], r["D"], r["max_cols"]) for r in rows[a] if r["status"] == "ok")
    bad = [r for r in rows[a] if r["status"] != "ok"]
    xs = [b + OFF[a] for b, _, _ in ok]
    axes[0].plot(xs, [d for _, d, _ in ok], color=COL[a], marker=MK[a], ms=7, mec="white", lw=1.5, label=a)
    axes[1].plot(xs, [c for _, _, c in ok], color=COL[a], marker=MK[a], ms=7, mec="white", lw=1.5, label=a)
    for r in bad:
        x = r["b"] + OFF[a]
        if r.get("D") is not None:  # timeout: partial values (lower bounds)
            axes[0].plot([x], [r["D"]], marker=MK[a], ms=9, mfc="none", mec=COL[a], mew=1.5, ls="none")
            axes[1].plot([x], [r["max_cols"]], marker=MK[a], ms=9, mfc="none", mec=COL[a], mew=1.5, ls="none")
        else:  # memory abort: no D recorded
            for ax, y in ((axes[0], 9.5), (axes[1], None)):
                if y is not None:
                    ax.plot([x], [y], marker="x", ms=8, color=COL[a], mew=1.8, ls="none")
                    ax.annotate(f"{a} abort (memory)", (x, y), xytext=(-8, 0), textcoords="offset points",
                                ha="right", va="center", fontsize=7, color="#555")
axes[0].set_ylabel("D = solving_degree_max (f4_fp, grevlex)")
axes[0].set_ylim(1.5, 10)
axes[1].set_ylabel("max_cols (widest F4 matrix; log2 scale)")
axes[1].set_yscale("log", base=2)
for ax in axes:
    ax.set_xlabel("b = log2|V| per summand (bits)")
    ax.grid(color="#e5e5e5", lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
fig.suptitle("m = 4 on p1 = 786433: F4 rows (t = 2, 3) and G1 rows (t >= 4, 3600 s, 15 GB), run-20261016T0558Z",
             fontsize=9.5)
fig.text(0.5, 0.005, "Filled: completed systems (both targets). Open: stopped at the budget, partial values. "
         "Cross: memory abort, no D.", ha="center", fontsize=7.5, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig(os.path.join(here, "m4_d_cols_vs_b.svg"))
fig.savefig(os.path.join(here, "m4_d_cols_vs_b.png"), dpi=150)
print("wrote m4_d_cols_vs_b.svg/png")
