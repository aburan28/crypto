"""Plot recorded values from raw/summary.json (written by the Rust `analyze`
subcommand).  Computes nothing new.   python3 -I plot.py"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
S = json.loads((HERE / "raw" / "summary.json").read_text())
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

ARMS = [("d4-ed", "D₄-ED (first-dependency stop), candidate A", "#2a6fdb"),
        ("d4-ed-all", "D₄-ED, all rows (exploratory)", "#74a9f0"),
        ("d4-adapt", "adaptive α (candidate B)", "#2b8a3e"),
        ("d4-2lp", "item-01 D₄ + 2LP, α = 0.95", "#d9480f"),
        ("d4-late110", "positive control (first LA at 1.1·cols)", "#999999")]


def save(fig, name):
    fig.tight_layout()
    for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
        fig.savefig(FIG / f"{name}.{ext}", **kw)


# Figure 1: paired ratios against D₄, per size, three phases.
fig, axs = plt.subplots(1, 3, figsize=(13, 4.3), sharex=True)
for ax, ph, title in zip(axs, ("whole", "relation", "la"),
                         ("whole charged cost", "relation phase", "linear algebra")):
    for k, (arm, label, c) in enumerate(ARMS):
        P = S["paired_vs_d4"][arm][ph]
        off = (k - 2) * 0.08
        for s in P["per_size"]:
            ax.scatter([s["log2_E"] + off] * len(s["per_curve"]), s["per_curve"], s=8, color=c, alpha=0.4)
        xs = [s["log2_E"] + off for s in P["per_size"]]
        ys = [s["ratio"] for s in P["per_size"]]
        lo = [s["ratio"] - s["ci95"][0] for s in P["per_size"]]
        hi = [s["ci95"][1] - s["ratio"] for s in P["per_size"]]
        pooled = P["pooled"]
        ax.errorbar(xs, ys, yerr=[lo, hi], color=c, lw=1.2, capsize=2,
                    label=f"{label}: {pooled['ratio']:.3f} [{pooled['ci95'][0]:.3f}, {pooled['ci95'][1]:.3f}]")
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("log₂ #E(F_{p⁴})")
    ax.set_title(title)
axs[0].set_ylabel("arm / D₄ (paired: same curve, target, stream)")
for ax in axs:
    ax.legend(fontsize=6.3, frameon=False, loc="best")
fig.suptitle("Item 04: paired ratios vs the D₄ incumbent, 27 held-out curves (p = 101…1601, 3 per size); "
             "bars = 95 % bootstrap over the 3 curves at each size; legend = pooled ratio [95 % CI]", fontsize=9.3)
save(fig, "paired_ratio")

# Figure 2: relations used per column, and the hindsight first dependency.
fig, ax = plt.subplots(figsize=(7.5, 4.3))
rows = S["rows"]
x = [r["log2_E"] for r in rows]
cols = [r["y_columns"] + 1 for r in rows]
ser = [("d4", "D₄ incumbent", "#d9480f", "o"), ("d4-ed", "D₄-ED", "#2a6fdb", "s"),
       ("d4-late110", "positive control", "#999999", "^")]
for arm, label, c, m in ser:
    ax.scatter(x, [r["arms"][arm]["full_relations"] / n for r, n in zip(rows, cols)], color=c, marker=m, s=18, label=label)
ax.scatter(x, [r["d4_first_dependency_rows_hindsight"] / n for r, n in zip(rows, cols)], facecolors="none",
           edgecolors="k", s=40, label="first linear dependency of D₄'s rows (hindsight, exact)")
ax.axhline(1.0, color="k", lw=0.6, ls=":")
ax.set_xlabel("log₂ #E(F_{p⁴})")
ax.set_ylabel("relations used / (columns + 1)")
ax.legend(fontsize=7.5, frameon=False)
ax.set_title("Relations each stop rule needs (same row sequence)", fontsize=10)
save(fig, "relations")

# Figure 3: phase costs and fitted exponents.
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
for ax, ph in zip(axs, ("relation", "la")):
    for arm, label, c in [("d4", "D₄ (incumbent)", "#d9480f"), ("d4-ed", "D₄-ED", "#2a6fdb"),
                          ("d4-2lp", "D₄ + 2LP α = 0.95", "#2b8a3e")]:
        ys = [r["arms"][arm][f"{ph}_gae"] for r in rows]
        f = S["fits"][arm][ph]
        ax.scatter(x, ys, s=10, color=c, alpha=0.7,
                   label=f"{label}: slope {f['slope']:.3f} [{f['ci95'][0]:.3f}, {f['ci95'][1]:.3f}]")
    ax.set_yscale("log", base=2)
    ax.set_xlabel("log₂ #E(F_{p⁴})")
    ax.set_ylabel(f"{'relation' if ph == 'relation' else 'linear algebra'} cost (GAE)")
    ax.legend(fontsize=7, frameon=False)
fig.suptitle("Phase costs by curve; slope = OLS exponent in #E with 95 % bootstrap CI", fontsize=9.5)
save(fig, "phases")
print("ok")
