"""Plot recorded values from results/run-*.jsonl. Layout only: no result is computed here.
Usage: python3 plot.py <results.jsonl> <outdir>"""
import json, math, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src, out = sys.argv[1], sys.argv[2]
rows = [json.loads(l) for l in open(src)]
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1a1a19", "#6b6a63", "#e4e3dc"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})

def save(fig, name):
    fig.savefig(f"{out}/{name}.svg", bbox_inches="tight")
    fig.savefig(f"{out}/{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

# Figure 1: E1, size of the factor-base span against the size of J(K)
e1 = [r for r in rows if r["exp"] == "E1"]
fields = sorted({(r["q"], r["n"]) for r in e1}, key=lambda t: (t[1], t[0]))
fig, ax = plt.subplots(figsize=(7.2, 3.6))
for i, (q, n) in enumerate(fields):
    for kind, col, mk, dx in (("weak", BLUE, "o", -0.18), ("random", ORANGE, "s", 0.18)):
        pts = [r for r in e1 if (r["q"], r["n"], r["kind"]) == (q, n, kind)]
        for j, r in enumerate(pts):
            ax.scatter(i + dx + 0.07 * (j - 1), math.log10(r["H"]) / math.log10(r["J_K_counts"]), s=42, color=col, marker=mk,
                       edgecolor="white", linewidth=1.0, zorder=3, label=None)
    sub = 2 if n % 2 == 0 else 1
    ax.hlines(sub / n, i - 0.4, i + 0.4, color=MUTED, linewidth=1, linestyle=(0, (3, 2)), zorder=2)
    ax.text(i + 0.42, sub / n, f"{sub}/{n}", va="center", ha="left", color=MUTED, fontsize=8)
ax.scatter([], [], s=42, color=BLUE, marker="o", label="weak curve, paper's factor base (3 per field)")
ax.scatter([], [], s=42, color=ORANGE, marker="s", label="random curve, factor base x in F_q (3 per field)")
ax.plot([], [], color=MUTED, linewidth=1, linestyle=(0, (3, 2)), label="log #J_C0(F_q^gcd(2,n)) / log #J(K), about gcd(2,n)/n")
ax.set_xticks(range(len(fields))); ax.set_xticklabels([f"q={q}, n={n}" for q, n in fields])
ax.set_ylabel("log |H| / log #J(K)"); ax.set_ylim(0, 1.08); ax.set_xlim(-0.6, len(fields) - 0.2)
ax.grid(axis="y", color=GRID, linewidth=0.8); ax.set_axisbelow(True)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
ax.set_title("E1: subgroup H spanned by the factor base, as a fraction of J(K) on a log scale", loc="left", color=INK, fontsize=10)
save(fig, "e1_span_fraction")

# Figure 2: E2, relation candidates and span per factor-base line (one weak curve per field)
e2 = [r for r in rows if r["exp"] == "E2"]
fields2 = sorted({(r["q"], r["n"]) for r in e2}, key=lambda t: (t[1], t[0]))
kinds = ["weak", "scaling", "shift", "affine", "moebius"]
fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.4), sharex=True)
for i, (q, n) in enumerate(fields2):
    for k, kind in enumerate(kinds):
        pts = [r for r in e2 if (r["q"], r["n"], r["line"]) == (q, n, kind)]
        col = BLUE if kind == "weak" else ORANGE
        for j, r in enumerate(pts):
            x = i * 6 + k + 0.2 * (j - (len(pts) - 1) / 2)
            axes[0].scatter(x, r["candidates"], s=36, color=col, marker="o" if kind == "weak" else "s", edgecolor="white", linewidth=1.0, zorder=3)
            axes[1].scatter(x, math.log10(r["H"]) / math.log10(r["J_K"]), s=36, color=col, marker="o" if kind == "weak" else "s", edgecolor="white", linewidth=1.0, zorder=3)
axes[0].set_yscale("symlog", linthresh=1); axes[0].set_ylabel("relation candidates (c, lambda)\nwith lambda^2 - c f_L in F_q[t]")
axes[1].set_ylabel("log |H| / log #J(K)"); axes[1].set_ylim(0, 1.08)
ticks = [i * 6 + k for i in range(len(fields2)) for k in range(5)]
axes[1].set_xticks(ticks); axes[1].set_xticklabels(kinds * len(fields2), rotation=60, ha="right", fontsize=7.5)
for i, (q, n) in enumerate(fields2):
    axes[0].text(i * 6 + 2, 1.02, f"q={q}, n={n}", transform=axes[0].get_xaxis_transform(), ha="center", color=INK, fontsize=9)
for a in axes: a.grid(axis="y", color=GRID, linewidth=0.8); a.set_axisbelow(True)
axes[0].scatter([], [], s=36, color=BLUE, marker="o", label="the weak line beta*F_q"); axes[0].scatter([], [], s=36, color=ORANGE, marker="s", label="other lines (3 each; Moebius uses the sextic model)")
axes[0].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.32), ncol=2)
save(fig, "e2_lines")
