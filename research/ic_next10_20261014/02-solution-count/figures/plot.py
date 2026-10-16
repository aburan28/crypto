"""Lay out recorded census values (data/ratios.csv, data/validation/*) as figures. No new computation
beyond averaging recorded per-curve values for display."""
import csv, math, sys, collections, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
rows = list(csv.DictReader(open(HERE.parent / "data" / "ratios.csv")))
ARMS = ["kummer", "subgroup", "dickson", "dickson-torus", "smallx", "random", "onlycurve", "dloginterval"]
COL = {"kummer": "#2a78d6", "subgroup": "#eb6834", "dickson": "#1baf7a", "dickson-torus": "#eda100",
       "smallx": "#e87ba4", "random": "#555555", "onlycurve": "#6250d6", "dloginterval": "#e34948"}
MK = {"kummer": "o", "subgroup": "s", "dickson": "^", "dickson-torus": "v", "smallx": "D",
      "random": "x", "onlycurve": "P", "dloginterval": "*"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e4e4e4", "grid.linewidth": 0.6})

# Figure 1: mean count per target vs B, mean over the 3 curves of each size, with reference lines.
fig, axes = plt.subplots(2, 3, figsize=(11, 6.6), sharex=True)
for i, m in enumerate([2, 3]):
    for j, bits in enumerate([20, 24, 28]):
        ax = axes[i][j]
        for arm in ARMS:
            if arm == "dloginterval":
                continue
            acc = collections.defaultdict(list)
            for r in rows:
                if r["arm"] == arm and int(r["m"]) == m and int(r["bits"]) == bits:
                    acc[int(r["B"])].append(float(r["mean"]))
            xs = sorted(b for b in acc if sum(acc[b]) > 0)
            ys = [sum(acc[b]) / len(acc[b]) for b in xs]
            if xs:
                ax.plot(xs, ys, marker=MK[arm], color=COL[arm], lw=1.2, ms=5 if arm != "random" else 7,
                        label=arm, mfc="none" if arm in ("subgroup", "dickson-torus") else COL[arm])
        p = 2 ** bits
        Bs = [2 ** t for t in range(2, 11)]
        ax.plot(Bs, [b ** m / p for b in Bs], color="#222222", lw=1, ls="--", label="product law $B^m/p$")
        ax.plot(Bs, [b ** (m - 1) for b in Bs], color="#888888", lw=1.6, ls="-", label="$B^{m-1}$ (§2.3 bar)")
        ax.plot(Bs, [b ** (m - 1) / 16 for b in Bs], color="#888888", lw=1, ls=":", label="$B^{m-1}/16$ (live-lead bar)")
        ax.set_xscale("log", base=2); ax.set_yscale("log")
        ax.set_title(f"m = {m}, p ≈ 2^{bits}")
        if j == 0: ax.set_ylabel("mean solutions per target  N̄")
        if i == 1: ax.set_xlabel("factor-base size B")
handles, labels = axes[0][0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False, fontsize=8)
fig.suptitle("Solutions of {S_{m+1}=0, L(x_i)=0} per target (F̄_p = F_p count): every base sits on B^m/p, far below B^{m-1}\n"
             "census 2026-10-14, 3 curves per size, 2000 (m=2) / 500 (m=3) targets per cell; dloginterval omitted (see fig. 2)", fontsize=9)
fig.tight_layout(rect=[0, 0.07, 1, 0.94])
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(HERE / f"fig1_count_vs_B.{ext}", **kw)

# Figure 2: (a) N̄ / ((2B_E)^m/n) pooled per arm; (b) dispersion Var/N̄ per arm.
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
for k, m in enumerate([2, 3]):
    for idx, arm in enumerate(ARMS):
        vals = [float(r["r_exact"]) for r in rows if r["arm"] == arm and int(r["m"]) == m
                and r["r_exact"] not in ("NA", "0") and float(r["mean"]) > 0
                and float(r["mean"]) / float(r["r_exact"]) >= 0.5]
        disp = [float(r["disp"]) for r in rows if r["arm"] == arm and int(r["m"]) == m
                and r["scored"] == "1" and r["disp"] != "NA"]
        x = idx + (k - 0.5) * 0.3
        if vals:
            mu = sum(vals) / len(vals)
            a1.errorbar(x, mu, yerr=[[mu - min(vals)], [max(vals) - mu]], fmt=MK[arm], color=COL[arm],
                        mfc=COL[arm] if m == 3 else "white", capsize=3, ms=6)
        if disp:
            a2.scatter([x] * len(disp), disp, marker=MK[arm], color=COL[arm], s=14,
                       facecolors=COL[arm] if m == 3 else "none")
a1.axhline(1, color="#222222", lw=1, ls="--")
a1.set_xticks(range(len(ARMS))); a1.set_xticklabels(ARMS, rotation=30, ha="right")
a1.set_ylabel("N̄ / ((2·B_E)^m / n)   (mean, min–max over cells)")
a1.set_title("(a) count per target ÷ counting identity\nhollow = m=2, filled = m=3; cells with prediction ≥ 0.5")
a2.set_yscale("log"); a2.axhline(2, color="#888", lw=1, ls=":"); a2.axhline(6, color="#888", lw=1, ls=":")
a2.set_xticks(range(len(ARMS))); a2.set_xticklabels(ARMS, rotation=30, ha="right")
a2.set_ylabel("index of dispersion Var(N)/N̄")
a2.set_title("(b) spread across targets (scored cells)\nreference lines: 2 (m=2) and 6 (m=3) = orderings of one decomposition")
fig.tight_layout()
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(HERE / f"fig2_ratio_dispersion.{ext}", **kw)
print("ok")
