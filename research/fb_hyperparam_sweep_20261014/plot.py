"""Lay out recorded values from data/sweep.json as SVG/PNG figures.

Plots only: every number drawn is read from the JSON written by
examples/fb_hyperparam_sweep.rs (measured values and its `analysis` /
`derived_cost_model` sections). Nothing here computes a result.

    python plot.py [data/sweep.json]
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

HERE = Path(__file__).parent
SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "data" / "sweep.json"
OUT = HERE / "figures"
OUT.mkdir(exist_ok=True)
doc = json.loads(SRC.read_text())
rows = doc["binary"]
TAG = f"study fb_hyperparam_sweep_20261014, seed {doc['seed']:#x}"

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
FAM = {  # fixed categorical order: blue, orange, aqua
    "coordinate": ("#2a78d6", "o"),
    "random": ("#eb6834", "s"),
    "subfield": ("#1baf7a", "^"),
}
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF,
    "axes.facecolor": SURF, "lines.linewidth": 2, "legend.frameon": False,
})


def save(fig, name):
    fig.savefig(OUT / f"{name}.svg")
    fig.savefig(OUT / f"{name}.png", dpi=150)
    plt.close(fig)


# 1. Measured yield vs heuristic |V|^2/2^n.
fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
for ax, cls in zip(axes, ["generic_b", "koblitz_b1"]):
    for fam, (col, mk) in FAM.items():
        sel = [r for r in rows if r["curve_class"] == cls and r["family"] == fam and r["yield"]["decomposed"] > 0]
        if not sel:
            continue
        h = [r["yield"]["heuristic_V2_over_2n"] for r in sel]
        m = [r["yield"]["rate"] for r in sel]
        lo = [r["yield"]["rate"] - r["yield"]["wilson95"][0] for r in sel]
        hi = [r["yield"]["wilson95"][1] - r["yield"]["rate"] for r in sel]
        ax.errorbar(h, m, yerr=[lo, hi], fmt=mk, ms=5, color=col, ecolor=col, elinewidth=1,
                    mec=SURF, mew=0.8, label=f"{fam} (n={len(sel)} cells)")
    lim = [1e-4, 1.5]
    ax.plot(lim, lim, color=INK2, lw=1, ls="--", label="measured = heuristic")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("heuristic yield |V|²/2ⁿ (derived)")
    ax.set_title(cls, color=INK, fontsize=10)
axes[0].set_ylabel("measured yield (4,096 targets/cell, Wilson 95%)")
axes[0].legend(loc="upper left", fontsize=8)
fig.suptitle("Relation yield: measured vs heuristic, m = 2 over F_{2ⁿ}, n = 8–16", fontsize=10)
fig.text(0.01, 0.005, TAG + "; cells with zero decompositions omitted (log scale)", fontsize=7, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1))
save(fig, "yield_vs_heuristic")

# 2. d_ff and mean D* vs n' (small multiples by family; generic_b).
ns = sorted({r["n"] for r in rows})
pick = [n for n in ns if n in (10, 12, 14, 16)] or ns[-4:]
fig, axes = plt.subplots(1, len(pick), figsize=(max(7.5, 2.6 * len(pick)), 3.8), sharey=True)
axes = axes if len(pick) > 1 else [axes]
for ax, n in zip(axes, pick):
    for fam, (col, mk) in FAM.items():
        sel = sorted([r for r in rows if r["curve_class"] == "generic_b" and r["family"] == fam and r["n"] == n],
                     key=lambda r: r["n_sub"])
        x = [r["n_sub"] for r in sel if r["degrees"]["dstar_mean"] is not None]
        ds = [r["degrees"]["dstar_mean"] for r in sel if r["degrees"]["dstar_mean"] is not None]
        sd = [r["degrees"]["dstar_sd"] or 0 for r in sel if r["degrees"]["dstar_mean"] is not None]
        if x:
            ax.errorbar(x, ds, yerr=sd, marker=mk, color=col, ms=6, mec=SURF, capsize=2, label=f"D* {fam}")
        xf = [r["n_sub"] for r in sel if r["degrees"]["first_fall_mean"] is not None]
        ff = [r["degrees"]["first_fall_mean"] for r in sel if r["degrees"]["first_fall_mean"] is not None]
        if xf:
            ax.plot(xf, ff, marker=mk, color=col, ls=":", lw=1.2, ms=4, mfc=SURF, label=f"d_ff {fam}")
    ax.set_title(f"n = {n}", fontsize=10, color=INK)
    ax.set_xlabel("factor-base dimension n'")
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
axes[0].set_ylabel("Macaulay degree (mean ± sd)")
axes[-1].legend(fontsize=7, loc="upper left")
fig.suptitle("First-fall degree d_ff (dotted) vs refutation degree D* (solid), generic b", fontsize=10)
fig.text(0.01, 0.005, TAG + f"; 16 non-decomposable targets per cell; d_max = {doc['budget']['d_max']}, censored D* excluded from the mean",
         fontsize=7, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1))
save(fig, "dff_vs_dstar")

# 3. Derived cost model vs n' for the coordinate family, against the rho reference.
fig, axes = plt.subplots(1, len(pick), figsize=(max(7.5, 2.6 * len(pick)), 3.8), sharey=True)
axes = axes if len(pick) > 1 else [axes]
for ax, n in zip(axes, pick):
    for cls, ls in (("generic_b", "-"), ("koblitz_b1", "--")):
        sel = sorted([r for r in rows if r["curve_class"] == cls and r["family"] == "coordinate" and r["n"] == n],
                     key=lambda r: r["n_sub"])
        for key, col, lab in (("total_seconds_root_find", "#2a78d6", "root finding"),
                              ("total_seconds_macaulay", "#eb6834", "Macaulay to D*")):
            pts = [(r["n_sub"], r["derived_cost_model"][key]) for r in sel if r["derived_cost_model"][key]]
            if pts:
                ax.plot(*zip(*pts), ls=ls, color=col, marker="o", ms=4, mec=SURF, label=f"{lab}, {cls}")
    rho = [r["derived_cost_model"]["rho_reference_seconds"] for r in rows if r["n"] == n][:1]
    if rho:
        ax.axhline(rho[0], color=INK2, lw=1, ls=":", label="rho reference")
    ax.set_yscale("log"); ax.set_title(f"n = {n}", fontsize=10, color=INK)
    ax.set_xlabel("factor-base dimension n'")
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
axes[0].set_ylabel("derived seconds (decomp. + LA)")
axes[0].legend(fontsize=6.5, loc="lower left")
fig.suptitle("Derived cost model, coordinate factor base (measured inputs, PROTOCOL.md model)", fontsize=10)
fig.text(0.01, 0.005, TAG + "; Macaulay curve is a lower bound where D* was censored", fontsize=7, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1))
save(fig, "cost_model")

# 4. Prime negative control: whole-solve seconds vs FB size, IC m=2/3 vs rho.
prime = doc["prime"]
bits = sorted({r["bits"] for r in prime})
fig, axes = plt.subplots(1, len(bits), figsize=(max(7.5, 3.2 * len(bits)), 3.8), sharey=True)
axes = axes if len(bits) > 1 else [axes]
for ax, b in zip(axes, bits):
    for meth, col, mk in (("ic_m2", "#2a78d6", "o"), ("ic_m3", "#eb6834", "s")):
        ok = [r for r in prime if r["bits"] == b and r["method"] == meth and r["solved"]]
        bad = [r for r in prime if r["bits"] == b and r["method"] == meth and not r["solved"]]
        if ok:
            ax.scatter([r["fb_size"] for r in ok], [r["seconds"] for r in ok], color=col, marker=mk, s=28,
                       edgecolors=SURF, linewidths=0.8, label=f"IC {meth[-2:]} solved", zorder=3)
        if bad:
            ax.scatter([r["fb_size"] for r in bad], [r["seconds"] for r in bad], facecolors="none",
                       edgecolors=col, marker=mk, s=28, label=f"IC {meth[-2:]} no solve", zorder=3)
    rho = [r["seconds"] for r in prime if r["bits"] == b and r["method"] == "rho_negation_1thread" and r["solved"]]
    if rho:
        ax.axhspan(min(rho), max(rho), color="#1baf7a", alpha=0.25, lw=0)
        ax.axhline(min(rho), color="#1baf7a", lw=1.5, label=f"rho, {len(rho)} solves (min–max band)")
        ax.axhline(max(rho), color="#1baf7a", lw=1.5)
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_title(f"{b}-bit a = −3 curve", fontsize=10, color=INK)
    ax.set_xlabel("factor-base size (points)")
axes[0].set_ylabel("whole solve, wall s (1 thread)")
axes[0].legend(fontsize=7, loc="upper left")
fig.suptitle("Prime-field negative control: small-x factor base vs rho, same curve and secret", fontsize=10)
fig.text(0.01, 0.005, TAG, fontsize=7, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1))
save(fig, "prime_control")
print("figures written to", OUT)
