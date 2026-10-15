"""Draw the 04b figures from analysis.json and tuning/kappa_choice.json.

Layout only: every plotted number is read from the Rust analysis output
(examples/ic_improve_04b_analyze.rs); nothing is computed here beyond unit
conversion for axes. Writes SVG and a 150 dpi PNG of each figure.

    python plots.py <analysis.json> <kappa_choice.json> <outdir>
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ARMS = {
    "BL": ("BL: κ = 2, tree-cycle elimination + Lanczos", "#2a78d6", "o", "-"),
    "A1": ("A1: κ retuned per size, as BL", "#eb6834", "s", "--"),
    "A2": ("A2: κ = 2, SGE + Lanczos", "#1baf7a", "^", "-."),
    "A3": ("A3 (secondary): κ retuned, SGE", "#eda100", "D", ":"),
}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})


def save(fig, out, name):
    fig.tight_layout()
    fig.savefig(f"{out}/{name}.svg")
    fig.savefig(f"{out}/{name}.png", dpi=150)
    plt.close(fig)


def series(an, arm, phase, what):
    d = an["point"][f"{arm}/{phase}"][what]
    ks = sorted(int(k) for k in d)
    return ks, [d[str(k)] for k in ks]


def ci(an, key):
    c = an["ci"].get(key)
    return (c["lo95"], c["hi95"]) if c else (None, None)


def main(apath, cpath, out):
    an = json.load(open(apath))
    ch = json.load(open(cpath))
    run = an.get("run_id", "")

    # 1. local exponents of total cost, with 4/3
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for i, (arm, (lab, col, mk, ls)) in enumerate(ARMS.items()):
        ks, es = series(an, arm, "total", "local")
        if not ks:
            continue
        x = [k + 0.5 + (i - 1.5) * 0.06 for k in ks]
        lo = [e - ci(an, f"{arm}/total/local/{k}")[0] for k, e in zip(ks, es)]
        hi = [ci(an, f"{arm}/total/local/{k}")[1] - e for k, e in zip(ks, es)]
        ax.errorbar(x, es, yerr=[lo, hi], color=col, marker=mk, ms=5, lw=1.6, ls=ls, capsize=2, label=lab)
    ax.axhline(4 / 3, color=INK, lw=1, ls=(0, (4, 3)))
    ax.text(0.99, 4 / 3 + 0.004, "GTTD 4/3", color=INK, fontsize=8, va="bottom", ha="right", transform=ax.get_yaxis_transform())
    ax.axhline(1.5, color=MUTED, lw=0.8, ls=(0, (1, 2)))
    ax.text(0.99, 1.5 + 0.004, "rho 3/2", color=MUTED, fontsize=8, va="bottom", ha="right", transform=ax.get_yaxis_transform())
    ax.set_xlabel("log₂ q (midpoint of the size step k → k+1)")
    ax.set_ylabel("local exponent of mean total cost")
    ax.set_title(f"Local cost exponents, 12 test instances per size, 95% bootstrap intervals {run}", fontsize=9)
    ax.legend(fontsize=7.5, frameon=False, loc="lower left")
    save(fig, out, "local_exponents")

    # 2. phases of BL and A2: local exponents of relation and linear algebra
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)
    for ax, arm in zip(axs, ["BL", "A2"]):
        lab, col, mk, _ = ARMS[arm]
        for phase, ls, m2, name in [("rel", "-", "o", "relation phase"), ("la", "--", "s", "linear algebra"), ("total", ":", "x", "total")]:
            ks, es = series(an, arm, phase, "local")
            lo = [e - ci(an, f"{arm}/{phase}/local/{k}")[0] for k, e in zip(ks, es)]
            hi = [ci(an, f"{arm}/{phase}/local/{k}")[1] - e for k, e in zip(ks, es)]
            ax.errorbar([k + 0.5 for k in ks], es, yerr=[lo, hi], color=col if phase != "total" else INK, marker=m2, ms=4, lw=1.4, ls=ls, capsize=2, label=name)
        ax.axhline(4 / 3, color=INK, lw=0.8, ls=(0, (4, 3)))
        ax.set_title(lab, fontsize=8.5)
        ax.set_xlabel("log₂ q (step midpoint)")
        ax.legend(fontsize=7.5, frameon=False)
    axs[0].set_ylabel("local exponent")
    save(fig, out, "phase_exponents")

    # 3. matrix weight entering Lanczos
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    for arm in ["BL", "A1", "A2", "A3"]:
        lab, col, mk, ls = ARMS[arm]
        ks = sorted(int(k.split("/")[1]) for k in an["weights"] if k.startswith(arm + "/"))
        if ks:
            ax.plot(ks, [an["weights"][f"{arm}/{k}"]["w_row"] for k in ks], color=col, marker=mk, ms=5, lw=1.6, ls=ls, label=lab)
    ax.set_xlabel("log₂ q (k)")
    ax.set_ylabel("nonzeros per relation (row) entering Lanczos")
    ax.set_title("Row weight of the Lanczos matrix, mean over 12 test instances", fontsize=9)
    ax.legend(fontsize=7.5, frameon=False)
    save(fig, out, "row_weight")

    # 4. paired ratios vs BL
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    for i, arm in enumerate(["A1", "A2", "A3"]):
        lab, col, mk, ls = ARMS[arm]
        ks = sorted(int(k.split("/")[2]) for k in an["ratios"] if k.startswith(f"{arm}/total/"))
        r = [an["ratios"][f"{arm}/total/{k}"]["geo_mean"] for k in ks]
        lo = [v - an["ratios"][f"{arm}/total/{k}"]["ci"]["lo95"] for k, v in zip(ks, r)]
        hi = [an["ratios"][f"{arm}/total/{k}"]["ci"]["hi95"] - v for k, v in zip(ks, r)]
        ax.errorbar([k + (i - 1) * 0.08 for k in ks], r, yerr=[lo, hi], color=col, marker=mk, ms=5, lw=1.6, ls=ls, capsize=2, label=lab)
    ax.axhline(1, color=INK, lw=0.8)
    ax.set_xlabel("log₂ q (k)")
    ax.set_ylabel("total cost / BL total cost (geometric mean)")
    ax.set_title("Paired cost ratios against the rebuilt baseline, 95% bootstrap intervals", fontsize=9)
    ax.legend(fontsize=7.5, frameon=False)
    save(fig, out, "paired_ratios")

    # 5. mean total cost against q, log-log, with a q^{4/3} guide through BL at k = 14
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    q = {int(k): v for k, v in an["q"].items()}
    for arm, (lab, col, mk, ls) in ARMS.items():
        ks, ms = series(an, arm, "total", "mean")
        if ks:
            ax.plot([q[k] for k in ks], ms, color=col, marker=mk, ms=5, lw=1.6, ls=ls, label=lab)
    ks, ms = series(an, "BL", "total", "mean")
    m14 = dict(zip(ks, ms)).get(14)
    if m14:
        xs = sorted(q.values())
        ax.plot(xs, [m14 * (x / q[14]) ** (4 / 3) for x in xs], color=INK, lw=0.9, ls=(0, (4, 3)), label="q^{4/3} through BL at 2^14 (guide)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("q")
    ax.set_ylabel("mean total multiplications")
    ax.set_title("Mean total cost per size (12 test instances)", fontsize=9)
    ax.legend(fontsize=7.5, frameon=False)
    save(fig, out, "total_cost")

    # 6. tuning curves
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    cmap = plt.get_cmap("Blues")
    ks = sorted(int(k) for k in ch["choice"])
    for j, k in enumerate(ks):
        tab = [t for t in ch["choice"][str(k)]["table"] if t["n"] > 0]
        if not any(abs(t["kappa"] - 2) < 1e-9 for t in tab):
            continue
        base =[t for t in tab if abs(t["kappa"] - 2) < 1e-9][0]["mean_total_bl_pipeline"]
        ax.plot([t["kappa"] for t in tab], [t["mean_total_bl_pipeline"] / base for t in tab], color=cmap(0.35 + 0.6 * j / max(1, len(ks) - 1)), marker="o", ms=4, lw=1.4, label=f"k = {k}, κ* = {ch['choice'][str(k)]['kappa_star']}")
    ax.axhline(1, color=MUTED, lw=0.8)
    ax.set_xscale("log", base=2)
    ax.set_xticks([0.5, 0.7, 1, 1.4, 2])
    ax.set_xticklabels(["0.5", "0.7", "1", "1.4", "2"])
    ax.set_xlabel("κ (factor base F = ceil(κ q^{2/3}))")
    ax.set_ylabel("mean total cost / cost at κ = 2 (BL pipeline)")
    ax.set_title("Tuning on training curves 5–7, targets 1–2 (6 runs per point)", fontsize=9)
    ax.legend(fontsize=7, frameon=False, ncol=2)
    save(fig, out, "tuning")


if __name__ == "__main__":
    main(*sys.argv[1:4])
