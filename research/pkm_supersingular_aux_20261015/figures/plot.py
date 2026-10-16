"""Plot the recorded values of research/pkm_supersingular_aux_20261015.

Reads only the committed TSVs (q1/pooled.tsv, q2/cells.tsv,
q3/tables/summary-arms.tsv, q3/tables/summary-ratios.tsv) and draws them; no
value is computed here beyond grouping rows for the axes.

    python3 figures/plot.py   # from the study directory
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
COL = {
    "ss-j1728": "#1f77b4", "ss-j1728-twist": "#17becf", "ss-j0": "#9467bd",
    "ss-j0-twist": "#c5b0d5", "ordinary": "#d62728", "torus": "#2ca02c",
    "split": "#98df8a", "null": "#7f7f7f",
}


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(fig, name):
    fig.savefig(HERE / f"{name}.svg")
    fig.savefig(HERE / f"{name}.png", dpi=150)
    plt.close(fig)


def q2():
    cells = rows(STUDY / "q2/cells.tsv")
    primes = [("786433", "p₀ = 3·2^18 + 1 (positive control)"),
              ("14680063", "T24 = 7·2^21 − 1"),
              ("2046820351", "T31 = 61·2^25 − 1"),
              ("18874367", "S25 = 9·2^21 − 1")]
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.2), sharey=True)
    for ax, (p, title) in zip(axes, primes):
        series = defaultdict(dict)
        for r in cells:
            if r["p"] != p or r["m"] != "2" or r["timed_out"] == "true" or r["oversize"] == "true":
                continue
            key = "null" if r["control"] == "null" else r["arm"]
            series[key].setdefault(int(r["N"]), set()).add(int(r["D"]))
        offs = {k: (i - len(series) / 2) * 0.12 for i, k in enumerate(sorted(series))}
        for k in sorted(series):
            ns = sorted(series[k])
            ds = [max(series[k][n]) for n in ns]
            mixed = [n for n in ns if len(series[k][n]) > 1]
            ax.plot([n + offs[k] for n in ns], [d + offs[k] * 0.3 for d in ds], marker="o",
                    ms=4, lw=1.4, color=COL.get(k, "k"), label=k)
            for n in mixed:
                ax.annotate("mixed", (n, max(series[k][n])), fontsize=7)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("N = m·t (m = 2)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="lower right")
    axes[0].set_ylabel("F4 solving degree D (both targets)")
    fig.suptitle("Q2: F4 solving degree of PKM m = 2 systems by tower family "
                 "(f4_fp_tower, 2 random targets per cell; small offsets separate equal lines)",
                 fontsize=10)
    fig.tight_layout()
    save(fig, "q2_d_vs_n")


def q1():
    pooled = rows(STUDY / "q1/pooled.tsv")
    fams = ["ss-j1728", "ss-j1728-twist", "ss-j0", "ss-j0-twist", "ordinary", "ss-j1728-iso16"]
    ts = sorted({int(r["t"]) for r in pooled})
    fig, ax = plt.subplots(figsize=(9, 4.2))
    width = 0.13
    for i, f in enumerate(fams):
        target = "ss-j1728" if f == "ss-j1728-iso16" else "torus"
        xs, ys, labels = [], [], []
        for t in ts:
            hit = [r for r in pooled if r["t"] == str(t) and r["family_a"] == f and r["family_b"] == target]
            if hit:
                xs.append(t + (i - 2.5) * width)
                ys.append(float(hit[0]["fraction"]))
                labels.append(f'{hit[0]["equivalent"]}/{hit[0]["sets_a"]}')
        lab = f"{f} → {target}" + (" (positive control)" if target != "torus" else "")
        bars = ax.bar(xs, ys, width, color=COL.get(f, "#ff7f0e" if "iso" in f else "k"), label=lab)
        for b, s in zip(bars, labels):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.01, s, rotation=90,
                    ha="center", va="bottom", fontsize=6)
    ax.set_xticks(ts)
    ax.set_xticklabels([f"t = {t}" for t in ts])
    ax.set_ylim(0, 1.6)
    ax.set_ylabel("fraction of sets Möbius-equivalent\nto some set of the target family")
    ax.set_title("Q1: tower sets up to x ↦ (αx+β)/(γx+δ) over F_p, pooled over primes "
                 "p ≡ 3 mod 4, 100 ≤ p ≤ 3000 (exact)", fontsize=9)
    ax.legend(fontsize=7, ncol=3, loc="upper center")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    save(fig, "q1_equivalence")


def q3():
    arms = rows(STUDY / "q3/tables/summary-arms.tsv")
    ratios = rows(STUDY / "q3/tables/summary-ratios.tsv")
    curves = ["s20", "s22", "s24"]
    order = ["ic-ss-rootfind", "ic-plain-ss", "ic-ord-rootfind", "ic-plain-ord", "ic-ss-mitm"]
    colors = {"ic-ss-rootfind": "#1f77b4", "ic-plain-ss": "#aec7e8", "ic-ord-rootfind": "#d62728",
              "ic-plain-ord": "#ff9896", "ic-ss-mitm": "#9467bd"}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.4))
    w = 0.15
    for i, arm in enumerate(order):
        xs, ys, sh = [], [], []
        for j, c in enumerate(curves):
            r = [x for x in arms if x["curve"] == c and x["arm"] == arm][0]
            xs.append(j + (i - 2) * w)
            ys.append(float(r["mean_S"]))
            sh.append(float(r["mean_aux_share"]) * float(r["mean_S"]))
        a1.bar(xs, ys, w, color=colors[arm], label=arm)
        a1.bar(xs, sh, w, color="none", edgecolor="k", hatch="///", lw=0.5)
    a1.set_xticks(range(3))
    a1.set_xticklabels(curves)
    a1.set_ylabel("mean S = GAE / √r (ecbench), 16 runs")
    a1.set_title("Q3: mean S by arm; hatched = auxiliary search + tower build share", fontsize=9)
    a1.legend(fontsize=7)
    a1.grid(axis="y", alpha=0.3)
    pairs = [("ic-ss-rootfind", "ic-plain-ss", "ss ÷ plain(ss size)"),
             ("ic-ord-rootfind", "ic-plain-ord", "ordinary ÷ plain(ord size)"),
             ("ic-ss-rootfind", "ic-ord-rootfind", "ss ÷ ordinary"),
             ("ic-plain-ss", "ic-plain-ord", "plain(ss size) ÷ plain(ord size)"),
             ("ic-ss-rootfind", "ic-ss-mitm", "rootfind ÷ mitm (positive control)")]
    for i, (b, a, lab) in enumerate(pairs):
        for j, c in enumerate(curves):
            r = [x for x in ratios if x["curve"] == c and x["B"] == b and x["A"] == a][0]
            y, lo, hi = float(r["ratio"]), float(r["lo"]), float(r["hi"])
            x = j + (i - 2) * 0.13
            a2.errorbar([x], [y], yerr=[[y - lo], [hi - y]], fmt="o", ms=4, capsize=2,
                        color=f"C{i}", label=lab if j == 0 else None)
    a2.axhline(1, color="k", lw=0.8)
    a2.axhspan(0.97, 1.03, color="#dddddd", alpha=0.6, lw=0)
    a2.set_xticks(range(3))
    a2.set_xticklabels(curves)
    a2.set_ylabel("paired ratio Σ S_B / Σ S_A, 95 % cluster bootstrap")
    a2.set_title("Q3: paired ratios (16 pairs on 8 workloads per curve; band = H3b's [0.97, 1.03])",
                 fontsize=9)
    a2.legend(fontsize=7)
    a2.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    save(fig, "q3_s_and_ratios")


if __name__ == "__main__":
    q3()
    if (STUDY / "q2/cells.tsv").exists():
        q2()
    if (STUDY / "q1/pooled.tsv").exists():
        q1()
