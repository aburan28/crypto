"""Lay out the recorded values written by `f4_generic_report` (AGENTS.md:
Python only plots; every number comes from the Rust tool's CSVs).

Writes, beside this file:
- fig_claims.svg/.png  — the five pre-registered speed ratios against their thresholds;
- fig_cost.svg/.png    — ns per multiply-add per backend (speed arms);
- fig_d_vs_n.svg/.png  — solving degree D against N per prime, m = 2 and m = 3, with the
                          recorded small-prime reference;
- fig_time_vs_n.svg/.png — F4 time against N per prime (m = 2).
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

here = Path(__file__).resolve().parent
speed = here / "runs" / "speed_report"
cells = here / "runs" / "cells_report"

INK, MUTED, GRID = "#1f2328", "#57606a", "#d0d7de"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "svg.fonttype": "none",
})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(here / f"{name}.svg")
    fig.savefig(here / f"{name}.png", dpi=150)
    plt.close(fig)


# ── claims ──
claims = list(csv.DictReader(open(speed / "claims.csv")))
fig, ax = plt.subplots(figsize=(6.4, 2.8))
ys = range(len(claims))
for y, c in zip(ys, claims):
    lo, hi, med, thr = (float(c[k]) for k in ("ratio_min", "ratio_max", "ratio_median_paired", "threshold"))
    ok = c["verdict"] == "PASS"
    ax.plot([lo, hi], [y, y], color=MUTED, lw=1.2)
    ax.plot(med, y, "o", color="#1a7f37" if ok else "#cf222e", ms=6)
    ax.plot(thr, y, "|", color=INK, ms=14, mew=1.6)
    ax.text(hi + 0.05, y, f"{med:.2f}× ({c['verdict']})", va="center", fontsize=8, color=INK)
ax.set_yticks(list(ys), [f"{c['id']} ({c['cell']}): {c['baseline']} → {c['candidate']}" for c in claims])
ax.invert_yaxis()
ax.set_xlabel("paired time ratio baseline / candidate\n(dot: median of 7, bar: min–max, tick: threshold)")
ax.set_xlim(0.8, 3.2)
save(fig, "fig_claims")

# ── cost per multiply-add ──
summ = list(csv.DictReader(open(speed / "summary.csv")))
label = {
    ("C1", "u32"): "p₁ 31-bit, f4_fp_tower", ("C1", "fp32"): "p₁ 31-bit, Fp32",
    ("C3", "u32"): "p₂ 32-bit, f4_fp_tower", ("C3", "fp32"): "p₂ 32-bit, Fp32",
    ("C5", "fp64"): "v40 45-bit, Fp64", ("C5", "fp52"): "v40 45-bit, Fp52 (IFMA)",
    ("C6", "fp64"): "v56 61-bit, Fp64", ("C8", "p192"): "P-192, FpP192",
    ("C7", "fp224"): "P-224, Fp224",
}
pts = [(label[(r["cell"], r["arm"])], float(r["ns_per_muladd"])) for r in summ if (r["cell"], r["arm"]) in label]
pts.sort(key=lambda t: t[1])
fig, ax = plt.subplots(figsize=(6.4, 3.0))
ax.barh([p[0] for p in pts], [p[1] for p in pts], color="#0969da", height=0.6)
for i, (_, v) in enumerate(pts):
    ax.text(v * 1.05, i, f"{v:.2f}", va="center", fontsize=8, color=INK)
ax.set_xscale("log")
ax.set_xlabel("ns per F_p multiply-add (median of 7 runs, F4 time / multiply-adds, one core)")
save(fig, "fig_cost")

# ── D and time against N ──
rows = list(csv.DictReader(open(cells / "plot_data.csv")))
series = defaultdict(list)
for r in rows:
    if r["control"] != "tower":
        continue
    series[(r["prime"], r["kind"], int(r["m"]))].append(r)
STYLE = {
    "29686813949953": ("Kummer v₂=40 (45-bit)", "#2e8b57", "o"),
    "4222124650659841": ("Kummer v₂=48 (52-bit)", "#0969da", "s"),
    "1945555039024054273": ("Kummer v₂=56 (61-bit)", "#8250df", "^"),
    "P-224": ("Kummer P-224 (v₂=96)", "#cf222e", "D"),
    "P-192": ("torus P-192 (v₂(p+1)=64)", "#bf8700", "v"),
    "P-256": ("torus P-256 (v₂(p+1)=96)", "#57606a", "x"),
}
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
for ax, m in zip(axes, (2, 3)):
    ref_done = False
    for (prime, kind, mm), rs in sorted(series.items()):
        if mm != m or prime not in STYLE:
            continue
        name, col, mk = STYLE[prime]
        rs.sort(key=lambda r: int(r["N"]))
        off = list(STYLE).index(prime) * 0.12 - 0.3
        ax.plot([int(r["N"]) + off for r in rs], [int(r["D"]) for r in rs], mk, color=col, label=name, ms=4)
        if not ref_done:
            refs = [(int(r["N"]), max(int(x) for x in r["ref_D"].split("/"))) for r in rs if r["ref_D"]]
            if refs:
                ax.plot([n for n, _ in refs], [d for _, d in refs], "-", color=INK, lw=0.8, label="recorded, 31–32-bit primes")
                ref_done = True
    ax.set_title(f"m = {m}", fontsize=9)
    ax.set_xlabel("N = m·t tower unknowns")
axes[0].set_ylabel("solving degree D")
axes[0].legend(fontsize=7, loc="upper left")
save(fig, "fig_d_vs_n")

fig, ax = plt.subplots(figsize=(6.4, 3.2))
for (prime, kind, m), rs in sorted(series.items()):
    if m != 2 or prime not in STYLE:
        continue
    name, col, mk = STYLE[prime]
    rs.sort(key=lambda r: int(r["N"]))
    ax.plot([int(r["N"]) for r in rs], [float(r["ms_median"]) / 1e3 for r in rs], "-" + mk, color=col, label=name, ms=4, lw=1)
ax.set_yscale("log")
ax.set_xlabel("N (m = 2)")
ax.set_ylabel("F4 time per system, s (median of 2, 2 threads)")
ax.legend(fontsize=7)
save(fig, "fig_time_vs_n")
