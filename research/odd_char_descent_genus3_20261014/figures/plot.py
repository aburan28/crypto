"""Lay out recorded values from the run's CSVs. Computes nothing new.

Usage: python plot.py <run-dir> <out-dir>
Reads min_genus.csv, genus_by_n.csv and counts.csv written by
examples/odd_char_ghs_genus.rs and draws genus_vs_n.svg (and a .png copy).
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"

run, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)


def rows(name):
    with open(run / name, newline="") as f:
        return list(csv.DictReader(f))


mins = rows("min_genus.csv")
byn = [r for r in rows("genus_by_n.csv") if r["period"] == r["n"]]
counts = rows("counts.csv")

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.2, 3.9), gridspec_kw={"width_ratios": [1.35, 1]})

# (a) genus against n
n_all = [int(r["n"]) for r in byn]
g_all = [int(r["genus"]) for r in byn]
ax.scatter(n_all, g_all, s=16, color="#b9b8b2", zorder=2, label="other genera, period-n types (n ≤ 8)")
nm = [int(r["n"]) for r in mins]
ax.plot(nm, [int(r["min_genus_period_n"]) for r in mins], "-o", color=BLUE, lw=2, ms=4.5,
        zorder=3, label="minimum genus, period-n types")
ax.plot(nm, [int(r["generic_genus"]) for r in mins], "--", color=ORANGE, lw=1.6, zorder=2,
        label="generic genus (n−1)·2ⁿ + 1")
ax.plot(nm, nm, ":", color=INK2, lw=1.2, zorder=1, label="g = n (dimension bound)")
g3 = [(int(r["n"]), int(r["genus"])) for r in byn if r["genus"] == "3"]
ax.scatter([a for a, _ in g3], [b for _, b in g3], s=70, facecolor="none", edgecolor=AQUA,
           lw=2, zorder=4, label="genus 3 (occurs only at n = 2, 3)")
mk = {int(r["n"]): r["min_genus_period_n"] for r in mins}
ax.text(0.98, 0.04, "minimum at n = 6, 10, 14: " + ", ".join(mk[n] for n in (6, 10, 14)),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=INK)
ax.set_yscale("log", base=2)
ax.set_xlabel("extension degree n  (E over F_{q^n}, q odd)")
ax.set_ylabel("genus of the GHS cover over F_q")
ax.set_xticks(range(1, 15))
ax.grid(True, color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.legend(fontsize=7, loc="upper left", frameon=False)
ax.set_title("(a) GHS genus, branch sets of exact period n", fontsize=9, color=INK, loc="left")

# (b) share of branch sets with genus 3 against q
for n, col in ((2, BLUE), (3, AQUA), (4, ORANGE)):
    pts = sorted((int(r["q"]), float(r["log10_share_genus3"])) for r in counts if r["n"] == str(n))
    bx.plot([q for q, _ in pts], [s for _, s in pts], "-o", color=col, lw=2, ms=4)
    q_last, s_last = pts[-1]
    bx.annotate(f"n = {n}", (q_last, s_last), xytext=(4, 0), textcoords="offset points",
                va="center", fontsize=8, color=INK)
bx.set_xscale("log")
bx.set_xlabel("q  (p ≤ 31, q = p and p²)")
bx.set_ylabel("log10(share of branch sets with genus 3)")
bx.grid(True, color=GRID, lw=0.6)
bx.set_axisbelow(True)
bx.set_xlim(4, 2500)
bx.set_title("(b) how rare genus 3 is among all x", fontsize=9, color=INK, loc="left")

fig.tight_layout()
fig.savefig(out / "genus_vs_n.svg")
fig.savefig(out / "genus_vs_n.png", dpi=150)
print("wrote", out / "genus_vs_n.svg")
