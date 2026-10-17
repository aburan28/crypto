"""Figure for item 3 (visualisation only; every number is read from the CSVs
next to this file, which the Rust tool and msolve produced).

    python3 -I plot.py   # from this directory
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(sys.argv[0]))


def rows(name):
    with open(os.path.join(here, name)) as f:
        return list(csv.DictReader(f))


f4 = rows("f4_committed.csv")
ms = rows("results.csv")
xl = rows("xl.csv")

colors = {"2": "#4C78A8", "3": "#E45756", "4": "#54A24B"}
fig, ax = plt.subplots(figsize=(8.6, 5.4))
for m in ("2", "3", "4"):
    pts = sorted({(int(r["N"]), int(r["D_f4"])) for r in f4 if r["m"] == m})
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "-o", color=colors[m], lw=1.6,
            ms=6, mfc="white", label=f"F4 (f4_fp_tower, committed), m = {m}")
    done = sorted({(int(r["N"]), int(r["D_msolve"])) for r in ms
                   if r["m"] == m and r["verdict"] == "no_solution"})
    if done:
        ax.plot([p[0] for p in done], [p[1] for p in done], "x", color=colors[m], ms=10,
                mew=2.2, label=f"msolve 0.6.5 (finished), m = {m}")
    lb = sorted({(int(r["N"]), int(r["D_last_productive"])) for r in ms
                 if r["m"] == m and r["verdict"] != "no_solution" and r["D_last_productive"] not in ("", "0")})
    if lb:
        ax.plot([p[0] for p in lb], [p[1] for p in lb], "^", color=colors[m], ms=9, mfc="none",
                mew=1.6, label=f"msolve stopped at the wall: D ≥ (last productive round), m = {m}")
    xs = sorted({(int(r["N"]), int(r["D_XL"]), r["exact"]) for r in xl if r["m"] == m})
    if xs:
        ax.plot([p[0] for p in xs], [p[1] for p in xs], "--s", color=colors[m], alpha=0.6, ms=5,
                label=f"Macaulay/XL degree d + δ, m = {m}")
nn = list(range(8, 28))
ax.plot(nn, [n + 1 for n in nn], ":", color="#888", label="top-part d_reg = N + 1 (m = 2, full)")
ax.axhline(7, color="#bbb", lw=0.8)
ax.set_xlabel("N = m·t (tower unknowns)")
ax.set_ylabel("solving degree")
ax.set_title("PKM Kummer tower at p₁ = 2013265921: F4 degree vs XL degree")
ax.set_xlim(7, 28)
ax.set_ylim(3, 28)
ax.grid(alpha=0.25)
ax.legend(fontsize=7.2, loc="upper left", frameon=False)
fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(os.path.join(here, f"degree_vs_N.{ext}"), dpi=160)
print("wrote degree_vs_N.svg/.png")
