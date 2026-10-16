"""Lay out recorded values from <run>/summary.json (no computation beyond
reading). Writes cost.svg / cost.png and shape.svg / shape.png into <outdir>.

    python plot.py <run-dir> <outdir>
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

run, outdir = sys.argv[1], sys.argv[2]
s = json.load(open(f"{run}/summary.json"))
arms = s["arms"]
cells = sorted({(a["n"], a["m"]) for a in arms}, key=lambda c: (c[1], c[0]))
get = {(a["n"], a["m"], a["arm"]): a for a in arms}
colors = {"sem": "#2a6fdb", "rr": "#d9822b"}
labels = {"sem": "Semaev, symmetrised (sem)", "rr": "Nagao / Riemann–Roch (rr)"}

fig, ax = plt.subplots(figsize=(8.5, 4.6))
w = 0.38
for k, arm in enumerate(["sem", "rr"]):
    for i, c in enumerate(cells):
        a = get.get((c[0], c[1], arm))
        x = i + (k - 0.5) * w
        if a and a["total_median"]:
            lo, hi = a["total_range"]
            med = a["total_median"]
            ax.bar(x, med, w, color=colors[arm], label=labels[arm] if i == 0 else None)
            ax.errorbar(x, med, yerr=[[med - lo], [hi - med]], fmt="none", ecolor="black", capsize=3, lw=1)
        else:
            fails = (a or {}).get("failures") or {}
            why = ",".join(f"{k2}×{v}" for k2, v in fails.items()) or "no data"
            ax.bar(x, 1, w, color="none", edgecolor=colors[arm], hatch="//", label=labels[arm] if i == 0 else None)
            ax.text(x, 3, why, rotation=90, ha="center", va="bottom", fontsize=8, color=colors[arm])
ax.set_yscale("log")
ax.set_xticks(range(len(cells)))
ax.set_xticklabels([f"(n={n}, m={m})" for n, m in cells])
ax.set_ylabel("F_p multiplications per target\n(F4 + extraction), median, min–max")
ax.set_title("Decomposition cost per target, E/F_{p^n}, p = 65521, one engine (f4_fp + solve_zero_dim)")
ax.legend(loc="upper left", fontsize=9)
ax.grid(axis="y", which="major", alpha=0.3)
fig.text(0.01, 0.01, "Run ic_nagao_semaev_fpn_20261016/run-20261016; hatched = no completed solve (reason × count)", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1))
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(f"{outdir}/cost.{ext}", **kw)

fig, axs = plt.subplots(1, 2, figsize=(9, 3.8))
for arm in ["sem", "rr"]:
    xs, sd, cols, xs2 = [], [], [], []
    for c in cells:
        a = get.get((c[0], c[1], arm))
        if a and a["solving_degree"]:
            xs.append(f"({c[0]},{c[1]})")
            sd.append(a["solving_degree"][1])
            cols.append(a["max_cols"][1])
    axs[0].plot(xs, sd, "o-", color=colors[arm], label=arm)
    axs[1].plot(xs, cols, "o-", color=colors[arm], label=arm)
axs[0].set_ylabel("solving degree (max over targets)")
axs[1].set_ylabel("widest F4 matrix (columns, max)")
axs[1].set_yscale("log")
for a in axs:
    a.set_xlabel("cell (n, m)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
fig.suptitle("System shape: rr has lower degree but many more unknowns", fontsize=10)
fig.tight_layout()
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(f"{outdir}/shape.{ext}", **kw)
print("wrote", outdir)
