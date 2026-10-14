#!/usr/bin/env python3
"""Figures for item 6: degree vs ell, XL (committed) against F4 (msolve).

python3 -I plot.py   (reads runs/*.jsonl, writes fig_degree_vs_ell.{svg,png}
and fig_g2_full_pass.{svg,png} when runs/g2_full_pass.jsonl exists)
"""
import json, os, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "runs")


def load(name):
    p = os.path.join(R, name)
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


f4 = [r for r in load("c2_calibration.jsonl") + load("m_main.jsonl") + load("m_main_b.jsonl")
      if "export" in r and r["export"].get("solutions", 1) <= 0 and r.get("gb_is_one")]
cells = {}
for r in f4:
    e = r["export"]
    if e.get("sat_ok"):
        continue
    cells.setdefault((e["m"], e["n"], e["ell"]), []).append(r["d_f4"])

# Committed XL refutation degrees (median reading; '>=' marks a lower bound).
# research/dreg_m4_chain_20261006/RESULTS.md, research/dreg_ell_grid_20260925,
# research/dreg_fixed_surplus_20260923, research/dreg_degree7_ell5_20260930.
xl = {4: [(1, 5, False), (2, 6, False), (3, 7, True)],
      3: [(3, 6, False), (4, 6, False), (5, 7, False), (5, 7, True)]}

C = {"xl": "#b4572e", "f4": "#2e6fb4", "rand": "#7a7a7a"}
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
for ax, m in zip(axes, (4, 3)):
    pts = [(ell, d, lb) for ell, d, lb in xl[m]]
    ax.plot([p[0] for p in pts if not p[2]], [p[1] for p in pts if not p[2]], "o",
            color=C["xl"], label="XL / Macaulay refutation degree (committed)")
    lbs = [p for p in pts if p[2]]
    if lbs:
        ax.plot([p[0] for p in lbs], [p[1] for p in lbs], "^", color=C["xl"], ms=9,
                label="XL lower bound (degree 6 built, not refuted)")
    for (mm, n, ell), ds in sorted(cells.items()):
        if mm != m:
            continue
        jitter = [ell + 0.06 * (i - 1.5) for i in range(len(ds))]
        ax.plot(jitter, ds, "s", color=C["f4"], ms=5, alpha=0.8)
        ax.annotate(f"n={n}", (ell, max(ds) + 0.15), ha="center", fontsize=7, color=C["f4"])
    ax.plot([], [], "s", color=C["f4"], label="F4 solving degree (msolve), each draw")
    ax.set_title(f"m = {m} chained system, unsatisfiable draws")
    ax.set_xlabel("ℓ (factor-base dimension)")
    ax.grid(alpha=0.3)
    ax.set_xticks(range(1, 7))
axes[0].set_ylabel("degree")
d = load("d_random.jsonl")
if d:
    r62 = [r["d_f4"] for r in d if r["label"].startswith("r62")]
    axes[0].axhline(min(r62), color=C["rand"], ls="--", lw=1)
    axes[0].text(1.0, min(r62) + 0.1, f"random systems, (6,2) shape: F4 reaches ≥{min(r62)}",
                 fontsize=7, color=C["rand"])
axes[0].legend(fontsize=7, loc="lower right")
axes[0].set_ylim(2.5, 8.5)
fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(os.path.join(HERE, f"fig_degree_vs_ell.{ext}"), dpi=150)

g = load("g2_full_pass.jsonl")
if g:
    g = sorted(g, key=lambda r: r["n"])
    fig, ax = plt.subplots(figsize=(6, 4))
    ns = [r["n"] for r in g]
    ax.semilogy(ns, [r["blocks_rref_cpu"] + r["inverse_cpu"] for r in g], "o-", color=C["f4"],
                label="equivariant: block RREF + inverse transform (CPU s)")
    fc = [(r["n"], r["full_rref_cpu"]) for r in g if r.get("full_rref_cpu") is not None]
    if fc:
        ax.semilogy([x for x, _ in fc], [y for _, y in fc], "s-", color=C["xl"],
                    label="full F₂ RREF of the symmetric matrix (CPU s)")
    ax.set_xlabel("n (barrel system, l = 3, degree 3)")
    ax.set_ylabel("CPU seconds")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(HERE, f"fig_g2_full_pass.{ext}"), dpi=150)
print({k: v for k, v in sorted(cells.items())})
