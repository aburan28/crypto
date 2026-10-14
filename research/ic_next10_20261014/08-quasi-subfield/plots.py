#!/usr/bin/env python3
"""Figures for item 8.  Usage: python3 -I plots.py <dir>"""
import csv
import glob
import json
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = sys.argv[1]
BETA_STAR = 1 / (2 * 4.876)


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def save(fig, name):
    fig.savefig(f"{D}/{name}.svg")
    fig.savefig(f"{D}/{name}.png", dpi=150)


# ---------------- Figure 1: beta frontier at prime n ----------------
a1 = rows(f"{D}/a1_stable_frontier.csv")
a3 = rows(f"{D}/a3_multiplicative.csv")
a2 = []
for p in sorted(glob.glob(f"{D}/a2_*.csv")):
    a2 += rows(p)
primes = sorted({int(r["n"]) for r in a1 if int(r["n"]) <= 61})

best1, best2, best3 = {}, {}, {}
for r in a1:
    n, n0, ell = int(r["n"]), int(r["n0"]), int(r["min_ell"])
    if n > 61 or n0 in (1, n - 1) or ell < 1:
        continue
    b = float(r["beta"])
    best1[n] = min(best1.get(n, 9e9), b)
for r in a2:
    if int(r["found"]) > 0 and int(r["j"]) >= 1:
        n = int(r["n"])
        best2[n] = min(best2.get(n, 9e9), float(r["beta_if_found"]))
for r in a3:
    n = int(r["n"])
    if r["subfield"] == "True" or n not in primes:
        continue
    best3[n] = min(best3.get(n, 9e9), float(r["beta"]))

fig, ax = plt.subplots(figsize=(8, 4.6))
ax.axhspan(0, BETA_STAR, color="#2a9d8f", alpha=0.15, label=f"HKPYY generic-beating region  β < {BETA_STAR:.4f}")
ax.axhline(1.0, color="#555", ls="--", lw=1, label="Definition 3.1 boundary  β = 1")
ax.axhline(0.75, color="#e76f51", ls=":", lw=1.2, label="Euler–Petit Thm 1 floor (linearised)  β = 3/4")
ax.plot(sorted(best1), [best1[n] for n in sorted(best1)], "o-", color="#264653", label="A1 Frobenius-stable linearised (all divisors of tⁿ−1)")
ax.plot(sorted(best2), [best2[n] for n in sorted(best2)], "s", mfc="none", ms=9, color="#e9c46a", mew=2, label="A2 general-coefficient linearised, j ≤ 1 (≤ 2 at n ≤ 13)")
ax.plot(sorted(best3), [best3[n] for n in sorted(best3)], "^", color="#9b2226", label="A3 multiplicative X^(2^n′) − X^a, root fraction ≥ 1/2, n′ < n")
ax.set_xlabel("prime extension degree n")
ax.set_ylabel("best β = n·log₂(deg λ) / n′²  (lower is better)")
ax.set_yscale("log")
ax.set_ylim(0.05, 40)
ax.set_title("Quasi-subfield polynomial frontier over GF(2^n), prime n ≤ 61")
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout()
save(fig, "fig_beta_frontier")

# ---------------- Figure 2: decomposition (if present) ----------------
dec = []
for p in sorted(glob.glob(f"{D}/decomp_*.jsonl")):
    with open(p) as fh:
        dec += [json.loads(l) for l in fh if l.strip()]
if dec:
    groups = {}
    for r in dec:
        groups.setdefault((r["n"], r["k"]), []).append(r)
    keys = sorted(groups)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for gi, key in enumerate(keys):
        for r in groups[key]:
            structured = not r["base"].startswith("random")
            y = r["yield"]
            nt = r["yield_targets"]
            se = math.sqrt(max(y * (1 - y), 1e-9) / nt)
            x = gi + (-0.15 if structured else 0.15)
            axes[0].errorbar(x, y, yerr=1.96 * se, fmt="o" if structured else "x",
                             color="#9b2226" if structured else "#264653", ms=7)
            degs = [t["refute_degree"] if t["refute_degree"] is not None else r["d_max"] + 1 for t in r["refute"]]
            if degs:
                axes[1].scatter([x] * len(degs), degs, marker="o" if structured else "x",
                                color="#9b2226" if structured else "#264653", alpha=0.6)
    for ax in axes:
        ax.set_xticks(range(len(keys)))
        ax.set_xticklabels([f"n={n}, k={k}" for n, k in keys], fontsize=8)
    axes[0].set_ylabel("relation yield (fraction of targets decomposed), 95% CI")
    axes[1].set_ylabel("refutation degree d (d_max+1 = not refuted by d_max)")
    axes[0].set_title("Yield: structured (●) vs random subspaces (×)")
    axes[1].set_title("Refutation degree on non-decomposable targets")
    fig.tight_layout()
    save(fig, "fig_decomposition")
print("ok")
