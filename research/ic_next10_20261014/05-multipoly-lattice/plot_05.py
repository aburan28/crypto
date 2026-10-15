"""Plots for item 05 from the recorded JSON (no new optimisation).

Usage: python3 -I plot_05.py <run_dir>
Writes reach_vs_family.svg/.png, finite_ell.svg/.png, lll.svg/.png into <run_dir>, and
extrapolation.json (least-squares fit δ_ℓ ≈ δ∞ − c/ℓ + e/ℓ² over the last five ℓ; derived, not exact).
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
run = sys.argv[1]
INK = "#1f2933"
C = {"single": "#3b6ea8", "pair": "#d07a2d", "triple": "#5a9e5a", "null": "#8a8a8a", "need": "#b03a48", "ceil": "#6b4fa0"}
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name):
    fig.savefig(os.path.join(run, name + ".svg"), bbox_inches="tight")
    fig.savefig(os.path.join(run, name + ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


lp = [r for r in json.load(open(os.path.join(run, "lp.json")))["records"] if r.get("part") == "lp"]
asym = json.load(open(os.path.join(HERE, "asymptotic_box.json")))

# ---------- 1. reach by family vs needed ----------
fig, ax = plt.subplots(figsize=(8.2, 4.6))
ms = [2, 3, 4]
labels = ["S₃ (m=2)", "S₄ (m=3)", "S₅ (m=4)"]
x = np.arange(3)
single = [next(a["delta_box"] for a in asym if a.get("system") == f"single S{m+1}") for m in ms]
bestmulti = []
bestlab = []
for m in ms:
    cand = [a for a in asym if a.get("m") == m and a.get("s") != m and "delta_box" in a]
    b = max(cand, key=lambda a: a["delta_box"])
    bestmulti.append(b["delta_box"])
    bestlab.append(b["system"])
ceil = [1 / (m * 2 ** (m - 1)) for m in ms]
need = [0.5, 0.25, 1 / 6]
w = 0.2
ax.bar(x - 1.5 * w, single, w, color=C["single"], label="single polynomial δ* (W10, reproduced)")
ax.bar(x - 0.5 * w, bestmulti, w, color=C["pair"], label="best shared-target δ* (box family, ℓ→∞)")
ax.bar(x + 0.5 * w, ceil, w, color=C["ceil"], label="ceiling 1/(md): every M, ℓ, sublattice")
ax.bar(x + 1.5 * w, need, w, color=C["need"], label="δ needed to beat rho (single target)")
for i in range(3):
    ax.text(x[i] - 0.5 * w, bestmulti[i] * 1.08, f"{bestmulti[i]:.4f}", ha="center", fontsize=7, rotation=90, va="bottom")
    ax.text(x[i] - 1.5 * w, single[i] * 1.08, f"{single[i]:.4f}", ha="center", fontsize=7, rotation=90, va="bottom")
ax.set_yscale("log")
ax.set_xticks(x, labels)
ax.set_ylabel("small-root bound exponent δ  (|x| < p^δ)")
ax.set_ylim(0.005, 1.5)
ax.legend(fontsize=8, loc="upper right", frameon=False)
ax.set_title("Coppersmith reach: multi-polynomial / shared-target vs needed δ (derived, exact LP)", fontsize=10)
save(fig, "reach_vs_family")

# ---------- 2. finite ℓ convergence for S3 systems ----------
fig, ax = plt.subplots(figsize=(7.4, 4.4))
series = [("C1 S3", "single S₃ (W10 control)", C["single"], "o"),
          ("N1 S3 x S3 disjoint", "null: two S₃, disjoint unknowns", C["null"], "x"),
          ("V3 S3 pair s=1", "shared target: 2 × S₃, 1 shared summand", C["pair"], "s"),
          ("V3 S3 triple s=1", "shared target: 3 × S₃, 1 shared summand", C["triple"], "^")]
extrap = {}
for key, lab, col, mk in series:
    pts = sorted((r["ell"], r["delta"]) for r in lp if r["system"] == key)
    if not pts:
        continue
    l, d = np.array(pts).T
    ax.plot(l, d, marker=mk, color=col, label=lab, lw=1.3, ms=5)
    if len(l) >= 5:
        L, D = l[-5:], d[-5:]
        A = np.vstack([np.ones_like(L), -1 / L, 1 / L ** 2]).T
        coef = np.linalg.lstsq(A, D, rcond=None)[0]
        extrap[key] = {"ells": L.tolist(), "delta_inf_fit": float(coef[0])}
for a, lab, col, dy in ((1 / 6, "1/6 (single, ℓ→∞)", C["single"], -0.004), (0.177886, "pair, box ℓ→∞: 0.1779", C["pair"], -0.003),
                    (0.182015, "triple, box ℓ→∞: 0.1820", C["triple"], 0.003), (0.25, "ceiling 1/4", C["ceil"], 0.003)):
    ax.axhline(a, color=col, ls=":", lw=1)
    ax.text(12.4, a + dy, lab, fontsize=7, va="center", color=col)
ax.set_xlabel("multiplicity ℓ")
ax.set_ylabel("exact optimum δ_ℓ (down-closed M, p → ∞)")
ax.set_xlim(0.5, 16.5)
ax.legend(fontsize=8, frameon=False, loc="lower right")
ax.set_title("S₃ systems: exact closure optimum per ℓ (Rust min-cut, exact rationals)", fontsize=10)
save(fig, "finite_ell")
json.dump(extrap, open(os.path.join(run, "extrapolation.json"), "w"), indent=1)

# ---------- 3. LLL ----------
fig, ax = plt.subplots(figsize=(7.4, 4.2))
any_ = False
for fn, bits, ls in (("lll40_l1.json", 40, "-"), ("lll40_l1_fine.json", 40, "-"), ("lll56_l1.json", 56, "--")):
    p = os.path.join(run, fn)
    if not os.path.exists(p):
        continue
    recs = [r for r in json.load(open(p))["records"] if r.get("part") == "lll"]
    for sysname, col, lab in (("C1 S3", C["single"], "single S₃ (dim 25)"), ("V3 S3 pair s=1", C["pair"], "S₃ pair (dim 45)")):
        rr = sorted((r["delta"], r["success"], r["trials"]) for r in recs if r["system"] == sysname)
        if not rr:
            continue
        dd = [a for a, _, _ in rr]
        fr = [s / t for _, s, t in rr]
        ax.plot(dd, fr, ls=ls, marker="o", color=col, ms=4,
                label=f"{lab}, {bits}-bit, n={rr[0][2]}/pt" + (" (fine grid)" if "fine" in fn else ""),
                alpha=0.55 if "fine" not in fn and bits == 40 else 1.0)
        any_ = True
if any_:
    ax.axvline(0.09, color=C["single"], ls=":", lw=1)
    ax.axvline(0.094737, color=C["pair"], ls=":", lw=1)
    ax.text(0.091, 1.03, "LP δ₁: 0.0900 / 0.0947", fontsize=7)
    ax.set_xlabel("δ (planted |x_i| < p^δ)")
    ax.set_ylabel("fraction of trials with ≥ n reduced rows vanishing over Z")
    ax.set_ylim(-0.03, 1.1)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("LLL at ℓ = 1, closure-optimal lattices (repo LLL, δ_LLL = 0.99)", fontsize=10)
    save(fig, "lll")
print(json.dumps(extrap, indent=1))
