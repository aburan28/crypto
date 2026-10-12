"""Plot recorded values only (rational.json, generic.json, elkies.json,
analysis.json written by examples/r4_16_3_vertical_margin.rs). No computation
beyond taking log2 of recorded numbers for axes and drawing recorded fits."""
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

D = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
rat = json.loads((D / "rational.json").read_text())["rows"]
gen = json.loads((D / "generic.json").read_text())["rows"]
elk = json.loads((D / "elkies.json").read_text())["rows"]
ana = json.loads((D / "analysis.json").read_text())

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "svg.fonttype": "none"})


def fit(series, rng):
    for f in ana["fits"]:
        if f["series"].startswith(series) and f["range"] == rng:
            return f
    raise KeyError(series)


def line(ax, f, x0, x1, color):
    xs = [x0, x1]
    ax.plot(xs, [f["log2_constant"] + f["exponent"] * x for x in xs], color=color,
            lw=1, ls="--", alpha=0.8)


# Figure 1: counted cost against ell.
fig, ax = plt.subplots(figsize=(7.2, 4.6))
L = lambda r: math.log2(r["ell"])  # noqa: E731
lg = lambda v: math.log2(v)  # noqa: E731
ax.plot([L(r) for r in rat], [lg(r["find"]["Meq"]) for r in rat], "o", ms=3, color=C[0],
        label="rational kernel: find (cofactor ladder)")
ax.plot([L(r) for r in rat if r["seed"] == 0], [lg(r["velu_1pt"]["Meq"]) for r in rat if r["seed"] == 0],
        "s", ms=3.5, color=C[1], label="rational kernel: Vélu, codomain + 1 image")
sq = [r for r in rat if r["seed"] == 0 and r["sqrtvelu_1pt"]]
ax.plot([L(r) for r in sq], [lg(r["sqrtvelu_1pt"]["Meq"]) for r in sq], "^", ms=3.5, color=C[2],
        label="rational kernel: √élu, codomain + 1 image")
ax.plot([L(r) for r in gen], [lg(r["find"]["Meq"]) for r in gen], "D", ms=2.5, color=C[3],
        label="generic kernel k=(ℓ−1)/2: find (ladder in F_{p^k} + minpoly)")
ax.plot([L(r) for r in gen], [lg(r["setup"]["Meq"]) for r in gen], "x", ms=3, color=C[4],
        label="generic kernel: F_{p^k} setup (irreducible search)")
ax.plot([L(r) for r in gen], [lg(r["eval_2pt_from_h"]["Meq"] / 2) for r in gen], "v", ms=2.5,
        color=C[5], label="generic kernel: evaluate 1 point from h (incl. 1 I)")
line(ax, fit("velu_1pt", "ell >= 2^10"), 10, 24, C[1])
line(ax, fit("sqrtvelu_1pt", "ell >= 2^16"), 16, 24, C[2])
line(ax, fit("generic_find_expected", "ell >= 30"), math.log2(30), math.log2(440), C[3])
r64 = [t for t in ana["crossover"] if t["g_mults_per_group_op"] == 6 and t["log2_r"] == 64][0]
ax.axhline(r64["log2_rho_mults"], color=INK2, lw=0.8, ls=":")
ax.text(13.0, r64["log2_rho_mults"] + 0.4,
        "rho on r=2^64 with 6 M/op (derived, 2^%.1f M)" % r64["log2_rho_mults"], color=INK2, fontsize=8)
ax.set_xlabel("log2 ℓ")
ax.set_ylabel("log2 F_p multiplication-equivalents (M + 91·I), 61-bit p")
ax.set_title("Vertical ℓ-isogeny on j=−3375 crater curves: counted cost per phase (measured)",
             fontsize=9.5, color=INK)
ax.grid(color="#e4e3df", lw=0.5)
ax.legend(fontsize=7, frameon=True, framealpha=0.9, edgecolor="none", loc="lower right")
fig.text(0.01, 0.005, "Dashed: OLS fits from analysis.json. Rational: 3 instances per ℓ (Vélu/√élu counts are"
         " identical across instances). Generic: 1 instance per ℓ.", fontsize=6.5, color=INK2)
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(D / "fig1_cost_vs_ell.svg")

# Figure 2: derived crossover ell* (g = 6), extrapolated.
fig, ax = plt.subplots(figsize=(7.2, 3.8))
tabs = [t for t in ana["crossover"] if t["g_mults_per_group_op"] == 6]
models = [r["model"] for r in tabs[0]["rows"]]
short = ["rational, Vélu", "rational, √élu (fit)", "rational, √élu (√ℓ log²ℓ model)",
         "generic (fit)", "generic (fit, CI low)", "generic (ℓ² log ℓ model)"]
mk = ["s", "^", "^", "D", "D", "D"]
for i, t in enumerate(tabs):
    x = t["log2_r"]
    for j, r in enumerate(t["rows"]):
        v = r["log2_ell_star"]
        if v is None:
            continue
        ax.plot(x + (j - 2.5) * 6, v, mk[j], color=C[j], ms=5, mfc=C[j] if j in (0, 1, 3) else "none",
                label=short[j] if i == 0 else None)
    ax.hlines(t["log2_ell_max"], x - 20, x + 20, color=INK, lw=1.2)
    ax.text(x + 21, t["log2_ell_max"], "ℓ_max=2√p", fontsize=7, color=INK, va="center")
ax.set_xticks([64, 128, 256], ["r = 2^64", "r = 2^128", "r = 2^256"])
ax.set_xlim(30, 300)
ax.set_ylabel("log2 ℓ* (find + evaluate = rho)")
ax.set_title("Derived crossover ℓ* (EXTRAPOLATED from 61-bit counts; 6 F_p mults per rho step)",
             fontsize=9.5, color=INK)
ax.grid(color="#e4e3df", lw=0.5)
ax.legend(fontsize=7, frameon=False, loc="upper left", ncol=2)
fig.text(0.01, 0.005, "Points above the ℓ_max bar: no conductor gap that large exists, so no crossover.",
         fontsize=6.5, color=INK2)
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(D / "fig2_crossover.svg")

# Figure 3: library Elkies route, contended wall time.
fig, ax = plt.subplots(figsize=(7.2, 3.0))
ax.plot([math.log2(r["ell"]) for r in elk], [math.log2(max(r["phi_wall_s_contended"], 1e-4)) for r in elk],
        "o", ms=3.5, color=C[0], label="library Φ_ℓ mod p (q-series)")
ax.plot([math.log2(r["ell"]) for r in elk],
        [math.log2(max(r["elkies_first_kernel_wall_s_contended"] or 1e-4, 1e-4)) for r in elk],
        "s", ms=3.5, color=C[1], label="library Elkies kernel from one root")
ax.plot([math.log2(r["ell"]) for r in gen], [math.log2(r["find_wall_s_contended"]) for r in gen],
        "D", ms=2.5, color=C[3], label="this example: F_{p^k} torsion route (incl. setup)")
ax.set_xlabel("log2 ℓ")
ax.set_ylabel("log2 seconds (contended)")
ax.set_title("Generic kernel finding, contended wall time (4 cores shared by ~6 agents)",
             fontsize=9.5, color=INK)
ax.grid(color="#e4e3df", lw=0.5)
ax.legend(fontsize=7, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(D / "fig3_wall_generic.svg")
print("ok")
