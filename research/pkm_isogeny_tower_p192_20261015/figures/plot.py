"""Lay out the recorded and derived values of this study (no computation).

Reads ../derived.json (written by examples/two_sylow_study_stats.rs) and
writes SVG + 150-dpi PNG for each figure next to this file.
"""
import json
import math
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
D = json.loads((HERE.parent / "derived.json").read_text())
plt.rcParams.update({"font.size": 9, "svg.fonttype": "none"})
C8, C2, CS, CG = "#1f6fb4", "#d9822b", "#2a9d5c", "#777777"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(HERE / f"{name}.svg")
    fig.savefig(HERE / f"{name}.png", dpi=150)
    plt.close(fig)


# 1. Draws per success against t, measured and extrapolated.
fig, ax = plt.subplots(figsize=(7.2, 4.4))
for fam, col, lab in (("z2xz8", C8, "Z/2×Z/8 family"), ("z2xz2", C2, "Z/2×Z/2 family")):
    # Counts below 3 are drawn separately (a single deep draw makes a flat run).
    rows = [r for r in D["distributions"][fam]["rows"] if r["count_ge"] >= 3]
    t = [r["t"] for r in rows]
    y = [math.log2(r["draws_per_success"]) for r in rows]
    lo = [y[i] - math.log2(r["draws_per_success_95"][0]) for i, r in enumerate(rows)]
    # A count whose 95% lower bound is 0 has an unbounded upper end (null in
    # derived.json); drawn as an arrow-length bar of 4 bits.
    hi = [(math.log2(r["draws_per_success_95"][1]) - y[i]) if r["draws_per_success_95"][1] else 4.0
          for i, r in enumerate(rows)]
    ax.errorbar(t, y, yerr=[lo, hi], fmt="o", ms=3.5, color=col, capsize=2,
                label=f"{lab}: measured 1/P(b≥t), 95% CI (N = {D['distributions'][fam]['draws']:,})")
deep = D["distributions"]["z2xz8"]["deepest"][0]
ax.annotate(f"deepest z2xz8 draw: b = {deep[1]} (seed {deep[0]}),\none draw in 2^28; counts < 3 not plotted",
            xy=(deep[1], 28), xytext=(36, 12), fontsize=7, color=C8,
            arrowprops=dict(arrowstyle="->", color=C8))
ax.plot([deep[1]], [28], "x", color=C8)
ts = list(range(2, 49))
ax.plot(ts[2:], [t - 3 for t in ts[2:]], color=C8, lw=1, ls="--",
        label="model 2^(t−3) (H2), extrapolated past the data")
ax.plot(ts, [t - math.log2(3) for t in ts], color=C2, lw=1, ls="--",
        label="model 2^t/3 (H1)")
ax.plot(ts, [t - 2 for t in ts], color=CG, lw=0.8, ls=":", label="card's earlier hypothesis 2^(t−2) (k = 2)")
ax.plot(ts, [t - 4 for t in ts], color=CG, lw=0.8, ls="-.", label="card's earlier hypothesis 2^(t−4) (k = 4)")
ax.scatter([48, 64], [0, 0], marker="*", s=90, color=CS, zorder=5,
           label="supersingular y² = x³ + x: 1–2 draws for any t ≤ 64 (C6)")
ax.axvline(48, color=CG, lw=0.6)
ax.text(48.3, 30, "t = 48\n(card)", fontsize=8, color=CG)
ax.set_xlabel("t (point of order 2^t wanted)")
ax.set_ylabel("log2(draws per curve with b ≥ t)")
ax.set_xlim(0, 66)
ax.set_title("P-192 auxiliary-curve search: draws needed vs t (C3, C4, C6)")
ax.legend(fontsize=6.5, loc="upper left", framealpha=0.95)
ax.set_ylim(-3, 62)
ax.grid(alpha=0.3)
save(fig, "draws_vs_t")

# 2. Per-draw cost by arm (C1).
fig, ax = plt.subplots(figsize=(6.4, 3.6))
arms = ["baseline", "fast-chi", "fast"]
for j, (fam, col) in enumerate((("z2xz8", C8), ("z2xz2", C2))):
    xs = [i + (j - 0.5) * 0.36 for i in range(3)]
    meds = [D["c1"][f"{fam}/{a}"]["median_us"] for a in arms]
    lo = [m - D["c1"][f"{fam}/{a}"]["min_us"] for m, a in zip(meds, arms)]
    hi = [D["c1"][f"{fam}/{a}"]["max_us"] - m for m, a in zip(meds, arms)]
    ax.bar(xs, meds, width=0.34, color=col, label=fam, yerr=[lo, hi], capsize=3)
    for x, m in zip(xs, meds):
        ax.text(x, m + 1.2, f"{m:.1f}", ha="center", fontsize=8)
ax.set_xticks(range(3), arms)
ax.set_ylabel("µs per draw (median of 7; bars min–max)")
ax.set_title("C1: cost per draw at P-192, 65,536 draws per rep, one core (shared host)")
ax.legend()
ax.grid(axis="y", alpha=0.3)
save(fig, "per_draw_cost")

# 3. Observed / model ratio of the depth counts.
fig, ax = plt.subplots(figsize=(6.4, 3.4))
for fam, col in (("z2xz8", C8), ("z2xz2", C2)):
    rows = [r for r in D["distributions"][fam]["rows"] if r["tested"]]
    t = [r["t"] for r in rows]
    ratio = [r["ratio_obs_model"] for r in rows]
    lo = [r["model_interval_999"][0] / r["model_count"] for r in rows]
    hi = [r["model_interval_999"][1] / r["model_count"] for r in rows]
    ax.plot(t, ratio, "o", ms=3.5, color=col, label=f"{fam}: observed / model count")
    ax.fill_between(t, lo, hi, color=col, alpha=0.15, label=f"{fam}: 99.9% Poisson band of the model")
ax.axhline(1, color=CG, lw=0.8)
ax.set_xlabel("t")
ax.set_ylabel("count(b ≥ t) / model count")
ax.set_title("H1/H2: depth counts against the Haar model (tested t: model count ≥ 20)")
ax.legend(fontsize=7)
ax.grid(alpha=0.3)
save(fig, "obs_vs_model")

# 4. Diagram of one draw (method, no data).
fig, ax = plt.subplots(figsize=(7.2, 3.4))
ax.axis("off")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
boxes = [
    (0.03, 0.78, "1. draw m (seed → SplitMix64)\nY² = X(X + B⁴)(X + A⁴)\nA = m² − 1, B = 2m, C = m² + 1"),
    (0.53, 0.78, "2. H₄ = (A²B², A²B²C²), order 4\nroots r = (AB, BC, AC): free\n⇒ b ≥ 3 for every member; η = χ(A² − B²)"),
    (0.03, 0.30, "3. half H of P_j (signs of r):\nx_H = x + r₁r₂ + r₁r₃ + r₂r₃\ni_k = (x_H − e_k)^((p−3)/4)  [one exponentiation]"),
    (0.53, 0.30, "4. i_k²·(x_H − e_k) ≠ 1 ⇒ stop, b = j + 1\nelse u = ((x_H − e_k)(x_H − e_o))^((p−3)/4)\nselects H or H + T_n; next roots r_k, r_o, r_m"),
]
for x, y, txt in boxes:
    ax.text(x, y, txt, fontsize=7.5, va="center", ha="left",
            bbox=dict(boxstyle="round", fc="#eef3fa", ec=C8))
ax.annotate("", xy=(0.52, 0.78), xytext=(0.45, 0.78), arrowprops=dict(arrowstyle="->"))
ax.annotate("", xy=(0.25, 0.42), xytext=(0.70, 0.64), arrowprops=dict(arrowstyle="->"))
ax.annotate("", xy=(0.52, 0.30), xytext=(0.47, 0.30), arrowprops=dict(arrowstyle="->"))
ax.annotate("", xy=(0.20, 0.16), xytext=(0.75, 0.16),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-0.2", color=CG))
ax.text(0.5, 0.02, "level j → j + 1 (P ≈ 1/2 per level, model H2); b ≥ t ⇒ Q = P_t has order exactly 2^t",
        fontsize=7.5, color=CG, ha="center")
ax.set_title("One draw of the fast arm (2-Sylow Z/2 × Z/2^b, p ≡ 3 mod 4)", fontsize=9)
save(fig, "halving_chain")
print("figures written")
