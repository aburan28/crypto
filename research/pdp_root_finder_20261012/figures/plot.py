# Lays out ANALYSIS.json (written by `pdp_rootfind analyze`) as the study's
# figures.  Plots recorded and derived values only; computes nothing new.
import json, math, os, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
A = json.load(open(os.path.join(here, "..", "ANALYSIS.json")))
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
     "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#e6e5e1", "surface": "#fcfcfb", "wash": "#9aa3ad"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"], "axes.titlecolor": C["ink"],
                     "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})
ARM_STYLE = {"baseline": (C["wash"], "o"), "group-scalar": (C["yellow"], "D"), "rf-pairs": (C["orange"], "s"), "rf-mitm": (C["blue"], "^"),
             "pairs3": (C["orange"], "s"), "mitm3": (C["blue"], "^"), "subtract2": (C["aqua"], "o"), "mitm4": (C["yellow"], "D")}
ARM_LABEL = {"baseline": "PS-0 baseline (pairs-and-solve, quartic + gcd)", "group-scalar": "PS-2 group-law kernel",
             "rf-pairs": "rf-pairs (this study, memoryless)", "rf-mitm": "rf-mitm (this study, 2-sum table)",
             "pairs3": "pairs3 (memoryless, m = 3)", "mitm3": "mitm3 (2-sum table, m = 3)", "subtract2": "subtract2 (m = 2)", "mitm4": "mitm4 (2-sum table, m = 4)"}

def style(ax):
    ax.grid(True, color=C["grid"], linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(here, name + ".svg"))
    fig.savefig(os.path.join(here, name + ".png"), dpi=150)

slopes = {(s["track"], s["arm"], s["quantity"]): s for s in A["slopes"]}

# ── Figure 1: binary per-target wall and counted multiplications vs ℓ ──
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
rows = [r for r in A["binary"] if r["n"] == 3 * r["l"]]
for arm in ["baseline", "group-scalar", "rf-pairs", "rf-mitm"]:
    pts = sorted([(r["l"], r["ns_median"], r["ns_min"]) for r in rows if r["arm"] == arm and r["ns_median"]])
    if pts:
        col, mk = ARM_STYLE[arm]
        sl = slopes.get(("binary", arm, "ns"))
        lab = ARM_LABEL[arm] + (f": {sl['slope_bits_per_rung']:.2f} bits/rung" if sl else "")
        ax1.plot([p[0] for p in pts], [p[1] / 1e3 for p in pts], marker=mk, markersize=7, linewidth=2, color=col,
                 markeredgecolor=C["surface"], markeredgewidth=1.5, label=lab)
        for p in pts:
            r0 = [r for r in rows if r["arm"] == arm and r["l"] == p[0]][0]
            ax1.scatter([p[0]] * len(r0["ns_all"]), [v / 1e3 for v in r0["ns_all"]], color="black", s=4, zorder=3)
ax1.set_yscale("log")
ax1.set_xlabel("factor-base dimension ℓ (n = 3ℓ, Koblitz-family curve y² + xy = x³ + x² + 1)")
ax1.set_ylabel("µs per refuted target (median of isolated rounds; dots: rounds)")
ax1.set_title("Binary: wall time per target, reserved core", fontsize=9)
ax1.legend(fontsize=7, loc="upper left", frameon=False)
style(ax1)
for arm in ["group-scalar", "rf-pairs", "rf-mitm"]:
    pts = sorted([(r["l"], r["mul"]) for r in rows if r["arm"] == arm and r["mul"]])
    if pts:
        col, mk = ARM_STYLE[arm]
        sl = slopes.get(("binary", arm, "mul"))
        lab = ARM_LABEL[arm] + (f": {sl['slope_bits_per_rung']:.2f} bits/rung" if sl else "")
        ax2.plot([p[0] for p in pts], [p[1] for p in pts], marker=mk, markersize=7, linewidth=2, color=col,
                 markeredgecolor=C["surface"], markeredgewidth=1.5, label=lab)
ax2.set_yscale("log")
ax2.set_xlabel("factor-base dimension ℓ")
ax2.set_ylabel("field multiplications per target (counted, deterministic)")
ax2.set_title("Binary: counted multiplications per target (enumeration slope 2 bits/rung, pair table 1)", fontsize=9)
ax2.legend(fontsize=7, loc="upper left", frameon=False)
style(ax2)
fig.suptitle("Root-finding oracles on the PS-0 rungs: per-target cost against the PS-0 and PS-2 references", fontsize=10)
save(fig, "binary_cost")

# ── Figure 2: prime ladders, multiplications and wall vs t ──
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
for m, arm in [(3, "pairs3"), (3, "mitm3"), (2, "subtract2"), (4, "mitm4")]:
    pts = sorted([(r["t"], r["mul"]) for r in A["prime"] if r["m"] == m and r["arm"] == arm and r["mul"]])
    if pts:
        col, mk = ARM_STYLE[arm]
        sl = slopes.get((f"prime m={m}", arm, "mul"))
        lab = ARM_LABEL[arm] + (f": {sl['slope_bits_per_rung']:.2f} bits/t" if sl else "")
        ax1.plot([p[0] for p in pts], [p[1] for p in pts], marker=mk, markersize=7, linewidth=2, color=col,
                 markeredgecolor=C["surface"], markeredgewidth=1.5, label=lab)
# Context: the committed F4 rows on the same kind of cell (finished runs only;
# the unit is F_p multiply-adds, compared one to one with a multiplication).
for m, col, mk in [(2, C["wash"], "o"), (3, C["wash"], "^")]:
    pts = sorted([(r["t"], r["f4_muladds"]) for r in A["f4_prime"] if r["kind"] == "kummer" and r["p"] == 2013265921 and r["m"] == m and r["arm"] in ("subtract2", "mitm3")])
    if pts:
        sl = slopes.get((f"prime m={m}", "f4_fp_tower kummer p=2013265921", "muladds"))
        lab = f"f4_fp_tower, m = {m} (committed rows, multiply-adds)" + (f": {sl['slope_bits_per_rung']:.2f} bits/t" if sl else "")
        ax1.plot([p[0] for p in pts], [p[1] for p in pts], marker=mk, markersize=6, linewidth=1.5, linestyle="--", color=col, label=lab)
ax1.set_yscale("log")
ax1.set_xlabel("tower height t (|V| = 2^t, Kummer towers at p₁ = 2013265921)")
ax1.set_ylabel("field multiplications per target (counted)")
ax1.set_title("Prime: counted multiplications per full scan", fontsize=9)
ax1.legend(fontsize=7, loc="upper left", frameon=False)
style(ax1)
for m, arm in [(3, "pairs3"), (3, "mitm3"), (2, "subtract2"), (4, "mitm4")]:
    pts = sorted([(r["t"], r["ns_median"]) for r in A["prime"] if r["m"] == m and r["arm"] == arm and r["ns_median"]])
    if pts:
        col, mk = ARM_STYLE[arm]
        sl = slopes.get((f"prime m={m}", arm, "ns"))
        lab = ARM_LABEL[arm] + (f": {sl['slope_bits_per_rung']:.2f} bits/t" if sl else "")
        ax2.plot([p[0] for p in pts], [p[1] / 1e3 for p in pts], marker=mk, markersize=7, linewidth=2, color=col,
                 markeredgecolor=C["surface"], markeredgewidth=1.5, label=lab)
        for p in pts:
            r0 = [r for r in A["prime"] if r["m"] == m and r["arm"] == arm and r["t"] == p[0]][0]
            ax2.scatter([p[0]] * len(r0["ns_all"]), [v / 1e3 for v in r0["ns_all"]], color="black", s=4, zorder=3)
ax2.set_yscale("log")
ax2.set_xlabel("tower height t")
ax2.set_ylabel("µs per full scan (median of isolated rounds; dots: rounds)")
ax2.set_title("Prime: wall time per full scan, reserved core", fontsize=9)
ax2.legend(fontsize=7, loc="upper left", frameon=False)
style(ax2)
fig.suptitle("Root-finding oracles on Kummer-tower bases over F_p: per-target cost against t", fontsize=10)
save(fig, "prime_cost")

# ── Figure 3: the committed F4 rows against the root-finders, same cells ──
fig, ax = plt.subplots(figsize=(11, 4.3))
rows = [r for r in A["f4_prime"] if r["kind"] == "kummer" and r["p"] == 2013265921]
cells = sorted({(r["m"], r["t"]) for r in rows})
xs = list(range(len(cells)))
w = 0.26
f4 = {(r["m"], r["t"]): r["f4_muladds"] for r in rows}
ax.bar([x - w for x in xs], [f4[c] for c in cells], width=w, color=C["wash"], label="f4_fp_tower: multiply-adds (committed rows, median)")
for j, (arm, col) in enumerate([("pairs3", C["orange"]), ("mitm3", C["blue"]), ("subtract2", C["aqua"]), ("mitm4", C["yellow"])]):
    vals = []
    for c in cells:
        v = [r["mul"] for r in rows if (r["m"], r["t"]) == c and r["arm"] == arm]
        vals.append(v[0] if v else float("nan"))
    if any(not math.isnan(v) for v in vals):
        ax.bar([x + (j % 2) * w for x in xs], vals, width=w, color=col, label=f"{arm}: field multiplications (same cell, same targets)")
ax.set_yscale("log")
ax.set_xticks(xs)
ax.set_xticklabels([f"m={c[0]}\nt={c[1]}\nN={c[0]*c[1]}" for c in cells], fontsize=8)
ax.set_ylabel("operations per target (log scale)")
ax.set_title("Committed F4 rows (Kummer towers, p₁) against the root-finders on the cells they measured: one bar group per cell", fontsize=9)
ax.legend(fontsize=7.5, loc="upper left", frameon=False)
style(ax)
save(fig, "f4_vs_rootfind")
print("ok")
