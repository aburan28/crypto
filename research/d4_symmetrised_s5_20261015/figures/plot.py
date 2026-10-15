"""Lay out the recorded values of research/d4_symmetrised_s5_20261015 as SVG
(and 150 dpi PNG) figures.  Plotting only: every number drawn is read from
raw/*.jsonl or summary.json (written by the Rust summariser); nothing is
computed here beyond placing points.

    python plot.py   (from this directory)
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
RAW = STUDY / "raw"

plt.rcParams.update({"font.size": 9, "svg.fonttype": "none"})
C_W, C_C8, C_PLAIN, C_GEN = "#1f6feb", "#d1495b", "#555555", "#2a9d8f"


def rows(name):
    p = RAW / name
    if not p.exists():
        return []
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()]


data = rows("tiny.jsonl") + rows("measured.jsonl") + rows("w_only.jsonl")
controls = rows("controls.jsonl")
summary = json.loads((STUDY / "summary.json").read_text())
ps = sorted({r["p"] for r in data})
xpos = {p: i for i, p in enumerate(ps)}

# Figure 1: cost of one solve, every target.
fig, ax = plt.subplots(figsize=(7.2, 4.0))
for r in data:
    x = xpos[r["p"]]
    if r["w_solved"]:
        rer = r["w_degree_used"] == 27
        ax.scatter(x + 0.12, r["w_solve"]["total_muls"], s=14, color=C_W,
                   marker="^" if rer else "o", alpha=0.75, zorder=3)
    if r.get("c8_solved"):
        ax.scatter(x - 0.12, r["c8_solve"]["total_muls"], s=14, color=C_C8,
                   alpha=0.75, zorder=3)
for c in controls:
    if c.get("mode") == "generic" and c["p"] in xpos:
        ax.scatter(xpos[c["p"]] + 0.3, c["total_muls"], s=14, color=C_GEN,
                   marker="s", alpha=0.8, zorder=3)
c4 = summary["c4_plain_quoted"]
ax.axhline(c4, color=C_PLAIN, lw=1.2)
ax.text(len(ps) - 0.5, c4 * 1.25, "plain S₄ form, quoted (§11.17): 1.213·10¹²  — not run, 3.42 GB > cap",
        ha="right", fontsize=7.5, color=C_PLAIN)
ax.axhline(c4 / 630, color=C_PLAIN, lw=1, ls="--")
ax.text(-0.45, c4 / 630 * 0.62, "projection C₄/630 (8^3.1)", fontsize=7.5, color=C_PLAIN)
ax.set_yscale("log")
ax.set_ylim(5e8, 4e12)
ax.set_xticks(range(len(ps)))
ax.set_xticklabels([str(p) for p in ps])
ax.set_xlabel("p   (curve over F_{p⁴}; 13–37: brute-force ground truth; 269–4093: meet in the middle)")
ax.set_ylabel("F_p multiplications per solve")
ax.scatter([], [], color=C_C8, s=14, label="C8: five unknowns, degree 14 (measured)")
ax.scatter([], [], color=C_W, s=14, label="w: four weighted unknowns, degree 26 (measured)")
ax.scatter([], [], color=C_W, s=14, marker="^", label="w, rerun at degree 27 (both runs charged)")
ax.scatter([], [], color=C_GEN, s=14, marker="s", label="w on planted generic systems (positive control)")
ax.legend(loc="center right", fontsize=7.5, frameon=False)
ax.set_title("One D₄-symmetrised S₅ solve, every target (1 thread; study d4_symmetrised_s5_20261015)",
             fontsize=9)
ax.grid(axis="y", which="major", lw=0.3, alpha=0.5)
fig.tight_layout()
fig.savefig(HERE / "c_per_solve.svg")
fig.savefig(HERE / "c_per_solve.png", dpi=150)

# Figure 2: the handover in E18's tau_T frame (derived).
fig, ax = plt.subplots(figsize=(7.2, 3.4))
hv = summary["handover_tau_frame"]
names = [h["solve"] for h in hv]
vals = [h["nstar_bits"] for h in hv]
cols = [C_PLAIN, C_C8, C_W, C_W]
y = list(range(len(hv)))[::-1]
ax.barh(y, vals, color=cols, alpha=0.85, height=0.6)
for yi, v, h in zip(y, vals, hv):
    ax.text(v + 0.3, yi, f"2^{v:.1f}   (C = {h['c']:.2e})", va="center", fontsize=7.5)
for ref, lab in [(151.1, "plain, prime order, unfolded walk (§11.19)"),
                 (156.1, "plain, prime order, negation walk (E18)"),
                 (146.0, "ψ-folded, j = 0 (E18)")]:
    ax.axvline(ref, color="#888", lw=0.8, ls=":")
    ax.text(ref + 0.3, -0.45, lab, rotation=90, fontsize=6, color="#555", ha="left", va="bottom")
ax.set_yticks(y)
ax.set_yticklabels(names, fontsize=8)
ax.set_xlim(110, 172)
ax.set_xlabel("log₂ n*  — τ_T-folded base, h = 4, y = 106 residuals/relation, r = 0.341 (E18); extrapolated")
ax.set_title("Where the k = 4 handover moves with the solve (derived, not measured)", fontsize=9)
fig.tight_layout()
fig.savefig(HERE / "handover.svg")
fig.savefig(HERE / "handover.png", dpi=150)

# Figure 3: phases of one paired target per measured size (recorded values).
fig, ax = plt.subplots(figsize=(7.2, 3.0))
phases = [("echelon_muls", "elimination"), ("nf_muls", "normal forms"),
          ("charpoly_muls", "char. poly."), ("eigenvector_muls", "eigenvectors")]
picked = [r for r in rows("measured.jsonl")
          if r["target"] == 1 and r.get("c8_solved") and r["w_solved"]]
width = 0.8 / max(1, 2 * len(picked))
k = 0
shades = {"C8": ["#d1495b", "#e07a87", "#eeb0b8"], "w": ["#1f6feb", "#6a9ff3", "#b0cdf9"]}
for i, r in enumerate(picked):
    for arm, key in [("C8", "c8_solve"), ("w", "w_solve")]:
        xs = [j + (k - len(picked) + 0.5) * width for j in range(len(phases))]
        ax.bar(xs, [max(r[key][ph], 1) for ph, _ in phases], width=width,
               color=shades[arm][i % 3], label=f"{arm}, p = {r['p']}")
        k += 1
ax.set_xticks(range(len(phases)))
ax.set_xticklabels([lab for _, lab in phases])
ax.set_yscale("log")
ax.set_ylim(1e6, 1e11)
ax.set_ylabel("F_p multiplications")
ax.legend(fontsize=7, frameon=False, ncol=3)
ax.set_title("Phases of one solve (target 1 of each measured cell; 0 shown at 1)", fontsize=9)
fig.tight_layout()
fig.savefig(HERE / "phases.svg")
fig.savefig(HERE / "phases.png", dpi=150)
print("wrote c_per_solve, handover, phases (.svg, .png)")
