# Lays out the recorded solves (raw/solves_*.jsonl, written by
# examples/gaudry_quartic_d4.rs, and the S4 solves of
# ../gaudry_quartic_replay_20261013/raw/s5_p269.jsonl) and draws §11.16's
# crossover model with the measured inputs.  Left: C4 by phase, S4 against D4.
# Right: S/rho against log2 n on the model, S4 (h = 1) and D4 (h = 1, 2).
import json, os, math, statistics as st
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
raw = os.path.join(here, "..", "raw")
s4 = [json.loads(l)["baseline"] for l in open(os.path.join(here, "..", "..", "gaudry_quartic_replay_20261013", "raw", "s5_p269.jsonl"))]
d4 = [json.loads(l)["solve"] for f in ("solves_planted.jsonl", "solves_uniform.jsonl") for l in open(os.path.join(raw, f))]
d15 = [json.loads(l)["solve"] for l in open(os.path.join(raw, "solves_planted_deg15.jsonl"))]
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "plum": "#8e5bd1", "ink": "#0b0b0b",
     "ink2": "#52514e", "grid": "#e6e5e1", "surface": "#fcfcfb", "wash": "#9aa3ad"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"],
                     "xtick.color": C["ink2"], "ytick.color": C["ink2"],
                     "figure.facecolor": C["surface"], "axes.facecolor": C["surface"]})
PH = [("echelon_muls", "elimination", C["blue"]), ("nf_muls", "normal forms", C["orange"]),
      ("charpoly_muls", "char. polynomial", C["aqua"]), ("eigenvector_muls", "eigenvectors", C["plum"])]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.6, 3.7), gridspec_kw={"width_ratios": [1, 1.25]})
arms = [("S₄, x-line\n(degree 29)", s4), ("D₄, Y-line\n(degree 14)", d4), ("D₄, Y-line\n(degree 15)", d15)]
for i, (lab, rows) in enumerate(arms):
    tot = st.median(r["total_muls"] for r in rows)
    a1.bar(i, tot, color=C["blue"] if i else C["wash"], width=0.6)
    a1.text(i, tot * 1.15, f"{tot:.3g}", ha="center", fontsize=8, color=C["ink"])
a1.set_yscale("log")
a1.set_ylim(1e9, 1e13)
a1.set_xticks(range(3), [a[0] for a in arms], fontsize=8)
a1.set_ylabel("F_p multiplications per decomposition (C₄)")
a1.grid(axis="y", color=C["grid"], lw=0.6)
a1.set_title("one decomposition, p = 269 (median)", fontsize=9)
# §11.16 model: S/rho = coef * C / p + r_inf, in rho units with cofactor h.
lg = [100 + 0.5 * k for k in range(141)]
c4s, c4d = st.median(r["total_muls"] for r in s4), st.median(r["total_muls"] for r in d4)
def ratio(log2n, coef, c, rinf, h):
    p = 2 ** (log2n / 4)
    return math.sqrt(h) * (coef * c / p + rinf)
for (lab, coef, c, rinf, h, col, ls) in [
        ("S₄ x-line, h = 1", 12 / 126, c4s, 0.34, 1, C["wash"], "-"),
        ("D₄ Y-line, h = 1", 181 * 0.25 / 126, c4d, 0.34 / 4, 1, C["blue"], "-"),
        ("D₄ Y-line, h = 2", 181 * 0.25 / 126, c4d, 0.34 / 4, 2, C["blue"], "--")]:
    a2.plot(lg, [ratio(x, coef, c, rinf, h) for x in lg], color=col, ls=ls, label=lab)
a2.axhline(1, color=C["ink2"], lw=0.8)
a2.text(101, 1.15, "rho", fontsize=8, color=C["ink2"])
a2.set_yscale("log")
a2.set_ylim(0.05, 1e4)
a2.set_xlabel("log₂ group order n")
a2.set_ylabel("index calculus ÷ rho (model, r∞ at φ = 0.73)")
a2.grid(color=C["grid"], lw=0.6)
a2.legend(frameon=False, fontsize=8)
a2.set_title("§11.16 crossover model with the measured C₄", fontsize=9)
fig.tight_layout()
for ext, kw in (("svg", {}), ("png", {"dpi": 150})):
    fig.savefig(os.path.join(here, f"c4_d4.{ext}"), **kw)
