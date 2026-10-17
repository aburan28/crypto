# Layout only: plots values recorded in ../summary.json (computed natively by
# examples/joux_vitse_summary.rs).  No result is computed here.
import json
import math
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
S = json.loads((HERE.parent / "summary.json").read_text())
rungs = S["rungs"]
P = [r["p"] for r in rungs]

plt.rcParams.update({"font.size": 9, "figure.dpi": 100})
BLUE, ORANGE, GREEN, GREY, RED = "#2a6fdb", "#e07b00", "#2c9d5a", "#777777", "#c0392b"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(HERE / f"{name}.svg")
    fig.savefig(HERE / f"{name}.png", dpi=150)
    plt.close(fig)


# 1. cost of one trial against p
fig, ax = plt.subplots(figsize=(6.2, 3.6))
mean = [r["c_trial_mean"] for r in rungs]
lo = [r["c_trial_mean"] - r["c_trial_min"] for r in rungs]
hi = [r["c_trial_max"] - r["c_trial_mean"] for r in rungs]
ax.errorbar(P, mean, yerr=[lo, hi], fmt="o-", color=BLUE, capsize=3,
            label="measured: F4 trial with trace replay (mean, min–max; 64 trials/rung)")
th = S["thesis_reference"]
ax.axhline(th["trial_muls_f4"], color=GREY, ls="--", lw=1, label="Vitse 2011, Rem. 7.3.6: F4 ≈ 10¹⁰")
ax.axhline(th["trial_muls_stated_7_3_3"], color=GREY, ls=":", lw=1.2, label="Vitse 2011, §7.3.3: ≈ 3·10⁹ per test")
ax.axhline(S["d4_n4_reference"]["c4"], color=GREEN, ls="-.", lw=1, label="n = 4 D₄ oracle C₄ (PR #1863): 6.9·10⁷")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("p   (curve over F_{p⁵})")
ax.set_ylabel("F_p multiplications per trial")
ax.set_title("One (n − 1)-decomposition trial at n = 5: cost against p")
ax.legend(fontsize=7, loc="center right")
save(fig, "trial_cost")

# 2. trials per relation against p
fig, ax = plt.subplots(figsize=(6.2, 3.6))
good = [r for r in rungs if r["rate_hits"] >= 1]
x = [r["p"] for r in good]
y = [r["trials_per_relation"] for r in good]
yl = [r["trials_per_relation"] - r["trials_per_relation_wilson95"][0] for r in good]
yh = [r["trials_per_relation_wilson95"][1] - r["trials_per_relation"] for r in good]
ax.errorbar(x, y, yerr=[yl, yh], fmt="o", color=BLUE, capsize=3,
            label="measured: trials per decomposition (Wilson 95 %)")
ax.plot(P, [24 * p for p in P], color=GREY, ls="--", label="thesis heuristic (n − 1)! q = 24 q")
ax.plot(P, [1 / r["expected_decompositions_per_trial_16C"] for r in rungs], color=ORANGE, ls=":",
        label="floor: 1 / E[#decompositions] = p⁵ / (16·C(|F|, 4))")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("p")
ax.set_ylabel("trials per relation")
ax.set_title("Decomposition rate (distinct base points; MITM-checked)")
ax.legend(fontsize=7)
save(fig, "trials_per_relation")

# 3. IC / rho: measured DLP arms and the model (extrapolation)
fig, ax = plt.subplots(figsize=(6.2, 3.6))
m = S["model_extrapolation"]
xs = [row["log2_q"] for row in m["rows"]]
ys = [row["ic_over_rho"] for row in m["rows"]]
ax.plot(xs, ys, color=GREY, ls="--", label="model (extrapolation): measured C_trial, trials = κq, rho with negation")
for d in S["dlp"]:
    if d["ic_complete"] and d["ic"]["verified"]:
        ax.plot(math.log2(int(d["p"])), d["ratio_online_ic_over_rho"], "o", color=RED)
ax.plot([], [], "o", color=RED, label="measured: complete DLP, IC online / rho (same target, F_p mults)")
ax.axhline(1, color="black", lw=0.8)
ax.axvline(m["log2_q_star"], color=ORANGE, ls=":", label=f"model crossover q* ≈ 2^{m['log2_q_star']:.1f}")
ax.axvline(th["crossover_base_field_bits_odd_p"], color=GREEN, ls="-.", lw=1, label="Vitse 2011 §7.3.3: q ≈ 2^60")
ax.set_yscale("log")
ax.set_xlabel("log₂ q   (group ≈ q⁵)")
ax.set_ylabel("IC cost / rho cost")
ax.set_title("(n − 1)-decomposition at n = 5 against rho")
ax.legend(fontsize=6.5, loc="upper right")
save(fig, "ic_over_rho")
print("ok", file=sys.stderr)
