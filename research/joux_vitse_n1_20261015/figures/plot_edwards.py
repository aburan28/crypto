# Layout only: plots values recorded in ../summary_edwards.json and
# ../summary.json (computed natively by examples/joux_vitse_edwards_summary.rs
# and examples/joux_vitse_summary.rs).  No result is computed here.
import json
import math
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
E = json.loads((HERE.parent / "summary_edwards.json").read_text())
W = json.loads((HERE.parent / "summary.json").read_text())
plt.rcParams.update({"font.size": 9})
BLUE, ORANGE, GREEN, GREY, RED, PURPLE = "#2a6fdb", "#e07b00", "#2c9d5a", "#777777", "#c0392b", "#7b4fd6"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(HERE / f"{name}.svg")
    fig.savefig(HERE / f"{name}.png", dpi=150)
    plt.close(fig)


# 1. one 2-torsion test: committed solve against trace replay, same residuals
fig, ax = plt.subplots(figsize=(6.2, 3.4))
labels = [f"p={c['p']}, s={c['seed']}\n|F|={c['base_points']}" for c in E["cost"]]
x = range(len(labels))
for off, key, col, lab in [(-0.12, "full", GREY, "committed test (f4_fp::solve), today's engine"),
                           (0.12, "replay", BLUE, "F4 trace replay (+ solve when consistent)")]:
    m = [c[f"{key}_mean"] for c in E["cost"]]
    lo = [c[f"{key}_mean"] - c[f"{key}_min"] for c in E["cost"]]
    hi = [c[f"{key}_max"] - c[f"{key}_mean"] for c in E["cost"]]
    ax.errorbar([i + off for i in x], m, yerr=[lo, hi], fmt="o", color=col, capsize=3, label=lab)
ax.axhline(E["committed_reference"]["c_second"], color=ORANGE, ls="--", lw=1,
           label="RESEARCH_K5 §3: C″ = 7.3·10⁷ (2026-09-29 engine)")
ax.set_xticks(list(x), labels)
ax.set_ylim(0, 8.5e7)
ax.set_ylabel("F_p multiplications per test")
ax.set_title("2-torsion Edwards test: replay against the committed solve (32 residuals each)")
ax.legend(fontsize=7, loc="lower right")
save(fig, "edwards_replay")

# 2. IC / rho: measured complete DLPs and the two models, against log2 of the group
fig, ax = plt.subplots(figsize=(6.2, 3.8))
me = E["model_extrapolation"]
ax.plot([r["log2_n"] for r in me["rows"]], [r["ic_over_rho"] for r in me["rows"]], color=PURPLE, ls="--",
        label=f"Edwards + replay, model (extrapolation): parity n ≈ 2^{me['log2_n_star']:.1f}")
mw = W["model_extrapolation"]
ax.plot([5 * r["log2_q"] for r in mw["rows"]], [r["ic_over_rho"] for r in mw["rows"]], color=GREY, ls="--",
        label=f"plain Weierstrass + replay, model: parity ≈ 2^{mw['log2_group_at_q_star']:.0f}")
for d in E["dlp"]:
    if d.get("complete") and d.get("verified"):
        ax.plot(d["log2_n"], d["ratio_online_ic_over_rho"], "o", color=PURPLE)
for d in W["dlp"]:
    if d["ic_complete"] and d["ic"]["verified"]:
        ax.plot(math.log2(int(d["ell"])), d["ratio_online_ic_over_rho"], "s", color=GREY)
ax.plot([], [], "o", color=PURPLE, label="Edwards, complete DLP (measured)")
ax.plot([], [], "s", color=GREY, label="plain, complete DLP (measured)")
ax.axhline(1, color="black", lw=0.8)
ax.axvline(233, color=ORANGE, ls=":", lw=1, label="RESEARCH_K5: 2^233 without replay")
ax.set_yscale("log")
ax.set_xlabel("log₂ (order of the group the logarithm lives in)")
ax.set_ylabel("IC / rho (F_p multiplications)")
ax.set_title("k = 5 Joux–Vitse routes against rho")
ax.legend(fontsize=6.5, loc="upper right")
save(fig, "edwards_ic_over_rho")
