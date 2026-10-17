"""Lay out the recorded values of the small-x yield runs (data/*.json).

Plotting only: every point is a field of a run's `report` object, written by
`ic prime`. Base-2 logarithms are the report's own `*_log2` fields.

    python3 graphs/plot_prime.py          # from this study directory
"""
import glob
import json
import math
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "legend.fontsize": 8.5})
HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
RUN = "run 2026-10-14, seed 7, 3M trials per toy run, 2M per P-192 run"


def load(pattern):
    return [json.loads(p.read_text())["report"] for p in sorted(DATA.glob(pattern))]


toy = load("toy*.json")
p192 = sorted(load("p192_*.json"), key=lambda r: r["base_abscissae"])


def save(fig, stem):
    fig.savefig(HERE / f"{stem}.svg", format="svg")
    fig.savefig(HERE / f"{stem}.png", format="png", dpi=150)
    plt.close(fig)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 4.2))
colors = {16: "#1f5fa8", 20: "#b45f06", 24: "#38761d"}
for bits in (16, 20, 24):
    rs = [r for r in toy if round(r["order_bits"]) == bits]
    ax1.errorbar(
        [r["exact_yield"] for r in rs],
        [r["mc_rate"] for r in rs],
        yerr=[
            [r["mc_rate"] - r["mc_ci95"][0] for r in rs],
            [r["mc_ci95"][1] - r["mc_rate"] for r in rs],
        ],
        fmt="o",
        ms=4,
        color=colors[bits],
        label=f"exact count, toy $N\\approx 2^{{{bits}}}$",
    )
    ax1.plot(
        [min(r["model_yield"], 2.0) for r in rs],
        [r["mc_rate"] for r in rs],
        "x",
        ms=5,
        color=colors[bits],
        alpha=0.6,
    )
ax1.plot([1e-3, 1], [1e-3, 1], "k:", lw=0.8)
ax1.plot([], [], "kx", alpha=0.6, label="collision-free model $2B^2/(N-1)$")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("predicted yield per trial")
ax1.set_ylabel("Monte Carlo rate (95% interval)")
ax1.set_title("Method check on toy primes: 17 of 17 exact counts inside the interval")
ax1.grid(alpha=0.25, which="both")
ax1.legend(loc="upper left")

xs = [math.log2(r["base_abscissae"]) for r in p192]
ys = [-r["exact_yield_log2"] for r in p192]
ax2.plot(xs, ys, "o-", color="#1f5fa8", label="exact: $-\\log_2$ yield per trial")
ax2.plot(xs, [-r["model_yield_log2"] for r in p192], "x", ms=8, color="black", label="collision-free model (identical)")
ax2.plot(xs, [r["trials_for_full_system_log2"] for r in p192], "s-", color="#b45f06", label="$\\log_2$ trials for $B+1$ relations")
ax2.axhline(p192[0]["rho_log2"], color="#38761d", ls="--", label="rho, $\\log_2\\sqrt{N}=96$")
mc = -math.log2(p192[0]["mc_ci95"][1])
ax2.axhline(mc, color="gray", ls=":", label=f"Monte Carlo bound, $2^{{-{mc:.1f}}}$ (2M trials)")
ax2.set_xlabel("$\\log_2$ of base abscissae $B$ (base points $2B$)")
ax2.set_ylabel("$\\log_2$ trials")
ax2.set_ylim(0, 195)
ax2.set_title("NIST P-192, small-x base, two-point sums")
ax2.grid(alpha=0.25)
ax2.legend(loc="center right", fontsize=8)
fig.suptitle(RUN, fontsize=8.5)
fig.tight_layout()
save(fig, "fig1_yield")
print("wrote fig1_yield")
