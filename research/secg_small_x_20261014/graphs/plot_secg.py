"""Lay out the recorded values of the SECG small-x yield runs (data/*.json).

Plotting only: every point is a field of a run's `report` object, written by
`ic prime`. Base-2 logarithms are the report's own `*_log2` fields.

    python3 graphs/plot_secg.py          # from this study directory
"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "legend.fontsize": 8.5})
HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
RUN = "run 2026-10-14, seed 7, 1M trials per SECG run, 3M per toy run"


def runs(pattern):
    out = []
    for p in sorted(DATA.glob(pattern)):
        out.append((p.stem, json.loads(p.read_text())["report"]))
    return out


def save(fig, stem):
    fig.savefig(HERE / f"{stem}.svg", format="svg")
    fig.savefig(HERE / f"{stem}.png", format="png", dpi=150)
    plt.close(fig)


# One row per SECG curve, at B = 4096.
secg = {}
for stem, r in runs("sec[pt]*_4096.json"):
    secg[stem.rsplit("_", 1)[0]] = r
order = sorted(secg, key=lambda n: (secg[n]["field_bits"], n))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.6), gridspec_kw={"width_ratios": [1.35, 1]})
xs = range(len(order))
colors = ["#1f5fa8" if secg[n]["field_kind"] == "binary" else "#b45f06" for n in order]
ax1.bar(xs, [secg[n]["trials_for_full_system_log2"] for n in order], color=colors, alpha=0.85)
ax1.plot(xs, [secg[n]["rho_log2"] for n in order], "k_", ms=14, mew=2, label="rho, $\\log_2\\sqrt{n}$")
ax1.set_xticks(list(xs))
ax1.set_xticklabels(order, rotation=70, ha="right", fontsize=8)
ax1.set_ylabel("$\\log_2$ trials")
ax1.set_title("$B+1$ relations from two-point sums, $B=4096$ (bars) against rho (dashes)")
ax1.bar([], [], color="#1f5fa8", label="binary field")
ax1.bar([], [], color="#b45f06", label="prime field")
ax1.grid(alpha=0.25, axis="y")
ax1.legend(loc="upper left")

# Method check on toy binary curves: Monte Carlo rate against the exact count.
toy = runs("toyb*.json")
cols = {13: "#1f5fa8", 17: "#b45f06", 19: "#38761d"}
inside = 0
for deg in (13, 17, 19):
    rs = [r for s, r in toy if s.startswith(f"toyb{deg}_")]
    ax2.errorbar(
        [r["exact_yield"] for r in rs],
        [r["mc_rate"] for r in rs],
        yerr=[
            [r["mc_rate"] - r["mc_ci95"][0] for r in rs],
            [r["mc_ci95"][1] - r["mc_rate"] for r in rs],
        ],
        fmt="o",
        ms=4,
        color=cols[deg],
        label=f"exact count, toy $F_{{2^{{{deg}}}}}$",
    )
    ax2.plot(
        [min(r["model_yield"], 2.0) for r in rs],
        [r["mc_rate"] for r in rs],
        "x",
        ms=5,
        color=cols[deg],
        alpha=0.6,
    )
for _, r in toy:
    lo, hi = r["mc_ci95"]
    inside += lo <= r["exact_yield"] <= hi
ax2.plot([1e-3, 1], [1e-3, 1], "k:", lw=0.8)
ax2.plot([], [], "kx", alpha=0.6, label="collision-free model $2B^2/(n-1)$")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("predicted yield per trial")
ax2.set_ylabel("Monte Carlo rate (95% interval)")
ax2.set_title(f"Method check on toy binary curves: {inside} of {len(toy)} inside")
ax2.grid(alpha=0.25, which="both")
ax2.legend(loc="upper left")
fig.suptitle(RUN, fontsize=8.5)
fig.tight_layout()
save(fig, "fig1_secg")
print("wrote fig1_secg:", len(order), "curves,", inside, "of", len(toy), "toy runs inside")
