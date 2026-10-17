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
from matplotlib.patches import Patch  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "legend.fontsize": 8.5})
HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
BIN, PRI = "#1f5fa8", "#b45f06"


def runs(pattern):
    return [(p.stem, json.loads(p.read_text())["report"]) for p in sorted(DATA.glob(pattern))]


def save(fig, stem):
    fig.savefig(HERE / f"{stem}.svg", format="svg")
    fig.savefig(HERE / f"{stem}.png", format="png", dpi=150)
    plt.close(fig)


def by_curve(pattern):
    out = {}
    for stem, r in runs(pattern):
        out[stem.rsplit("_", 1)[0].removeprefix("s3_")] = r
    return out


def order(d):
    return sorted(d, key=lambda n: (d[n]["field_bits"], n))


def color(r):
    return BIN if r["field_kind"] == "binary" else PRI


def legend_patches(ax, extra=()):
    h = [Patch(color=BIN, alpha=0.85, label="binary field"), Patch(color=PRI, alpha=0.85, label="prime field")]
    ax.legend(handles=list(extra) + h, loc="upper left")


# ── Figure 1: yield against rho, and the toy method check ─────────────
secg = by_curve("sec[pt]*_4096.json")
names = order(secg)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.8), gridspec_kw={"width_ratios": [1.5, 1]})
xs = range(len(names))
ax1.bar(xs, [secg[n]["trials_for_full_system_log2"] for n in names], color=[color(secg[n]) for n in names], alpha=0.85)
ax1.plot(xs, [secg[n]["rho_log2"] for n in names], "k_", ms=12, mew=2)
ax1.set_xticks(list(xs))
ax1.set_xticklabels(names, rotation=70, ha="right", fontsize=7.5)
ax1.set_ylabel("$\\log_2$ trials")
ax1.set_title(f"$B+1$ relations from two-point sums, $B=4096$ (bars) against rho (dashes), {len(names)} curves")
ax1.grid(alpha=0.25, axis="y")
legend_patches(ax1, [plt.Line2D([], [], color="k", marker="_", ls="", ms=10, mew=2, label="rho, $\\log_2\\sqrt{n}$")])

toy = runs("toy*.json")
style = {  # prefix -> (marker, color, label)
    "toyb": ("o", BIN, "binary"),
    "toyp": ("s", PRI, "prime"),
}
inside = 0
for stem, r in toy:
    lo, hi = r["mc_ci95"]
    inside += lo <= r["exact_yield"] <= hi
for prefix, (marker, col, label) in style.items():
    for h in (1, 2, 4):
        rs = [r for s, r in toy if s.startswith(prefix) and r["cofactor"] == h]
        if not rs:
            continue
        ax2.errorbar(
            [r["exact_yield"] for r in rs],
            [r["mc_rate"] for r in rs],
            yerr=[[r["mc_rate"] - r["mc_ci95"][0] for r in rs], [r["mc_ci95"][1] - r["mc_rate"] for r in rs]],
            fmt=marker,
            ms=4,
            color=col,
            alpha={1: 0.5, 2: 0.8, 4: 1.0}[h],
            mfc="none" if h == 4 else col,
            label=f"toy {label}, $h={h}$",
        )
ax2.plot([1e-3, 1], [1e-3, 1], "k:", lw=0.8)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("exact yield per trial")
ax2.set_ylabel("Monte Carlo rate (95% interval)")
ax2.set_title(f"Method check on toy curves: {inside} of {len(toy)} inside")
ax2.grid(alpha=0.25, which="both")
ax2.legend(loc="upper left")
fig.tight_layout()
save(fig, "fig1_secg")
print("fig1_secg:", len(names), "curves,", inside, "of", len(toy), "toy runs inside")

# ── Figure 2: what a trial costs by a Semaev S3 solve ─────────────────
s3 = by_curve("s3_sec[pt]*_4096.json")
s3n = order(s3)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.8))
xs = range(len(s3n))
ax1.bar(xs, [s3[n]["s3"]["s3_us_per_trial"] / 1000 for n in s3n], color=[color(s3[n]) for n in s3n], alpha=0.85)
ax1.set_yscale("log")
ax1.set_xticks(list(xs))
ax1.set_xticklabels(s3n, rotation=70, ha="right", fontsize=7.5)
ax1.set_ylabel("ms per target, one thread")
ax1.set_title("Deciding one target by a Semaev $S_3$ solve, $B=4096$")
ax1.grid(alpha=0.25, axis="y", which="both")
legend_patches(ax1)
ax2.bar(xs, [s3[n]["s3"]["bits_slower_than_rho"] for n in s3n], color=[color(s3[n]) for n in s3n], alpha=0.85)
ax2.set_xticks(list(xs))
ax2.set_xticklabels(s3n, rotation=70, ha="right", fontsize=7.5)
ax2.set_ylabel("bits (log$_2$ of the time ratio)")
ax2.set_title("$B+1$ relations by $S_3$ solves against rho, same machine")
ax2.grid(alpha=0.25, axis="y")
fig.tight_layout()
save(fig, "fig2_s3_cost")
print("fig2_s3_cost:", len(s3n), "curves")

# ── Figure 3: the raw small-x base ────────────────────────────────────
raw = {s.removeprefix("raw_").rsplit("_", 1)[0]: r for s, r in runs("raw_*.json")}
rn = sorted(raw, key=lambda n: (raw[n]["cofactor"], raw[n]["field_bits"]))
fig, ax = plt.subplots(figsize=(8.2, 4.2))
for i, n in enumerate(rn):
    r = raw[n]["raw_base"]
    ax.errorbar(
        i,
        r["fraction"],
        yerr=[[r["fraction"] - r["ci95"][0]], [r["ci95"][1] - r["fraction"]]],
        fmt="o",
        color=color(raw[n]),
        ms=5,
    )
    ax.plot([i - 0.3, i + 0.3], [r["one_over_h"]] * 2, "k-", lw=1.5)
ax.set_xticks(range(len(rn)))
ax.set_xticklabels([f"{n}\n$h={raw[n]['cofactor']}$" for n in rn], fontsize=7.5)
ax.set_ylabel("share of raw small-x pair sums in the subgroup")
ax.set_title("Raw small-x base: sampled share (dots, 95%) against $1/h$ (bars)")
ax.set_ylim(0, 0.6)
ax.grid(alpha=0.25, axis="y")
legend_patches(ax)
fig.tight_layout()
save(fig, "fig3_raw_base")
print("fig3_raw_base:", len(rn), "curves")
