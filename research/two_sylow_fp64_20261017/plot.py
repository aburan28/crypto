"""Lay out the recorded results of the two_sylow_search timing (SVG + 150 dpi PNG).

Plotting only: every value drawn is read from the files `hwcost measure` and
`hwcost compare` wrote; nothing is computed here (AGENTS.md "Implementation
language").

    python plot.py RUNS_DIR OUT_DIR
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RUNS = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 9, "svg.fonttype": "none"})
RUN_NAMES = ["base-1", "cand-1", "base-2", "cand-2"]
STYLE = {"base-1": ("#1b6ca8", "o"), "base-2": ("#1b6ca8", "s"), "cand-1": ("#d1495b", "o"), "cand-2": ("#d1495b", "s")}
CLASS_COLOR = {"candidate faster": "#2e8540", "candidate slower": "#d1495b", "no detected effect": "#888888"}


def save(fig, name):
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def read_csv(p):
    with open(p) as f:
        return list(csv.DictReader(f))


def started(run):
    p = RUNS / run / "started_utc.txt"
    return p.read_text().strip() if p.exists() else "?"


runs = [r for r in RUN_NAMES if (RUNS / r / "samples.csv").exists()]

# Every tsy sample, per arm and run.
if runs:
    samples = {r: [s for s in read_csv(RUNS / r / "samples.csv") if s["group"] == "tsy"] for r in runs}
    arms = sorted({s["arm"] for r in runs for s in samples[r]})
    fig, ax = plt.subplots(figsize=(11, 4.2))
    width = 0.18
    for k, r in enumerate(runs):
        c, mk = STYLE[r]
        for i, a in enumerate(arms):
            vs = [float(s["value"]) / 1000 for s in samples[r] if s["arm"] == a]
            x = i + (k - 1.5) * width
            ax.plot([x] * len(vs), vs, mk, color=c, mfc="none", ms=4, label=r if i == 0 else None)
            if vs:
                med = sorted(vs)[len(vs) // 2]
                ax.plot([x - 0.07, x + 0.07], [med, med], color=c, lw=2)
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels([a.replace("tsy.", "") for a in arms], fontsize=8)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("µs per draw (markers: samples; bar: median)")
    ax.legend(fontsize=7, ncol=4)
    ax.grid(axis="y", alpha=0.3)
    ax.set_title("two_sylow_search Halver::run_draw on fp_backend::Fp64 (p = 2^63 − 25, the largest prime ≡ 3 mod 4 below 2^63), baseline (blue) vs candidate (red)\nruns " + ", ".join(f"{r} {started(r)}" for r in runs), fontsize=8)
    save(fig, "fig_tsy")

# Ratios with drift margins: tsy arms and the two kernel controls at B ≥ 4096.
cmp_path = RUNS / "compare" / "compare.csv"
if cmp_path.exists():
    rows = [r for r in read_csv(cmp_path) if r["arm"].startswith("tsy.") or (r["arm"].startswith("kernel.") and int(r["arm"].split(".")[2]) >= 4096)]
    fig, ax = plt.subplots(figsize=(11, 4))
    for i, r in enumerate(rows):
        m = float(r["margin"])
        ax.plot([i, i], [1 - m, 1 + m], color="#dddddd", lw=4, zorder=0)
        ax.plot(i, float(r["r1"]), "o", color=CLASS_COLOR[r["class"]], ms=5)
        ax.plot(i, float(r["r2"]), "s", mfc="none", color=CLASS_COLOR[r["class"]], ms=5)
    ax.axhline(1, color="black", lw=0.6)
    ax.axhline(1.10, color="#d1495b", lw=0.8, ls=":", label="S1 bound 1.10 (tsy arms)")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([r["arm"].replace("kernel.", "k.") for r in rows], rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("candidate / baseline")
    ax.legend(fontsize=7)
    ax.grid(axis="y", alpha=0.3)
    ax.set_title("r1 (filled), r2 (open); grey = ±max(2·drift, 2%); green faster, red slower, grey no detected effect.\nk.fp64.* is the positive control (known gain), k.w1.* the negative control (unchanged code)", fontsize=8)
    save(fig, "fig_ratios")

print("wrote", sorted(q.name for q in OUT.glob("fig_*")))
