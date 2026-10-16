"""Lay out the recorded results of the Fp64 branch-free study (SVG + 150 dpi PNG).

Plotting only: every value drawn is read from the files `hwcost measure`,
`fit`, `compare`, `pair` and `disasm.sh` wrote; nothing is computed here
beyond unit conversions for display (AGENTS.md "Implementation language").

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
STYLE = {
    "base-1": ("#1b6ca8", "o", "-"),
    "base-2": ("#1b6ca8", "s", "--"),
    "cand-1": ("#d1495b", "o", "-"),
    "cand-2": ("#d1495b", "s", "--"),
}
FIELDS = [
    ("fp64", "fp_backend::Fp64 (2^63-25), under test"),
    ("pf64", "prime_fast::Fp64 (2^62-57), under test"),
    ("fp64cm", "hwcost Fp64Cm (2^63-25), control"),
    ("w1", "FpWide<1> (2^64-59), control"),
]
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


runs = [r for r in RUN_NAMES if (RUNS / r / "validation.csv").exists()]

# Kernel time per field, baseline vs candidate, all four runs.
if runs:
    val = {r: read_csv(RUNS / r / "validation.csv") for r in runs}
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, (f, title) in zip(axes.flatten(), FIELDS):
        for r in runs:
            rs = sorted((x for x in val[r] if x["field"] == f), key=lambda x: int(x["batch"]))
            c, mk, ls = STYLE[r]
            ax.plot([int(x["batch"]) for x in rs], [float(x["measured_ns"]) for x in rs], ls, marker=mk, ms=3, color=c, label=r)
        ax.axvspan(4096, 2**22 * 1.3, color="#f0f0f0", zorder=0)
        ax.set_xscale("log", base=2)
        ax.set_ylim(bottom=0)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("walks per batch B")
        ax.set_ylabel("ns per rho iteration (median of 7)")
        ax.grid(alpha=0.3)
    axes[0][0].legend(fontsize=7)
    fig.suptitle("hwcost rho kernel, baseline (blue) vs candidate (red) binaries, runs " + ", ".join(f"{r} {started(r)}" for r in runs) + "\ngrey: B ≥ 4096 (gate cells)", fontsize=8)
    fig.tight_layout()
    save(fig, "fig_kernel")

# Candidate / baseline ratio per arm (r1, r2) with the drift margin.
cmp_path = RUNS / "compare" / "compare.csv"
if cmp_path.exists():
    rows = read_csv(cmp_path)
    groups = [
        ("field operations", [r for r in rows if r["arm"].split(".")[0] in ("fp64", "pf64", "fp64cm", "w1")]),
        ("rho kernel", [r for r in rows if r["arm"].startswith("kernel.")]),
        ("thread scaling and production rho (pfrho)", [r for r in rows if r["arm"].startswith(("scaling.", "pfrho."))]),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(13, 9.5))
    for ax, (title, rs) in zip(axes, groups):
        for i, r in enumerate(rs):
            m = float(r["margin"])
            ax.plot([i, i], [1 - m, 1 + m], color="#dddddd", lw=4, zorder=0)
            ax.plot(i, float(r["r1"]), "o", color=CLASS_COLOR[r["class"]], ms=4)
            ax.plot(i, float(r["r2"]), "s", mfc="none", color=CLASS_COLOR[r["class"]], ms=4)
        ax.axhline(1, color="black", lw=0.6)
        ax.set_xticks(range(len(rs)))
        ax.set_xticklabels([r["arm"].replace("kernel.", "") for r in rs], rotation=90, fontsize=5.5)
        ax.set_ylabel("candidate / baseline")
        ax.set_title(f"{title}: r1 (filled), r2 (open); grey = ±max(2·drift, 2%); green faster, red slower, grey no detected effect", fontsize=8)
        ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    save(fig, "fig_ratios")

# Production rho: every recorded sample of pfrho, per run.
samples = {r: [s for s in read_csv(RUNS / r / "samples.csv") if s["group"] == "pfrho"] for r in runs if (RUNS / r / "samples.csv").exists()}
if any(samples.values()):
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    for i, r in enumerate(runs):
        vs = [float(s["value"]) for s in samples[r]]
        c, mk, _ = STYLE[r]
        ax.plot([i] * len(vs), vs, mk, color=c, mfc="none", ms=5)
        if vs:
            med = sorted(vs)[len(vs) // 2]
            ax.plot([i - 0.25, i + 0.25], [med, med], color=c, lw=2)
    arm = samples[runs[0]][0]["arm"] if samples[runs[0]] else "pfrho"
    ax.set_xticks(range(len(runs)))
    ax.set_xticklabels(runs)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("ns per walk step (one solve per sample)")
    ax.set_title(f"{arm}: prime_fast::rho_parallel, one thread, known-answer solves\n(markers: samples; bar: median of the run)", fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    save(fig, "fig_pfrho")

# Disassembly counts.
dis = RUNS / "disasm" / "counts.csv"
if dis.exists():
    cnt = read_csv(dis)
    fig, ax = plt.subplots(figsize=(8, 3.6))
    w = 0.25
    xs = range(len(cnt))
    for k, (col, c) in enumerate([("jcc", "#d1495b"), ("cmov", "#2e8540"), ("mulx", "#e1a100")]):
        ax.bar([i + (k - 1) * w for i in xs], [int(x[col]) for x in cnt], w, color=c, label=col)
    for i, x in enumerate(cnt):
        ax.text(i, max(int(x["jcc"]), int(x["cmov"])) + 1, f"{x['instructions']} ins", ha="center", fontsize=6.5)
    ax.set_xticks(list(xs))
    ax.set_xticklabels([f"{x['label']}\n{x['function']}" for x in cnt], fontsize=7)
    ax.set_ylim(0, 1.2 * max(int(x[k]) for x in cnt for k in ("jcc", "cmov")) + 3)
    ax.set_ylabel("instructions in the function")
    ax.set_title("objdump -d of the non-inlined kernel step (runs/disasm)", fontsize=8)
    ax.legend(fontsize=7)
    ax.grid(axis="y", alpha=0.3)
    save(fig, "fig_disasm")

print("wrote", sorted(q.name for q in OUT.glob("fig_*")))
