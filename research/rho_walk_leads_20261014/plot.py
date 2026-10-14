"""Lay out recorded values as SVG/PNG graphs. Reads the summary.jsonl files
written by `rho_walk_leads summary`; places those numbers and computes nothing.

    python plot.py runs/<main> runs/<split> runs/<batch> reports
"""
import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

main_run, split_run, batch_run, out = map(pathlib.Path, sys.argv[1:5])
out.mkdir(parents=True, exist_ok=True)


def load(d):
    return [json.loads(l) for l in (d / "summary.jsonl").read_text().splitlines() if l.strip()]


rows, split, batch = load(main_run), load(split_run), load(batch_run)
plt.rcParams.update({"font.size": 8.5, "svg.fonttype": "none"})
BLUE, RED, GREY, ORANGE = "#2f6db5", "#c0392b", "#888888", "#d68910"

LABEL = {
    "pm_coin": "P±R, sign from hash of P (control)",
    "pm_score": "best of P+R, P−R",
    "pm_split": "best of P+R_i, P−R_j (split control)",
    "mm_score": "best of P−R_i, P−R_j",
    "two_score": "best of P+R_i, P+R_j",
    "best4_score": "best of P±R_i, P±R_j",
    "k2": "random function, best of 2",
    "k4": "random function, best of 4",
    "lazy_dp": "DP test on lazy x",
    "lazy_index": "table index on lazy x",
    "pair_table": "pair table R_i+R_j, one addition",
    "shared_index": "index shared across walks (lower bound)",
    "neg_canon": "negation map, fixed rule",
    "neg_hashrule": "negation map, rule per point",
    "neg_epoch": "negation map, rule per epoch",
}
PANELS = [
    ("e1", 36, rows, "Lead 1, one walk, all points stored, doubling escape (ref: plain walk)"),
    ("e1", 36, split, "Lead 1, split-sign control, same setup (separate run)"),
    ("e1shift", 36, rows, "Lead 1, one walk, shift escape (ref: plain walk)"),
    ("e1dp", 36, rows, "Lead 1, distinguished points, 64 walks (ref: plain walk)"),
    ("rf", 36, rows, "Positive control: abstract random function (ref: best of 1)"),
    ("e2", 36, rows, "Leads 2 and 3, distinguished points, 64 walks (ref: canonical)"),
    ("e4", 36, rows, "Lead 4, distinguished points, 64 walks, r = 1024 (ref: no negation)"),
]
heights = [len([r for r in src if r.get("of") == e and r["bits"] == b and "speedup_ref_over_arm" in r]) for e, b, src, _ in PANELS]
fig, axes = plt.subplots(len(PANELS), 1, figsize=(7.4, 1.0 + 0.36 * sum(heights) + 0.45 * len(PANELS)),
                         sharex=True, gridspec_kw={"height_ratios": heights})
for ax, (exp, bits, src, title) in zip(axes, PANELS):
    sel = [r for r in src if r.get("of") == exp and r["bits"] == bits and "speedup_ref_over_arm" in r]
    for y, r in enumerate(sel):
        lo, hi = r["ci95"]
        v = r["speedup_ref_over_arm"]
        ax.errorbar(v, y, xerr=[[v - lo], [hi - v]], fmt="o", color=BLUE, ms=3.5, capsize=2, lw=1)
        if "predicted" in r:
            ax.plot(r["predicted"], y, marker="|", color=RED, ms=11, mew=2)
    ax.set_yticks(range(len(sel)), [f'{LABEL.get(r["arm"].removeprefix("sel_"), r["arm"])}  n={r["paired_n"]}' for r in sel])
    ax.set_ylim(-0.6, len(sel) - 0.4)
    ax.axvline(1.0, color=GREY, lw=0.8)
    ax.set_title(title, loc="left", fontsize=8.5)
    ax.set_xscale("log")
    ax.set_xticks([0.07, 0.1, 0.2, 0.5, 0.7, 1, 1.15, 1.41, 2], ["0.07", "0.1", "0.2", "0.5", "0.7", "1", "1.15", "1.41", "2"])
axes[-1].set_xlabel("steps speedup = reference mean / arm mean, paired by trial, 95% bootstrap\n"
                    "blue: measured; red tick: derived prediction; N ≈ 2^36")
fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(out / f"speedups.{ext}", dpi=160)

# cost against batch size
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 2.9))
for bits, name, col in ((256, "P-256 (4 limbs)", BLUE), (64, "2^62−57 (1 limb)", ORANGE)):
    pts = sorted((r["batch_walks"], r["cost_ratio_vs_base"]) for r in batch if r["bits"] == bits and r["arm"] == "pm_score")
    a1.plot([m for m, _ in pts], [c for _, c in pts], "o-", color=col, ms=3.5, label=name)
    for arm, ls in (("base", "-"), ("pm_score", "--")):
        q = sorted((r["batch_walks"], r["median_ns_per_step"]) for r in batch if r["bits"] == bits and r["arm"] == arm)
        a2.plot([m for m, _ in q], [v for _, v in q], ls, marker="o", ms=3, color=col,
                label=f'{name}, {"plain step" if arm == "base" else "best of P±R"}')
a1.axhline(1.1547, color=RED, lw=1, ls=":")
a1.text(9, 1.1547 + 0.006, "break-even = predicted steps speedup 1.155", color=RED, fontsize=7)
a1.set_xscale("log", base=2)
a1.set_xlabel("walks per batched inversion")
a1.set_ylabel("ns per step, best of P±R / plain")
a1.legend(fontsize=7, frameon=False, loc="lower right")
a1.set_title("Cost ratio (median of 21 reps)", loc="left", fontsize=8.5)
a2.set_xscale("log", base=2)
a2.set_yscale("log")
a2.set_xlabel("walks per batched inversion")
a2.set_ylabel("ns per step (median)")
a2.legend(fontsize=6, frameon=False)
a2.set_title("Absolute cost per step", loc="left", fontsize=8.5)
fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(out / f"cost_vs_batch.{ext}", dpi=160)
print("wrote graphs to", out)
