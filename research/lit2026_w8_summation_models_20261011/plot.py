"""Plot recorded W8 values (layout only: every number is read from the raw JSON
written by examples/lit2026_w8_summation_models.rs; nothing is recomputed
except medians already reported there and their arrangement).

usage: python -I plot.py <study dir>
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt

root = Path(sys.argv[1])
run = root / "reports" / "run-20261011T2017Z" / "raw"
speed = root / "reports" / "run-20261011T2131Z"
figs = run.parent  # diagnostic run report directory
sfigs = speed  # speed-session run report directory
plt.rcParams.update({"font.size": 8, "svg.fonttype": "none"})
COL = {"sym": "#4c72b0", "eg": "#dd8452", "u-plain": "#55a868", "x-direct": "#8172b3",
       "x-frob-aux": "#c44e52", "e-split-repo": "#4c72b0", "e-split-convolution": "#c44e52",
       "weierstrass-x": "#4c72b0", "wu-kummer-ratio": "#dd8452", "wu-u": "#55a868", "binary-edwards-w": "#c44e52"}


def load(p):
    return json.loads(Path(p).read_text())


def bars(ax, cells, arms_of, value, title, ylabel, log=True):
    names = [c for c, _ in cells]
    arms = []
    for _, c in cells:
        for a in arms_of(c):
            if a["arm"] not in arms:
                arms.append(a["arm"])
    w = 0.8 / max(1, len(arms))
    for j, arm in enumerate(arms):
        xs, ys = [], []
        for i, (_, c) in enumerate(cells):
            for a in arms_of(c):
                if a["arm"] == arm and value(a) is not None:
                    xs.append(i + (j - (len(arms) - 1) / 2) * w)
                    ys.append(value(a))
        ax.bar(xs, ys, w, label=arm, color=COL.get(arm))
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=20, ha="right")
    if log:
        ax.set_yscale("log")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=7)


# 1. E/G check: median word ops per arm, every eg cell.
if (run / "eg.json").exists():
    eg = load(run / "eg.json")
else:  # the last cell did not finish; plot the cells recorded by the checkpoint
    eg = {"result": load(run / "eg.json.partial")}
cells = [(c["cell"], c) for c in eg["result"]["cells"] if "arms" in c]
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
bars(axes[0], cells, lambda c: c["arms"], lambda a: max(a["median_word_ops"], 1),
     "E/G check: median F4 word ops per decomposition solve", "word XORs (log)")
bars(axes[1], cells, lambda c: c["arms"], lambda a: a["system"]["n_vars"],
     "Boolean unknowns of the descended system", "variables", log=False)
fig.text(0.01, 0.01, "Run run-20261011T2017Z (eg.json): 8 targets x 5 interleaved rounds per cell, MatrixF4 cap 4, "
         "KIC_F4_MAX_DEGREE_ONLY=0. Same factor base F_u for sym/eg/u-plain; x-direct on the Frobenius x-base.", fontsize=6)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(figs / "eg_arms.svg")

# 2. Interleaved speed session: original / previous / candidate.
cmp = load(speed / "compare.json")["result"]["cells"]
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
names = list(cmp.keys())
for k, (num, den, lab, col) in enumerate([("orig/sym", "cand/eg", "original sym / candidate eg", "#4c72b0"),
                                          ("cand/sym", "cand/eg", "previous sym / candidate eg", "#dd8452")]):
    for i, name in enumerate(names):
        rs = [r for r in cmp[name][0]["ratios"] if r["numerator"] == num and r["denominator"] == den]
        if not rs:
            continue
        r = rs[0]
        x = i + (k - 0.5) * 0.25
        axes[0].scatter([x] * len(r["per_round_ratio"]), r["per_round_ratio"], s=10, color=col, alpha=0.6,
                        label=lab if i == 0 else None)
        axes[0].scatter([x], [r["median_ratio_of_median_ms"]], marker="_", s=200, color=col)
        axes[1].bar(x, r["word_ops_ratio"], 0.25, color=col, label=lab if i == 0 else None)
for ax in axes:
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=20, ha="right")
    ax.axhline(1.0, color="#999", lw=0.6)
    ax.legend(fontsize=7)
axes[0].set_title("Wall time ratio per round (dots) and ratio of medians (bar), L0")
axes[0].set_ylabel("ratio of per-round total ms")
axes[1].set_title("Deterministic F4 word-op ratio")
axes[1].set_ylabel("ratio")
fig.text(0.01, 0.01, "Run run-20261011T2131Z (compare.json): original = harness on the 4b6e09b78 library, sym arm; "
         "candidate binary runs sym (previous) and eg in-process; 5 rounds ABAB, 8 targets per cell.", fontsize=6)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(sfigs / "eg_speed_session.svg")

# 3. Wu / binary Edwards models.
wu = load(run / "wu.json")
wcells = [(c["curve"], c) for c in wu["result"]["cells"] if "arms" in c]
fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
bars(axes[0], wcells[:1], lambda c: c["arms"], lambda a: a["poly"]["terms"],
     "S3 terms in each coordinate (" + wcells[0][0] + ")", "terms", log=False)
bars(axes[1], wcells[:1], lambda c: c["arms"], lambda a: a["poly"]["total_degree"],
     "S3 total degree", "degree", log=False)
bars(axes[2], wcells, lambda c: c["arms"], lambda a: max(a["median_word_ops"], 1),
     "median F4 word ops (descended S3, m = 2)", "word XORs (log)")
fig.text(0.01, 0.01, "Run run-20261011T2017Z (wu.json): factor base c(P) in <1, z, ..., z^(l-1)> for each coordinate c; "
         "8 targets x 5 rounds; MatrixF4 cap 4.", fontsize=6)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(figs / "wu_models.svg")

# 4. Frobenius auxiliary variable and convolution splitting.
fr = load(run / "frob.json")
fri = load(run / "frob-inherited.json")
lu = load(run / "lucas.json")
fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
bars(axes[0], [(c["cell"], c) for c in fr["result"]["cells"]], lambda c: c["arms"],
     lambda a: max(a["median_word_ops"], 1), "Frobenius aux: MatrixF4", "word XORs (log)")
bars(axes[1], [(c["cell"], c) for c in fri["result"]["cells"]], lambda c: c["arms"],
     lambda a: max(a["median_word_ops"], 1), "Frobenius aux: InheritedF4 (linear elimination)", "word XORs (log)")
bars(axes[2], [(c["cell"].replace("S4 symmetrised, b=1, ", ""), c) for c in lu["result"]["cells"]], lambda c: c["arms"],
     lambda a: max(a["median_word_ops"], 1), "e_i splitting: repo vs convolution", "word XORs (log)")
fig.text(0.01, 0.01, "Run run-20261011T2017Z (frob.json, frob-inherited.json, lucas.json); 8 targets x 5 rounds; "
         "root sets identical between arms in every cell.", fontsize=6)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(figs / "frob_lucas.svg")
print("wrote", sorted(str(p) for p in list(figs.glob("*.svg")) + list(sfigs.glob("*.svg"))))
