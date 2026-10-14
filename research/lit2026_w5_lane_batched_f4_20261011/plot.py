"""Plot recorded values of a W5 run (layout only; computes nothing new).

Usage: python plot.py <run dir>
Reads <run dir>/summary.json (written by `f4_lane_probe summarize`) and the
raw diagnostic records raw/diag-*.jsonl, and writes arms.svg, lanes.svg,
costmodel.svg and lockstep.svg into <run dir>.
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

run = Path(sys.argv[1])
summary = json.loads((run / "summary.json").read_text())
rid = run.name


def records(name):
    out = []
    for line in (run / "raw" / name).read_text().splitlines():
        v = json.loads(line)
        if "uptime" not in v:
            out.append(v)
    return out


CELLS = ["K1_17", "K1_23", "K0_23", "K0_31"]
LABEL = {"K1_17": "K1/2^17 d6", "K1_23": "K1/2^23 d7", "K0_23": "K0/2^23 d8", "K0_31": "K0/2^31 d9"}
ARMS = ["prev", "gate", "lanes", "rref"]
COLORS = {"prev": "#7f7f7f", "gate": "#1f77b4", "lanes": "#d62728", "rref": "#ff7f0e", "gaterref": "#9467bd"}

# 1. arms: stage-unit ratio and wall ratio, original / arm
ARMS = [a for a in ["prev", "gate", "lanes", "rref", "gaterref"]
        if any(a in c["arms"] for c in summary["cells"])]
fig, axes = plt.subplots(3, 2, figsize=(11, 10.5))
for col, engine in enumerate(["f4", "f6p"]):
    rows = {c["cell"]: c for c in summary["cells"] if c["engine"] == engine}
    for row, (key, base, title) in enumerate(
        [("word_xors", "orig", "stage unit (word_xors), original / arm"),
         ("wall_round_ratios", "orig", "wall per round (L0), original / arm: median, min-max"),
         ("cpu_round_ratios", "prev", "process CPU per round (L0), previous / arm: median, min-max")]
    ):
        ax = axes[row][col]
        width = 0.8 / len(ARMS)
        for i, arm in enumerate(ARMS):
            if arm == base:
                continue
            xs, ys, lo, hi = [], [], [], []
            for j, cell in enumerate(CELLS):
                r = rows.get(cell, {}).get("ratios", {}).get(f"{base}/{arm}")
                if not r:
                    continue
                if key == "word_xors":
                    if r["word_xors"] is None:
                        continue
                    xs.append(j + (i - (len(ARMS) - 1) / 2) * width)
                    ys.append(r["word_xors"])
                    lo.append(0)
                    hi.append(0)
                else:
                    w = r[key]
                    if not w:
                        continue
                    xs.append(j + (i - (len(ARMS) - 1) / 2) * width)
                    ys.append(w["median"])
                    lo.append(w["median"] - w["min"])
                    hi.append(w["max"] - w["median"])
            ax.bar(xs, ys, width, color=COLORS[arm], label=arm,
                   yerr=[lo, hi] if key != "word_xors" else None, capsize=2)
        ax.axhline(1.0, color="black", lw=0.8)
        if key == "wall_round_ratios":
            ax.set_yscale("log")
        ax.set_xticks(range(len(CELLS)))
        ax.set_xticklabels([LABEL[c] for c in CELLS], fontsize=8)
        ax.set_title(f"{engine}: {title}", fontsize=9)
        ax.set_ylabel("ratio (>1: arm does less)")
        if row == 0 and col == 0:
            ax.legend(fontsize=8)
fig.suptitle(f"W5 arms, {rid}, 5 interleaved rounds, reps 3 per round; host load "
             f"{summary['load1']['min']:.0f}-{summary['load1']['max']:.0f}", fontsize=10)
fig.tight_layout()
fig.savefig(run / "arms.svg")

if not (run / "raw" / "diag-lanes-echelon.jsonl").exists():
    print("arms plot written; no diagnostics in this run")
    sys.exit(0)

# 2. lanes: root phase, sequential root stage units vs lane-batched word equivalents
lanes = records("diag-lanes-echelon.jsonl")
lanes_rref = records("diag-lanes-rref.jsonl")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
for cell in sorted({r["cell"] for r in lanes}):
    rs = [r for r in lanes if r["cell"] == cell]
    rr = [r for r in lanes_rref if r["cell"] == cell]
    ks = [r["k"] for r in rs]
    axes[0].plot(ks, [r["lane_word_equivalent"] / r["seq_root_xors"] for r in rs], "o-", label=f"{cell} vs echelon")
    if rr:
        axes[0].plot([r["k"] for r in rr], [r["lane_word_equivalent"] / r["seq_root_xors"] for r in rr], "x--", label=f"{cell} vs RREF")
    axes[1].plot(ks, [r["lane_ms"][0] / r["seq_ms"][0] for r in rs], "o-", label=cell)
    axes[2].plot(ks, [r["shared_columns"] / r["pivot_columns"] for r in rs], "o-", label=f"{cell} shared")
    axes[2].plot(ks, [r["pivot_passes"] / r["pivot_columns"] / r["k"] for r in rs], "x:", label=f"{cell} passes/col/k")
for ax in axes:
    ax.set_xscale("log", base=2)
    ax.set_xlabel("lanes k (targets [1]G..[k]G)")
axes[0].set_yscale("log")
axes[0].set_title("root elimination: lane word-equivalents / sequential", fontsize=9)
axes[1].set_yscale("log")
axes[1].set_title("root phase wall (L0, median of 5): lanes / sequential echelon", fontsize=9)
axes[2].set_title("pivot pattern: columns shared by all lanes;\npivot passes per column per lane", fontsize=9)
for ax in axes:
    ax.axhline(1.0, color="black", lw=0.8)
    ax.legend(fontsize=6)
fig.suptitle(f"W5 lane-batched roots, {rid} (diag-lanes-*.jsonl)", fontsize=10)
fig.tight_layout()
fig.savefig(run / "lanes.svg")

# 3. cost model: measured / generic Boolean counts at degree 3
cm = [r for r in records("diag-costmodel.jsonl")]
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
cells = [r["cell"] for r in cm if r["degree"] == 3]
for ax, idx, what in [(axes[0], 0, "rows"), (axes[1], 1, "columns")]:
    width = 0.2
    for i, (key, label) in enumerate([
        ("descended_support_local", "descended, support-local (production)"),
        ("descended_occurring", "descended, occurring multipliers"),
        ("dense_support_local", "dense random, same shape"),
    ]):
        ys = []
        for r in cm:
            if r["degree"] != 3:
                continue
            v = r[key]
            ys.append(v[idx] / r["generic"][idx] if v else 0)
        ax.bar([j + (i - 1) * width for j in range(len(cells))], ys, width, label=label)
    ax.axhline(1.0, color="black", lw=0.8)
    ax.set_xticks(range(len(cells)))
    ax.set_xticklabels(cells, fontsize=8)
    ax.set_title(f"degree-3 root Macaulay {what}: measured / generic Boolean count", fontsize=9)
    ax.set_ylim(0, 1.4)
    ax.legend(fontsize=7, loc="upper left")
fig.suptitle(f"W5 F4 cost-model check, {rid} (diag-costmodel.jsonl)", fontsize=10)
fig.tight_layout()
fig.savefig(run / "costmodel.svg")

# 4. lockstep: share of tree work at nodes reached by two or more targets
ls = records("diag-lockstep.jsonl") + records("diag-lockstep-k0-31-4targets.jsonl")
ls = [r for r in ls if r["targets"] > 1]
fig, ax = plt.subplots(figsize=(8, 3.6))
labels = [f"{r['cell']} {r['arm']} ({r['targets']} t)" for r in ls]
ax.barh(range(len(ls)), [r["work_share_shared"] for r in ls], color="#2ca02c")
for i, r in enumerate(ls):
    ax.text(r["work_share_shared"] + 0.01, i, f"{r['mean_lanes_per_node']:.2f} lanes/node", va="center", fontsize=7)
ax.set_yticks(range(len(ls)))
ax.set_yticklabels(labels, fontsize=7)
ax.set_xlim(0, 1.25)
ax.set_xlabel("share of stage unit spent at nodes reached by >= 2 targets")
ax.set_title(f"W5 lockstep upper bound, {rid} (diag-lockstep*.jsonl)", fontsize=9)
fig.tight_layout()
fig.savefig(run / "lockstep.svg")
print("plots written")
