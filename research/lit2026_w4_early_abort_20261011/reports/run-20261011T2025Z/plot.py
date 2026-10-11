"""Plot recorded values of run-20261011T2025Z (lit2026 W4). Reads summary.json
(written by `f4_refutation_degree_probe --summarize raw`) and raw/time_*.jsonl;
computes nothing beyond picking values out. Writes degrees.svg and arms.svg."""
import json
import pathlib

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).parent
S = json.loads((HERE / "summary.json").read_text())
CELLS = [
    ("k1n17_dim6", "K1/2^17 dim 6\n8 attempts"),
    ("k1n23_dim7", "K1/2^23 dim 7\n4 attempts"),
    ("k0n23_dim8", "K0/2^23 dim 8\n4 attempts"),
    ("k0n31_dim9", "K0/2^31 dim 9\n1 attempt, budget 2000"),
]
plt.rcParams.update({"font.size": 9, "svg.fonttype": "none"})

# --- degrees.svg: where nodes are decided (f4 arm, mode 0) and stage-unit share
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
cats = [
    ("refuted, 1 first at degree 2", "#1f77b4"),
    ("refuted, 1 first at degree 3", "#9ecae1"),
    ("forced, first at degree 3", "#ff7f0e"),
    ("stalled through degree 3 (split)", "#bbbbbb"),
]
x = range(len(CELLS))
bottoms = [0] * len(CELLS)
for ci, (label, colour) in enumerate(cats):
    vals = []
    for key, _ in CELLS:
        t = S["traces"][f"trace_{key}_mode0"]["f4"]
        rf, ri = t["refuted_first_degree_fresh"], t["refuted_first_degree_inherited"]
        vals.append(
            [rf[2] + ri[2], rf[3] + ri[3], t["forced_first_degree"][3], t["stalled_from_scratch"]][ci]
        )
    ax1.bar(x, vals, bottom=bottoms, color=colour, label=label)
    bottoms = [b + v for b, v in zip(bottoms, vals)]
ax1.set_xticks(list(x), [c[1] for c in CELLS])
ax1.set_ylabel("F4 reductions (nodes), summed over attempts")
ax1.set_title("First from-scratch degree deciding each node\n(inherited F4 arm, KIC_F4_EARLY_ABORT=0)")
ax1.legend(fontsize=7, loc="upper left")

ref = []
spl = []
for key, _ in CELLS:
    t = S["traces"][f"trace_{key}_mode0"]["f4"]
    ref.append(t["word_ops_refuted_nodes"] / 1e6)
    spl.append(t["word_ops_split_nodes"] / 1e6)
ax2.bar(x, ref, color="#1f77b4", label="nodes that end refuted")
ax2.bar(x, spl, bottom=ref, color="#bbbbbb", label="nodes that end in a split")
for i, (key, _) in enumerate(CELLS):
    share = S["traces"][f"trace_{key}_mode0"]["f4"]["refuted_share_of_node_word_ops"]
    ax2.text(i, ref[i] + spl[i], f"refuted {100 * share:.1f}%", ha="center", va="bottom", fontsize=7)
ax2.set_xticks(list(x), [c[1] for c in CELLS])
ax2.set_ylabel("stage-unit word operations (millions)")
ax2.set_yscale("log")
ax2.set_title("Where the stage unit goes, by node outcome\n(same runs; children excluded)")
ax2.legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(HERE / "degrees.svg")

# --- arms.svg: word_xors ratios (deterministic) and wall per round
TCELLS = [("k1n17", "K1/2^17 dim 6"), ("k0n23", "K0/2^23 dim 8"), ("k0n31", "K0/2^31 dim 9")]
ARMS = [("orig", "original 4b6e09b78", "#555555"), ("prev", "previous (EA=0)", "#2ca02c"),
        ("cand1", "candidate EA=1", "#1f77b4"), ("cand2", "candidate EA=2", "#d62728")]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
w = 0.18
labels = []
pos = 0
for ci, (cell, cname) in enumerate(TCELLS):
    for ei, eng in enumerate(["f4", "f6p"]):
        arms = S["timing"][cell][eng]
        for ai, (arm, alabel, colour) in enumerate(ARMS):
            if arm == "orig":
                val = 1.0
            else:
                val = arms[arm]["word_xors_ratio_orig_over_this"]
            ax.bar(pos + (ai - 1.5) * w, val, w, color=colour, label=alabel if pos == 0 else None)
        labels.append((pos, f"{cname}\n{eng}"))
        pos += 1
ax.set_xticks([p for p, _ in labels], [l for _, l in labels], fontsize=7, rotation=35, ha="right")
ax.set_ylim(0.995, 1.005)
ax.axhline(1.0, color="black", lw=0.6)
ax.set_ylabel("word_xors(original) / word_xors(arm)")
ax.set_title("Stage unit, original over arm (deterministic,\nidentical in all 5 rounds)")
ax.legend(fontsize=7, loc="lower left")

ax = axes[1]
pos = 0
for ci, (cell, cname) in enumerate(TCELLS):
    for eng in ["f4", "f6p"]:
        arms = S["timing"][cell][eng]
        for ai, (arm, alabel, colour) in enumerate(ARMS):
            ms = arms[arm]["total_ms"]
            xs = [pos + (ai - 1.5) * w] * len(ms)
            ax.scatter(xs, ms, s=8, color=colour, label=alabel if pos == 0 else None)
            ax.plot([xs[0] - w / 3, xs[0] + w / 3], [arms[arm]["median_ms"]] * 2, color=colour, lw=2)
        pos += 1
ax.set_xticks([p for p, _ in labels], [l for _, l in labels], fontsize=7, rotation=35, ha="right")
ax.set_yscale("log")
ax.set_ylabel("probe total_ms (dots: 5 rounds; bar: median)")
ax.set_title("Wall, L0, load average 63.6-110.9\n(shared M4 Pro; spread dominates)")
ax.legend(fontsize=7, loc="lower left")
fig.tight_layout()
fig.savefig(HERE / "arms.svg")
