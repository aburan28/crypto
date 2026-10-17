"""Layout only: draws the outcomes recorded in ../../runs/<cell>.d<K>.<arm>.jsonl
as a grid (rows = draws, columns = arms): filled = the degree-4 closure refutes,
hollow = degree 4 saturates without refuting."""
import glob, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
arms = ["A0", "A1", "B", "C"]
labels = ["A0\nfix summand coord = 0\n(6,7,7)", "A1\nfix summand coord = 1\n(6 coset,7,7)",
          "B\nfix one coord per summand\n(6,6,6) positive control", "C\nfix one intermediate coord\n(7,7,7) null control"]
rows = {}
for p in sorted(glob.glob(os.path.join(here, "../../runs/m3-*.jsonl"))):
    cell, d, arm = os.path.basename(p)[:-6].split(".")
    v = json.loads(open(p).read().splitlines()[0])
    rows.setdefault((cell, int(d[1:])), {})[arm] = v["closure_degree"] is not None
keys = sorted(rows)
fig, ax = plt.subplots(figsize=(9, 5.2))
for i, k in enumerate(keys):
    for j, a in enumerate(arms):
        ref = rows[k][a]
        ax.scatter(j, i, s=260, marker="o", facecolors="tab:blue" if ref else "white",
                   edgecolors="tab:blue", linewidths=2)
ax.set_xticks(range(4)); ax.set_xticklabels(labels, fontsize=7.5)
ax.set_yticks(range(len(keys)))
ax.set_yticklabels([f"{c.replace('m3-','m3 (').replace('-7',',7)').replace('-',',')} draw {d}" for c, d in keys], fontsize=8)
ax.invert_yaxis(); ax.set_xlim(-0.6, 3.6)
counts = [sum(rows[k][a] for k in keys) for a in arms]
for j, c in enumerate(counts):
    ax.annotate(f"{c}/8 refute", (j, len(keys) - 0.35), ha="center", fontsize=8)
ax.scatter([], [], s=90, facecolors="tab:blue", edgecolors="tab:blue", label="degree-4 closure refutes")
ax.scatter([], [], s=90, facecolors="white", edgecolors="tab:blue", label="degree 4 saturates without 1")
ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
ax.set_title("One fixed variable on the m = 3, ℓ = 7 draws (unfixed: none refute at 4)\n"
             "research/dreg_closure_specialise_20261017, run 2026-10-17T02:11Z", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(here, "specialise.svg")); fig.savefig(os.path.join(here, "specialise.png"), dpi=150)
