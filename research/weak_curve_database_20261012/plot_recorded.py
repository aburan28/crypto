"""Layout only: plot the Rust-generated fixed CLASS_SUMMARY.tsv for the PDF."""
import csv
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

with (ROOT / "CLASS_SUMMARY.tsv").open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
assert len(rows) == 10

fig, axes = plt.subplots(2, 1, figsize=(10, 8.0),
                         gridspec_kw={"height_ratios": [7, 3]})
groups = [
    (axes[0], rows[:7], "Exact two-branch ordinary trace classes"),
    (axes[1], rows[7:], "Cubic norm-one diagnostic; randomized counter"),
]
for ax, group, title in groups:
    primes = [int(r["p"]) for r in group]
    ordinary = [int(r["ordinary"]) for r in group]
    positive = [int(r["positive_ordinary"]) for r in group]
    zero = [int(r["zero_ordinary"]) for r in group]
    assert all(a + b == c for a, b, c in zip(positive, zero, ordinary))
    yes = [a / n for a, n in zip(positive, ordinary)]
    no = [b / n for b, n in zip(zero, ordinary)]
    y = list(range(len(group)))
    ax.barh(y, yes, color="#075985", label="support recorded")
    ax.barh(y, no, left=yes, color="#bfcbd6", label="zero row in named mechanism")
    ax.set_yticks(y, [f"p={p}" for p in primes])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Share of ordinary rows in this mechanism")
    ax.set_title(title, loc="left", fontsize=13)
    for index, (a, n) in enumerate(zip(positive, ordinary)):
        ax.text(1.01, index, f"{a:,}/{n:,}", va="center", fontsize=9)
    ax.grid(axis="x", alpha=0.2)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.948),
           ncol=2, fontsize=9)
fig.suptitle("Recorded class support by base prime", y=0.985,
             fontsize=16, fontweight="bold")
fig.subplots_adjust(left=0.11, right=0.82, top=0.88, bottom=0.12, hspace=0.45)
fig.text(0.02, 0.018,
         "Source: CLASS_SUMMARY.tsv. The top and bottom panels have different branch scopes and evidence grades.",
         fontsize=8)
fig.savefig(ROOT / "class_counts_plot.pdf", bbox_inches="tight")
fig.savefig(ROOT / "class_counts.png", dpi=150, bbox_inches="tight")
