import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

base = [83.96, 76.38, 41.94]
cand = [40.76, 23.40, 22.42]

fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(["baseline (main 9f7fb66e6b)", "candidate (u64 RootBatch)"],
       [sorted(base)[1], sorted(cand)[1]], color=["#888", "#2a7"])
ax.scatter([0] * 3, base, color="k", zorder=3, s=14, label="reps")
ax.scatter([1] * 3, cand, color="k", zorder=3, s=14)
ax.set_title("Koblitz orbit-DLP fast path, u64 pipeline, n=41\nwhole-process wall (median of 3, s; shared host)")
ax.set_ylabel("seconds")
ax.text(1, sorted(cand)[1] + 3, "3.26x median\n(1.87x min-to-min)", ha="center", fontweight="bold")
ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig("benchmarks/20261015-u64-rootbatch/ab-comparison.svg")
fig.savefig("benchmarks/20261015-u64-rootbatch/ab-comparison.png", dpi=150)

with PdfPages("benchmarks/20261015-u64-rootbatch/run-report.pdf") as pdf:
    pdf.savefig(fig)
    f2 = plt.figure(figsize=(8.5, 11))
    f2.text(0.05, 0.95,
            "u64 S3 root-solve batching run report (2026-10-15)\n"
            "See research/u64_root_batch_20261015/README.md.\n"
            "Key result: n=41 whole-process median 76.4 s -> 23.4 s (3.26x; 1.87x min-to-min),\n"
            "identical recovered scalars, relations and probe counts on every rep (3/3).",
            va="top", fontsize=10)
    pdf.savefig(f2)
print("done")
