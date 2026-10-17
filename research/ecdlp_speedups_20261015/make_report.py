import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

base = [65.82, 60.93, 57.02]
cand = [37.46, 32.56, 67.59]
j0_base = [14.795, 13.992, 11.588]
j0_cand = [13.813, 11.897, 11.941]

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

ax = axes[0]
ax.bar(["baseline", "candidate"], [sorted(base)[1], sorted(cand)[1]], color=["#888", "#2a7"])
for vals in (base, cand):
    ax.scatter([0 if vals is base else 1] * 3, vals, color="k", zorder=3, s=12)
ax.set_title("Koblitz orbit-DLP fast path, n=41\nwhole-process wall (median of 3, s)")
ax.set_ylabel("seconds")
ax.text(1, sorted(cand)[1] + 2, "1.63x", ha="center", fontweight="bold")

ax = axes[1]
ax.bar(["baseline", "candidate"], [sorted(j0_base)[1], sorted(j0_cand)[1]], color=["#888", "#2a7"])
for vals in (j0_base, j0_cand):
    ax.scatter([0 if vals is j0_base else 1] * 3, vals, color="k", zorder=3, s=12)
ax.set_title("Prime j=0 staged IC, 16-bit\nrelation-stage wall (median of 3, s)")
ax.set_ylabel("seconds")
ax.text(1, sorted(j0_cand)[1] + 0.4, "1.17x", ha="center", fontweight="bold")

fig.suptitle("ECDLP engineering speedups, 2026-10-15 (identical outputs, matched A/B)")
fig.tight_layout()
fig.savefig("benchmarks/20261015-ecdlp-speedups/ab-comparison.svg")
fig.savefig("benchmarks/20261015-ecdlp-speedups/ab-comparison.png", dpi=150)

with PdfPages("benchmarks/20261015-ecdlp-speedups/run-report.pdf") as pdf:
    pdf.savefig(fig)
    f2 = plt.figure(figsize=(8.5, 11))
    f2.suptitle("ECDLP speedups across approaches — run report (2026-10-15)")
    f2.clf()
    txt = open("research/ecdlp_speedups_20261015/README.md").read()
    f2 = plt.figure(figsize=(8.5, 11))
    f2.text(0.05, 0.97, "ECDLP speedups across approaches — run report (2026-10-15)\nSee research/ecdlp_speedups_20261015/README.md for the full report.\nKey results: Koblitz orbit fast path n=41 whole-process median 60.9 s -> 37.5 s (1.63x);\nj=0 relation stage 16-bit median 13.99 s -> 11.90 s (1.17x). Identical outputs, both arms.",
             va="top", fontsize=10)
    pdf.savefig(f2)
print("done")
