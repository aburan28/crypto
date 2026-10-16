"""Lay out the classification recorded in papers.csv (counts only; no research computation)."""
import csv, collections, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open(sys.argv[1])))
impacts = ["Solving-degree analysis tool", "Constant factor", "Feasibility check (S5 oracle)",
           "Lower bound (negative)", "None for ECDLP", "Unsupported or unverified claim"]
colors = ["#2a6fdb", "#3aa76d", "#c98a1b", "#7a4cc2", "#9aa0a6", "#d0453b"]
comps = sorted({r["component"] for r in rows}, key=lambda c: -sum(r["component"] == c for r in rows))
cnt = collections.Counter((r["component"], r["impact"]) for r in rows)

fig, ax = plt.subplots(figsize=(9.5, 5.2))
left = [0] * len(comps)
for imp, col in zip(impacts, colors):
    vals = [cnt[(c, imp)] for c in comps]
    ax.barh(comps, vals, left=left, color=col, label=imp, edgecolor="white")
    left = [a + b for a, b in zip(left, vals)]
ax.invert_yaxis()
ax.set_xlabel("papers (count, after merging companion versions)")
ax.set_title(f"ECDLP-relevant preprints, 2024-10 to 2026-10, by repo component and expected impact (N={len(rows)})\n"
             "Zero items classified as exponent-level improvements for elliptic-curve DLP", fontsize=10)
ax.legend(fontsize=8, loc="lower right", frameon=False)
ax.xaxis.get_major_locator().set_params(integer=True)
fig.text(0.01, 0.01, "Source: research/literature_search_20261011/papers.csv; classification from abstracts (search round 2026-10-11).",
         fontsize=7, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(sys.argv[2])
