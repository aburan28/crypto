"""A15's graph: each pair's cold-time ratio v8/a15 per size (open: a
contended run), with the clean pairs' geometric mean and 95% t interval
from summary.tsv. Lays out pairs.tsv and summary.tsv; computes nothing new."""
import csv, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

pairs = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
summ = {r["curve"]: r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
out = sys.argv[3]
order = [("icv1-f2m53-tm56619371-dac20a85", "44.3"), ("icv1-f2m59-tm943548413-98844ecc", "44.5"), ("icv1-f2m61-t158598901-ab42b6c5", "47.2")]
fig, ax = plt.subplots(figsize=(7.5, 4.2))
for i, (c, lr) in enumerate(order):
    for k, p in enumerate(r for r in pairs if r["curve"] == c):
        clean = p["v8_contended"] == "false" and p["a15_contended"] == "false"
        x = i - 0.25 + 0.04 * (k % 12)
        ax.plot(x, float(p["v8_over_a15"]), "o", ms=3.5, mfc="#7f7f7f" if clean else "none", mec="#7f7f7f")
    s = summ[c]
    g, lo, hi = float(s["v8_over_a15"]), float(s["lo95"]), float(s["hi95"])
    ax.errorbar(i + 0.32, g, yerr=[[g - lo], [hi - g]], fmt="D", color="#1f5fa8", capsize=4, ms=5,
                label="clean pairs: geometric mean, 95% t" if i == 0 else None)
ax.axhline(1.0, color="black", lw=0.8)
ax.set_xticks(range(3))
ax.set_xticklabels([f"{lr}\n{c}" for c, lr in order], fontsize=7)
ax.set_xlabel("log₂ r and ICV1 slug ")
ax.set_ylabel("cold time, v8 / v8 + slots hashed once")
ax.set_title("A15: the scan's filter slots hashed once a block, 6 rounds × 2 rows a size", fontsize=10)
ax.legend(loc="lower right", fontsize=8)
fig.text(0.01, 0.01, "research/ic_tool_program/explorations/A15-filter-slots-20261017/pairs.tsv; 72 isolated processes; open circles: a pair with a contended run", fontsize=6.5)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=150)
