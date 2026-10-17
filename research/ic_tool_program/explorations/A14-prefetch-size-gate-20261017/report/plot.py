"""A14's graph: each pair's cold-time ratio v8/gate per size (open: a
contended run), with the clean pairs' geometric mean and 95% t interval
from summary.tsv. Lays out pairs.tsv and summary.tsv; computes nothing new."""
import csv, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

pairs = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
summ = {r["curve"]: r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
out = sys.argv[3]
order = [("icv1-f2m37-tm534059-32aad96b", "27.8"), ("icv1-f2m57-tm747311035-c1f545af", "38.0"), ("icv1-f2m53-tm56619371-dac20a85", "44.3")]
fig, ax = plt.subplots(figsize=(7.5, 4.2))
for i, (c, lr) in enumerate(order):
    for k, p in enumerate(r for r in pairs if r["curve"] == c):
        clean = p["v8_contended"] == "false" and p["gate_contended"] == "false"
        x = i - 0.25 + 0.04 * (k % 12)
        ax.plot(x, float(p["v8_over_gate"]), "o", ms=3.5, mfc="#7f7f7f" if clean else "none", mec="#7f7f7f")
    s = summ[c]
    g, lo, hi = float(s["v8_over_gate"]), float(s["lo95"]), float(s["hi95"])
    ax.errorbar(i + 0.32, g, yerr=[[g - lo], [hi - g]], fmt="D", color="#1f5fa8", capsize=4, ms=5,
                label="clean pairs: geometric mean, 95% t" if i == 0 else None)
ax.axhline(1.0, color="black", lw=0.8)
ax.axvline(1.5, color="#c2571a", lw=0.8, ls=":")
ax.text(1.52, ax.get_ylim()[1], " gate: 2^20 stored pairs", color="#c2571a", fontsize=7, va="top")
ax.set_xticks(range(3))
ax.set_xticklabels([f"{lr}\n{c}" for c, lr in order], fontsize=7)
ax.set_xlabel("log₂ r and ICV1 slug (left of the dotted line: the gated arm skips the prefetch)")
ax.set_ylabel("cold time, v8 / v8 + size gate")
ax.set_title("A14: R12's run prefetch gated on table size, 6 rounds × 2 rows a size", fontsize=10)
ax.legend(loc="lower right", fontsize=8)
fig.text(0.01, 0.01, "research/ic_tool_program/explorations/A14-prefetch-size-gate-20261017/pairs.tsv; 72 isolated processes; open circles: a pair with a contended run", fontsize=6.5)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=150)
