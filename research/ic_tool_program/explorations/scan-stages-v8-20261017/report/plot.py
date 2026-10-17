"""The scan's stages a scanned summand, v3 beside v8, per target size,
stacked. Lays out the two stages.tsv files; computes nothing new."""
import csv, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

v3 = {r["curve"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
v8 = {r["curve"]: r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
out = sys.argv[3]
order = [("icv1-f2m53-tm56619371-dac20a85", "44.3"), ("icv1-f2m59-tm943548413-98844ecc", "44.5"), ("icv1-f2m61-t158598901-ab42b6c5", "47.2")]
stages = [("subtract_ns", "subtract", "#1f5fa8"), ("key_ns", "key", "#c2571a"), ("filter_ns", "filter", "#4c9a2a"), ("admitted_ns", "admitted", "#8e5aa8")]
fig, ax = plt.subplots(figsize=(7.5, 4.2))
ticks, labels = [], []
for i, (c, lr) in enumerate(order):
    for j, (name, tab) in enumerate((("v3", v3), ("v8", v8))):
        x = i * 3 + j
        bottom = 0.0
        for key, label, colour in stages:
            v = float(tab[c][key])
            ax.bar(x, v, bottom=bottom, color=colour, width=0.8, label=label if (i, j) == (0, 0) else None)
            bottom += v
        ax.text(x, bottom + 0.4, f"{bottom:.1f}", ha="center", fontsize=7)
        ticks.append(x)
        labels.append(f"{name}\n{lr}")
ax.set_xticks(ticks)
ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel("ns a scanned summand (mean over clean processes)")
ax.set_xlabel("baseline and log₂ r (ICV1 slugs in the report)")
ax.set_title("The collection scan's stages, R04's probes: v3 (2026-10-06) and v8 (2026-10-17)", fontsize=9.5)
ax.legend(fontsize=8, loc="upper right")
fig.text(0.01, 0.01, "explorations/scan-stages-v3-20261006/stages.tsv (2-4 processes a size); scan-stages-v8-20261017/stages.tsv (7-8); v8 counts R12's run prefetch in 'admitted'", fontsize=6.3)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=150)
