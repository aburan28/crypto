"""A18's graph: at 2^47.2, each arm's paired cold-time ratio against v8 on
the suite's recipe (per pair, and the geometric mean with its 95% t
interval), and each arm's geometric-mean phases. Lays out pairs-*.tsv,
summary-*.tsv and phases.tsv; computes nothing new."""
import csv, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

d, out = sys.argv[1], sys.argv[2]
arms = ["v8-c400", "v8-c470", "a17-c320", "a17-c400", "a17-c470"]
fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.5, 4.0), gridspec_kw={"width_ratios": [1.2, 1]})
for i, a in enumerate(arms):
    pairs = list(csv.DictReader(open(f"{d}/pairs-{a}.tsv"), delimiter="\t"))
    key = [k for k in pairs[0] if "_over_" in k][0]
    for k, p in enumerate(pairs):
        ax.plot(i - 0.2 + 0.05 * k, float(p[key]), "o", ms=3, color="#7f7f7f")
    s = next(csv.DictReader(open(f"{d}/summary-{a}.tsv"), delimiter="\t"))
    skey = [k for k in s if "_over_" in k][0]
    g, lo, hi = float(s[skey]), float(s["lo95"]), float(s["hi95"])
    ax.errorbar(i + 0.3, g, yerr=[[g - lo], [hi - g]], fmt="D", color="#1f5fa8", capsize=4, ms=5)
ax.axhline(1.0, color="black", lw=0.8)
ax.axhline(1.03, color="#666666", lw=0.8, ls="--")
ax.set_xticks(range(len(arms)))
ax.set_xticklabels(arms, fontsize=8)
ax.set_ylabel("cold time, v8-c320 / arm (above 1: arm faster)")
ax.set_title("paired against v8 on the suite's recipe", fontsize=9)
ph = list(csv.DictReader(open(f"{d}/phases.tsv"), delimiter="\t"))
names = [r["arm"] for r in ph]
bottom = [0.0] * len(ph)
for key, colour in (("collect_ms", "#1f5fa8"), ("build_ms", "#c2571a"), ("rest_ms", "#bbbbbb")):
    vals = [float(r[key]) for r in ph]
    bx.bar(range(len(ph)), vals, bottom=bottom, color=colour, label=key[:-3])
    bottom = [b + v for b, v in zip(bottom, vals)]
bx.set_xticks(range(len(ph)))
bx.set_xticklabels(names, fontsize=7, rotation=30)
bx.set_ylabel("ms, geometric mean over 8 processes")
bx.set_title("phases per arm", fontsize=9)
bx.legend(fontsize=7)
fig.suptitle("A18: the factor base widened at 2^47.2 (icv1-f2m61-t158598901-ab42b6c5), with and without A17's filter", fontsize=9.5)
fig.text(0.01, 0.01, "explorations/A18-wider-base-2p47-20261017; 48 isolated processes, none contended; dashed: 1.03", fontsize=6.5)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=150)
