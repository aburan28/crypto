"""Cold-time ratio per target size for several arm comparisons: each pair
(open: a contended run) and the clean pairs' geometric mean with its 95% t
interval from the summary file. Lays out pairs-*.tsv and summary-*.tsv;
computes nothing new.
Usage: arms_plot.py OUT TITLE YLABEL FOOTER (LABEL PAIRS SUMMARY)..."""
import csv, sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

out, title, ylabel, footer = sys.argv[1:5]
specs = [sys.argv[i:i + 3] for i in range(5, len(sys.argv), 3)]
order = [("icv1-f2m53-tm56619371-dac20a85", "44.3"), ("icv1-f2m59-tm943548413-98844ecc", "44.5"), ("icv1-f2m61-t158598901-ab42b6c5", "47.2")]
colours = ["#1f5fa8", "#c2571a", "#4c9a2a"]
fig, ax = plt.subplots(figsize=(7.8, 4.4))
width = 0.8 / len(specs)
for j, (label, pf, sf) in enumerate(specs):
    pairs = list(csv.DictReader(open(pf), delimiter="\t"))
    summ = {r["curve"]: r for r in csv.DictReader(open(sf), delimiter="\t")}
    ratio_key = [k for k in pairs[0] if "_over_" in k][0]
    flags = [k for k in pairs[0] if k.endswith("_contended")]
    skey = [k for k in next(iter(summ.values())) if "_over_" in k][0]
    for i, (c, _) in enumerate(order):
        x0 = i - 0.4 + width * (j + 0.5)
        rows = [p for p in pairs if p["curve"] == c]
        for k, p in enumerate(rows):
            clean = all(p[f] == "false" for f in flags)
            ax.plot(x0 - width * 0.3 + width * 0.4 * k / max(1, len(rows) - 1), float(p[ratio_key]), "o", ms=2.5,
                    mfc=colours[j] if clean else "none", mec=colours[j], alpha=0.45)
        s = summ[c]
        g, lo, hi = float(s[skey]), float(s["lo95"]), float(s["hi95"])
        ax.errorbar(x0 + width * 0.3, g, yerr=[[g - lo], [hi - g]], fmt="D", color=colours[j], capsize=3, ms=5,
                    label=label if i == 0 else None)
ax.axhline(1.0, color="black", lw=0.8)
ax.axhline(1.03, color="#666666", lw=0.8, ls="--")
ax.set_xticks(range(3))
ax.set_xticklabels([f"{lr}\n{c}" for c, lr in order], fontsize=7)
ax.set_xlabel("log₂ r and ICV1 slug")
ax.set_ylabel(ylabel)
ax.set_title(title, fontsize=9.5)
ax.legend(fontsize=7.5, loc="lower left")
fig.text(0.01, 0.01, footer, fontsize=6.3)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=150)
